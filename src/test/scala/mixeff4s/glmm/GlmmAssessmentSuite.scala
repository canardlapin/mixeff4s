package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.fixtures.GlmmFrames
import mixeff4s.formula.Formula
import mixeff4s.model.Family

class GlmmAssessmentSuite extends munit.FunSuite:
  private def context(
      formula: String,
      frame: ModelFrame,
      family: Family,
      options: GlmmOptions = GlmmOptions.jointLaplace
  ): GlmmFitContext =
    val design = Glmm.compile(formula, frame).toOption.get
    val response = GlmmResponse.prepare(family, Glmm.responseValues(design), options).toOption.get
    new GlmmFitContext(design, response, options, family, family.canonicalLink)

  test("mixed grouped-binomial rows and positive Poisson counts provide explicit sufficient witnesses"):
    val cbpp = GlmmFrames.cbpp
    val sizes = cbpp.numeric("size").get.map(_.toInt)
    val ys = cbpp.numeric("incidence").get.zip(sizes).map((k, n) => k / n.toDouble)
    val frame = cbpp.addNumeric("proportion", ys).toOption.get
    val c = context(
      "proportion ~ 1 + period + (1 | herd)",
      frame,
      Family.Binomial,
      GlmmOptions.jointAgq(5).copy(trials = Some(sizes))
    )
    GlmmAssessment.fixedEffectEvidence(c) match
      case FixedEffectRecession.Excluded("mixed_binomial_rows_span_retained_design", 4) => ()
      case other => fail(s"expected sufficient mixed-row rank witness, got $other")
    val data = GlmmFrames.synthetic_poisson_offset
    val formula = Formula.parse("count ~ 1 + x + (1 | group)").toOption.get
    val options = GlmmOptions.jointAgq(5).copy(offset = Some(data.numeric("exposure").get.map(math.log)))
    val fit = Glmm.fit(formula, data, Family.Poisson, options = options).fold(e => fail(e.message), identity)
    val assessment = GlmmAssessment.assess(fit)
    assert(assessment.eligibleForBootstrap, clues(assessment))
    assert(assessment.outer.exists(_.firstOrderPassed))
    assertEquals(
      assessment.fixedEffectRecession,
      FixedEffectRecession.Excluded("positive_count_rows_span_retained_design", 2)
    )
    assert(assessment.riskFlags.contains("finite_joint_mle_not_assessed"))

  test("pure groups are retained risk flags while their Gaussian-penalized modes are finite"):
    val successes = Vector(0, 6, 1, 5, 2, 4)
    val y = successes.flatMap(k => Vector.tabulate(6)(j => if j < k then 1.0 else 0.0))
    val frame = ModelFrame
      .of("y" -> ModelFrame.numeric(y), "g" -> ModelFrame.factor(Vector.tabulate(y.length)(i => s"g${i / 6}")))
      .toOption
      .get
    val fit = Glmm
      .fit(Formula.parse("y ~ 1 + (1 | g)").toOption.get, frame, Family.Bernoulli, options = GlmmOptions.jointAgq(9))
      .fold(e => fail(e.message), identity)
    val assessment = GlmmAssessment.assess(fit)
    assert(assessment.conditional.isRight)
    assert(assessment.riskFlags.exists(_.startsWith("pure_outcome_groups:g:2")))
    assessment.fixedEffectRecession match
      case FixedEffectRecession.Excluded(_, 1) => ()
      case other                               => fail(s"global overlap is not disproved by pure groups: $other")
    assert(assessment.eligibleForBootstrap, clues(assessment))

  test("separation and unavailable high-dimensional evidence are not mistaken for an overlap certificate"):
    val x = Vector(-3.0, -2.0, -1.0, 1.0, 2.0, 3.0)
    val y = x.map(v => if v > 0 then 1.0 else 0.0)
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(y),
        "x" -> ModelFrame.numeric(x),
        "z" -> ModelFrame.numeric(x.map(v => v * v)),
        "g" -> ModelFrame.factor(Vector("a", "b", "c", "a", "b", "c"))
      )
      .toOption
      .get
    GlmmAssessment.fixedEffectEvidence(context("y ~ 1 + x + (1 | g)", frame, Family.Bernoulli)) match
      case FixedEffectRecession.Detected(_, direction) =>
        x.zip(y)
          .foreach: (value, outcome) =>
            val signed = (if outcome == 1.0 then 1.0 else -1.0) * (direction(0) + direction(1) * value)
            assert(signed > 0.0)
      case other => fail(s"expected separating direction, got $other")
    GlmmAssessment.fixedEffectEvidence(context("y ~ 1 + x + z + (1 | g)", frame, Family.Bernoulli)) match
      case FixedEffectRecession.Unavailable("no_sufficient_fixed_effect_witness") => ()
      case other => fail(s"higher-dimensional absence of a finding must remain unavailable: $other")

  test("strict one-dimensional outcome hull overlap includes an interior singleton class"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(0.0, 1.0, 0.0, 0.0)),
        "x" -> ModelFrame.numeric(Vector(-1.0, 0.0, 1.0, 2.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "b", "b"))
      )
      .toOption
      .get
    assertEquals(
      GlmmAssessment.fixedEffectEvidence(context("y ~ 1 + x + (1 | g)", frame, Family.Bernoulli)),
      FixedEffectRecession.Excluded("strict_outcome_hull_overlap", 2)
    )

  test("numerical boundary fits and fast fits remain ineligible rather than gaining a certificate"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(0.0, 1.0, 1.0, 0.0, 0.0, 1.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c"))
      )
      .toOption
      .get
    val formula = Formula.parse("y ~ 1 + (1 | g)").toOption.get
    val fit = Glmm.fit(formula, frame, Family.Bernoulli, options = GlmmOptions.jointLaplace).toOption.get
    val assessment = GlmmAssessment.assess(fit)
    assert(assessment.numericalTerminationAccepted)
    assert(!assessment.eligibleForBootstrap)
    assert(assessment.refusalCodes.contains("bootstrap_covariance_boundary"))
    val fast = Glmm.fit(formula, frame, Family.Bernoulli).toOption.get
    assert(GlmmAssessment.assess(fast).refusalCodes.contains("bootstrap_requires_joint"))
    assert(
      GlmmAssessment
        .assess(fit, GlmmAssessmentOptions(gradientTolerance = Double.NaN))
        .refusalCodes
        .contains("glmm_assessment_invalid_options")
    )

  test("touching outcome hulls refuse and near-dependent mixed-row witnesses do not pass rank"):
    val touching = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(0.0, 0.0, 1.0, 1.0)),
        "x" -> ModelFrame.numeric(Vector(-1.0, 0.0, 0.0, 1.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "b", "b"))
      )
      .toOption
      .get
    GlmmAssessment.fixedEffectEvidence(context("y ~ 1 + x + (1 | g)", touching, Family.Bernoulli)) match
      case FixedEffectRecession.Detected(_, _) | FixedEffectRecession.Unavailable(_) => ()
      case other => fail(s"touching hulls cannot establish overlap: $other")
    val data = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(0.5, 0.5, 0.0, 1.0)),
        "x" -> ModelFrame.numeric(Vector(1.0, 1.0 + 1e-12, -1.0, 2.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "b", "b"))
      )
      .toOption
      .get
    val evidence = GlmmAssessment.fixedEffectEvidence(
      context(
        "y ~ 1 + x + (1 | g)",
        data,
        Family.Binomial,
        GlmmOptions.jointLaplace.copy(trials = Some(Vector.fill(4)(2)))
      )
    )
    evidence match
      case FixedEffectRecession.Excluded("mixed_binomial_rows_span_retained_design", _) =>
        fail("nearly collinear mixed rows must not be promoted as the spanning witness")
      case _ => ()
