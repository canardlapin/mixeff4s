package mixeff4s.stats

import mixeff4s.fixtures.Sleepstudy
import mixeff4s.comparison.{BootstrapFrozen, JsonValue}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.lmm.*
import mixeff4s.optimizer.TrustBqOptions

class BootstrapSuite extends munit.FunSuite:
  private def right[A](value: FitResult[A]): A = value.fold(e => fail(s"${e.code}: ${e.message}"), identity)
  private val formula = "reaction ~ 1 + days + (1 | subj)"
  private val strict = FitOptions.ml.copy(
    optimizer = LmmOptimizer.AnalyticGradient,
    optimizerControls = Some(
      TrustBqOptions(
        finalRadius = 1e-10,
        ftolAbs = 0,
        ftolRel = 0,
        stallFtolAbs = 0,
        stallFtolRel = 0,
        stallIterations = 20,
        gradientTolerance = 1e-9
      )
    )
  )
  private lazy val base = right(
    Lmm.fit(formula, Sleepstudy.frame, FitOptions.ml.copy(optimizer = LmmOptimizer.AnalyticGradient))
  )
  private val targets =
    Vector(BootstrapTarget.Coefficient("days"), BootstrapTarget.Contrast("two days", Map("days" -> 2.0)))
  private def options(n: Int) = BootstrapOptions(n, 20261026L, targets, fitOptions = Some(strict))

  test("explicit tight controls produce 39 independently assessed bootstrap refits and intervals"):
    val result = right(Bootstrap.parametric(base, options(39)))
    assertEquals(result.attempts.length, 39)
    assertEquals(
      result.qualifiedCount,
      39,
      clues(result.attempts.filterNot(_.qualified).map(a => (a.index, a.assessment, a.refit.left.toOption)))
    )
    assertEquals(result.failedOrUnqualifiedIndices, Vector.empty)
    result.attempts.foreach: attempt =>
      assertEquals(attempt.refit.toOption.get.options, strict)
      assertEquals(attempt.response.toOption.get, Simulation.response(base, attempt.seed).toOption.get)
      assert(attempt.assessment.get.gradient.isInstanceOf[GradientCheck.Pass])
    val first = right(result.intervals.head)
    val second = right(result.intervals(1))
    assert(first.lower < first.upper)
    assertEqualsDouble(second.lower, first.lower * 2, 1e-12)
    assertEqualsDouble(second.upper, first.upper * 2, 1e-12)
    println(
      s"BOOTSTRAP_WITNESS qualified=${result.qualifiedCount} lower=${first.lower} upper=${first.upper} level=${first.level}"
    )

  test("ML and REML originals can request tight controls without changing assessment thresholds"):
    for criterion <- Vector(Criterion.ML, Criterion.REML) do
      val fit = right(Lmm.fit(formula, Sleepstudy.frame, strict.copy(criterion = criterion)))
      assert(LmmAssessment.assess(fit).eligibleForBootstrap, clues(LmmAssessment.assess(fit)))
      val result =
        right(Bootstrap.parametric(fit, options(39).copy(fitOptions = Some(strict.copy(criterion = criterion)))))
      assertEquals(result.qualifiedCount, 39, clues(criterion, result.failedOrUnqualifiedIndices))
      assert(result.intervals.forall(_.isRight))

  test("actual generated responses reproduce independent lme4 refits and R percentile intervals"):
    val frozen = JsonValue.parse(BootstrapFrozen.json)
    val result = right(Bootstrap.parametric(base, options(39)))
    val rows = frozen.req("results").asArray
    assertEquals(rows.length, result.attempts.length)
    result.attempts
      .zip(rows)
      .foreach: (attempt, reference) =>
        assertEquals(attempt.index, reference.req("index").asInt)
        assertEquals(attempt.seed.toString, reference.req("seed").asString)
        val fitted = right(attempt.refit)
        fitted.beta.zip(reference.req("beta").asArray.map(_.asDouble)).foreach((a, b) => assertEqualsDouble(a, b, 1e-7))
        assertEqualsDouble(fitted.objective, reference.req("objective").asDouble, 1e-5)
        assertEqualsDouble(fitted.sigma, reference.req("sigma").asDouble, 1e-4)
    val interval = right(result.intervals.head)
    val expected = frozen.req("percentile").req("days").asArray.map(_.asDouble)
    assertEqualsDouble(interval.lower, expected(0), 1e-7)
    assertEqualsDouble(interval.upper, expected(1), 1e-7)

  test("injected middle failure retains later responses, seeds and cold-fit statistics"):
    val regular = right(Bootstrap.parametric(base, options(5)))
    val session = right(base.validatedContext.toOption.get.compiled.session(strict))
    val broken = right(Bootstrap.withRefitter(base, options(5)): (index, y) =>
      if index == 2 then Left(MixedModelError.Interrupted("injected replicate failure"))
      else session.refit(y, RefitStart.Cold))
    assertEquals(broken.attempts.length, 5)
    assertEquals(broken.attempts(2).refit.left.toOption.get.code, "interrupted")
    assert(broken.intervals.forall(_.left.toOption.exists(_.code == "bootstrap_replicate_unqualified")))
    for index <- Vector(0, 1, 3, 4) do
      assertEquals(broken.attempts(index).seed, regular.attempts(index).seed)
      assertEquals(broken.attempts(index).response, regular.attempts(index).response)
      assertEquals(broken.attempts(index).statistics, regular.attempts(index).statistics)
    val prefix = right(Bootstrap.parametric(base, options(3)))
    assertEquals(prefix.attempts.map(_.response), regular.attempts.take(3).map(_.response))

  test("unqualified replicas remain visible and authoritative fit budget overrides low-level controls"):
    val tiny =
      strict.copy(maxEvaluations = 1, optimizerControls = strict.optimizerControls.map(_.copy(maxEvaluations = 100000)))
    val result = right(Bootstrap.parametric(base, options(4).copy(fitOptions = Some(tiny))))
    assertEquals(result.attempts.length, 4)
    assertEquals(result.qualifiedCount, 0)
    assert(result.attempts.forall(_.refit.left.toOption.exists(_.code == "optimizer_budget_exhausted")))
    assert(result.intervals.forall(_.isLeft))

  test("percentile convention is R type 7 and confidence level is explicit"):
    assertEqualsDouble(Bootstrap.type7(Vector.tabulate(10)(_.toDouble), 0.025), 0.225, 1e-14)
    assertEqualsDouble(Bootstrap.type7(Vector.tabulate(10)(_.toDouble), 0.975), 8.775, 1e-14)
    assertEqualsDouble(Bootstrap.type7(Vector(-1e308, 1e308), 0.5), 0.0, 0.0)
    val narrow = right(Bootstrap.parametric(base, options(5).copy(level = 0.8)))
    assertEquals(right(narrow.intervals.head).level, 0.8)
    assert(Bootstrap.parametric(base, options(5).copy(level = Double.NaN)).isLeft)

  test("target validation and unqualified original fit refusals happen before sampling"):
    assertEquals(
      Bootstrap
        .parametric(base, options(2).copy(targets = Vector(BootstrapTarget.Coefficient("missing"))))
        .left
        .toOption
        .get
        .code,
      "invalid_argument"
    )
    val coarse =
      right(Lmm.fit(formula, Sleepstudy.frame, FitOptions.reml.copy(optimizer = LmmOptimizer.AnalyticGradient)))
    assertEquals(Bootstrap.parametric(coarse, options(2)).left.toOption.get.code, "bootstrap_unqualified_original")
    assertEquals(
      Bootstrap.parametric(base.copy(sigma = base.sigma + 1), options(2)).left.toOption.get.code,
      "fit_snapshot_mismatch"
    )
    assertEquals(Bootstrap.parametric(base, 10).left.toOption.get.code, "bootstrap_seed_required")

  test("invalid exposed numerical controls are refused by fits and sessions"):
    val plan = right(Lmm.prepare(formula, Sleepstudy.frame))
    val invalid = Vector(
      TrustBqOptions(ftolAbs = Double.NaN),
      TrustBqOptions(ftolRel = -0.1),
      TrustBqOptions(finalRadius = Double.PositiveInfinity),
      TrustBqOptions(stallIterations = 0)
    )
    invalid.foreach: controls =>
      val bad = strict.copy(optimizerControls = Some(controls))
      assertEquals(plan.fit(bad).left.toOption.get.code, "invalid_argument")
      assertEquals(plan.session(bad).left.toOption.get.code, "invalid_argument")
