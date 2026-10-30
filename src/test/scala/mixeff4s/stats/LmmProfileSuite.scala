package mixeff4s.stats

import mixeff4s.error.MixedModelError
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.lmm.{FitOptions, Lmm}
import mixeff4s.optimizer.TrustBqOptions

class LmmProfileSuite extends munit.FunSuite:
  private val tight = FitOptions.ml.copy(
    optimizerControls = Some(
      TrustBqOptions(
        finalRadius = 1e-9,
        ftolAbs = 0.0,
        ftolRel = 0.0,
        stallFtolAbs = 0.0,
        stallFtolRel = 0.0,
        stallIterations = 20
      )
    )
  )

  test("REML profile refuses unless an ML refit is explicit"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.reml).getOrElse(fail("fit"))
    LmmProfile.confint(fit, LmmProfileTarget.Coefficient("days")) match
      case Left(error: MixedModelError.InferenceUnavailable) =>
        assertEquals(error.code, "profile_reml_requires_ml_refit")
      case other => fail(s"expected explicit-ML refusal, got $other")

  test("unknown coefficient is a stable refusal after an ML fit"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml).getOrElse(fail("fit"))
    LmmProfile.confint(fit, LmmProfileTarget.Coefficient("dropped")) match
      case Left(error: MixedModelError.InferenceUnavailable) =>
        assertEquals(error.code, "profile_unknown_or_dropped_coefficient")
      case other => fail(s"expected retained-basis refusal, got $other")

  test("qualified ML coefficient profile retains finite likelihood endpoints"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, tight).getOrElse(fail("fit"))
    val options = LmmProfileOptions(
      maxEvaluations = 2000,
      maxBracketSteps = 10,
      maxBisectionSteps = 20,
      optimizerControls = tight.optimizerControls.get
    )
    val profile = LmmProfile.confint(fit, LmmProfileTarget.Coefficient("days"), options).getOrElse(fail("profile"))
    assert(profile.lower.value.exists(_ < profile.estimate), clues(profile))
    assert(profile.upper.value.exists(_ > profile.estimate), clues(profile))
    assert(
      profile.evaluations.forall(_.objective.toOption.forall(_ >= profile.baseline.objective - 1e-6)),
      clues(profile)
    )

  test("residual and scalar random-effect targets use the likelihood profiler"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, tight).getOrElse(fail("fit"))
    val options = LmmProfileOptions(
      maxEvaluations = 2000,
      maxBracketSteps = 10,
      maxBisectionSteps = 20,
      optimizerControls = tight.optimizerControls.get
    )
    val residual = LmmProfile.confint(fit, LmmProfileTarget.ResidualStdDev, options).getOrElse(fail("residual"))
    val sd = LmmProfile.confint(fit, LmmProfileTarget.RandomEffectStdDev("subj"), options).getOrElse(fail("SD"))
    val variance =
      LmmProfile.confint(fit, LmmProfileTarget.RandomEffectVariance("subj"), options).getOrElse(fail("variance"))
    assert(residual.estimate > 0.0, clues(residual))
    assertEqualsDouble(variance.estimate, sd.estimate * sd.estimate, 1e-10)
    assert(sd.evaluations.nonEmpty)
    assert(variance.evaluations.nonEmpty)

  test("explicit ML refit preserves a response-only REML snapshot rather than training y"):
    val plan = Lmm.prepare("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame).toOption.get
    val changed = Sleepstudy.reaction.map(_ * 2.0)
    val reml = plan.fit(changed, tight.copy(criterion = mixeff4s.lmm.Criterion.REML)).toOption.get
    val profile = LmmProfile
      .confint(
        reml,
        LmmProfileTarget.Coefficient("days"),
        LmmProfileOptions(refitWithML = true, maxEvaluations = 2000, optimizerControls = tight.optimizerControls.get)
      )
      .fold(e => fail(s"${e.code}: ${e.message}"), identity)
    assertEquals(profile.originalCriterion, mixeff4s.lmm.Criterion.REML)
    assertEquals(profile.baseline.options.criterion, mixeff4s.lmm.Criterion.ML)
    assertEquals(profile.baseline.y, changed)
    assert(profile.lower.value.nonEmpty && profile.upper.value.nonEmpty)

  test("budget failures retain optimizer evidence and unresolved brackets never fabricate endpoints"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, tight).toOption.get
    val limited = LmmProfile
      .confint(
        fit,
        LmmProfileTarget.Coefficient("days"),
        LmmProfileOptions(maxEvaluations = 1, optimizerControls = tight.optimizerControls.get)
      )
      .toOption
      .get
    assert(limited.lower.value.isEmpty && limited.upper.value.isEmpty)
    assert(limited.evaluations.exists(_.optimizer.exists(_.requireConvergence.isLeft)))
    val unresolved = LmmProfile
      .confint(
        fit,
        LmmProfileTarget.Coefficient("days"),
        LmmProfileOptions(maxEvaluations = 2000, maxBisectionSteps = 1, optimizerControls = tight.optimizerControls.get)
      )
      .toOption
      .get
    assertEquals(unresolved.lower.reason, Some("profile_root_unresolved"))
    assertEquals(unresolved.upper.reason, Some("profile_root_unresolved"))
