package mixeff4s.lmm

import mixeff4s.comparison.Datasets
import mixeff4s.fixtures.{Pastes, Penicillin, Pw2Nested, Sleepstudy}

class GradientFitSuite extends munit.FunSuite:
  for reml <- Vector(false, true) do
    test(s"gradient strategy preserves classic nested/crossed/slope fits: REML=$reml"):
      val options = FitOptions(if reml then Criterion.REML else Criterion.ML)
      val cases = Vector(
        ("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame),
        ("strength ~ 1 + (1 | batch / cask)", Pastes.frame),
        ("diameter ~ 1 + (1 | plate) + (1 | sample)", Penicillin.frame)
      )
      cases.foreach: (formula, frame) =>
        val old = Lmm.fit(formula, frame, options).fold(e => fail(e.message), identity)
        val fit = Lmm
          .fit(formula, frame, options.copy(optimizer = LmmOptimizer.AnalyticGradient))
          .fold(e => fail(e.message), identity)
        assertEqualsDouble(fit.objective, old.objective, 0.002, clues(formula))
        fit.beta.zip(old.beta).foreach((a, b) => assertEqualsDouble(a, b, 0.002))
        assert(fit.optimizer.gradientEvaluations > 0)
        assert(fit.optimizer.projectedGradientNorm.exists(_.isFinite))
        println(
          s"GRADIENT_FIT formula=$formula reml=$reml objective=${fit.objective} gradient_evals=${fit.optimizer.fevals} interpolation_evals=${old.optimizer.fevals}"
        )

  test("pw2_nested default-budget REML matches frozen lme4 adoption reference"):
    val fit = Lmm
      .fit(Pw2Nested.formula, Pw2Nested.frame, FitOptions.reml.copy(optimizer = LmmOptimizer.AnalyticGradient))
      .fold(e => fail(e.message), identity)
    assertEqualsDouble(fit.beta.head, 49.7641455983855, 1e-4)
    assertEqualsDouble(fit.sigma, 0.602065789420922, 1e-3)
    assertEqualsDouble(fit.loglikelihood, -219.865761241694, 1e-3)
    assertEquals(fit.theta.length, 2)
    println(s"GRADIENT_PW2 objective=${fit.objective} evaluations=${fit.optimizer.fevals}")

  test("maximal singular REML meets the Rust gradient regression bound within default budget"):
    val frame = Datasets.load("singular").frame
    val fit = Lmm
      .fit(
        "y ~ 1 + A * B * C + (A * B * C | group)",
        frame,
        FitOptions.reml.copy(optimizer = LmmOptimizer.AnalyticGradient)
      )
      .fold(e => fail(e.message), identity)
    assertEquals(fit.theta.length, 36)
    assert(fit.objective < 770.0, clues(fit.objective))
    assert(fit.optimizer.fevals <= 1000)
    val assessed = LmmAssessment.assess(fit)
    assert(assessed.numericalTerminationAccepted)
    assessed.gradient match
      case GradientCheck.Failed(norm, tolerance) =>
        assert(norm > 0.1)
        assert(norm > tolerance)
      case other => fail(s"singular diagnostic must expose its failed first-order check: $other")
    assert(!assessed.eligibleForBootstrap)
    println(
      s"GRADIENT_MAXIMAL objective=${fit.objective} evaluations=${fit.optimizer.fevals} stop=${fit.optimizer.stopReason} projected_gradient=${fit.optimizer.projectedGradientNorm}"
    )
