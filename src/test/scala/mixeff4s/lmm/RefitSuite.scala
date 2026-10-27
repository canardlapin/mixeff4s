package mixeff4s.lmm

import mixeff4s.data.ModelFrame
import mixeff4s.error.MixedModelError
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.stats.Lrt

class RefitSuite extends munit.FunSuite:
  private val formula = "reaction ~ 1 + days + (1 + days | subj)"
  private val yA = Sleepstudy.reaction
  private val yB = yA.zip(Sleepstudy.days).map((y, x) => y + 0.1 * x * x)
  private def frame(y: Vector[Double]): ModelFrame = ModelFrame
    .of(
      "reaction" -> ModelFrame.numeric(y),
      "days" -> ModelFrame.numeric(Sleepstudy.days),
      "subj" -> ModelFrame.factor(Sleepstudy.subj)
    )
    .toOption
    .get
  private def equivalent(got: LmmFit, expected: LmmFit, tolerance: Double = 1e-7): Unit =
    assertEquals(got.y, expected.y)
    assertEquals(got.options, expected.options)
    assertEquals(got.feNames, expected.feNames)
    assertEqualsDouble(got.objective, expected.objective, tolerance)
    got.beta.zip(expected.beta).foreach((a, b) => assertEqualsDouble(a, b, tolerance))
    assertEqualsDouble(got.sigma, expected.sigma, tolerance)
    got.theta.zip(expected.theta).foreach((a, b) => assertEqualsDouble(a, b, tolerance))

  for reml <- Vector(false, true); strategy <- LmmOptimizer.values do
    test(s"compiled fit and cold A-B-A refits equal rebuilds: REML=$reml $strategy"):
      val options = FitOptions(if reml then Criterion.REML else Criterion.ML, optimizer = strategy)
      val plan = Lmm.prepare(formula, Sleepstudy.frame).toOption.get
      val session = plan.session(options).toOption.get
      val freshA = Lmm.fit(formula, frame(yA), options).fold(e => fail(e.message), identity)
      equivalent(plan.fit(options).fold(e => fail(e.message), identity), freshA)
      val old = session.refit(yA).fold(e => fail(e.message), identity)
      equivalent(old, freshA)
      val freshB = Lmm.fit(formula, frame(yB), options).fold(e => fail(e.message), identity)
      equivalent(session.refit(yB).fold(e => fail(e.message), identity), freshB)
      equivalent(session.refit(yA).fold(e => fail(e.message), identity), freshA)
      equivalent(old, freshA)
      // New response is bound to new evidence; older snapshots retain their original response.
      Lrt.statistic(old, old) match
        case Left(e: MixedModelError.InferenceUnavailable) => assert(e.details.contains("more parameters"))
        case other => fail(s"expected valid evidence followed by count refusal: $other")

  test("previous start requires success and explicit warm refits preserve options and objective"):
    val options = FitOptions.ml.copy(optimizer = LmmOptimizer.AnalyticGradient)
    val plan = Lmm.prepare(formula, Sleepstudy.frame).toOption.get
    val session = plan.session(options).toOption.get
    assert(session.refit(yA, RefitStart.Previous).isLeft)
    val cold = session.refit(yA).fold(e => fail(e.message), identity)
    val warm = session.refit(yA, RefitStart.Previous).fold(e => fail(e.message), identity)
    assertEquals(warm.options, options)
    assertEqualsDouble(warm.objective, cold.objective, 1e-6)
    assert(session.refit(yA.dropRight(1), RefitStart.Previous).isLeft)
    assert(session.refit(yA.updated(0, Double.NaN), RefitStart.Previous).isLeft)
    val again = session.refit(yA, RefitStart.Previous).fold(e => fail(e.message), identity)
    assertEqualsDouble(again.objective, cold.objective, 1e-6)
    assertEquals(again.y, yA)

  test("failed numerical fit does not poison a later cold refit"):
    val session = Lmm.prepare(formula, Sleepstudy.frame).toOption.get.session(FitOptions.ml).toOption.get
    val first = session.refit(yA).fold(e => fail(e.message), identity)
    assert(session.refit(Vector.fill(180)(0.0)).isLeft)
    equivalent(session.refit(yA).fold(e => fail(e.message), identity), first)

  test("independent batches reset starts and separate sessions do not share scratch"):
    val plan = Lmm.prepare(formula, Sleepstudy.frame).toOption.get
    val options = FitOptions.ml
    val batch = plan.fitResponses(Vector(yA, yB, yA), options).fold(e => fail(e.message), identity)
    equivalent(batch.head, batch.last)
    val reversed = plan.fitResponses(Vector(yB, yA), options).fold(e => fail(e.message), identity)
    equivalent(batch(1), reversed.head)
    equivalent(batch.head, reversed.last)
    val warm =
      plan.fitResponses(Vector(yA, yB), options, BatchStart.WarmSequential).fold(e => fail(e.message), identity)
    assertEqualsDouble(warm.last.objective, batch(1).objective, 1e-4)
    val budgetSession = plan.session(options.copy(maxEvaluations = 1)).toOption.get
    assert(budgetSession.refit(yA).isLeft)
    assert(budgetSession.refit(yA, RefitStart.Previous).isLeft)
