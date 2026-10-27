package mixeff4s.stats

import mixeff4s.data.ModelFrame
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.examples.GlmmAnalysis
import mixeff4s.fixtures.GlmmAnalysisFrames
import mixeff4s.glmm.*
import mixeff4s.lmm.{Simulation, SimulationMode}
import mixeff4s.model.Family
import scala.concurrent.duration.*

class GlmmBootstrapSuite extends munit.FunSuite:
  override val munitTimeout = 120.seconds
  private def right[A](value: FitResult[A]): A = value.fold(e => fail(s"${e.code}: ${e.message}"), identity)
  private lazy val plan = right(Glmm.prepare(GlmmAnalysis.formula, GlmmAnalysisFrames.binomial, Family.Binomial))
  private lazy val base = right(plan.fit(GlmmAnalysis.fitOptions))
  private def options(n: Int) = GlmmAnalysis.bootstrapOptions.copy(nsim = n)

  test("injected middle failure retains later seeds, responses and independent cold-fit statistics"):
    val regular = right(GlmmBootstrap.parametric(base, options(5)))
    val session = right(base.context.compiled.session(base.options))
    val broken = right(GlmmBootstrap.withRefitter(base, options(5)): (index, response) =>
      if index == 2 then Left(MixedModelError.Interrupted("injected middle failure"))
      else session.refit(response))
    assertEquals(broken.attempts.length, 5)
    assertEquals(broken.attempts(2).refit.left.toOption.get.code, "interrupted")
    assertEquals(broken.failedOrUnqualifiedIndices, Vector(2))
    assert(broken.intervals.forall(_.left.toOption.exists(_.code == "bootstrap_replicate_unqualified")))
    Vector(0, 1, 3, 4).foreach: index =>
      val a = regular.attempts(index); val b = broken.attempts(index)
      assertEquals(a.seed, b.seed)
      assertEquals(a.response, b.response)
      right(a.statistics).zip(right(b.statistics)).foreach((x, y) => assertEqualsDouble(x, y, 1e-6))
      assert(b.qualified, clues(b.assessment))

  test("every budget failure is recorded without replacement responses"):
    val capped = base.options.copy(maxEvaluations = 1)
    val result = right(GlmmBootstrap.parametric(base, options(3).copy(fitOptions = Some(capped))))
    assertEquals(result.attempts.length, 3)
    assert(result.attempts.forall(_.refit.left.toOption.exists(_.code == "optimizer_budget_exhausted")))
    result.attempts.foreach: attempt =>
      assertEquals(right(attempt.response), right(GlmmSimulation.response(base, attempt.seed)))
      assertEquals(attempt.seed, Simulation.seedFor(options(3).seed, attempt.index.toLong))
    assertEquals(result.qualifiedCount, 0)
    assert(result.intervals.forall(_.isLeft))

  test("foreign designs and estimator/options mismatches cannot count as successful replicates"):
    val foreign = right(Glmm.prepare(GlmmAnalysis.formula, GlmmAnalysisFrames.binomial, Family.Binomial))
    val foreignSession = right(foreign.session(base.options))
    val wrongDesign = right(GlmmBootstrap.withRefitter(base, options(2))((_, y) => foreignSession.refit(y)))
    assert(wrongDesign.attempts.forall(_.statistics.left.toOption.exists(_.code == "bootstrap_refit_mismatch")))
    val wrongOptions = base.options.copy(nAgq = 1)
    val session = right(base.context.compiled.session(wrongOptions))
    val wrongEstimator = right(GlmmBootstrap.withRefitter(base, options(2))((_, y) => session.refit(y)))
    assert(wrongEstimator.attempts.forall(_.statistics.left.toOption.exists(_.code == "bootstrap_refit_mismatch")))
    assert(
      GlmmBootstrap
        .parametric(base, options(2).copy(fitOptions = Some(wrongOptions)))
        .left
        .toOption
        .exists(_.code == "bootstrap_estimator_mismatch")
    )

  test("invalid targets, levels and missing seeds remain precise refusals"):
    assert(GlmmBootstrap.parametric(base, options(2).copy(level = 1.0)).isLeft)
    assert(
      GlmmBootstrap.parametric(base, options(2).copy(targets = Vector(BootstrapTarget.Coefficient("missing")))).isLeft
    )
    assert(
      GlmmBootstrap
        .parametric(base, options(2).copy(targets = Vector(BootstrapTarget.Contrast("bad", Map("missing" -> 1.0)))))
        .isLeft
    )
    assertEquals(Bootstrap.parametric(base, 100).left.toOption.get.code, "bootstrap_seed_required")
    val one = right(GlmmBootstrap.parametric(base, options(1)))
    assert(one.intervals.forall(_.left.toOption.exists(_.code == "bootstrap_insufficient_replicates")))

  test("rare constant responses fail visibly and are never silently replaced"):
    val session = right(base.context.compiled.session(base.options))
    val result = right(GlmmBootstrap.withRefitter(base, options(3)): (index, y) =>
      if index == 1 then session.refit(Vector.fill(y.length)(0.0)) else session.refit(y))
    assertEquals(result.attempts(1).refit.left.toOption.get.code, "constant_response")
    assertEquals(result.attempts.length, 3)
    assert(result.attempts(2).refit.isRight)
    assert(result.intervals.forall(_.isLeft))

  test("unqualified boundary originals refuse before any bootstrap attempts"):
    val data = right(
      ModelFrame.of(
        "y" -> ModelFrame.numeric(Vector(0.0, 1.0, 1.0, 0.0, 0.0, 1.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c"))
      )
    )
    val fit = right(Glmm.fit("y ~ 1 + (1 | g)", data, Family.Bernoulli, None, GlmmOptions.jointLaplace))
    val request = GlmmBootstrapOptions(3, 42L, Vector(BootstrapTarget.Coefficient("(Intercept)")))
    assert(GlmmBootstrap.parametric(fit, request).left.toOption.exists(_.code == "bootstrap_unqualified_original"))
    val fast = right(Glmm.fit("y ~ 1 + (1 | g)", data, Family.Bernoulli))
    assert(GlmmBootstrap.parametric(fast, request).isLeft)

  test("conditional-mode experiments remain explicitly labelled in the result"):
    val result =
      right(GlmmBootstrap.parametric(base, options(2).copy(simulationMode = SimulationMode.ConditionOnModes)))
    assertEquals(result.options.simulationMode, SimulationMode.ConditionOnModes)
    assertEquals(result.attempts.length, 2)
