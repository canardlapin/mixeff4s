package mixeff4s.stats

import mixeff4s.data.ModelFrame
import mixeff4s.error.FitResult
import mixeff4s.glmm.*
import mixeff4s.model.Family

class LikelihoodBoundarySuite extends munit.FunSuite:
  private def right[A](x: FitResult[A]): A = x.fold(e => fail(s"${e.code}: ${e.message}"), identity)
  private val x = Vector.fill(12)(Vector(-1.0, 0.0, 1.0)).flatten
  private val frame = right(
    ModelFrame.of(
      "x" -> ModelFrame.numeric(x),
      "g" -> ModelFrame.factor(Vector.tabulate(36)(i => s"g${i / 3}")),
      "y" -> ModelFrame.numeric(x.map(v => if v < 0 then .25 else if v == 0 then .4 else .55))
    )
  )
  private val options = GlmmOptions.jointLaplace.copy(trials = Some(Vector.fill(36)(20)))
  private lazy val large = right(Glmm.fit("y ~ 1 + x + (1 | g)", frame, Family.Binomial, None, options))
  private lazy val small = right(Glmm.fit("y ~ 1 + (1 | g)", frame, Family.Binomial, None, options))
  test("strict scalar boundary likelihood evidence differs from interior coefficient eligibility"):
    println(
      s"BOUNDARY_FIT small=${small.theta}/${small.beta} large=${large.theta}/${large.beta} assessment=${GlmmAssessment.assess(large)}"
    )
    assert(LikelihoodQualification.glmmBoundary(large))
    assert(!GlmmAssessment.assess(large).eligibleForBootstrap)
    assertEquals(LikelihoodQualification.glmm(large), Right(()))
    assertEquals(LikelihoodQualification.glmm(small), Right(()))
    assert(right(GlmmLrt.statistic(small, large)).likelihoodRatio > 1)
    assertEquals(GlmmLrt.compare(small, large).left.toOption.get.code, "lrt_reference_unavailable")
    assertEquals(GlmmInference.coefficients(large).left.toOption.get.code, "coefficient_inference_unqualified")
  test("feasible variance score detects zero-theta saddle even though theta derivative vanishes"):
    val data = right(
      ModelFrame.of(
        "g" -> ModelFrame.factor(Vector.tabulate(36)(i => s"g${i / 3}")),
        "x" -> ModelFrame.numeric(x),
        "y" -> ModelFrame.numeric(Vector.tabulate(36)(i => if i / 3 < 6 then .1 else .9))
      )
    )
    val formula = mixeff4s.formula.Formula.parse("y ~ 1 + x + (1 | g)").toOption.get
    val zero =
      right(Glmm.jointDeviance(formula, data, Family.Binomial, Vector(0.0, 0.0), Vector(0.0), options = options))
    val positive =
      right(Glmm.jointDeviance(formula, data, Family.Binomial, Vector(0.0, 0.0), Vector(.01), options = options))
    val derivativeInVariance = (positive.objective - zero.objective) / 1e-4
    assert(derivativeInVariance < -1)
  test("crossed scalar boundary attempt refuses failed tangent stationarity"):
    val expanded = right(
      ModelFrame.of((frame.columns :+ ("batch" -> ModelFrame.factor(Vector.tabulate(36)(i => s"b${(i / 3) % 4}"))))*)
    )
    val pairOptions = options.copy(optimizerControls =
      Some(mixeff4s.optimizer.TrustBqOptions(finalRadius = 1e-6, ftolRequiresLocalRadius = true, stallIterations = 10))
    )
    val fit = right(Glmm.fit("y ~ 1 + x + (1 | g) + (1 | batch)", expanded, Family.Binomial, None, pairOptions))
    assertEquals(GlmmLrt.compatible(large, fit), Right(()))
    val comparison = GlmmLrt.statistic(large, fit)
    assertEquals(comparison.left.toOption.get.code, "likelihood_tangent_gradient")
    assertEquals(GlmmLrt.compare(large, fit).left.toOption.get.code, "likelihood_tangent_gradient")
    println(
      s"CROSSED_BOUNDARY_REFUSAL theta=${fit.theta} beta=${fit.beta} refusal=${comparison.left.toOption.map(_.code)}"
    )
  test("scalar boundary null comparison retains all seven generated model pairs"):
    val bootstrap = right(LikelihoodBootstrap.compare(small, large, LikelihoodBootstrapOptions(7, 2026102903L)))
    println(
      s"BOUNDARY_BOOTSTRAP attempts=${bootstrap.attempts.length} failures=${bootstrap.failedOrUnqualifiedIndices} calibration=${bootstrap.calibration}"
    )
    assertEquals(bootstrap.attempts.length, 7)
    bootstrap.attempts.foreach(a =>
      println(s"BOUNDARY_PAIR index=${a.index} seed=${a.seed} statistic=${a.statistic.map(_.likelihoodRatio)}")
    )
    if bootstrap.failedOrUnqualifiedIndices.nonEmpty then assert(bootstrap.calibration.isLeft)
    else assert(bootstrap.calibration.toOption.get.pvalue > 0)
