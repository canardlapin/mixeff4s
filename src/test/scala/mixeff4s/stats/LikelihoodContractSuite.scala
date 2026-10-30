package mixeff4s.stats

import mixeff4s.examples.LikelihoodInference
import mixeff4s.fixtures.{GlmmAnalysisFrames, Sleepstudy}
import mixeff4s.glmm.*
import mixeff4s.lmm.Lmm
import mixeff4s.model.Family
import mixeff4s.error.FitResult

class LikelihoodContractSuite extends munit.FunSuite:
  private def right[A](x: FitResult[A]): A = x.fold(e => fail(s"${e.code}: ${e.message}"), identity)
  private val options = GlmmOptions.jointAgq(5).copy(trials = Some(GlmmAnalysisFrames.trials))
  private lazy val small = right(
    Glmm.fit("y ~ 1 + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options)
  )
  private lazy val large = right(
    Glmm.fit("y ~ 1 + x + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options)
  )
  test("rowresponse, trials, offsets and power likelihood remain explicit comparison contracts"):
    val otherOffset = right(
      Glmm.fit(
        "y ~ 1 + x + (1 | g)",
        GlmmAnalysisFrames.binomial,
        Family.Binomial,
        None,
        options.copy(offset = Some(Vector.fill(90)(.2)))
      )
    )
    assertEquals(GlmmLrt.compatible(small, otherOffset).left.toOption.get.code, "lrt_response_mismatch")
    val power = right(
      Glmm.fit(
        "y ~ 1 + x + (1 | g)",
        GlmmAnalysisFrames.binomial,
        Family.Binomial,
        None,
        options.copy(weights = Some(Vector.fill(90)(2.0)))
      )
    )
    assertEquals(GlmmInference.coefficients(power).left.toOption.get.code, "coefficient_inference_unqualified")
    assertEquals(GlmmLrt.compatible(power, power).left.toOption.get.code, "lrt_case_weights")
    val plan = right(Glmm.prepare("y ~ 1 + x + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial))
    val changed = right(plan.fit(large.y.updated(0, .5), options))
    assertEquals(GlmmLrt.compatible(small, changed).left.toOption.get.code, "lrt_response_mismatch")
  test("invalid targets and information budgets refuse without covariance substitution"):
    val unknown =
      GlmmInference.coefficients(large, GlmmInferenceOptions(Vector(BootstrapTarget.Coefficient("unknown"))))
    assertEquals(unknown.left.toOption.get.code, "invalid_argument")
    val limited =
      GlmmInference.coefficients(large, GlmmInferenceOptions(curvature = CurvatureOptions(maxParameters = 1)))
    assertEquals(limited.left.toOption.get.code, "curvature_dimension")
    val zero =
      GlmmInference.coefficients(large, GlmmInferenceOptions(Vector(BootstrapTarget.Contrast("zero", Map("x" -> 0.0)))))
    assertEquals(zero.left.toOption.get.code, "contrast_information_unavailable")
  test("Gaussian null bootstrap attempts both models and retains qualification failures"):
    val controls = LikelihoodInference.tightML.copy(
      optimizer = mixeff4s.lmm.LmmOptimizer.AnalyticGradient,
      optimizerControls = LikelihoodInference.tightML.optimizerControls.map(_.copy(gradientTolerance = 1e-10))
    )
    val large = right(Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, controls))
    val small = right(Lmm.fit("reaction ~ 1 + (1 | subj)", Sleepstudy.frame, controls))
    println(
      s"LMM_NULL_ORIGINAL small=${mixeff4s.lmm.LmmAssessment.assess(small)} large=${mixeff4s.lmm.LmmAssessment.assess(large)}"
    )
    val result = right(LikelihoodBootstrap.compare(small, large, LikelihoodBootstrapOptions(7, 2026102904L)))
    assertEquals(result.attempts.length, 7)
    result.attempts.foreach: a =>
      assert(a.response.isRight)
      println(s"LMM_NULL_PAIR index=${a.index} statistic=${a.statistic.map(_.likelihoodRatio)}")
    println(s"LMM_NULL_BOOTSTRAP qualified=${result.attempts.count(_.qualified)} calibration=${result.calibration}")
    if result.failedOrUnqualifiedIndices.nonEmpty then assert(result.calibration.isLeft)
    val capped = right(
      LikelihoodBootstrap.compare(small, large, LikelihoodBootstrapOptions(3, 2026102904L, maxEvaluations = Some(1)))
    )
    assertEquals(capped.attempts.length, 3)
    assert(capped.calibration.isLeft)
