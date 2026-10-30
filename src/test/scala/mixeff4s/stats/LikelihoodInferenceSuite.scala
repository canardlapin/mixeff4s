package mixeff4s.stats

import mixeff4s.glmm.*
import mixeff4s.fixtures.GlmmAnalysisFrames
import mixeff4s.model.Family
import mixeff4s.error.{FitResult, MixedModelError}

class LikelihoodInferenceSuite extends munit.FunSuite:
  private def right[A](x: FitResult[A]): A = x.fold(e => fail(s"${e.code}: ${e.message}"), identity)
  private val options = GlmmOptions.jointAgq(5).copy(trials = Some(GlmmAnalysisFrames.trials))
  private lazy val large = right(
    Glmm.fit("y ~ 1 + x + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options)
  )
  private lazy val small = right(
    Glmm.fit("y ~ 1 + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options)
  )
  test("selected-estimator observed covariance and linear contrasts retain nuisance uncertainty"):
    val output = right(
      GlmmInference.coefficients(
        large,
        GlmmInferenceOptions(targets =
          Vector(BootstrapTarget.Coefficient("x"), BootstrapTarget.Contrast("twox", Map("x" -> 2.0)))
        )
      )
    )
    assertEquals(output.nAgq, 5)
    val a = output.intervals(0); val b = output.intervals(1)
    assertEqualsDouble(b.estimate, 2 * a.estimate, 1e-12)
    assertEqualsDouble(b.standardError, 2 * a.standardError, 1e-12)
    assert(a.lower < a.estimate && a.upper > a.estimate)
    println(
      s"INFERENCE_COV beta=${large.beta} theta=${large.theta} covariance=${output.covariance} rows=${output.intervals} stability=${output.curvature.hessianRelativeDifference}/${output.curvature.covarianceRelativeDifference}"
    )
  test("nested GLMM statistic matches explicit objectives and incompatible estimators refuse"):
    val result = right(Lrt.statistic(small, large))
    assertEqualsDouble(result.likelihoodRatio, small.deviance - large.deviance, 1e-10)
    assertEquals(result.parameterDifference, 1)
    assert(right(Lrt.compare(small, large)).pvalue < .001)
    val laplace =
      right(Glmm.fit("y ~ 1 + x + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options.copy(nAgq = 1)))
    assertEquals(Lrt.statistic(small, laplace).left.toOption.get.code, "lrt_estimator_mismatch")
    assertEquals(Lrt.statistic(large, small).left.toOption.get.code, "lrt_nesting_unavailable")
  test("null simulation preserves complete failed-pair ledger and later indexed seeds"):
    val run =
      right(LikelihoodBootstrap.compare(small, large, LikelihoodBootstrapOptions(3, 123L, maxEvaluations = Some(1))))
    assertEquals(run.attempts.length, 3)
    assert(run.attempts.forall(a => a.nullFit.isLeft && a.alternativeFit.isLeft))
    assertEquals(run.calibration.left.toOption.get.code, "lrt_bootstrap_unqualified")
    assertEquals(run.failedOrUnqualifiedIndices, Vector(0, 1, 2))
  test("tail calibration includes ties, plus-one correction and Monte Carlo interval at all-zero extremes"):
    val observed = LrtStatistic(1, Vector("small", "large"), Vector(1, 2), Vector(0, 0), Vector(0, 0), 2, 1)
    val attempted = scala.collection.mutable.ArrayBuffer.empty[Int]
    val run = LikelihoodBootstrap.execute(
      observed,
      LikelihoodBootstrapOptions(4, 7L),
      _ => Right(Vector(0.0)),
      (i, _) => { attempted += i; Right(i) },
      (i, _) => Right(i),
      (i, _, _) =>
        if i == 1 then Left(MixedModelError.Optimization("injected")) else Right(observed.copy(likelihoodRatio = 2))
    )
    assertEquals(attempted.toVector, Vector(0, 1, 2, 3))
    assertEquals(run.failedOrUnqualifiedIndices, Vector(1))
    assert(run.calibration.isLeft)
    val zero = LikelihoodBootstrap.tail(0, 39, .95)
    assertEqualsDouble(zero.pvalue, 1.0 / 40, 1e-15)
    assert(zero.upper > 0 && zero.lower == 0)
    val all = LikelihoodBootstrap.tail(39, 39, .95)
    assertEquals(all.pvalue, 1.0); assert(all.lower < 1)

  test("numerical zero ties have a recorded tolerance and can be counted exactly when requested"):
    val observed = LrtStatistic(1, Vector("small", "large"), Vector(1, 2), Vector(0, 0), Vector(0, 0), 2e-13, 1)
    def run(tolerance: Double) = LikelihoodBootstrap.execute(
      observed,
      LikelihoodBootstrapOptions(3, 7L, statisticTolerance = tolerance),
      _ => Right(Vector(0.0)),
      (_, _) => Right(0),
      (_, _) => Right(0),
      (_, _, _) => Right(observed.copy(likelihoodRatio = 0))
    )
    val tolerant = run(1e-7)
    assertEquals(tolerant.calibration.toOption.get.pvalue, 1.0)
    assertEquals(tolerant.options.statisticTolerance, 1e-7)
    assertEqualsDouble(run(0).calibration.toOption.get.pvalue, .25, 1e-15)
