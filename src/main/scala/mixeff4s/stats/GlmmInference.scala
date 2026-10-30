package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.{GlmmFit, GlmmAssessment, JointLikelihood}

final case class GlmmInferenceOptions(
    targets: Vector[BootstrapTarget] = Vector.empty,
    level: Double = 0.95,
    curvature: CurvatureOptions = CurvatureOptions()
)

/** Approximate normal interval from selected-estimator observed information; no p-value is supplied. */
final case class GlmmCoefficientInterval(
    label: String,
    estimate: Double,
    standardError: Double,
    z: Double,
    lower: Double,
    upper: Double,
    level: Double
)

final class GlmmInferenceResult private[stats] (
    val coefficientNames: Vector[String],
    val covariance: Vector[Vector[Double]],
    val curvature: CurvatureEvidence,
    val intervals: Vector[GlmmCoefficientInterval],
    val estimator: String,
    val nAgq: Int
)

object GlmmInference:
  def coefficients(
      fit: GlmmFit,
      options: GlmmInferenceOptions = GlmmInferenceOptions()
  ): FitResult[GlmmInferenceResult] =
    for
      _ <-
        if options.level.isFinite && options.level > 0 && options.level < 1 then Right(())
        else Left(MixedModelError.InvalidArgument("normal interval level must be in (0,1)"))
      assessment = GlmmAssessment.assess(fit)
      _ <-
        if assessment.eligibleForBootstrap then Right(())
        else
          Left(
            MixedModelError
              .InferenceUnavailable("coefficient_inference_unqualified", assessment.refusalCodes.mkString(", "))
          )
      targets <- Bootstrap.compileTargets(
        fit.fixedEffects,
        if options.targets.nonEmpty then options.targets else fit.feNames.map(BootstrapTarget.Coefficient.apply)
      )
      evaluator <- JointLikelihood(fit)
      evidence <- ObservedInformation.compute(
        fit.beta ++ fit.theta,
        Vector.fill(fit.beta.length)(Double.NegativeInfinity) ++ fit.context.design.parmap.map { case (_, row, col) =>
          if row == col then 0.0 else Double.NegativeInfinity
        },
        options.curvature
      )(evaluator.evaluate)
      covariance = evidence.covariance.take(fit.beta.length).map(_.take(fit.beta.length))
      rows <- targets.foldLeft[FitResult[Vector[GlmmCoefficientInterval]]](Right(Vector.empty)): (acc, target) =>
        acc.flatMap: prior =>
          val estimate = target.weights.zip(fit.beta).map(_ * _).sum
          val variance = target.weights.indices
            .map(i => target.weights.indices.map(j => target.weights(i) * covariance(i)(j) * target.weights(j)).sum)
            .sum
          val se = math.sqrt(variance)
          val z = math.sqrt(InferenceQuantiles.chiSquare1(options.level))
          if !se.isFinite || se <= 0 || !estimate.isFinite || !(estimate - z * se).isFinite || !(estimate + z * se).isFinite
          then
            Left(
              MixedModelError.InferenceUnavailable(
                "contrast_information_unavailable",
                "Contrast has zero or invalid observed variance."
              )
            )
          else
            Right(
              prior :+ GlmmCoefficientInterval(
                target.label,
                estimate,
                se,
                estimate / se,
                estimate - z * se,
                estimate + z * se,
                options.level
              )
            )
    yield new GlmmInferenceResult(fit.feNames, covariance, evidence, rows, fit.algorithmLabel, fit.nAgq)

private[stats] object InferenceQuantiles:
  def chiSquare1(level: Double): Double =
    var lo = 0.0; var hi = 1.0
    while ChiSquare.sf(hi, 1) > 1 - level do hi *= 2
    var i = 0
    while i < 80 do
      val mid = (lo + hi) / 2
      if ChiSquare.sf(mid, 1) > 1 - level then lo = mid else hi = mid
      i += 1
    (lo + hi) / 2
