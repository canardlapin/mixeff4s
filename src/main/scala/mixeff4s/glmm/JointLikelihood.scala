package mixeff4s.glmm

import mixeff4s.error.{FitResult, MixedModelError}

/** Package bridge for inference: owns fresh scratch, evaluates the selected full-density objective. It never fits,
  * changes approximation, or mutates the captured fit.
  */
private[mixeff4s] final class JointLikelihood private (private val workspace: JointWorkspace):
  def evaluate(parameters: Vector[Double]): FitResult[Double] = synchronized:
    workspace.evaluateParameters(parameters).map(_.objective)

private[mixeff4s] object JointLikelihood:
  def apply(fit: GlmmFit, scoreTolerance: Double = 1e-10): FitResult[JointLikelihood] =
    if fit.approximation == Approximation.FastPirls then
      Left(
        MixedModelError.InferenceUnavailable(
          "inference_requires_joint",
          "Inference requires an explicitly selected joint estimator."
        )
      )
    else
      JointWorkspace
        .create(
          fit.context.design,
          fit.family,
          fit.link,
          fit.context.response,
          fit.options.copy(conditionalScoreTolerance = scoreTolerance)
        )
        .map(new JointLikelihood(_))
