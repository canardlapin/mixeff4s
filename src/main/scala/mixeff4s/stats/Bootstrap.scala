package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.GlmmFit
import mixeff4s.lmm.LmmFit

/** Parametric bootstrap. Unavailable until a certified simulator exists. */
object Bootstrap:
  def parametric(fit: LmmFit, nsim: Int): FitResult[Nothing] =
    val _ = (fit, nsim)
    Left(
      MixedModelError.InferenceUnavailable(
        "bootstrap_unavailable",
        "Parametric bootstrap is not implemented."
      )
    )

  def parametric(fit: GlmmFit, nsim: Int): FitResult[Nothing] =
    val _ = (fit, nsim)
    Left(
      MixedModelError.InferenceUnavailable(
        "bootstrap_unavailable",
        "GLMM parametric bootstrap is refused until a certified response simulator exists."
      )
    )
