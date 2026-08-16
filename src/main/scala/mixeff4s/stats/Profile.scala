package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.GlmmFit
import mixeff4s.lmm.LmmFit

/** Profile-likelihood intervals. Unavailable until the profiler is ported. */
object Profile:
  def confint(fit: LmmFit, level: Double = 0.95): FitResult[Nothing] =
    val _ = (fit, level)
    Left(
      MixedModelError.InferenceUnavailable(
        "profile_unavailable",
        "Profile-likelihood intervals are not implemented."
      )
    )

  def confint(fit: GlmmFit): FitResult[Nothing] =
    confint(fit, 0.95)

  def confint(fit: GlmmFit, level: Double): FitResult[Nothing] =
    val _ = (fit, level)
    Left(
      MixedModelError.InferenceUnavailable(
        "profile_unavailable",
        "GLMM profile-likelihood intervals are refused."
      )
    )
