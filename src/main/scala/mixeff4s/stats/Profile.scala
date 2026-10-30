package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.GlmmFit
import mixeff4s.lmm.LmmFit

/** LMM profiles require an explicit scalar target; generic and GLMM profile requests remain typed refusals. */
object Profile:
  def confint(fit: LmmFit, target: LmmProfileTarget, options: LmmProfileOptions): FitResult[LmmProfileResult] =
    LmmProfile.confint(fit, target, options)
  def confint(fit: LmmFit, target: LmmProfileTarget): FitResult[LmmProfileResult] = LmmProfile.confint(fit, target)
  def confint(fit: LmmFit, level: Double = 0.95): FitResult[Nothing] =
    val _ = (fit, level)
    Left(
      MixedModelError.InferenceUnavailable(
        "profile_unavailable",
        "Specify an LMM profile target; generic profile-likelihood intervals are unavailable."
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
