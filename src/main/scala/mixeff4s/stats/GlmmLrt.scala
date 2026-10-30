package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.{Approximation, GlmmFit}

/** Comparable approximate likelihoods, not an automatic asymptotic reference law. */
object GlmmLrt:
  def statistic(smaller: GlmmFit, larger: GlmmFit): FitResult[LrtStatistic] =
    for
      _ <- compatible(smaller, larger)
      _ <- LikelihoodQualification.glmm(smaller)
      _ <- LikelihoodQualification.glmm(larger)
      value <- summarize(smaller, larger)
    yield value

  def compare(smaller: GlmmFit, larger: GlmmFit): FitResult[Lrt] =
    for
      value <- statistic(smaller, larger)
      _ <-
        if LikelihoodNesting.sameRandom(
            LikelihoodNesting.capture(smaller.context.design),
            LikelihoodNesting.capture(larger.context.design)
          ) &&
          !LikelihoodQualification.glmmBoundary(smaller) && !LikelihoodQualification.glmmBoundary(larger)
        then Right(())
        else
          unavailable("lrt_reference_unavailable", "Covariance changes or boundary fits require explicit calibration.")
    yield Lrt(
      value.nobs,
      value.formulas,
      value.dof,
      value.loglik,
      value.deviance,
      value.likelihoodRatio,
      value.parameterDifference,
      ChiSquare.sf(value.likelihoodRatio, value.parameterDifference)
    )

  private[stats] def compatible(a: GlmmFit, b: GlmmFit): FitResult[Unit] =
    val ea = LikelihoodNesting.capture(a.context.design); val eb = LikelihoodNesting.capture(b.context.design)
    if a.approximation == Approximation.FastPirls || b.approximation == Approximation.FastPirls then
      unavailable("lrt_requires_joint", "Fast-PIRLS is not a joint likelihood comparison.")
    else if a.family != b.family || a.link != b.link || a.options.nAgq != b.options.nAgq || a.approximation != b.approximation
    then unavailable("lrt_estimator_mismatch", "Family, link, approximation and quadrature must match.")
    else if a.y != b.y || a.context.response.trials != b.context.response.trials ||
      a.context.response.offset != b.context.response.offset || a.context.response.caseWeights != b.context.response.caseWeights
    then unavailable("lrt_response_mismatch", "Responses, trials, case weights and offsets must match row for row.")
    else if a.context.response.caseWeights.exists(_ != 1) then
      unavailable("lrt_case_weights", "Calibration of power likelihood is unavailable.")
    else if !LikelihoodNesting.nested(ea, eb) then
      unavailable("lrt_nesting_unavailable", "Exact fitted-column/grouping/covariance inclusion is not witnessed.")
    else if b.beta.length + b.theta.length <= a.beta.length + a.theta.length then
      unavailable("lrt_parameter_order", "Larger model must have more parameters.")
    else Right(())

  private[stats] def summarize(a: GlmmFit, b: GlmmFit): FitResult[LrtStatistic] =
    val difference = a.deviance - b.deviance
    if difference < -1e-7 then
      unavailable(
        "lrt_negative_statistic",
        s"Larger model objective is worse by ${-difference}; no failed pair is discarded."
      )
    else
      Right(
        LrtStatistic(
          a.n,
          Vector(a.formula.toString, b.formula.toString),
          Vector(a.beta.length + a.theta.length, b.beta.length + b.theta.length),
          Vector(-a.deviance / 2, -b.deviance / 2),
          Vector(a.deviance, b.deviance),
          math.max(0, difference),
          b.beta.length + b.theta.length - a.beta.length - a.theta.length
        )
      )
  private def unavailable(code: String, message: String): FitResult[Nothing] = Left(
    MixedModelError.InferenceUnavailable(code, message)
  )
