package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.GlmmFit
import mixeff4s.lmm.{
  FitOptions,
  LmmAssessment,
  LmmAssessmentOptions,
  LmmFit,
  LmmFitContext,
  RefitStart,
  Simulation,
  SimulationMode
}

enum BootstrapTarget:
  case Coefficient(name: String)
  case Contrast(label: String, weightsByName: Map[String, Double])

object BootstrapTarget:
  def retainedBasisCoefficient(name: String): BootstrapTarget = BootstrapTarget.Coefficient(name)
  def retainedBasisContrast(label: String, weightsByName: Map[String, Double]): BootstrapTarget =
    BootstrapTarget.Contrast(label, weightsByName)

final case class BootstrapOptions(
    nsim: Int,
    seed: Long,
    targets: Vector[BootstrapTarget],
    simulationMode: SimulationMode = SimulationMode.DrawRandomEffects,
    fitOptions: Option[FitOptions] = None,
    assessmentOptions: LmmAssessmentOptions = LmmAssessmentOptions(),
    level: Double = 0.95
)
final case class BootstrapAttempt(
    index: Int,
    seed: Long,
    response: FitResult[Vector[Double]],
    refit: FitResult[LmmFit],
    assessment: Option[LmmAssessment],
    statistics: FitResult[Vector[Double]]
):
  def qualified: Boolean = refit.isRight && statistics.isRight && assessment.exists(_.eligibleForBootstrap)

final case class BootstrapInterval(target: String, lower: Double, upper: Double, level: Double)

/** Immutable record of every requested attempt. Failed or unqualified replicates are never dropped for intervals. */
final class BootstrapResult private[stats] (
    val options: BootstrapOptions,
    val originalAssessment: LmmAssessment,
    val attempts: Vector[BootstrapAttempt],
    val intervals: Vector[FitResult[BootstrapInterval]]
):
  def qualifiedCount: Int = attempts.count(_.qualified)
  def failedOrUnqualifiedIndices: Vector[Int] = attempts.filterNot(_.qualified).map(_.index)

object Bootstrap:
  def parametric(fit: LmmFit, options: BootstrapOptions): FitResult[BootstrapResult] =
    for
      prepared <- prepare(fit, options)
      session <- prepared.context.compiled.session(options.fitOptions.getOrElse(fit.options))
    yield execute(fit, options, prepared, (_, response) => session.refit(response, RefitStart.Cold))

  /** Test seam for failure accounting; the original fit and every returned replicate remain subject to assessment. */
  private[stats] def withRefitter(fit: LmmFit, options: BootstrapOptions)(
      refit: (Int, Vector[Double]) => FitResult[LmmFit]
  ): FitResult[BootstrapResult] =
    prepare(fit, options).map(prepared => execute(fit, options, prepared, refit))

  private final case class Target(label: String, weights: Vector[Double])
  private final case class Prepared(context: LmmFitContext, targets: Vector[Target], original: LmmAssessment)

  private def prepare(fit: LmmFit, options: BootstrapOptions): FitResult[Prepared] =
    for
      context <- fit.validatedContext
      _ <-
        if options.nsim > 0 && options.targets.nonEmpty && options.level > 0 && options.level < 1 && options.level.isFinite
        then Right(())
        else Left(MixedModelError.InvalidArgument("bootstrap requires positive nsim, targets and 0 < level < 1"))
      _ <- options.assessmentOptions.validate
      targets <- compileTargets(fit, options.targets)
      original = LmmAssessment.assess(fit, options.assessmentOptions)
      _ <-
        if original.eligibleForBootstrap then Right(())
        else
          Left(
            MixedModelError.InferenceUnavailable(
              "bootstrap_unqualified_original",
              original.bootstrapIneligibilityReason.getOrElse("original fit is not bootstrap eligible")
            )
          )
    yield Prepared(context, targets, original)

  private def execute(
      fit: LmmFit,
      options: BootstrapOptions,
      prepared: Prepared,
      refit: (Int, Vector[Double]) => FitResult[LmmFit]
  ): BootstrapResult =
    val attempts = Vector.tabulate(options.nsim): index =>
      val seed = Simulation.seedFor(options.seed, index.toLong)
      val response = Simulation.response(fit, seed, options.simulationMode)
      val fitted = response.flatMap(y => refit(index, y))
      val assessment = fitted.toOption.map(LmmAssessment.assess(_, options.assessmentOptions))
      val statistics = fitted.flatMap: result =>
        if result.fixedEffects != fit.fixedEffects || result.y != response.toOption.get then
          Left(
            MixedModelError.InferenceUnavailable(
              "bootstrap_refit_mismatch",
              "Refit does not match the simulated response or retained basis."
            )
          )
        else valuesFor(result, prepared.targets)
      BootstrapAttempt(index, seed, response, fitted, assessment, statistics)
    val intervals =
      prepared.targets.zipWithIndex.map((target, index) => strictInterval(target, index, attempts, options.level))
    new BootstrapResult(options, prepared.original, attempts, intervals)

  private def compileTargets(fit: LmmFit, requested: Vector[BootstrapTarget]): FitResult[Vector[Target]] =
    requested.foldLeft[FitResult[Vector[Target]]](Right(Vector.empty)): (acc, target) =>
      acc.flatMap: current =>
        val next = target match
          case BootstrapTarget.Coefficient(name) =>
            fit.feNames.indexOf(name) match
              case i if i >= 0 =>
                Right(
                  Target(
                    s"retained-basis coefficient: $name",
                    Vector.tabulate(fit.beta.length)(j => if i == j then 1.0 else 0.0)
                  )
                )
              case _ if fit.fixedEffects.droppedNames.contains(name) =>
                Left(
                  MixedModelError.InferenceUnavailable(
                    "bootstrap_dropped_coefficient",
                    s"`$name` is dropped from the original fixed-effect basis"
                  )
                )
              case _ => Left(MixedModelError.InvalidArgument(s"unknown retained coefficient `$name`"))
          case BootstrapTarget.Contrast(label, weights) =>
            val unknown = weights.keys.filterNot(fit.feNames.contains).toVector
            if label.isBlank || weights.isEmpty || !weights.values.forall(_.isFinite) then
              Left(MixedModelError.InvalidArgument("contrast needs a label and finite weights"))
            else if unknown.nonEmpty then
              Left(
                MixedModelError.InferenceUnavailable(
                  "bootstrap_nonretained_contrast",
                  s"contrast `$label` references non-retained coefficients: ${unknown.mkString(", ")}"
                )
              )
            else
              Right(Target(s"retained-basis contrast: $label", fit.feNames.map(name => weights.getOrElse(name, 0.0))))
        next.map(current :+ _)

  private def valuesFor(fit: LmmFit, targets: Vector[Target]): FitResult[Vector[Double]] =
    val values = targets.map(t => t.weights.zip(fit.beta).map(_ * _).sum)
    if values.forall(_.isFinite) then Right(values)
    else
      Left(MixedModelError.InferenceUnavailable("bootstrap_nonfinite_statistic", "bootstrap statistic is non-finite"))

  private def strictInterval(
      target: Target,
      targetIndex: Int,
      attempts: Vector[BootstrapAttempt],
      level: Double
  ): FitResult[BootstrapInterval] =
    if attempts.length < 2 then
      Left(
        MixedModelError.InferenceUnavailable(
          "bootstrap_insufficient_replicates",
          "Percentile intervals require at least two replicates."
        )
      )
    else if attempts.exists(!_.qualified) then
      Left(
        MixedModelError.InferenceUnavailable(
          "bootstrap_replicate_unqualified",
          s"At least one replicate failed or is unqualified for ${target.label}; all attempts are retained."
        )
      )
    else
      val values = attempts.map(_.statistics.toOption.get(targetIndex)).sorted
      val tail = (1.0 - level) / 2.0
      val lower = type7(values, tail)
      val upper = type7(values, 1.0 - tail)
      if !lower.isFinite || !upper.isFinite then
        Left(
          MixedModelError.InferenceUnavailable("bootstrap_nonfinite_interval", "Percentile endpoints are non-finite.")
        )
      else Right(BootstrapInterval(target.label, lower, upper, level))

  /** R's default type-7 empirical quantile: linear interpolation between adjacent order statistics. */
  private[stats] def type7(sorted: Vector[Double], probability: Double): Double =
    val h = (sorted.length - 1) * probability
    val lo = math.floor(h).toInt
    val hi = math.ceil(h).toInt
    if lo == hi || sorted(lo) == sorted(hi) then sorted(lo)
    else
      val fraction = h - lo
      val low = sorted(lo)
      val high = sorted(hi)
      val value =
        if low < 0 && high > 0 then (1.0 - fraction) * low + fraction * high else low + fraction * (high - low)
      math.max(low, math.min(high, value))

  def parametric(fit: LmmFit, nsim: Int): FitResult[Nothing] =
    val _ = (fit, nsim)
    Left(
      MixedModelError.InferenceUnavailable(
        "bootstrap_seed_required",
        "Use BootstrapOptions with an explicit seed and retained-basis targets."
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
