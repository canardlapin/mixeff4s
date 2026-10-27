package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.{GlmmAssessment, GlmmAssessmentOptions, GlmmFit, GlmmOptions, GlmmSimulation}
import mixeff4s.lmm.{Simulation, SimulationMode}

final case class GlmmBootstrapOptions(
    nsim: Int,
    seed: Long,
    targets: Vector[BootstrapTarget],
    simulationMode: SimulationMode = SimulationMode.DrawRandomEffects,
    fitOptions: Option[GlmmOptions] = None,
    assessmentOptions: GlmmAssessmentOptions = GlmmAssessmentOptions(),
    level: Double = 0.95
)

final case class GlmmBootstrapAttempt(
    index: Int,
    seed: Long,
    response: FitResult[Vector[Double]],
    refit: FitResult[GlmmFit],
    assessment: Option[GlmmAssessment],
    statistics: FitResult[Vector[Double]]
):
  def qualified: Boolean = refit.isRight && statistics.isRight && assessment.exists(_.eligibleForBootstrap)

/** Complete immutable replicate ledger. Intervals never omit or replace failed attempts. */
final class GlmmBootstrapResult private[stats] (
    val options: GlmmBootstrapOptions,
    val originalAssessment: GlmmAssessment,
    val attempts: Vector[GlmmBootstrapAttempt],
    val intervals: Vector[FitResult[BootstrapInterval]]
):
  def qualifiedCount: Int = attempts.count(_.qualified)
  def failedOrUnqualifiedIndices: Vector[Int] = attempts.filterNot(_.qualified).map(_.index)

object GlmmBootstrap:
  def parametric(fit: GlmmFit, options: GlmmBootstrapOptions): FitResult[GlmmBootstrapResult] =
    for
      prepared <- prepare(fit, options)
      session <- fit.context.compiled.session(prepared.controls)
    yield execute(fit, options, prepared, (_, response) => session.refit(response))

  private[stats] def withRefitter(fit: GlmmFit, options: GlmmBootstrapOptions)(
      refit: (Int, Vector[Double]) => FitResult[GlmmFit]
  ): FitResult[GlmmBootstrapResult] = prepare(fit, options).map(execute(fit, options, _, refit))

  private final case class Prepared(
      targets: Vector[Bootstrap.Target],
      controls: GlmmOptions,
      original: GlmmAssessment
  )

  private def prepare(fit: GlmmFit, options: GlmmBootstrapOptions): FitResult[Prepared] =
    val controls = options.fitOptions.getOrElse(fit.options)
    for
      _ <-
        if options.nsim > 0 && options.targets.nonEmpty && options.level.isFinite && options.level > 0.0 && options.level < 1.0
        then Right(())
        else Left(MixedModelError.InvalidArgument("GLMM bootstrap requires positive nsim, targets and 0 < level < 1"))
      _ <- options.assessmentOptions.validate
      _ <-
        if sameResponseContract(fit.options, controls) then Right(())
        else
          Left(
            MixedModelError.InferenceUnavailable(
              "bootstrap_estimator_mismatch",
              "Bootstrap controls must preserve estimator, quadrature, weights, trials and offsets."
            )
          )
      targets <- Bootstrap.compileTargets(fit.fixedEffects, options.targets)
      _ <-
        if targets.map(_.label).distinct.length == targets.length then Right(())
        else Left(MixedModelError.InvalidArgument("Bootstrap target labels must be distinct"))
      original = GlmmAssessment.assess(fit, options.assessmentOptions)
      _ <-
        if original.eligibleForBootstrap then Right(())
        else
          Left(
            MixedModelError.InferenceUnavailable("bootstrap_unqualified_original", original.refusalCodes.mkString(", "))
          )
    yield Prepared(targets, controls, original)

  private def sameResponseContract(a: GlmmOptions, b: GlmmOptions): Boolean =
    a.fast == b.fast && a.nAgq == b.nAgq && a.weights == b.weights && a.trials == b.trials && a.offset == b.offset

  private def execute(
      fit: GlmmFit,
      options: GlmmBootstrapOptions,
      prepared: Prepared,
      refit: (Int, Vector[Double]) => FitResult[GlmmFit]
  ): GlmmBootstrapResult =
    val attempts = Vector.tabulate(options.nsim): index =>
      val seed = Simulation.seedFor(options.seed, index.toLong)
      val response = GlmmSimulation.response(fit, seed, options.simulationMode)
      val fitted = response.flatMap(y => refit(index, y))
      val assessment = fitted.toOption.map(GlmmAssessment.assess(_, options.assessmentOptions))
      val statistics = fitted.flatMap: result =>
        response.flatMap: simulated =>
          if result.y != simulated || result.options != prepared.controls || result.family != fit.family ||
            result.link != fit.link || result.fixedEffects != fit.fixedEffects ||
            !(result.context.design eq fit.context.design)
          then
            Left(
              MixedModelError.InferenceUnavailable(
                "bootstrap_refit_mismatch",
                "Refit does not match the simulated response, estimator, controls or frozen design."
              )
            )
          else Bootstrap.valuesFor(result.beta, prepared.targets)
      GlmmBootstrapAttempt(index, seed, response, fitted, assessment, statistics)
    val intervals = prepared.targets.zipWithIndex.map: (target, j) =>
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
            "At least one GLMM replicate failed or is unqualified; all attempts are retained."
          )
        )
      else
        val sorted = attempts.map(_.statistics.toOption.get(j)).sorted
        val tail = (1.0 - options.level) / 2.0
        val lower = Bootstrap.type7(sorted, tail)
        val upper = Bootstrap.type7(sorted, 1.0 - tail)
        if lower.isFinite && upper.isFinite then Right(BootstrapInterval(target.label, lower, upper, options.level))
        else
          Left(
            MixedModelError.InferenceUnavailable(
              "bootstrap_nonfinite_interval",
              "GLMM percentile endpoints are non-finite."
            )
          )
    new GlmmBootstrapResult(options, prepared.original, attempts, intervals)
