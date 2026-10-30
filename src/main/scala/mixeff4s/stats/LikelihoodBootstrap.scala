package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.{GlmmFit, GlmmSimulation}
import mixeff4s.lmm.{LmmFit, Simulation, Criterion}

final case class LikelihoodBootstrapOptions(
    nsim: Int,
    seed: Long,
    maxEvaluations: Option[Int] = None,
    monteCarloLevel: Double = .95,
    statisticTolerance: Double = 1e-7
)
final case class LikelihoodBootstrapAttempt[A](
    index: Int,
    seed: Long,
    response: FitResult[Vector[Double]],
    nullFit: FitResult[A],
    alternativeFit: FitResult[A],
    statistic: FitResult[LrtStatistic]
):
  def qualified: Boolean = statistic.isRight

/** Plus-one Monte Carlo tail estimate and Wilson uncertainty for the underlying empirical exceedance probability. */
final case class MonteCarloTail(
    exceedances: Int,
    replicates: Int,
    pvalue: Double,
    empiricalTail: Double,
    monteCarloStandardError: Double,
    lower: Double,
    upper: Double,
    level: Double
)
final class LikelihoodBootstrapResult[A] private[stats] (
    val observed: LrtStatistic,
    val options: LikelihoodBootstrapOptions,
    val attempts: Vector[LikelihoodBootstrapAttempt[A]],
    val calibration: FitResult[MonteCarloTail]
):
  def failedOrUnqualifiedIndices: Vector[Int] = attempts.filterNot(_.qualified).map(_.index)

object LikelihoodBootstrap:
  def compare(
      smaller: GlmmFit,
      larger: GlmmFit,
      options: LikelihoodBootstrapOptions
  ): FitResult[LikelihoodBootstrapResult[GlmmFit]] =
    for
      _ <- validate(options)
      observed <- GlmmLrt.statistic(smaller, larger)
      smallOptions = smaller.options.copy(maxEvaluations =
        options.maxEvaluations.getOrElse(smaller.options.maxEvaluations)
      )
      largeOptions = larger.options.copy(maxEvaluations =
        options.maxEvaluations.getOrElse(larger.options.maxEvaluations)
      )
      small <- smaller.context.compiled.session(smallOptions)
      large <- larger.context.compiled.session(largeOptions)
    yield execute(
      observed,
      options,
      seed => GlmmSimulation.response(smaller, seed),
      (_, y) => small.refit(y),
      (_, y) => large.refit(y),
      (a, b, y) =>
        if a.y != y || b.y != y || a.options != smallOptions || b.options != largeOptions ||
          !(a.context.design eq smaller.context.design) || !(b.context.design eq larger.context.design)
        then unavailable("lrt_refit_mismatch", "Replicate does not match response, controls or frozen design.")
        else GlmmLrt.statistic(a, b)
    )

  def compare(
      smaller: LmmFit,
      larger: LmmFit,
      options: LikelihoodBootstrapOptions
  ): FitResult[LikelihoodBootstrapResult[LmmFit]] =
    for
      _ <- validate(options)
      _ <-
        if smaller.options.criterion == Criterion.ML && larger.options.criterion == Criterion.ML then Right(())
        else unavailable("lrt_requires_ml", "Explicit ML refits are required before null bootstrap comparison.")
      a <- smaller.validatedContext
      b <- larger.validatedContext
      _ <- LikelihoodQualification.lmm(smaller)
      _ <- LikelihoodQualification.lmm(larger)
      observed <- Lrt.statistic(smaller, larger)
      smallOptions = smaller.options.copy(maxEvaluations =
        options.maxEvaluations.getOrElse(smaller.options.maxEvaluations)
      )
      largeOptions = larger.options.copy(maxEvaluations =
        options.maxEvaluations.getOrElse(larger.options.maxEvaluations)
      )
      small <- a.compiled.session(smallOptions)
      large <- b.compiled.session(largeOptions)
    yield execute(
      observed,
      options,
      seed => Simulation.response(smaller, seed),
      (_, y) => small.refit(y),
      (_, y) => large.refit(y),
      (x, z, y) =>
        if x.y != y || z.y != y || x.options != smallOptions || z.options != largeOptions ||
          !x.validatedContext.toOption.exists(_.design eq a.design) || !z.validatedContext.toOption
            .exists(_.design eq b.design)
        then unavailable("lrt_refit_mismatch", "Replicate does not match response, controls or frozen design.")
        else
          for
            _ <- LikelihoodQualification.lmm(x)
            _ <- LikelihoodQualification.lmm(z)
            value <- Lrt.statistic(x, z)
          yield value
    )

  /** Test seam also used by independent witnesses; every response/pair is retained. */
  private[stats] def execute[A](
      observed: LrtStatistic,
      options: LikelihoodBootstrapOptions,
      simulate: Long => FitResult[Vector[Double]],
      small: (Int, Vector[Double]) => FitResult[A],
      large: (Int, Vector[Double]) => FitResult[A],
      statistic: (A, A, Vector[Double]) => FitResult[LrtStatistic]
  ): LikelihoodBootstrapResult[A] =
    val attempts = Vector.tabulate(options.nsim): index =>
      val seed = Simulation.seedFor(options.seed, index.toLong)
      val response = simulate(seed)
      // Both models are attempted even when the other fit fails.
      val a = response.flatMap(y => small(index, y))
      val b = response.flatMap(y => large(index, y))
      val value = for
        y <- response
        x <- a
        z <- b
        v <- statistic(x, z, y)
      yield v
      LikelihoodBootstrapAttempt(index, seed, response, a, b, value)
    val calibration =
      if attempts.exists(!_.qualified) then
        unavailable(
          "lrt_bootstrap_unqualified",
          "At least one simulated model pair failed qualification; all pairs are retained."
        )
      else
        Right(
          tail(
            attempts.count(
              _.statistic.toOption.get.likelihoodRatio >= observed.likelihoodRatio - options.statisticTolerance
            ),
            options.nsim,
            options.monteCarloLevel
          )
        )
    new LikelihoodBootstrapResult(observed, options, attempts, calibration)

  private[stats] def tail(k: Int, n: Int, level: Double): MonteCarloTail =
    val p = k.toDouble / n
    val z = math.sqrt(InferenceQuantiles.chiSquare1(level)); val zz = z * z
    val center = (p + zz / (2 * n)) / (1 + zz / n)
    val radius = z * math.sqrt(p * (1 - p) / n + zz / (4.0 * n * n)) / (1 + zz / n)
    MonteCarloTail(
      k,
      n,
      (k + 1.0) / (n + 1.0),
      p,
      math.sqrt(p * (1 - p) / n),
      math.max(0, center - radius),
      math.min(1, center + radius),
      level
    )

  private def validate(options: LikelihoodBootstrapOptions): FitResult[Unit] =
    if options.nsim < 1 || !options.monteCarloLevel.isFinite || options.monteCarloLevel <= 0 || options.monteCarloLevel >= 1 ||
      options.maxEvaluations.exists(_ < 1) || !options.statisticTolerance.isFinite || options.statisticTolerance < 0
    then
      Left(
        MixedModelError.InvalidArgument(
          "Invalid null-bootstrap count, budget, Monte Carlo level or statistic tolerance."
        )
      )
    else Right(())
  private def unavailable(code: String, message: String): FitResult[Nothing] = Left(
    MixedModelError.InferenceUnavailable(code, message)
  )
