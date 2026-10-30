package mixeff4s.examples

import mixeff4s.error.FitResult
import mixeff4s.fixtures.{GlmmAnalysisFrames, Sleepstudy}
import mixeff4s.glmm.*
import mixeff4s.lmm.{Lmm, FitOptions}
import mixeff4s.model.Family
import mixeff4s.optimizer.TrustBqOptions
import mixeff4s.stats.*

/** Chosen-estimator uncertainty, ML profiles and a complete fitted-null comparison ledger. */
object LikelihoodInference:
  final case class Output(
      coefficients: GlmmInferenceResult,
      comparison: LikelihoodBootstrapResult[GlmmFit],
      profiles: Vector[LmmProfileResult]
  )
  val tightML: FitOptions = FitOptions.ml.copy(optimizerControls =
    Some(
      TrustBqOptions(
        finalRadius = 1e-9,
        ftolAbs = 0,
        ftolRel = 0,
        stallFtolAbs = 0,
        stallFtolRel = 0,
        stallIterations = 20
      )
    )
  )
  def run(): FitResult[Output] =
    val options = GlmmOptions.jointAgq(5).copy(trials = Some(GlmmAnalysisFrames.trials))
    for
      large <- Glmm.fit("y ~ 1 + x + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options)
      small <- Glmm.fit("y ~ 1 + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options)
      coefficients <- GlmmInference.coefficients(large)
      comparison <- LikelihoodBootstrap.compare(small, large, LikelihoodBootstrapOptions(39, 20261029L))
      gaussian <- Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, tightML)
      profiles <- Vector(
        LmmProfileTarget.Coefficient("days"),
        LmmProfileTarget.RandomEffectStdDev("subj"),
        LmmProfileTarget.RandomEffectVariance("subj"),
        LmmProfileTarget.ResidualStdDev
      )
        .foldLeft[FitResult[Vector[LmmProfileResult]]](Right(Vector.empty)): (acc, target) =>
          acc.flatMap(prior =>
            LmmProfile
              .confint(
                gaussian,
                target,
                LmmProfileOptions(maxEvaluations = 2000, optimizerControls = tightML.optimizerControls.get)
              )
              .map(prior :+ _)
          )
    yield Output(coefficients, comparison, profiles)
  def main(args: Array[String]): Unit =
    val _ = args
    run() match
      case Left(error)   => throw IllegalStateException(s"${error.code}: ${error.message}")
      case Right(output) =>
        println(
          s"Observed-information intervals (${output.coefficients.estimator}, nAgq=${output.coefficients.nAgq}): ${output.coefficients.intervals}"
        )
        println(
          s"Null-bootstrap observed LR=${output.comparison.observed.likelihoodRatio}; attempts=${output.comparison.attempts.length}; unqualified=${output.comparison.failedOrUnqualifiedIndices}; calibration=${output.comparison.calibration}"
        )
        output.profiles.foreach(p =>
          println(s"ML profile ${p.target}: ${p.lower} / ${p.upper}; evaluations=${p.evaluations.length}")
        )
