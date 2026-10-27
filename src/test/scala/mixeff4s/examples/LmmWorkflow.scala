package mixeff4s.examples

import mixeff4s.prelude.*

/** Runnable public-API example: sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.LmmWorkflow'. */
object LmmWorkflow:
  final case class Output(
      fit: LmmFit,
      assessment: LmmAssessment,
      predictions: Vector[Double],
      residuals: Vector[Double],
      simulated: Vector[Double],
      bootstrap: BootstrapResult,
      interval: FitResult[BootstrapInterval]
  )

  def run(frame: ModelFrame = mixeff4s.fixtures.Sleepstudy.frame): FitResult[Output] =
    val controls = TrustBqOptions(
      finalRadius = 1e-10,
      ftolAbs = 0,
      ftolRel = 0,
      stallFtolAbs = 0,
      stallFtolRel = 0,
      stallIterations = 20,
      gradientTolerance = 1e-9
    )
    val options = FitOptions.ml.copy(optimizer = LmmOptimizer.AnalyticGradient, optimizerControls = Some(controls))
    for
      plan <- Lmm.prepare("reaction ~ 1 + days + (1 | subj)", frame)
      fit <- plan.fit(options)
      assessment = LmmAssessment.assess(fit)
      newData <- ModelFrame.of("days" -> numeric(Vector(0.0, 5.0, 10.0)))
      predictions <- Prediction.predict(fit, newData, PredictionMode.Population)
      residuals <- Prediction.residuals(fit)
      simulated <- Simulation.response(fit, seed = 42L, mode = SimulationMode.DrawRandomEffects)
      bootstrap <- Bootstrap.parametric(
        fit,
        BootstrapOptions(39, seed = 20261026L, targets = Vector(BootstrapTarget.retainedBasisCoefficient("days")))
      )
    yield Output(fit, assessment, predictions, residuals, simulated, bootstrap, bootstrap.intervals.head)

  def main(args: Array[String]): Unit =
    run() match
      case Left(error)   => throw IllegalStateException(s"${error.code}: ${error.message}")
      case Right(result) =>
        println(s"Gradient assessment: ${result.assessment.gradient}")
        println(s"Population predictions: ${result.predictions.mkString(", ")}")
        println(s"Qualified bootstrap refits: ${result.bootstrap.qualifiedCount}/${result.bootstrap.attempts.length}")
        println(s"Failed or unqualified indices: ${result.bootstrap.failedOrUnqualifiedIndices.mkString(", ")}")
        result.interval match
          case Right(interval) =>
            println(
              s"95% percentile interval for retained-basis days coefficient: [${interval.lower}, ${interval.upper}]"
            )
          case Left(error) => println(s"Interval refused: ${error.code}: ${error.message}")
