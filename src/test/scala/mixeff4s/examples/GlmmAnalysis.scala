package mixeff4s.examples

import mixeff4s.fixtures.GlmmAnalysisFrames
import mixeff4s.prelude.*

/** Complete chosen-estimator GLMM workflow, with all bootstrap outcomes retained. */
object GlmmAnalysis:
  final case class Output(
      fit: GlmmFit,
      assessment: GlmmAssessment,
      prediction: Vector[Double],
      residuals: Vector[Double],
      simulated: Vector[Double],
      bootstrap: GlmmBootstrapResult,
      interval: FitResult[BootstrapInterval]
  )
  val formula: String = "y ~ 1 + x + (1 | g)"
  def fitOptions: GlmmOptions = GlmmOptions.jointAgq(5).copy(trials = Some(GlmmAnalysisFrames.trials))
  def bootstrapOptions: GlmmBootstrapOptions = GlmmBootstrapOptions(
    39,
    20261027L,
    Vector(BootstrapTarget.Coefficient("x"), BootstrapTarget.Contrast("two units of x", Map("x" -> 2.0)))
  )

  def run(): FitResult[Output] =
    for
      compiled <- Glmm.prepare(formula, GlmmAnalysisFrames.binomial, Family.Binomial)
      fit <- compiled.fit(fitOptions)
      assessment = GlmmAssessment.assess(fit)
      newData <- ModelFrame.of("x" -> ModelFrame.numeric(Vector(-1.0, 0.0, 1.0)))
      predictions <- GlmmPrediction.predict(fit, newData, PredictionMode.Population)
      residuals <- GlmmPrediction.residuals(fit)
      simulated <- GlmmSimulation.response(fit, 42L)
      bootstrap <- Bootstrap.parametric(fit, bootstrapOptions)
    yield Output(fit, assessment, predictions, residuals, simulated, bootstrap, bootstrap.intervals.head)

  def main(args: Array[String]): Unit =
    val _ = args
    run() match
      case Left(error)   => throw IllegalStateException(s"${error.code}: ${error.message}")
      case Right(output) =>
        println(s"${output.fit.algorithmLabel}, nAgq=${output.fit.nAgq}, objective=${output.fit.deviance}")
        println(s"Assessment: ${output.assessment}")
        println(s"Predictions at zero random effects: ${output.prediction.mkString(", ")}")
        println(s"Qualified bootstrap refits: ${output.bootstrap.qualifiedCount}/${output.bootstrap.attempts.length}")
        println(s"Failed or unqualified indices: ${output.bootstrap.failedOrUnqualifiedIndices.mkString(", ")}")
        output.bootstrap.attempts
          .filterNot(_.qualified)
          .foreach(a =>
            println(s"Attempt ${a.index}: ${a.refit.left.toOption.map(_.code)}; ${a.assessment.map(_.refusalCodes)}")
          )
        output.interval match
          case Right(interval) =>
            println(s"First-order-qualified percentile interval: [${interval.lower}, ${interval.upper}]")
          case Left(error) => println(s"Interval refusal ${error.code}: ${error.message}")
