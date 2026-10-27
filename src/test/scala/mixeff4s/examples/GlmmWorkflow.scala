package mixeff4s.examples

import mixeff4s.fixtures.GlmmFrames
import mixeff4s.prelude.*

/** Runnable grouped-binomial example with separate inner and outer evidence. */
object GlmmWorkflow:
  def run(): FitResult[GlmmFit] =
    val raw = GlmmFrames.cbpp
    val successes = raw.numeric("incidence").get
    val sizes = raw.numeric("size").get
    for
      frame <- raw.addNumeric("proportion", successes.zip(sizes).map(_ / _))
      formula <- Formula.parse("proportion ~ 1 + period + (1 | herd)").left.map(MixedModelError.Formula.apply)
      fit <- Glmm.fit(
        formula,
        frame,
        Family.Binomial,
        options = GlmmOptions.jointLaplace.copy(trials = Some(sizes.map(_.toInt)))
      )
    yield fit

  def main(args: Array[String]): Unit =
    val _ = args
    run() match
      case Left(error) => throw IllegalStateException(s"${error.code}: ${error.message}")
      case Right(fit)  =>
        println(s"${fit.algorithmLabel}: ${fit.objectiveDefinition} = ${fit.deviance}")
        println(s"Numerical termination: ${fit.optimizer.stopReason}; evaluations=${fit.optimizer.fevals}")
        println(s"Conditional evidence: ${fit.inner}")
        println(s"Outer assessment: ${fit.jointAssessment}")
        println(s"Coefficients: ${fit.feNames.zip(fit.beta).mkString(", ")}")
        println(s"Conditional means (first four): ${fit.conditionalMean.take(4).mkString(", ")}")
