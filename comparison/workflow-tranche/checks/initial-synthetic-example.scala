package mixeff4s.examples

import mixeff4s.prelude.*

/** Runnable public-API example: sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.LmmWorkflow'. */
object LmmWorkflow:
  final case class Output(fit:LmmFit,assessment:LmmAssessment,predictions:Vector[Double],residuals:Vector[Double],simulated:Vector[Double],bootstrap:BootstrapResult,interval:BootstrapInterval)

  def run():FitResult[Output]=
    val n=180
    val days=Vector.tabulate(n)(i=>(i%10).toDouble)
    val groups=Vector.tabulate(n)(i=>s"subject${i/10}")
    val response=Vector.tabulate(n)(i=>240.0+10.0*days(i)+15.0*math.sin((i/10)*1.3)+5.0*math.cos(i*2.1))
    val controls=TrustBqOptions(finalRadius=1e-10,ftolAbs=0,ftolRel=0,stallFtolAbs=0,stallFtolRel=0,stallIterations=20,gradientTolerance=1e-9)
    val options=FitOptions.ml.copy(optimizer=LmmOptimizer.AnalyticGradient,optimizerControls=Some(controls))
    for
      frame<-ModelFrame.of("reaction"->numeric(response),"days"->numeric(days),"subj"->factorCol(groups))
      plan<-Lmm.prepare("reaction ~ 1 + days + (1 | subj)",frame)
      fit<-plan.fit(options)
      assessment=LmmAssessment.assess(fit)
      newData<-ModelFrame.of("days"->numeric(Vector(0.0,5.0,10.0)))
      predictions<-Prediction.predict(fit,newData,PredictionMode.Population)
      residuals<-Prediction.residuals(fit)
      simulated<-Simulation.response(fit,seed=42L,mode=SimulationMode.DrawRandomEffects)
      bootstrap<-Bootstrap.parametric(fit,BootstrapOptions(39,seed=20261026L,targets=Vector(BootstrapTarget.retainedBasisCoefficient("days"))))
      interval<-bootstrap.intervals.head
    yield Output(fit,assessment,predictions,residuals,simulated,bootstrap,interval)

  def main(args:Array[String]):Unit=
    run() match
      case Left(error)=>throw IllegalStateException(s"${error.code}: ${error.message}")
      case Right(result)=>
        println(s"Gradient assessment: ${result.assessment.gradient}")
        println(s"Population predictions: ${result.predictions.mkString(", ")}")
        println(s"Qualified bootstrap refits: ${result.bootstrap.qualifiedCount}/${result.bootstrap.attempts.length}")
        println(s"95% percentile interval for retained-basis days coefficient: [${result.interval.lower}, ${result.interval.upper}]")
