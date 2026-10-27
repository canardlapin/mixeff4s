package mixeff4s.lmm
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.optimizer.{TrustBq,ObjectiveGradient}
class WorkflowPrecisionProbe extends munit.FunSuite:
  test("measure explicit tight gradient controls"):
    val formula = "reaction ~ 1 + days + (1 | subj)"
    val opts = FitOptions.ml.copy(optimizer=LmmOptimizer.AnalyticGradient)
    val base = Lmm.fit(formula,Sleepstudy.frame,opts).toOption.get
    val design = base.validatedContext.toOption.get.design
    val controls = Pls.gradientFamilyOptions(1).copy(finalRadius=1e-10,ftolAbs=0,ftolRel=0,stallFtolAbs=0,stallFtolRel=0,stallIterations=20,stallRequiresStableX=true)
    var failed = 0
    for i <- 0 until 39 do
      val y = Simulation.response(base,Simulation.seedFor(20261026L,i.toLong)).toOption.get
      val ws = PlsWorkspace(design,false).toOption.get
      ws.updateResponse(y).toOption.get
      val optimized = TrustBq.minimizeWithGradient(design.theta,ws.lowerBounds,ws.upperBounds,controls)(point=>ws.objectiveAndGradientAt(point).map(v=>ObjectiveGradient(v.objective,v.gradient))).toOption.get
      ws.objectiveAt(optimized.x).toOption.get
      val fit = Lmm.snapshot(PlsFit(ws,optimized),opts)
      val assessed=LmmAssessment.assess(fit)
      if !assessed.eligibleForBootstrap then failed+=1
      println(s"PRECISION index=$i evals=${optimized.fevals} stop=${optimized.stopReason} gradient=${assessed.gradient} eligible=${assessed.eligibleForBootstrap}")
    println(s"PRECISION unqualified=$failed of39")
