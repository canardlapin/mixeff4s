package mixeff4s.stats
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.lmm.*
class WorkflowQualificationProbe extends munit.FunSuite:
  test("measure strict bootstrap endpoint qualification"):
    for criterion <- Vector(Criterion.ML, Criterion.REML) do
      val options = FitOptions(criterion, optimizer = LmmOptimizer.AnalyticGradient)
      val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, options).toOption.get
      val assessed = LmmAssessment.assess(fit)
      println(s"QUALIFY base=$criterion gradient=${assessed.gradient} eligible=${assessed.eligibleForBootstrap}")
      val result = Bootstrap.parametric(fit, BootstrapOptions(39, 20261026L, Vector(BootstrapTarget.Coefficient("days"))))
      result match
        case Left(error) => println(s"QUALIFY refused=${error.code}: ${error.message}")
        case Right(batch) =>
          println(s"QUALIFY attempts=${batch.attempts.length} qualified=${batch.qualifiedCount} intervals=${batch.intervals}")
          batch.attempts.filterNot(_.qualified).foreach(a => println(s"QUALIFY failed=${a.index} refit=${a.refit.left.toOption.map(_.code)} gradient=${a.assessment.map(_.gradient)} reasons=${a.assessment.map(_.reasonCodes)}"))
