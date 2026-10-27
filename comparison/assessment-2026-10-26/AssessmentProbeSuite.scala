package mixeff4s.assessment

import mixeff4s.compiler.Compiler
import mixeff4s.data.ModelFrame
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.formula.Formula
import mixeff4s.glmm.Glmm
import mixeff4s.lmm.{FitOptions, Lmm}
import mixeff4s.model.Family
import mixeff4s.pathology.FitStatus
import mixeff4s.stats.Lrt

/** Observational probes: passing confirms the recorded defect, not correct behavior. */
class AssessmentProbeSuite extends munit.FunSuite:
  test("an unfitted theta is labelled converged"):
    val design = Lmm.compile("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame).toOption.get
    val assessed = Compiler.assess(Compiler.artifact(design), Vector(12345.0))
    println(s"ASSESSMENT unfitted theta status=${assessed.pathology.fitStatus}")
    assertEquals(assessed.pathology.fitStatus, FitStatus.ConvergedInterior)

  test("invalid GLMM theta becomes successful finite penalty"):
    val frame = ModelFrame.of(
      "y" -> ModelFrame.numeric(Vector(0.0, 1.0, 1.0, 0.0, 0.0, 1.0)),
      "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c"))
    ).toOption.get
    val value = Glmm.profiledDeviance(Formula.parse("y ~ 1 + (1 | g)").toOption.get,
      frame, Family.Bernoulli, Vector.empty)
    println(s"ASSESSMENT wrong theta dimension=$value")
    assertEquals(value, Right(1e300))

  test("dropped middle column misaligns standard errors"):
    val frame = Sleepstudy.frame.addNumeric("zero", Vector.fill(Sleepstudy.frame.nRows)(0.0)).toOption.get
    val fit = Lmm.fit("reaction ~ 1 + zero + days + (1 | subj)", frame, FitOptions.ml).toOption.get
    println(s"ASSESSMENT names=${fit.feNames} beta=${fit.beta} se=${fit.stderror} z=${fit.coefTable.zValues}")
    assertEquals(fit.feNames.length, 2)
    assertEquals(fit.stderror.length, 3)
    assert(fit.coefTable.zValues(1).isNaN)

  test("LRT accepts different nonnested predictor values under the same name"):
    val changed = ModelFrame.of(Sleepstudy.frame.columns.map { case (name, col) =>
      name -> (if name == "days" then ModelFrame.numeric(Sleepstudy.frame.numeric("days").get.map(x => x*x)) else col)
    }*).toOption.get
    val smaller = Lmm.fit("reaction ~ 1 + days + (1 | subj)", changed, FitOptions.ml).toOption.get
    val expanded = Sleepstudy.frame.addNumeric("cube", Sleepstudy.frame.numeric("days").get.map(x => x*x*x)).toOption.get
    val larger = Lmm.fit("reaction ~ 1 + days + cube + (1 | subj)", expanded, FitOptions.ml).toOption.get
    val result = Lrt.compare(smaller, larger)
    println(s"ASSESSMENT nonnested same-name LRT=$result")
    assert(result.isRight)

  test("LRT refuses the genuinely nested random intercept to random slope comparison"):
    val smaller = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml).toOption.get
    val larger = Lmm.fit("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame, FitOptions.ml).toOption.get
    val result = Lrt.compare(smaller, larger)
    println(s"ASSESSMENT valid nested RE comparison=$result")
    assert(result.left.toOption.exists(_.message.contains("random-effect term structures are not nested")))
