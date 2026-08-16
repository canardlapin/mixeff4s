package mixeff4s.lmm

import mixeff4s.fixtures.Sleepstudy
import mixeff4s.formula.Formula

class SleepstudyFitSuite extends munit.FunSuite:
  private val formula =
    Formula.parse("reaction ~ 1 + days + (1 + days | subj)").getOrElse(fail("parse"))

  test("objective at Julia ML theta matches MixedModels.jl"):
    val design = Lmm.compile(formula, Sleepstudy.frame).getOrElse(fail("compile"))
    val workspace = PlsWorkspace(design, reml = false).getOrElse(fail("workspace"))
    val juliaTheta = Vector(0.9292145007480422, 0.018171743703656894, 0.222646742690827)
    val obj = workspace.objectiveAt(juliaTheta).getOrElse(fail("objective"))
    assertEqualsDouble(obj, 1751.939344474405, 1e-4)
    assertEqualsDouble(workspace.beta(0), 251.40510484848537, 1e-3)
    assertEqualsDouble(workspace.beta(1), 10.4672859595958, 1e-3)
    assertEqualsDouble(workspace.sigma, 25.591813564885108, 1e-3)

  test("sleepstudy ML fit matches MixedModels.jl"):
    val fit = Lmm.fit(formula, Sleepstudy.frame, FitOptions.ml).getOrElse(fail("fit"))
    assertEqualsDouble(fit.objective, 1751.9393444636682, 1e-2)
    assertEquals(fit.theta.length, 3)
    assertEqualsDouble(fit.theta(0), 0.9292297167514472, 1e-3)
    assertEqualsDouble(fit.theta(1), 0.01816466496782548, 1e-3)
    assertEqualsDouble(fit.theta(2), 0.22264601131030412, 1e-3)
    assertEqualsDouble(fit.beta(0), 251.40510484848454, 1e-3)
    assertEqualsDouble(fit.beta(1), 10.467285959596126, 1e-3)
    assertEqualsDouble(fit.sigma, 25.591813564885108, 1e-3)

  test("sleepstudy REML fit matches MixedModels.jl"):
    val fit = Lmm.fit(formula, Sleepstudy.frame, FitOptions.reml).getOrElse(fail("fit"))
    assertEqualsDouble(fit.objective, 1743.6282719599442, 1e-2)
    assertEqualsDouble(fit.theta(0), 0.9667417690560796, 1e-3)
    assertEqualsDouble(fit.theta(1), 0.015169059384716037, 1e-3)
    assertEqualsDouble(fit.theta(2), 0.2309099529619309, 1e-3)
    assertEqualsDouble(fit.beta(0), 251.40510484848528, 1e-3)
    assertEqualsDouble(fit.beta(1), 10.467285959595493, 1e-3)
    assertEqualsDouble(fit.sigma, 25.591795732317802, 1e-3)
