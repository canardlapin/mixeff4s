package mixeff4s.lmm

import mixeff4s.comparison.FitCompare
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

  test("sleepstudy VarCorr and Wald SE match MixedModels.jl"):
    val fit = Lmm.fit(formula, Sleepstudy.frame, FitOptions.ml).getOrElse(fail("fit"))
    assertEquals(fit.varcorr.components.length, 1)
    val re = fit.varcorr.components.head
    assertEquals(re.group, "subj")
    assertEquals(re.names, Vector("(Intercept)", "days"))
    assertEqualsDouble(re.stdDev(0), 23.78066438213187, 0.1)
    assertEqualsDouble(re.stdDev(1), 5.7168446983832775, 0.1)
    assertEquals(re.correlations.length, 1)
    assertEqualsDouble(re.correlations(0), 0.0813, 0.01)
    assertEqualsDouble(fit.varcorr.residualSd, fit.sigma, 1e-12)
    val se = fit.stderror
    assertEqualsDouble(se(0), 6.632295312722272, 0.01)
    assertEqualsDouble(se(1), 1.5022387911441102, 0.01)
    val feCorr = fit.vcov(0)(1) / (se(0) * se(1))
    assertEqualsDouble(feCorr, -0.13755599049585931, 0.01)
    val table = fit.coefTable
    assertEquals(table.names, Vector("(Intercept)", "days"))
    assertEquals(table.pValueCode, "p_value_unavailable")
    assertEqualsDouble(table.zValues(0), fit.beta(0) / se(0), 1e-12)
    assertEqualsDouble(fit.logdetRe, 73.90350673367566, 0.1)

  test("sleepstudy zerocorr ML optimizes two diagonal theta slots"):
    val formula = "reaction ~ 1 + days + (1 + days || subj)"
    val design = Lmm.compile(formula, Sleepstudy.frame).getOrElse(fail("compile"))
    assertEquals(design.nTheta, 2)
    assertEquals(design.parmap, Vector((0, 0, 0), (0, 1, 1)))
    assertEquals(design.theta, Vector(1.0, 1.0))
    val workspace = PlsWorkspace(design, reml = false).getOrElse(fail("workspace"))
    val row = FitCompare.row("sleepstudy", "ML", formula)
    val ref = FitCompare.claimed(row)
    val juliaTheta = ref.theta.getOrElse(fail("frozen theta"))
    val atJulia = workspace.objectiveAt(juliaTheta).getOrElse(fail("objective"))
    ref.objective.foreach(obj => assertEqualsDouble(atJulia, obj, 1e-4))
    val fit = Lmm.fit(formula, Sleepstudy.frame, FitOptions.ml).getOrElse(fail("fit"))
    assertEquals(fit.theta.length, 2)
    ref.objective.foreach(obj => assertEqualsDouble(fit.objective, obj, FitCompare.objectiveTol(row, ref)))
    FitCompare.alignedTheta(fit.theta, juliaTheta, scalarTerms = false).zipWithIndex.foreach: (pair, i) =>
      val (got, expected) = pair
      assertEqualsDouble(got, expected, FitCompare.thetaTol(row, ref), clues(i))
    ref.beta.foreach: beta =>
      val names = ref.coefNames.getOrElse(fail("coef_names"))
      FitCompare.alignedBeta(fit.feNames, fit.beta, names, beta).zipWithIndex.foreach: (triple, i) =>
        val (name, got, expected) = triple
        assertEqualsDouble(got, expected, FitCompare.betaTol(row, ref), clues(name, i))
    ref.sigma.foreach(sigma => assertEqualsDouble(fit.sigma, sigma, FitCompare.sigmaTol(row, ref)))
    val re = fit.varcorr.components.head
    assertEquals(re.names, Vector("(Intercept)", "days"))
    assertEquals(re.correlations.length, 1)
    assertEqualsDouble(re.correlations(0), 0.0, 1e-12)
    assertEqualsDouble(re.stdDev(0), 24.171269957611873, 0.1)
    assertEqualsDouble(re.stdDev(1), 5.79939919963132, 0.1)
    val se = fit.stderror
    assertEqualsDouble(se(0), 6.707646513654387, 0.01)
    assertEqualsDouble(se(1), 1.5193112497954953, 0.01)
    assertEqualsDouble(fit.logdetRe, 74.4694698615524, 0.1)
    val diag = Lmm
      .fit("reaction ~ 1 + days + diag(1 + days | subj)", Sleepstudy.frame, FitOptions.ml)
      .getOrElse(fail("diag fit"))
    assertEquals(diag.theta.length, 2)
    assertEqualsDouble(diag.objective, fit.objective, 1e-8)
    assertEquals(diag.theta, fit.theta)
