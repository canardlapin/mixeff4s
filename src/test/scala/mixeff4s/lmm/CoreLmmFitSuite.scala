package mixeff4s.lmm

import mixeff4s.comparison.{Datasets, FitCompare}
import mixeff4s.pathology.{FitStatus, Pathology}

class CoreLmmFitSuite extends munit.FunSuite:
  test("oats REML matches frozen lme4"):
    assertMatches("oats", "yield ~ 1 + Variety * nitro + (1 | Block) + (1 | Block:Variety)")

  test("orthodont REML matches frozen lme4"):
    assertMatches("orthodont", "distance ~ 1 + age * Sex + (1 + age | Subject)")

  test("oxide nested intercepts match MixedModels.jl"):
    assertMatches("oxide", "Thickness ~ 1 + (1 | Lot) + (1 | Lot:Wafer)")

  test("oxide nesting sugar matches the Lot + Lot:Wafer fit"):
    val explicit = Lmm
      .fit(
        "Thickness ~ 1 + (1 | Lot) + (1 | Lot:Wafer)",
        Datasets.load("oxide").frame,
        FitOptions.reml
      )
      .fold(err => fail(err.message), identity)
    val sugar = Lmm
      .fit("Thickness ~ 1 + (1 | Lot/Wafer)", Datasets.load("oxide").frame, FitOptions.reml)
      .fold(err => fail(err.message), identity)
    assertEqualsDouble(sugar.objective, explicit.objective, 1e-8)
    assertEquals(sugar.theta.sorted, explicit.theta.sorted)
    assertMatches("oxide", "Thickness ~ 1 + (1 | Lot/Wafer)")

  test("machines worker-machine intercepts match frozen lme4"):
    assertMatches("machines", "score ~ 1 + Machine + (1 | Worker) + (1 | Worker:Machine)")

  test("machines categorical random slope matches frozen lme4"):
    assertMatches("machines", "score ~ 1 + Machine + (1 + Machine | Worker)")

  private def assertMatches(dataset: String, formula: String): Unit =
    val row = FitCompare.row(dataset, "REML", formula)
    val ref = FitCompare.claimed(row)
    val frame = Datasets.load(dataset).frame
    val design = Lmm.compile(formula, frame).getOrElse(fail(s"$dataset compile"))
    val fit = Lmm.fit(formula, frame, FitOptions.reml).fold(err => fail(err.message), identity)
    val assessed = Pathology.assessFit(Pathology.certify(design), fit.theta, design.parmap)
    assertEquals(assessed.fitStatus, FitStatus.ConvergedInterior, clues(dataset, formula))
    ref.objective.foreach: obj =>
      assertEqualsDouble(fit.objective, obj, FitCompare.objectiveTol(row, ref))
    ref.sigma.foreach: sigma =>
      assertEqualsDouble(fit.sigma, sigma, FitCompare.sigmaTol(row, ref))
    ref.theta.foreach: theta =>
      val scalarTerms =
        fit.varcorr.components.forall(_.stdDev.length == 1) &&
          fit.varcorr.components.length == fit.theta.length
      FitCompare.alignedTheta(fit.theta, theta, scalarTerms).zipWithIndex.foreach: (pair, i) =>
        val (got, expected) = pair
        assertEqualsDouble(got, expected, FitCompare.thetaTol(row, ref), clues(dataset, formula, i))
    ref.beta.foreach: beta =>
      val names = ref.coefNames.getOrElse(fail(s"$dataset missing coef_names"))
      FitCompare.alignedBeta(fit.feNames, fit.beta, names, beta).zipWithIndex.foreach: (triple, i) =>
        val (name, got, expected) = triple
        assertEqualsDouble(got, expected, FitCompare.betaTol(row, ref), clues(dataset, name, i))
