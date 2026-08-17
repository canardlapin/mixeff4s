package mixeff4s.lmm

import mixeff4s.comparison.{Datasets, FitCompare}
import mixeff4s.pathology.{FitStatus, Pathology}

class DyestuffFitSuite extends munit.FunSuite:
  private val formula = "Yield ~ 1 + (1 | Batch)"

  test("dyestuff ML matches frozen lme4"):
    assertMatches("dyestuff", FitOptions.ml, boundary = false)

  test("dyestuff REML matches frozen lme4"):
    assertMatches("dyestuff", FitOptions.reml, boundary = false)

  test("dyestuff2 ML is a certified boundary against frozen lme4"):
    assertMatches("dyestuff2", FitOptions.ml, boundary = true)

  test("dyestuff2 REML is a certified boundary against frozen lme4"):
    assertMatches("dyestuff2", FitOptions.reml, boundary = true)

  private def assertMatches(dataset: String, options: FitOptions, boundary: Boolean): Unit =
    val estimator = if options.criterion == Criterion.ML then "ML" else "REML"
    val row = FitCompare.row(dataset, estimator)
    val ref = FitCompare.claimed(row)
    val frame = Datasets.load(dataset).frame
    val design = Lmm.compile(formula, frame).getOrElse(fail(s"$dataset compile"))
    val fit = Lmm.fit(formula, frame, options).fold(err => fail(err.message), identity)
    val assessed = Pathology.assessFit(Pathology.certify(design), fit.theta, design.parmap)
    assertEquals(fit.beta.length, 1)
    assertEquals(fit.theta.length, 1)
    ref.objective.foreach: obj =>
      assertEqualsDouble(fit.objective, obj, FitCompare.objectiveTol(row, ref))
    ref.beta.foreach: beta =>
      assertEqualsDouble(fit.beta(0), beta.head, FitCompare.betaTol(row, ref))
    ref.theta.foreach: theta =>
      assertEqualsDouble(fit.theta(0), theta.head, FitCompare.thetaTol(row, ref))
    ref.sigma.foreach: sigma =>
      assertEqualsDouble(fit.sigma, sigma, FitCompare.sigmaTol(row, ref))
    if boundary then
      assertEquals(assessed.fitStatus, FitStatus.ConvergedBoundary)
      assertEqualsDouble(fit.theta(0), 0.0, 1e-8)
    else assertEquals(assessed.fitStatus, FitStatus.ConvergedInterior)
