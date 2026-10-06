package mixeff4s.lmm

import mixeff4s.comparison.{Datasets, FitCompare}
import mixeff4s.pathology.{FitStatus, Pathology}

class SimpleLmmFitSuite extends munit.FunSuite:
  test("rail REML matches frozen lme4"):
    assertMatches("rail")

  test("ergostool REML matches frozen lme4"):
    assertMatches("ergostool")

  test("cake REML matches frozen lme4 after aligning contrast names"):
    assertMatches("cake")

  private def assertMatches(dataset: String): Unit =
    val row = FitCompare.row(dataset, "REML")
    val ref = FitCompare.claimed(row)
    val frame = Datasets.load(dataset).frame
    val design = Lmm.compile(row.key.formula, frame).getOrElse(fail(s"$dataset compile"))
    val fit = Lmm.fit(row.key.formula, frame, FitOptions.reml).fold(err => fail(err.message), identity)
    val assessed = Pathology.assessFit(Pathology.certify(design), fit.theta, design.parmap)
    assertEquals(assessed.fitStatus, FitStatus.ConvergedInterior, clues(dataset))
    ref.objective.foreach: obj =>
      assertEqualsDouble(fit.objective, obj, FitCompare.objectiveTol(row, ref))
    ref.sigma.foreach: sigma =>
      assertEqualsDouble(fit.sigma, sigma, FitCompare.sigmaTol(row, ref))
    ref.theta.foreach: theta =>
      assertEquals(fit.theta.length, theta.length, clues(dataset))
      fit.theta.zip(theta).zipWithIndex.foreach: (pair, i) =>
        val (got, expected) = pair
        assertEqualsDouble(got, expected, FitCompare.thetaTol(row, ref), clues(dataset, i))
    ref.beta.foreach: beta =>
      val names = ref.coefNames.getOrElse(fail(s"$dataset missing coef_names"))
      val aligned = FitCompare.alignedBeta(fit.feNames, fit.beta, names, beta)
      aligned.zipWithIndex.foreach: (triple, i) =>
        val (name, got, expected) = triple
        assertEqualsDouble(got, expected, FitCompare.betaTol(row, ref), clues(dataset, name, i))
