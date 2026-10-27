package mixeff4s.lmm

import mixeff4s.comparison.{Datasets, FitCompare}
import mixeff4s.error.MixedModelError
import mixeff4s.pathology.{FitStatus, Pathology}

class StationFitSuite extends munit.FunSuite:
  test("station intercepts match frozen lme4"):
    assertMatches(
      "effect ~ 1 + duration + (1 | sites) + (1 | season)",
      boundary = false
    )

  test("station duration slopes have reduced covariance rank at the boundary against frozen lme4"):
    assertMatches(
      "effect ~ 1 + duration + (1 + duration | sites) + (1 + duration | season)",
      boundary = true
    )

  test("station fixed-only formula is refused"):
    val formula = "effect ~ 1 + duration * season * sites"
    val frame = Datasets.load("station_season_duration").frame
    Lmm.compile(formula, frame) match
      case Left(MixedModelError.NoRandomEffects) => ()
      case other                                 => fail(s"expected NoRandomEffects compile, got $other")
    Lmm.fit(formula, frame, FitOptions.ml) match
      case Left(MixedModelError.NoRandomEffects) => ()
      case other                                 => fail(s"expected NoRandomEffects fit, got $other")

  private def assertMatches(formula: String, boundary: Boolean): Unit =
    val row = FitCompare.row("station_season_duration", "ML", formula)
    val ref = FitCompare.claimed(row)
    val frame = Datasets.load("station_season_duration").frame
    val design = Lmm.compile(formula, frame).getOrElse(fail("compile"))
    val fit = Lmm.fit(formula, frame, FitOptions.ml).fold(err => fail(err.message), identity)
    val assessed = Pathology.assessFit(Pathology.certify(design), fit, design.parmap)
    if boundary then assertEquals(assessed.fitStatus, FitStatus.ConvergedReducedRank, clues(formula, fit.theta))
    else assertEquals(assessed.fitStatus, FitStatus.ConvergedInterior, clues(formula))
    ref.objective.foreach: obj =>
      assertEqualsDouble(fit.objective, obj, FitCompare.objectiveTol(row, ref))
    ref.sigma.foreach: sigma =>
      assertEqualsDouble(fit.sigma, sigma, FitCompare.sigmaTol(row, ref))
    ref.theta.foreach: theta =>
      val scalarTerms =
        fit.varcorr.components.forall(_.stdDev.length == 1) &&
          fit.varcorr.components.length == fit.theta.length
      FitCompare
        .alignedTheta(fit.theta, theta, scalarTerms)
        .zipWithIndex
        .foreach: (pair, i) =>
          val (got, expected) = pair
          assertEqualsDouble(got, expected, FitCompare.thetaTol(row, ref), clues(formula, i))
    ref.beta.foreach: beta =>
      val names = ref.coefNames.getOrElse(fail("missing coef_names"))
      FitCompare
        .alignedBeta(fit.feNames, fit.beta, names, beta)
        .zipWithIndex
        .foreach: (triple, i) =>
          val (name, got, expected) = triple
          assertEqualsDouble(got, expected, FitCompare.betaTol(row, ref), clues(name, i))
