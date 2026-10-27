package mixeff4s.lmm

import mixeff4s.fixtures.Sleepstudy

class LmmAssessmentSuite extends munit.FunSuite:
  test("a zero-theta scalar variance maximum fails the feasible boundary probe"):
    val result = LmmAssessment.assessBoundary(
      Vector(0.0),
      Vector((0, 0, 0)),
      baseline = 0.0,
      tolerance = 1e-12,
      LmmAssessmentOptions(boundaryProbe = 0.1)
    )(theta => Right(-theta.head * theta.head))
    result match
      case BoundaryAdequacy.Failed(0, step, objective) =>
        assertEqualsDouble(step, 0.1, 1e-15)
        assertEqualsDouble(objective, -0.01, 1e-15)
      case other => fail(s"expected a feasible boundary descent, got $other")

  test("relative covariance rank retains nonzero off-diagonal degeneracy"):
    val design = Lmm.compile("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame).toOption.get
    val ranks = LmmAssessment.covarianceRanks(
      design.reterms,
      Vector(1.0, 1.0, 0.0),
      LmmAssessmentOptions(covarianceRankRelativeTolerance = 1e-10)
    )
    assertEquals(ranks.length, 1)
    assertEquals(ranks.head.rank, Some(1))
    assert(ranks.head.eigenvalues.exists(_ > 1.0))
    assert(ranks.head.eigenvalues.exists(math.abs(_) <= ranks.head.tolerance))

  test("tampered fit evidence is unavailable rather than trusting stored optimizer norms"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml).toOption.get
    val tampered = fit.copy(optimizer = fit.optimizer.copy(projectedGradientNorm = Some(0.0)))
    val assessed = LmmAssessment.assess(tampered)
    assertEquals(assessed.gradient, GradientCheck.Unavailable("fit_snapshot_mismatch"))
    assertEquals(assessed.bootstrapIneligibilityReason, Some("fit_snapshot_mismatch"))

  test("weak full-rank covariance uses its own relative scale"):
    val design = Lmm.compile("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame).toOption.get
    val ranks = LmmAssessment.covarianceRanks(design.reterms, Vector(1e-8, 0.0, 1e-8), LmmAssessmentOptions())
    assertEquals(ranks.head.rank, Some(2))
    assert(ranks.head.tolerance < 1e-16)

  test("gradient scratch refusal preserves independent covariance rank evidence"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml).toOption.get
    val assessed = LmmAssessment.assess(fit, LmmAssessmentOptions(gradientOptions = ProfiledGradientOptions(1)))
    assert(assessed.gradient.isInstanceOf[GradientCheck.Unavailable])
    assertEquals(assessed.covarianceRanks.map(_.rank), Vector(Some(1)))
    assert(!assessed.eligibleForBootstrap)
