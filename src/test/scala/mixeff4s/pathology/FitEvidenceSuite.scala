package mixeff4s.pathology

import mixeff4s.compiler.Compiler
import mixeff4s.error.MixedModelError
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.lmm.{FitOptions, Lmm, LmmAssessment, LmmAssessmentOptions}
import mixeff4s.optimizer.TrustBqStopReason

class FitEvidenceSuite extends munit.FunSuite:
  private val formula = "reaction ~ 1 + days + (1 | subj)"
  private lazy val design = Lmm.compile(formula, Sleepstudy.frame).toOption.get
  private lazy val fit = Lmm.fit(formula, Sleepstudy.frame, FitOptions.ml).toOption.get

  test("an arbitrary raw theta does not certify optimization"):
    val cert = Pathology.certify(design)
    val assessed = Pathology.assessFit(cert, Vector(12345.0), design.parmap)
    assertEquals(assessed.fitStatus, FitStatus.NotAssessed)
    assertEquals(assessed.thetaGeometry, ThetaGeometry.Interior)
    assertEquals(Pathology.assessOutcome(cert, Right(Vector(12345.0)), design.parmap).fitStatus, FitStatus.NotAssessed)

  test("raw reassessment clears previous convergence evidence and stale geometry"):
    val fitted = Pathology.assessFit(Pathology.certify(design), fit, design.parmap)
    assertEquals(fitted.fitStatus, FitStatus.ConvergedInterior)
    val raw = Pathology.assessFit(fitted, Vector(0.0), design.parmap)
    assertEquals(raw.fitStatus, FitStatus.NotAssessed)
    assertEquals(raw.thetaGeometry, ThetaGeometry.Boundary)
    val invalid = Pathology.assessFit(fitted, Vector(Double.NaN), design.parmap)
    assertEquals(invalid.fitStatus, FitStatus.NotAssessed)
    assertEquals(invalid.thetaGeometry, ThetaGeometry.NotAssessed)

  test("a mismatched or exhausted optimizer result cannot certify a fit"):
    val cert = Pathology.certify(design)
    val changedTheta = fit.copy(theta = Vector(12345.0))
    assertEquals(Pathology.assessFit(cert, changedTheta, design.parmap).fitStatus, FitStatus.NotAssessed)
    val changedBeta = fit.copy(beta = fit.beta.updated(0, fit.beta.head + 1.0))
    val invalid = Pathology.assessFit(cert, changedBeta, design.parmap)
    assertEquals(invalid.fitStatus, FitStatus.NotAssessed)
    assert(invalid.notes.contains("fit_snapshot_mismatch"))
    val exhausted = fit.copy(optimizer = fit.optimizer.copy(stopReason = TrustBqStopReason.MaxEvaluations))
    assertEquals(Pathology.assessFit(cert, exhausted, design.parmap).fitStatus, FitStatus.NotOptimized)

  test("fitted covariance rank recognizes a nonzero off-diagonal reduced-rank Lambda"):
    val slopeDesign = Lmm.compile("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame).toOption.get
    val rank = LmmAssessment
      .covarianceRanks(
        slopeDesign.reterms,
        Vector(1.0, 1.0, 0.0),
        LmmAssessmentOptions(covarianceRankRelativeTolerance = 1e-10)
      )
      .head
    assertEquals(rank.rank, Some(1))
    assert(rank.eigenvalues.exists(_ > 1.0))
    assert(rank.eigenvalues.exists(math.abs(_) <= rank.tolerance))

  test("failed fit outcomes clear geometry and retain the stable reason code"):
    val cert = Pathology.assessFit(Pathology.certify(design), fit, design.parmap)
    val result = Pathology.assessFitOutcome(cert, Left(MixedModelError.PirlsNotConverged(1)), design.parmap)
    assertEquals(result.fitStatus, FitStatus.NotOptimized)
    assertEquals(result.thetaGeometry, ThetaGeometry.NotAssessed)
    assert(result.notes.contains("pirls_not_converged"))

  test("compiler serializes geometry separately and requires a matching fit"):
    val artifact = Compiler.artifact(design)
    val assessed = Compiler.assess(artifact, fit)
    assertEquals(assessed.schema.version, 2)
    assertEquals(assessed.pathology.fitStatus, FitStatus.ConvergedInterior)
    assert(assessed.toJson.contains("\"theta_geometry\": \"interior\""))
    val raw = Compiler.assess(assessed, fit.theta)
    assertEquals(raw.pathology.fitStatus, FitStatus.NotAssessed)
    val other = Compiler.compile("reaction ~ 1 + (1 | subj)", Sleepstudy.frame).toOption.get
    assertEquals(Compiler.assess(other, fit).pathology.fitStatus, FitStatus.NotAssessed)
