package mixeff4s.pathology

import mixeff4s.lmm.{FitOptions, Lmm}

class CorpusSuite extends munit.FunSuite:
  test("easy full-rank truth expects an interior fit"):
    val spec = GeneratorSpec.lmm(
      "easy",
      Vector.fill(30)(6),
      nFePredictors = 1,
      nReSlopes = 1,
      reCovTruth = Vector(Vector(4.0, 0.5), Vector(0.5, 1.0))
    )
    val cert = Pathology.certify(spec)
    assertEquals(cert.stratum, Stratum.Easy)
    assertEquals(cert.reRankTruth, 2)
    assertEquals(cert.reRankRequested, 2)
    assertEquals(cert.boundaries, Vector.empty)
    assertEquals(cert.structuralIssue, None)
    assertEquals(cert.expectedStatuses, Vector(FitStatus.ConvergedInterior))
    assertEquals(cert.n, 180)
    assertEquals(cert.minGroupSize, 6)

  test("zero slope variance is a contract boundary"):
    val spec = GeneratorSpec.lmm(
      "boundary_zero_slope",
      Vector.fill(30)(6),
      nFePredictors = 1,
      nReSlopes = 1,
      reCovTruth = Vector(Vector(4.0, 0.0), Vector(0.0, 0.0))
    )
    val cert = Pathology.certify(spec)
    assert(cert.boundaries.contains(BoundaryKind.ZeroVariance(1)), clues(cert.boundaries))
    assert(cert.expectedStatuses.contains(FitStatus.ConvergedBoundary), clues(cert.expectedStatuses))

  test("unit correlation is reduced-rank truth"):
    val spec = GeneratorSpec.lmm(
      "reduced_rank",
      Vector.fill(30)(6),
      nFePredictors = 1,
      nReSlopes = 1,
      reCovTruth = Vector(Vector(4.0, 2.0), Vector(2.0, 1.0))
    )
    val cert = Pathology.certify(spec)
    assertEquals(cert.stratum, Stratum.ReducedRank)
    assertEquals(cert.reRankTruth, 1)
    assertEquals(cert.reRankRequested, 2)
    assert(cert.boundaries.contains(BoundaryKind.UnitCorrelation(0, 1)), clues(cert.boundaries))
    assertEquals(
      cert.expectedStatuses,
      Vector(FitStatus.ConvergedReducedRank, FitStatus.ConvergedBoundary)
    )

  test("singletons with a slope are a refusal spec"):
    val spec = GeneratorSpec.lmm(
      "refusal_singletons",
      Vector.fill(6)(1),
      nFePredictors = 1,
      nReSlopes = 1,
      reCovTruth = Vector(Vector(4.0, 0.5), Vector(0.5, 1.0))
    )
    val cert = Pathology.certify(spec)
    assertEquals(cert.stratum, Stratum.Refusal)
    assertEquals(cert.structuralIssue.map(_.code), Some("singletons_with_slope"))
    assert(cert.expectedStatuses.contains(FitStatus.NotIdentifiable), clues(cert.expectedStatuses))

  test("a misshapen truth covariance is a malformed spec"):
    val spec = GeneratorSpec.lmm(
      "malformed_dim",
      Vector.fill(10)(6),
      nFePredictors = 1,
      nReSlopes = 1,
      reCovTruth = Vector(Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0))
    )
    val cert = Pathology.certify(spec)
    assertEquals(cert.structuralIssue.map(_.code), Some("malformed_spec"))
    assertEquals(cert.stratum, Stratum.Refusal)

  test("generate draws a deterministic Gaussian frame for an easy spec"):
    val spec = easySpec
    val first = Pathology.generate(spec).getOrElse(fail("generate"))
    val second = Pathology.generate(spec).getOrElse(fail("generate again"))
    assertEquals(first.formula, "y ~ 1 + x1 + (1 + x1 | g)")
    assertEquals(first.frame.nRows, 180)
    assertEquals(first.frame.numeric("y"), second.frame.numeric("y"))
    assertEquals(first.frame.numeric("x1"), second.frame.numeric("x1"))

  test("generate refuses a misshapen truth covariance"):
    val spec = GeneratorSpec.lmm(
      "malformed_dim",
      Vector.fill(10)(6),
      nFePredictors = 1,
      nReSlopes = 1,
      reCovTruth = Vector(Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0))
    )
    Pathology.generate(spec) match
      case Left(err) =>
        assertEquals(err.code, "invalid_argument")
      case Right(other) =>
        fail(s"expected generate refusal, got $other")

  test("an easy generated LMM fits to an interior status"):
    val spec = easySpec
    val cert = Pathology.certify(spec)
    val generated = Pathology.generate(spec).getOrElse(fail("generate"))
    val design = Lmm.compile(generated.formula, generated.frame).getOrElse(fail("compile"))
    val fit = Lmm.fit(generated.formula, generated.frame, FitOptions.ml).getOrElse(fail("fit"))
    val assessed = Pathology.assessFit(cert, fit.theta, design.parmap)
    assert(
      cert.expectedStatuses.contains(assessed.fitStatus),
      clues(assessed.fitStatus, fit.theta, cert.expectedStatuses)
    )
    assertEquals(assessed.fitStatus, FitStatus.ConvergedInterior)

  private def easySpec: GeneratorSpec =
    GeneratorSpec.lmm(
      "easy",
      Vector.fill(30)(6),
      nFePredictors = 1,
      nReSlopes = 1,
      reCovTruth = Vector(Vector(4.0, 0.5), Vector(0.5, 1.0)),
      seed = 42L,
      feTruth = Vector(1.0, 2.0)
    )
