package mixeff4s.pathology

import mixeff4s.model.{Family, Link}

class SeparationSuite extends munit.FunSuite:
  test("complete FE separation is a strict sign split"):
    val x = design(Vector(-2.0, -1.5, -1.0, 1.0, 1.5, 2.0))
    val y = Vector(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
    val (kind, beta) = Separation.detectFe(x, y)
    assertEquals(kind, Some(FeSeparationKind.Complete))
    val slope = beta.getOrElse(fail("beta"))(1)
    assert(slope > Separation.MarginTol, clues(beta))

  test("quasi-complete FE separation has a tie at zero"):
    val x = design(Vector(-2.0, -1.0, 0.0, 0.0, 1.0, 2.0))
    val y = Vector(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
    val (kind, _) = Separation.detectFe(x, y)
    assertEquals(kind, Some(FeSeparationKind.QuasiComplete))

  test("overlapping classes are not FE-separated"):
    val x = design(Vector(-1.0, -0.5, 0.0, 0.5, 1.0, 1.5))
    val y = Vector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    val (kind, beta) = Separation.detectFe(x, y)
    assertEquals(kind, None)
    assertEquals(beta, None)

  test("conditional scan finds all-zero and all-one groups"):
    val y = Vector(0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0)
    val groups = Vector(0, 0, 0, 1, 1, 2, 2, 2)
    assertEquals(Separation.detectConditional(y, groups), Vector(0, 2))

  test("every singleton group is conditionally separated"):
    val y = Vector(0.0, 1.0, 0.0)
    val groups = Vector(0, 1, 2)
    assertEquals(Separation.detectConditional(y, groups), Vector(0, 1, 2))

  test("empty inputs are not conditionally separated"):
    assertEquals(Separation.detectConditional(Vector.empty, Vector.empty), Vector.empty)

  test("a Gaussian spec has an empty separation report"):
    val spec = GeneratorSpec.lmm(
      "easy",
      Vector.fill(10)(4),
      nFePredictors = 1,
      nReSlopes = 0,
      reCovTruth = Vector(Vector(1.0))
    )
    assertEquals(Pathology.detectSeparation(spec), SeparationReport.empty)

  test("Bernoulli generate draws a binary response"):
    val spec = GeneratorSpec.extremePrevalence(
      GeneratorSpec.lmm(
        "binary",
        Vector(6, 6),
        nFePredictors = 1,
        nReSlopes = 0,
        reCovTruth = Vector(Vector(1.0)),
        seed = 3L,
        feTruth = Vector(0.0, 0.5)
      ),
      interceptShift = 0.0
    )
    val generated = Pathology.generate(spec).getOrElse(fail("generate"))
    val y = generated.frame.numeric("y").getOrElse(fail("y"))
    assert(y.forall(v => v == 0.0 || v == 1.0), clues(y))
    assertEquals(y.length, 12)

  test("Bernoulli generate refuses an unsupported link"):
    val spec = GeneratorSpec
      .lmm(
        "probit",
        Vector(4, 4),
        nFePredictors = 1,
        nReSlopes = 0,
        reCovTruth = Vector(Vector(1.0))
      )
      .copy(family = Family.Bernoulli, link = Link.Probit)
    Pathology.generate(spec) match
      case Left(err) =>
        assertEquals(err.code, "unsupported_family_link")
      case Right(other) =>
        fail(s"expected family/link refusal, got $other")

  private def design(xs: Vector[Double]): Vector[Vector[Double]] =
    xs.map(x => Vector(1.0, x))
