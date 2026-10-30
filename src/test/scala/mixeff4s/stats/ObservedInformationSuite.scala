package mixeff4s.stats

import mixeff4s.error.MixedModelError

class ObservedInformationSuite extends munit.FunSuite:
  test("full nuisance-adjusted inverse has correct factor two and differs from fixed-nuisance curvature"):
    // D=2b²+2bt+3t²: H=[[4,2],[2,6]], 2H^-1=[[.6,-.2],[-.2,.4]].
    val result = ObservedInformation
      .compute(Vector(0.0, 0.0), Vector.fill(2)(Double.NegativeInfinity), CurvatureOptions()) { x =>
        Right(2 * x(0) * x(0) + 2 * x(0) * x(1) + 3 * x(1) * x(1))
      }
      .toOption
      .get
    assertEqualsDouble(result.covariance(0)(0), .6, 1e-12)
    assertEqualsDouble(result.covariance(0)(1), -.2, 1e-12)
    assertEqualsDouble(result.covariance(1)(1), .4, 1e-12)
    assert(math.abs(result.covariance(0)(0) - .5) > .09)
    assertEquals(result.evaluations, 17)
  test("boundary, indefinite information and failed probes never obtain a fallback covariance"):
    val boundary =
      ObservedInformation.compute(Vector(0.0), Vector(0.0), CurvatureOptions())(x => Right(x.head * x.head))
    assertEquals(boundary.left.toOption.get.code, "curvature_boundary")
    val indefinite =
      ObservedInformation.compute(Vector(0.0, 0.0), Vector.fill(2)(Double.NegativeInfinity), CurvatureOptions())(x =>
        Right(x(0) * x(0) - x(1) * x(1))
      )
    assertEquals(indefinite.left.toOption.get.code, "curvature_not_positive_definite")
    val error = ObservedInformation.compute(Vector(1.0), Vector(0.0), CurvatureOptions())(_ =>
      Left(MixedModelError.Interrupted("probe"))
    )
    assertEquals(error.left.toOption.get.code, "interrupted")
  test("step-unstable curvature refuses and scalar chi-square threshold matches independent R"):
    val x = ObservedInformation.compute(
      Vector(0.0),
      Vector(Double.NegativeInfinity),
      CurvatureOptions(relativeStep = .1, stabilityTolerance = 1e-6)
    )(x => Right(x.head * x.head + 100 * math.pow(x.head, 4)))
    assertEquals(x.left.toOption.get.code, "curvature_step_unstable")
    assertEqualsDouble(InferenceQuantiles.chiSquare1(.95), 3.841458820694124, 1e-5)

  test("equilibrated checks detect unstable middle coordinate hidden by raw matrix scale"):
    val result = ObservedInformation.compute(
      Vector(0.0, 0.0, 0.0),
      Vector.fill(3)(Double.NegativeInfinity),
      CurvatureOptions(relativeStep = .1, stabilityTolerance = .02)
    ) { x =>
      Right(5e11 * x(0) * x(0) + .5 * x(1) * x(1) + 100 * math.pow(x(1), 4) + 5e-13 * x(2) * x(2))
    }
    assertEquals(result.left.toOption.get.code, "curvature_step_unstable")

  test("extreme admitted confidence levels expand the chi-square bracket"):
    val level = 1.0 - 1e-16
    val q = InferenceQuantiles.chiSquare1(level)
    assert(q > 64 && q < 80, clues(q))
    assertEqualsDouble(ChiSquare.sf(q, 1), 1 - level, 1e-25)

  test("diagonal equilibration avoids overflow and underflow in products of admissible diagonals"):
    Vector(1e-200, 1e200).foreach: scale =>
      val result = ObservedInformation
        .compute(Vector(0.0), Vector(Double.NegativeInfinity), CurvatureOptions())(x => Right(scale * x.head * x.head))
        .toOption
        .get
      val relative = math.abs(result.covariance.head.head * scale - 1)
      assert(relative < 1e-12, clues(scale, result))
