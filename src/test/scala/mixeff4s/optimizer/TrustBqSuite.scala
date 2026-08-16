package mixeff4s.optimizer

class TrustBqSuite extends munit.FunSuite:
  test("recovers a shifted quadratic"):
    val result = TrustBq
      .minimize(
        Vector(0.0, 0.0),
        Vector(Double.NegativeInfinity, Double.NegativeInfinity),
        Vector(Double.PositiveInfinity, Double.PositiveInfinity),
        TrustBqOptions(finalRadius = 1e-8, ftolAbs = 1e-12, ftolRel = 1e-12)
      )(x => Right((x(0) - 3.0) * (x(0) - 3.0) + (x(1) + 1.0) * (x(1) + 1.0)))
      .getOrElse(fail("minimize"))
    assertEqualsDouble(result.x(0), 3.0, 1e-4)
    assertEqualsDouble(result.x(1), -1.0, 1e-4)
    assertEqualsDouble(result.fmin, 0.0, 1e-8)
    assert(result.stopReason.isAcceptableConvergence)
