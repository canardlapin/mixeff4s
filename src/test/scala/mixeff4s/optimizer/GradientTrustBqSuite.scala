package mixeff4s.optimizer

class GradientTrustBqSuite extends munit.FunSuite:
  private val unbounded = (Vector.fill(2)(Double.NegativeInfinity), Vector.fill(2)(Double.PositiveInfinity))

  test("gradient oracle solves bounded quadratic and records first-order evidence"):
    val (lower, upper) = unbounded
    val result = TrustBq.minimizeWithGradient(
      Vector(4.0, -5.0),
      lower,
      upper,
      TrustBqOptions(initialRadius = 1.0, finalRadius = 1e-8, maxEvaluations = 100, gradientTolerance = 1e-9)
    ): x =>
      val f = (x(0) - 1.25) * (x(0) - 1.25) + 2.0 * (x(1) + 2.0) * (x(1) + 2.0)
      Right(ObjectiveGradient(f, Vector(2.0 * (x(0) - 1.25), 4.0 * (x(1) + 2.0))))
    val fitted = result.getOrElse(fail("gradient minimization"))
    assert(fitted.stopReason.isAcceptableConvergence)
    assertEqualsDouble(fitted.x(0), 1.25, 1e-5)
    assertEqualsDouble(fitted.x(1), -2.0, 1e-5)
    assert(fitted.projectedGradientNorm.exists(_ <= 1e-8))
    assert(fitted.gradientEvaluations == fitted.fevals)

  test("multidimensional box-boundary quadratic respects the active bound"):
    val result = TrustBq.minimizeWithGradient(
      Vector(3.0, 2.0),
      Vector(Double.NegativeInfinity, 0.0),
      Vector(Double.PositiveInfinity, Double.PositiveInfinity),
      TrustBqOptions(initialRadius = 1.0, maxEvaluations = 100, gradientTolerance = 1e-9)
    ): x =>
      Right(
        ObjectiveGradient(
          (x(0) - 1.25) * (x(0) - 1.25) + 2.0 * (x(1) + 2.0) * (x(1) + 2.0),
          Vector(2.0 * (x(0) - 1.25), 4.0 * (x(1) + 2.0))
        )
      )
    val fitted = result.getOrElse(fail("bounded quadratic"))
    assertEqualsDouble(fitted.x(0), 1.25, 1e-5)
    assertEquals(fitted.x(1), 0.0)

  test("zero derivative at a bounded maximum is probed before first-order stop"):
    val result = TrustBq.minimizeWithGradient(
      Vector(0.0),
      Vector(0.0),
      Vector(1.0),
      TrustBqOptions(initialRadius = 0.25, finalRadius = 1e-8, maxEvaluations = 100, gradientTolerance = 1e-9)
    ): x =>
      val t = x.head
      Right(ObjectiveGradient(-t * t + t * t * t * t, Vector(-2.0 * t + 4.0 * t * t * t)))
    val fitted = result.getOrElse(fail("gradient minimization"))
    assert(fitted.x.head > 0.5, fitted.toString)
    assert(fitted.fmin < -0.2, fitted.toString)

  test("budget includes initial gradient oracle call and seed probes"):
    val result =
      TrustBq.minimizeWithGradient(Vector(0.0), Vector(-1.0), Vector(1.0), TrustBqOptions(maxEvaluations = 1)) { x =>
        Right(ObjectiveGradient(x.head * x.head, Vector(2.0 * x.head)))
      }
    val fitted = result.getOrElse(fail("budget result"))
    assertEquals(fitted.stopReason, TrustBqStopReason.MaxEvaluations)
    assertEquals(fitted.fevals, 1)
    assert(fitted.requireConvergence.isLeft)

  test("shifted Newton direction minimizes the original quadratic along its feasible ray"):
    val result = TrustBq.minimizeWithGradient(
      Vector(0.0),
      Vector(Double.NegativeInfinity),
      Vector(0.4),
      TrustBqOptions(initialRadius = 1.0, finalRadius = 1e-8, maxEvaluations = 20, gradientTolerance = 0.0)
    ): x =>
      val t = x.head
      Right(ObjectiveGradient(-t - 0.5 * t * t, Vector(-1.0 - t)))
    val fitted = result.getOrElse(fail("indefinite ray"))
    assertEqualsDouble(fitted.x.head, 0.4, 1e-12)

  test("disabled gradient tolerance cannot report a gradient stop"):
    val result = TrustBq.minimizeWithGradient(
      Vector(0.0),
      Vector(-1.0),
      Vector(1.0),
      TrustBqOptions(finalRadius = 1e-5, maxEvaluations = 100, gradientTolerance = 0.0)
    ): x =>
      Right(ObjectiveGradient(x.head * x.head, Vector(2.0 * x.head)))
    val fitted = result.getOrElse(fail("disabled gradient stop"))
    assertNotEquals(fitted.stopReason, TrustBqStopReason.GradientBelowTolerance)

  test("first-order stop validates a scalar feasible escape at a zero-gradient maximum"):
    def polynomial(x: Vector[Double]) =
      val t = x.head
      ObjectiveGradient(-t * t + 2.5 * t * t * t - t * t * t * t, Vector(-2.0 * t + 7.5 * t * t - 4.0 * t * t * t))
    val result = TrustBq.minimizeWithGradient(
      Vector(1.0),
      Vector(0.0),
      Vector(1.0),
      TrustBqOptions(initialRadius = 1.0, maxEvaluations = 100, gradientTolerance = 1e-8)
    )(x => Right(polynomial(x)))
    val fitted = result.getOrElse(fail("scalar escape"))
    assert(fitted.fmin < 0.0, fitted.toString)

  test("first-order stop validates a feasible mixed negative-curvature escape"):
    val result = TrustBq.minimizeWithGradient(
      Vector(0.0, 0.0),
      Vector(0.0, 0.0),
      Vector(1.0, 1.0),
      TrustBqOptions(initialRadius = 1.0, maxEvaluations = 100, gradientTolerance = 1e-8)
    ): x =>
      Right(
        ObjectiveGradient(
          x(0) * x(0) + x(1) * x(1) - 4.0 * x(0) * x(1),
          Vector(2.0 * x(0) - 4.0 * x(1), 2.0 * x(1) - 4.0 * x(0))
        )
      )
    val fitted = result.getOrElse(fail("mixed escape"))
    assert(fitted.fmin < 0.0, fitted.toString)

  test("returned projected-gradient evidence belongs to returned best point"):
    def polynomial(x: Vector[Double]) =
      val t = x.head
      ObjectiveGradient(-t * t + 2.5 * t * t * t - t * t * t * t, Vector(-2.0 * t + 7.5 * t * t - 4.0 * t * t * t))
    val result = TrustBq.minimizeWithGradient(
      Vector(1.0),
      Vector(0.0),
      Vector(1.0),
      TrustBqOptions(
        initialRadius = 1.0,
        maxEvaluations = 3,
        etaAccept = 0.9,
        etaExpand = 0.95,
        gradientTolerance = 1e-8
      )
    )(x => Right(polynomial(x)))
    val fitted = result.fold(error => fail(error.toString), identity)
    assertEquals(fitted.stopReason, TrustBqStopReason.MaxEvaluations)
    assertEqualsDouble(fitted.x.head, 0.0, 1e-12)
    assert(fitted.projectedGradientNorm.contains(0.0), fitted.toString)

  test("negative curvature hidden from all axial and pair directions is escaped"):
    val result = TrustBq.minimizeWithGradient(
      Vector(0.0, 0.0, 0.0),
      Vector.fill(3)(0.0),
      Vector.fill(3)(1.0),
      TrustBqOptions(initialRadius = 1.0, maxEvaluations = 100, gradientTolerance = 1e-8)
    ): x =>
      val f = x.map(v => v * v).sum - 1.2 * (x(0) * x(1) + x(0) * x(2) + x(1) * x(2))
      Right(
        ObjectiveGradient(
          f,
          Vector(
            2.0 * x(0) - 1.2 * x(1) - 1.2 * x(2),
            2.0 * x(1) - 1.2 * x(0) - 1.2 * x(2),
            2.0 * x(2) - 1.2 * x(0) - 1.2 * x(1)
          )
        )
      )
    val fitted = result.getOrElse(fail("three way saddle"))
    assert(fitted.fmin < 0.0, fitted.toString)

  test("Cauchy fallback projects outward components onto the feasible tangent cone"):
    val result = TrustBq.minimizeWithGradient(
      Vector(0.0, 0.5),
      Vector(0.0, 0.0),
      Vector(0.001, 1.0),
      TrustBqOptions(initialRadius = 1.0, maxEvaluations = 100, gradientTolerance = 0.0)
    ): x =>
      Right(ObjectiveGradient(1000.0 * x(0) - 500000.0 * x(0) * x(0) + x(1), Vector(1000.0 - 1000000.0 * x(0), 1.0)))
    val fitted = result.getOrElse(fail("tangent Cauchy"))
    assert(fitted.x(1) < 0.5, fitted.toString)

  test("invalid trust-region constants are refused before oracle evaluation"):
    val result =
      TrustBq.minimizeWithGradient(Vector(0.0), Vector(-1.0), Vector(1.0), TrustBqOptions(shrinkFactor = 1.0))(_ =>
        fail("oracle should not run")
      )
    assert(result.isLeft)

  test("invalid gradient values and dimensions are refused"):
    val nonFinite = TrustBq.minimizeWithGradient(Vector(0.0), Vector(-1.0), Vector(1.0), TrustBqOptions())(_ =>
      Right(ObjectiveGradient(0.0, Vector(Double.NaN)))
    )
    assert(nonFinite.isLeft)
    val wrongLength = TrustBq.minimizeWithGradient(Vector(0.0), Vector(-1.0), Vector(1.0), TrustBqOptions())(_ =>
      Right(ObjectiveGradient(0.0, Vector.empty))
    )
    assert(wrongLength.isLeft)

  test("a fresh escape probe remains the returned best point when the next probe exhausts budget"):
    val result = TrustBq
      .minimizeWithGradient(
        Vector(0.0, 0.0),
        Vector(0.0, 0.0),
        Vector(1.0, 1.0),
        TrustBqOptions(maxEvaluations = 4, gradientTolerance = 2.0)
      )(x => Right(ObjectiveGradient(-x.sum, Vector(-1.0, -1.0))))
      .toOption
      .get
    assertEquals(result.fevals, 4)
    assertEquals(result.stopReason, TrustBqStopReason.MaxEvaluations)
    assertEqualsDouble(result.fmin, -0.0002, 1e-12)
    assertEqualsDouble(result.x.sum, 0.0002, 1e-12)

  test("projected Cauchy ray uses its own directional derivative after blocking a coordinate"):
    val result = TrustBq.minimizeWithGradient(
      Vector(0.0, 0.75),
      Vector(0.0, 0.0),
      Vector(0.001, 1.0),
      TrustBqOptions(maxEvaluations = 4, initialRadius = 0.75, gradientTolerance = 0.0)
    ): x =>
      Right(
        ObjectiveGradient(
          1000 * x(0) - 500000 * x(0) * x(0) + (x(1) - 0.5) * (x(1) - 0.5),
          Vector(1000 - 1000000 * x(0), 2 * (x(1) - 0.5))
        )
      )
    val best = result.toOption.get
    assertEqualsDouble(best.x(1), 0.5, 1e-8)
    assertEqualsDouble(best.fmin, 0.0, 1e-12)
