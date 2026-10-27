package mixeff4s.glmm

class CountSamplerSuite extends munit.FunSuite:
  private def right[A](value: Either[?, A]): A = value.fold(error => fail(error.toString), identity)

  test("Bernoulli and binomial boundaries are exact and deterministic"):
    assertEquals(right(CountSamplers.bernoulli(0.0, 1L)), 0.0)
    assertEquals(right(CountSamplers.bernoulli(1.0, 1L)), 1.0)
    assertEquals(right(CountSamplers.binomial(1, 0.0, 2L)), 0.0)
    assertEquals(right(CountSamplers.binomial(1, 1.0, 2L)), 1.0)
    assertEquals(CountSamplers.binomial(10, 0.2, 812L), CountSamplers.binomial(10, 0.2, 812L))

  test("binomial samples stay integral and bounded across inversion and BTPE regimes"):
    Vector((Int.MaxValue, 1e-9), (Int.MaxValue, 0.5), (100, 0.2), (10000, 0.7)).zipWithIndex.foreach:
      (parameters, index) =>
        val (n, p) = parameters
        val draw = right(CountSamplers.binomial(n, p, 1203L + index))
        assert(draw >= 0.0 && draw <= n.toDouble)
        assertEquals(draw, math.rint(draw))
    // This seed enters BTPE's p2 branch with adjusted v > 1; it must reject
    // the candidate rather than evaluating a ratio loop from y = 0 to its mode.
    val p2Rejected = right(CountSamplers.binomial(Int.MaxValue, 0.5, 116L))
    assert(p2Rejected >= 0.0 && p2Rejected <= Int.MaxValue.toDouble)

  test("count moments and tail bins remain plausible without cross-row stream coupling"):
    val binomial = Vector.tabulate(12000)(i => right(CountSamplers.binomial(40, 0.35, i.toLong + 91L)))
    val poisson = Vector.tabulate(12000)(i => right(CountSamplers.poisson(37.0, i.toLong + 523L)))
    assertEqualsDouble(binomial.sum / binomial.length, 14.0, 0.12)
    assertEqualsDouble(poisson.sum / poisson.length, 37.0, 0.25)
    assert(binomial.count(_ <= 8) > 0)
    assert(poisson.count(_ >= 50) > 0)
    val second = right(CountSamplers.poisson(1000.0, 52L))
    assertEquals(second, right(CountSamplers.poisson(1000.0, 52L)))

  test("invalid and unsupported sampler inputs refuse with stable codes"):
    assertEquals(CountSamplers.binomial(-1, 0.5, 1L).left.toOption.get.code, "invalid_argument")
    assertEquals(CountSamplers.poisson(1000000.1, 1L).left.toOption.get.code, "simulation_lambda_limit")
    assertEquals(CountSamplers.poisson(Double.NaN, 1L).left.toOption.get.code, "invalid_argument")

  test("stable Poisson PMF agrees with independent frozen density values"):
    val lambdaTen = Vector(-6.087976994571854, -4.884004190245918, -3.967713458371761, -3.2745662778118163,
      -2.763740654045826, -2.4070657101070942, -2.1839221587928836, -2.0785616431350533, -2.0785616431350533,
      -2.1738718229393825, -2.356193379733334, -2.6185576442008234, -2.9550298808220425)
    (2 to 14)
      .zip(lambdaTen)
      .foreach: (k, expected) =>
        assertEqualsDouble(CountSamplers.poissonLogPmf(k.toDouble, 10.0), expected, 1e-12)
    Vector(
      (999000, -8.326360395178199),
      (999997, -7.826696896925569),
      (1000000, -7.826693898066878),
      (1000003, -7.826699895784259),
      (1001000, -8.327027063816786)
    ).foreach: (k, expected) =>
      assertEqualsDouble(CountSamplers.poissonLogPmf(k.toDouble, 1e6), expected, 1e-8)
