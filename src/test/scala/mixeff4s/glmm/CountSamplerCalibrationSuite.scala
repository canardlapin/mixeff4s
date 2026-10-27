package mixeff4s.glmm

/** Distribution-law checks against R 4.5.1 frozen q/p{binom,pois} values. They exercise sampler output only, never the
  * sampler's acceptance formula.
  */
class CountSamplerCalibrationSuite extends munit.FunSuite:
  private val n = 20000
  private final case class Cdf(q: Vector[Double], p: Vector[Double])
  private def draws(f: Long => Double): Vector[Double] = Vector.tabulate(n)(i => f(900000L + i))
  private def check(xs: Vector[Double], mean: Double, variance: Double, fourth: Double, cdf: Cdf): Unit =
    val sampleMean = xs.sum / n
    assert(math.abs(sampleMean - mean) <= 7.0 * math.sqrt(variance / n), clues(sampleMean, mean))
    val central2 = xs.iterator.map(x => (x - mean) * (x - mean)).sum / n
    val se2 = math.sqrt((fourth - variance * variance) / n)
    assert(math.abs(central2 - variance) <= 7.0 * se2, clues(central2, variance, se2))
    cdf.q
      .zip(cdf.p)
      .foreach: (q, probability) =>
        val observed = xs.count(_ <= q).toDouble / n
        val tolerance = 7.0 * math.sqrt(probability * (1.0 - probability) / n) + 1.0 / n
        assert(math.abs(observed - probability) <= tolerance, clues(q, observed, probability, tolerance))

  private def binomial(nTrials: Int, p: Double, cdf: Cdf): Unit =
    val mean = nTrials.toDouble * p; val variance = mean * (1.0 - p)
    val fourth = 3.0 * variance * variance + variance * (1.0 - 6.0 * p * (1.0 - p))
    val fixed = Vector.tabulate(8)(i => CountSamplers.binomial(nTrials, p, 77123L + i).toOption.get)
    println(s"COUNT_VECTOR binomial n=$nTrials p=$p values=$fixed")
    check(draws(seed => CountSamplers.binomial(nTrials, p, seed).toOption.get), mean, variance, fourth, cdf)

  private def poisson(lambda: Double, cdf: Cdf): Unit =
    println(s"COUNT_VECTOR poisson lambda=$lambda values=${Vector
        .tabulate(8)(i => CountSamplers.poisson(lambda, 77123L + i).toOption.get)}")
    check(
      draws(seed => CountSamplers.poisson(lambda, seed).toOption.get),
      lambda,
      lambda,
      3.0 * lambda * lambda + lambda,
      cdf
    )

  test("independent frozen binomial calibration"):
    binomial(100, .5, Cdf(Vector(38, 50, 62), Vector(.010489367838925786, .539794618693589134, .993983512137318281)))
    binomial(
      10000,
      .7,
      Cdf(Vector(6893, 7000, 7106), Vector(.010231127066842514, .503772377545773864, .990112044710585559))
    )
    binomial(
      Int.MaxValue,
      5e-9,
      Cdf(Vector(4, 11, 19), Vector(.018015322109273146, .610580152605755266, .992714221405822883))
    )
    binomial(
      Int.MaxValue,
      .5,
      Cdf(Vector(1073687921, 1073741823, 1073795726), Vector(.010000673209828904, .5, .990000477063108786))
    )

  test("independent frozen Poisson calibration"):
    poisson(9.9, Cdf(Vector(3, 10, 18), Vector(.011119788299312656, .595548658425270716, .993494782037312896)))
    poisson(10, Cdf(Vector(3, 10, 18), Vector(.010336050675925716, .583039750192985373, .992813495396145718)))
    poisson(10.1, Cdf(Vector(4, 10, 18), Vector(.027416776038824555, .570530821115144082, .992075473029838339)))
    poisson(1000, Cdf(Vector(927, 1000, 1074), Vector(.010282157148410293, .508409367168506043, .990168693070978101)))
    poisson(
      1e6,
      Cdf(Vector(997674, 1000000, 1002327), Vector(.010002996249147574, .500265961486283661, .990011101271484173))
    )
