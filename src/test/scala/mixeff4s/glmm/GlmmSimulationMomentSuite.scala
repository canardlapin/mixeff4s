package mixeff4s.glmm

import mixeff4s.fixtures.{GlmmFrames, GlmmAnalysisFrames}
import mixeff4s.formula.Formula
import mixeff4s.lmm.{Simulation, SimulationMode}
import mixeff4s.model.Family
import scala.concurrent.duration.*

/** Independent Gaussian-mixture moment witnesses for the public simulator. */
class GlmmSimulationMomentSuite extends munit.FunSuite:
  override val munitTimeout = 120.seconds
  private val draws = 4000
  private def right[A](value: Either[?, A]): A = value.fold(e => fail(e.toString), identity)
  private def mean(xs: Vector[Double]): Double = xs.sum / xs.length
  private def gaussianIntegral(f: Double => Double): Double =
    val intervals = 20000
    val step = 20.0 / intervals
    var sum = 0.0
    var i = 0
    while i <= intervals do
      val z = -10.0 + i * step
      val weight = if i == 0 || i == intervals then 1.0 else if i % 2 == 1 then 4.0 else 2.0
      sum += weight * f(z) * math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.Pi)
      i += 1
    sum * step / 3.0

  test("drawn correlated and crossed effects reproduce multivariate Poisson-lognormal moments"):
    val raw = GlmmFrames.synthetic_poisson_crossed
    val options = GlmmOptions.jointLaplace.copy(offset = Some(raw.numeric("exposure").get.map(math.log)))
    val fit = right(
      Glmm.fit(
        Formula.parse("count ~ 1 + x + (1 + x | site) + (1 | batch)").toOption.get,
        raw,
        Family.Poisson,
        options = options
      )
    )
    assertEquals(fit.randomEffects.map(_.groupingName), Vector("site", "batch"))
    val selected = Vector(3, 8, 28, 33)
    val values = Vector.tabulate(draws): index =>
      val y = right(GlmmSimulation.response(fit, Simulation.seedFor(1987112L, index.toLong)))
      selected.map(y(_))
    val x = raw.numeric("x").get
    val exposure = raw.numeric("exposure").get
    val sites = raw.factor("site").get.values
    val batches = raw.factor("batch").get.values
    val s00 = fit.theta(0) * fit.theta(0)
    val s01 = fit.theta(0) * fit.theta(1)
    val s11 = fit.theta(1) * fit.theta(1) + fit.theta(2) * fit.theta(2)
    val batchVariance = fit.theta(3) * fit.theta(3)
    def fixed(i: Int) = fit.beta(0) + fit.beta(1) * x(i) + math.log(exposure(i))
    def covariance(i: Int, j: Int): Double =
      (if sites(i) == sites(j) then s00 + (x(i) + x(j)) * s01 + x(i) * x(j) * s11 else 0.0) +
        (if batches(i) == batches(j) then batchVariance else 0.0)
    def stirling(power: Int): Vector[(Int, Double)] = power match
      case 0 => Vector(0 -> 1.0)
      case 1 => Vector(1 -> 1.0)
      case 2 => Vector(1 -> 1.0, 2 -> 1.0)
      case 3 => Vector(1 -> 1.0, 2 -> 3.0, 3 -> 1.0)
      case 4 => Vector(1 -> 1.0, 2 -> 7.0, 3 -> 6.0, 4 -> 1.0)
      case _ => fail("unsupported reference moment")
    def moment(i: Int, j: Int, r: Int, s: Int): Double =
      (for (k, a) <- stirling(r); (l, b) <- stirling(s)
      yield a * b * math.exp(
        k * fixed(i) + l * fixed(j) + 0.5 *
          (k * k * covariance(i, i) + l * l * covariance(j, j) + 2 * k * l * covariance(i, j))
      )).sum
    selected.zipWithIndex.foreach: (row, column) =>
      val expected = math.exp(fixed(row) + 0.5 * covariance(row, row))
      val variance = expected + expected * expected * math.expm1(covariance(row, row))
      val observed = mean(values.map(_(column)))
      assert(math.abs(observed - expected) <= 6.0 * math.sqrt(variance / draws), clues(row, observed, expected))
      val second = mean(values.map(v => math.pow(v(column) - expected, 2)))
      val fourth = moment(row, row, 4, 0) - 4 * expected * moment(row, row, 3, 0) +
        6 * expected * expected * moment(row, row, 2, 0) - 3 * math.pow(expected, 4)
      val secondSe = math.sqrt(math.max(0.0, fourth - variance * variance) / draws)
      assert(math.abs(second - variance) <= 6.0 * secondSe, clues(row, second, variance, secondSe))
    Vector((0, 1), (0, 2), (0, 3)).foreach: (a, b) =>
      val i = selected(a); val j = selected(b)
      val mi = moment(i, j, 1, 0); val mj = moment(i, j, 0, 1)
      val expected = mi * mj * math.expm1(covariance(i, j))
      val ca = Vector(mi * mi, -2.0 * mi, 1.0)
      val cb = Vector(mj * mj, -2.0 * mj, 1.0)
      val productSecond = (for r <- 0 to 2; s <- 0 to 2 yield ca(r) * cb(s) * moment(i, j, r, s)).sum
      val se = math.sqrt(math.max(0.0, productSecond - expected * expected) / draws)
      val observed = mean(values.map(v => (v(a) - mi) * (v(b) - mj)))
      assert(math.abs(observed - expected) <= 6.0 * se, clues(i, j, observed, expected, se))

  test("grouped binomial simulation agrees with independently integrated logistic-normal moments"):
    val options = GlmmOptions.jointAgq(5).copy(trials = Some(GlmmAnalysisFrames.trials))
    val fit = right(
      Glmm.fit(
        Formula.parse("y ~ 1 + x + (1 | g)").toOption.get,
        GlmmAnalysisFrames.binomial,
        Family.Binomial,
        options = options
      )
    )
    val x = Vector(-1.0, 0.0)
    def probability(i: Int, z: Double): Double =
      1.0 / (1.0 + math.exp(-fit.beta(0) - fit.beta(1) * x(i) - fit.theta.head * z))
    val mu = Vector.tabulate(2)(i => gaussianIntegral(z => probability(i, z)))
    val variance = Vector.tabulate(2): i =>
      gaussianIntegral: z =>
        val p = probability(i, z)
        p * (1.0 - p) / 20.0 + math.pow(p - mu(i), 2)
    val expectedCovariance = gaussianIntegral(z => probability(0, z) * probability(1, z)) - mu(0) * mu(1)
    val productSecond = gaussianIntegral: z =>
      val a = probability(0, z); val b = probability(1, z)
      (a * (1.0 - a) / 20.0 + math.pow(a - mu(0), 2)) *
        (b * (1.0 - b) / 20.0 + math.pow(b - mu(1), 2))
    val samples = Vector.tabulate(draws): index =>
      val response = right(GlmmSimulation.response(fit, Simulation.seedFor(443220L, index.toLong)))
      assert(response.forall(y => y >= 0.0 && y <= 1.0 && math.abs(y * 20.0 - math.rint(y * 20.0)) < 1e-10))
      response.take(2)
    mu.indices.foreach: i =>
      val observed = mean(samples.map(_(i)))
      assert(math.abs(observed - mu(i)) <= 6.0 * math.sqrt(variance(i) / draws), clues(i, observed, mu(i)))
    val observedCovariance = mean(samples.map(v => (v(0) - mu(0)) * (v(1) - mu(1))))
    val covarianceSe = math.sqrt((productSecond - expectedCovariance * expectedCovariance) / draws)
    assert(
      math.abs(observedCovariance - expectedCovariance) <= 6.0 * covarianceSe,
      clues(observedCovariance, expectedCovariance, covarianceSe)
    )

  test("conditional modes produce independent Poisson observations with fitted means and no RE resampling"):
    val raw = GlmmFrames.synthetic_poisson_offset
    val fit = right(
      Glmm.fit(
        Formula.parse("count ~ 1 + x + (1 | group)").toOption.get,
        raw,
        Family.Poisson,
        options = GlmmOptions.jointLaplace.copy(offset = Some(raw.numeric("exposure").get.map(math.log)))
      )
    )
    val samples = Vector.tabulate(draws): index =>
      right(GlmmSimulation.response(fit, Simulation.seedFor(19294L, index.toLong), SimulationMode.ConditionOnModes))
        .take(2)
    val mu = fit.conditionalMean.take(2)
    mu.indices.foreach: i =>
      val observed = mean(samples.map(_(i)))
      assert(math.abs(observed - mu(i)) <= 6.0 * math.sqrt(mu(i) / draws), clues(i, observed, mu(i)))
    val observed = mean(samples.map(v => (v(0) - mu(0)) * (v(1) - mu(1))))
    assert(math.abs(observed) <= 6.0 * math.sqrt(mu(0) * mu(1) / draws), clues(observed))

  test("Bernoulli mixtures reproduce independently integrated probability and within-group covariance"):
    val successes = Vector(0, 6, 1, 5, 2, 4)
    val ys = successes.flatMap(k => Vector.tabulate(6)(j => if j < k then 1.0 else 0.0))
    val frame = mixeff4s.data.ModelFrame
      .of(
        "y" -> mixeff4s.data.ModelFrame.numeric(ys),
        "g" -> mixeff4s.data.ModelFrame.factor(Vector.tabulate(ys.length)(i => s"g${i / 6}"))
      )
      .toOption
      .get
    val fit = right(
      Glmm.fit(
        Formula.parse("y ~ 1 + (1 | g)").toOption.get,
        frame,
        Family.Bernoulli,
        options = GlmmOptions.jointAgq(9)
      )
    )
    def p(z: Double) = 1.0 / (1.0 + math.exp(-fit.beta.head - fit.theta.head * z))
    val mu = gaussianIntegral(p)
    val covariance = gaussianIntegral(z => p(z) * p(z)) - mu * mu
    val productSecond = gaussianIntegral: z =>
      val conditionalCentralSecond = p(z) * (1.0 - p(z)) + math.pow(p(z) - mu, 2)
      conditionalCentralSecond * conditionalCentralSecond
    val samples = Vector.tabulate(draws): index =>
      right(GlmmSimulation.response(fit, Simulation.seedFor(78112L, index.toLong))).take(2)
    val observedMean = mean(samples.map(_.head))
    assert(math.abs(observedMean - mu) <= 6.0 * math.sqrt(mu * (1.0 - mu) / draws))
    val observedCovariance = mean(samples.map(v => (v(0) - mu) * (v(1) - mu)))
    val se = math.sqrt((productSecond - covariance * covariance) / draws)
    assert(math.abs(observedCovariance - covariance) <= 6.0 * se, clues(observedCovariance, covariance, se))
