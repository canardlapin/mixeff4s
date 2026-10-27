package mixeff4s.lmm

import mixeff4s.comparison.Datasets
import mixeff4s.fixtures.{Sleepstudy, Penicillin}

class SimulationSuite extends munit.FunSuite:
  private lazy val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml).toOption.get

  test("seeded response simulation is deterministic and index streams differ"):
    val first = Simulation.response(fit, 42L).toOption.get
    assertEquals(first, Simulation.response(fit, 42L).toOption.get)
    assertNotEquals(Simulation.seedFor(42L, 0L), Simulation.seedFor(42L, 1L))
    assertNotEquals(first, Simulation.response(fit, Simulation.seedFor(42L, 1L)).toOption.get)

  test("a copied fit without its sealed snapshot is refused"):
    assertEquals(
      Simulation.response(fit.copy(workflow = None), 1L).left.toOption.map(_.code),
      Some("fit_snapshot_mismatch")
    )
    assertEquals(
      Simulation.response(fit.copy(sigma = fit.sigma * 2), 1L).left.toOption.map(_.code),
      Some("fit_snapshot_mismatch")
    )

  test("uniform conversion stays strictly open at long extremes"):
    for bits <- Vector(Long.MinValue, Long.MaxValue, -1L, 0L) do
      assert(Simulation.uniform(bits) > 0.0)
      assert(Simulation.uniform(bits) < 1.0)

  /** Independent Gaussian moment identities. Covariances are centered at the known model mean, so the standard error is
    * sqrt((var_i*var_j+cov_ij^2)/B), including diagonal entries.
    */
  private def moments(
      name: String,
      fitted: LmmFit,
      mode: SimulationMode,
      rows: Vector[Int],
      covariance: Vector[Vector[Double]],
      count: Int
  ): Unit =
    val predictionMode =
      if mode == SimulationMode.DrawRandomEffects then PredictionMode.Population else PredictionMode.Conditional
    val mean = Prediction.fitted(fitted, predictionMode).toOption.get
    val sums = Array.fill(rows.length)(0.0)
    val products = Array.fill(rows.length, rows.length)(0.0)
    var replicate = 0
    while replicate < count do
      val draw = Simulation.response(fitted, Simulation.seedFor(728139L, replicate.toLong), mode).toOption.get
      var i = 0
      while i < rows.length do
        val x = draw(rows(i)) - mean(rows(i))
        sums(i) += x
        var j = 0
        while j < rows.length do
          products(i)(j) += x * (draw(rows(j)) - mean(rows(j)))
          j += 1
        i += 1
      replicate += 1
    for i <- rows.indices do
      assertEqualsDouble(sums(i) / count, 0.0, 6 * math.sqrt(covariance(i)(i) / count) + 1e-9)
      for j <- rows.indices do
        val expected = covariance(i)(j)
        val observed = products(i)(j) / count
        val se = math.sqrt((covariance(i)(i) * covariance(j)(j) + expected * expected) / count)
        assertEqualsDouble(observed, expected, 6 * se + 1e-9, clues(name, mode, i, j))
        if j <= i then
          println(s"SIM_MOMENT name=$name mode=$mode pair=$i,$j n=$count expected=$expected observed=$observed se=$se")

  test("random-intercept simulation has analytic within-group covariance and variance"):
    val tau2 = math.pow(fit.varcorr.components.head.stdDev.head, 2)
    val sigma2 = fit.sigma * fit.sigma
    val expected =
      Vector(Vector(tau2 + sigma2, tau2, 0.0), Vector(tau2, tau2 + sigma2, 0.0), Vector(0.0, 0.0, tau2 + sigma2))
    moments("intercept", fit, SimulationMode.DrawRandomEffects, Vector(0, 1, 10), expected, 4000)
    val conditional = Vector.tabulate(3, 3)((i, j) => if i == j then sigma2 else 0.0)
    moments("intercept", fit, SimulationMode.ConditionOnModes, Vector(0, 1, 10), conditional, 4000)

  test("random-slope covariance includes intercept/slope correlation and is shared only within a group"):
    val fitted = Lmm.fit("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame, FitOptions.ml).toOption.get
    val component = fitted.varcorr.components.head
    val a = component.stdDev(0) * component.stdDev(0)
    val b = component.stdDev(0) * component.stdDev(1) * component.correlations.head
    val c = component.stdDev(1) * component.stdDev(1)
    val noise = fitted.sigma * fitted.sigma
    val expected =
      Vector(Vector(a + noise, a + b, 0.0), Vector(a + b, a + 2 * b + c + noise, 0.0), Vector(0.0, 0.0, a + noise))
    moments("slope", fitted, SimulationMode.DrawRandomEffects, Vector(0, 1, 10), expected, 2000)

  test("crossed intercept simulation adds the covariance of each shared grouping factor"):
    val fitted = Lmm.fit("diameter ~ 1 + (1 | plate) + (1 | sample)", Penicillin.frame, FitOptions.ml).toOption.get
    val plate = math.pow(fitted.varcorr.components.find(_.group == "plate").get.stdDev.head, 2)
    val sample = math.pow(fitted.varcorr.components.find(_.group == "sample").get.stdDev.head, 2)
    val total = plate + sample + fitted.sigma * fitted.sigma
    val expected = Vector(
      Vector(total, plate, sample, 0.0),
      Vector(plate, total, 0.0, sample),
      Vector(sample, 0.0, total, plate),
      Vector(0.0, sample, plate, total)
    )
    moments("crossed", fitted, SimulationMode.DrawRandomEffects, Vector(0, 1, 6, 7), expected, 2000)

  test("zero fitted random variance removes its covariance without removing residual noise"):
    val frame = Datasets.load("dyestuff2").frame
    val zero = Lmm.fit("Yield ~ 1 + (1 | Batch)", frame, FitOptions.ml).toOption.get
    assert(zero.theta.forall(_ == 0.0))
    val a = Simulation.response(zero, 829L, SimulationMode.DrawRandomEffects).toOption.get
    val b = Simulation.response(zero, 829L, SimulationMode.ConditionOnModes).toOption.get
    assertEquals(a, b)
    val noise = zero.sigma * zero.sigma
    moments(
      "zero",
      zero,
      SimulationMode.DrawRandomEffects,
      Vector(0, 1),
      Vector(Vector(noise, 0.0), Vector(0.0, noise)),
      2000
    )
