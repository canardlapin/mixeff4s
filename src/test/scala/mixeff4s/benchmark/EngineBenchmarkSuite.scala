package mixeff4s.benchmark

import mixeff4s.comparison.Datasets
import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.fixtures.{Pastes, Penicillin, Sleepstudy}
import mixeff4s.lmm.{BatchStart, FitOptions, Lmm, LmmFit, LmmOptimizer, PlsWorkspace}

/** Opt in with MIXEFF4S_BENCHMARK=1. Raw samples are evidence, never a performance assertion. */
class EngineBenchmarkSuite extends munit.FunSuite:
  override def munitTests(): Seq[munit.Test] =
    if BenchmarkPlatform.enabled then super.munitTests() else Seq.empty

  private var sink = 0.0
  private final case class Workload(name: String, formula: String, frame: ModelFrame, response: String)
  private def synthetic(crossed: Boolean): Workload =
    val n = 1200
    val x = Vector.tabulate(n)(i => (i % 5 - 2).toDouble)
    val g = Vector.tabulate(n)(i => i / 10)
    val h = Vector.tabulate(n)(i => if crossed then i % 12 else i / 100)
    val y = Vector.tabulate(n)(i =>
      4.0 + 0.7 * x(i) + math.sin(g(i) * 1.7) + 0.8 * math.cos(h(i) * 0.9) + 0.2 * math.sin(i * 2.3)
    )
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(y),
        "x" -> ModelFrame.numeric(x),
        "g" -> ModelFrame.factor(g.map(_.toString)),
        "h" -> ModelFrame.factor(h.map(_.toString))
      )
      .fold(err => fail(err.message), identity)
    Workload(if crossed then "crossed1200" else "nested1200", "y ~ 1 + x + (1 | g) + (1 | h)", frame, "y")

  private def measure[A](workload: String, phase: String, repetitions: Int)(op: => A)(consume: A => Unit): Unit =
    var i = 0
    while i < repetitions do
      val heap0 = BenchmarkPlatform.heapBytes
      val alloc0 = BenchmarkPlatform.allocatedBytes
      val t0 = System.nanoTime()
      val result = op
      val elapsed = System.nanoTime() - t0
      val alloc1 = BenchmarkPlatform.allocatedBytes
      val heap1 = BenchmarkPlatform.heapBytes
      consume(result)
      val allocation = if alloc0 < 0 || alloc1 < 0 then "null" else (alloc1 - alloc0).toString
      println(
        s"BENCH {\"platform\":\"${BenchmarkPlatform.name}\",\"runtime\":\"${BenchmarkPlatform.runtime}\",\"workload\":\"$workload\",\"phase\":\"$phase\",\"sample\":$i,\"ns\":$elapsed,\"allocated_bytes\":$allocation,\"heap_delta_bytes\":${heap1 - heap0}}"
      )
      i += 1

  test("matched engine measurements"):
    val workloads = Vector(
      Workload("sleepstudy", "reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame, "reaction"),
      Workload("pastes", "strength ~ 1 + (1 | batch / cask)", Pastes.frame, "strength"),
      Workload("penicillin", "diameter ~ 1 + (1 | plate) + (1 | sample)", Penicillin.frame, "diameter"),
      Workload(
        "station_boundary",
        "effect ~ 1 + duration + (1 + duration | sites) + (1 + duration | season)",
        Datasets.load("station_season_duration").frame,
        "effect"
      ),
      synthetic(false),
      synthetic(true)
    )
    workloads.foreach: w =>
      def fit(data: ModelFrame): LmmFit =
        Lmm.fit(w.formula, data, FitOptions.ml).fold(err => fail(s"${w.name}: ${err.message}"), identity)
      var reference: Option[LmmFit] = None
      measure(w.name, "first_fit", 1)(fit(w.frame)): result =>
        reference = Some(result)
        sink += result.objective
        println(
          s"BENCH_FIT {\"platform\":\"${BenchmarkPlatform.name}\",\"workload\":\"${w.name}\",\"objective\":${result.objective},\"evaluations\":${result.optimizer.fevals},\"n\":${result.n},\"theta_count\":${result.theta.length}}"
        )
      val expected = reference.get
      // Untimed warmup: distinguish first use from steady repeated work.
      (0 until 3).foreach(_ => sink += fit(w.frame).objective)
      measure(w.name, "compile", 7)(Lmm.compile(w.formula, w.frame).fold(err => fail(err.message), identity))(d =>
        sink += d.nTheta
      )
      val design = Lmm.compile(w.formula, w.frame).toOption.get
      measure(w.name, "workspace", 7)(PlsWorkspace(design, false).fold(err => fail(err.message), identity))(ws =>
        sink += ws.nTheta
      )
      val workspace = PlsWorkspace(design, false).toOption.get
      measure(w.name, "objective", 20)(workspace.objectiveAt(expected.theta).fold(err => fail(err.message), identity)):
        value =>
          assertEqualsDouble(value, expected.objective, 1e-6)
          sink += value
      measure(w.name, "warm_fit", 7)(fit(w.frame)): result =>
        assertEqualsDouble(result.objective, expected.objective, 1e-6)
        sink += result.objective
      val plan = Lmm.prepare(w.formula, w.frame).fold(err => fail(err.message), identity)
      measure(w.name, "compiled_fit", 7)(plan.fit(FitOptions.ml).fold(err => fail(err.message), identity)): result =>
        assertEqualsDouble(result.objective, expected.objective, 1e-6)
        sink += result.objective
      val gradientOptions = FitOptions.ml.copy(optimizer = LmmOptimizer.AnalyticGradient)
      val gradientFirst = plan.fit(gradientOptions).fold(err => fail(err.message), identity)
      println(
        s"BENCH_FIT {\"platform\":\"${BenchmarkPlatform.name}\",\"workload\":\"${w.name}\",\"strategy\":\"analytic_gradient\",\"objective\":${gradientFirst.objective},\"evaluations\":${gradientFirst.optimizer.fevals},\"n\":${gradientFirst.n},\"theta_count\":${gradientFirst.theta.length}}"
      )
      assertEqualsDouble(gradientFirst.objective, expected.objective, 0.002)
      (0 until 3).foreach(_ => sink += plan.fit(gradientOptions).fold(err => fail(err.message), identity).objective)
      measure(w.name, "compiled_fit_gradient", 7)(plan.fit(gradientOptions).fold(err => fail(err.message), identity)):
        result =>
          assertEqualsDouble(result.objective, gradientFirst.objective, 1e-6)
          sink += result.objective
      measure(w.name, "objective_gradient", 20)(
        workspace.objectiveAndGradientAt(expected.theta).fold(err => fail(err.message), identity)
      ): result =>
        assertEqualsDouble(result.objective, expected.objective, 1e-6)
        sink += result.gradient.sum
      measure(w.name, "inspection", 20)((expected.coefTable, expected.varcorr)): result =>
        sink += result._1.estimates.sum + result._2.components.length
      val responses = (0 until 5).map: i =>
        val cols = w.frame.columns.map:
          case (name, Column.Numeric(values)) if name == w.response =>
            name -> ModelFrame.numeric(values.map(_ + 0.01 * i))
          case other => other
        ModelFrame.of(cols*).fold(err => fail(err.message), identity)
      measure(w.name, "response_batch_rebuild_5", 3)(responses.map(fit)): results =>
        results.foreach(r => assertEqualsDouble(r.objective, expected.objective, 1e-4))
        sink += results.map(_.objective).sum
      val dfSession = plan.session(FitOptions.ml).toOption.get
      val gradientSession = plan.session(gradientOptions).toOption.get
      val responseValues = responses.map(_.numeric(w.response).get).toVector
      measure(w.name, "response_batch_refit_independent_5", 3)(
        dfSession.refitAll(responseValues).fold(err => fail(err.message), identity)
      ): results =>
        results.foreach(r => assertEqualsDouble(r.objective, expected.objective, 1e-4))
        sink += results.map(_.objective).sum
      measure(w.name, "response_batch_refit_gradient_independent_5", 3)(
        gradientSession.refitAll(responseValues).fold(err => fail(err.message), identity)
      ): results =>
        results.foreach(r => assertEqualsDouble(r.objective, gradientFirst.objective, 1e-4))
        sink += results.map(_.objective).sum
      measure(w.name, "response_batch_refit_gradient_warm_5", 3)(
        gradientSession.refitAll(responseValues, BatchStart.WarmSequential).fold(err => fail(err.message), identity)
      ): results =>
        results.foreach(r => assertEqualsDouble(r.objective, gradientFirst.objective, 0.002))
        sink += results.map(_.objective).sum
    assert(sink.isFinite)
