package mixeff4s.benchmark

import mixeff4s.data.ModelFrame
import mixeff4s.fixtures.{Sleepstudy, Pastes, Penicillin}
import mixeff4s.lmm.*
import mixeff4s.optimizer.TrustBqOptions
import mixeff4s.error.FitResult

/** Opt-in full-rebuild/response-refit experiment on fresh simulated responses. */
class WorkflowBenchmarkSuite extends munit.FunSuite:
  override def munitTests(): Seq[munit.Test] = if BenchmarkPlatform.enabled then super.munitTests() else Seq.empty
  private val controls = TrustBqOptions(
    finalRadius = 1e-10,
    ftolAbs = 0,
    ftolRel = 0,
    stallFtolAbs = 0,
    stallFtolRel = 0,
    stallIterations = 20,
    gradientTolerance = 1e-9
  )
  private val accurate =
    FitOptions.ml.copy(optimizer = LmmOptimizer.AnalyticGradient, optimizerControls = Some(controls))
  private def right[A](value: FitResult[A]): A = value.fold(e => fail(e.message), identity)
  private case class Workload(name: String, formula: String, frame: ModelFrame, response: String)
  private def replace(data: ModelFrame, name: String, values: Vector[Double]): ModelFrame =
    right(ModelFrame.of(data.columns.map { case (key, column) =>
      key -> (if key == name then ModelFrame.numeric(values) else column)
    }*))

  private def measured(workload: String, strategy: LmmOptimizer, index: Int, phase: String)(
      operation: => FitResult[LmmFit]
  ): FitResult[LmmFit] =
    val allocation0 = BenchmarkPlatform.allocatedBytes
    val start = System.nanoTime()
    val result = operation
    val elapsed = System.nanoTime() - start
    val allocation1 = BenchmarkPlatform.allocatedBytes
    val allocated = if allocation0 < 0 || allocation1 < 0 then "null" else (allocation1 - allocation0).toString
    val assess0 = BenchmarkPlatform.allocatedBytes
    val assessStart = System.nanoTime()
    val assessed = result.toOption.map(LmmAssessment.assess(_))
    val assessElapsed = System.nanoTime() - assessStart
    val assess1 = BenchmarkPlatform.allocatedBytes
    val assessAllocated = if assess0 < 0 || assess1 < 0 then "null" else (assess1 - assess0).toString
    val objective = result.toOption.map(_.objective.toString).getOrElse("null")
    val evaluations = result.toOption.map(_.optimizer.fevals.toString).getOrElse("null")
    val error = result.left.toOption.map(e => s"\"${e.code}\"").getOrElse("null")
    val projected = assessed
      .map(_.gradient)
      .collect { case GradientCheck.Pass(n, _) => n; case GradientCheck.Failed(n, _) => n }
      .map(_.toString)
      .getOrElse("null")
    val reasons =
      assessed.map(_.reasonCodes).getOrElse(Vector.empty).map(code => "\"" + code + "\"").mkString("[", ",", "]")
    println(
      s"FLOW_BENCH {\"platform\":\"${BenchmarkPlatform.name}\",\"runtime\":\"${BenchmarkPlatform.runtime}\",\"workload\":\"$workload\",\"strategy\":\"$strategy\",\"index\":$index,\"phase\":\"$phase\",\"ns\":$elapsed,\"allocated_bytes\":$allocated,\"assessment_ns\":$assessElapsed,\"assessment_allocated_bytes\":$assessAllocated,\"objective\":$objective,\"evaluations\":$evaluations,\"error\":$error,\"qualified\":${assessed.exists(_.eligibleForBootstrap)},\"projected_norm\":$projected,\"reasons\":$reasons}"
    )
    result

  test("realistic simulated responses match full rebuilds and report all outcomes"):
    val workloads = Vector(
      Workload("sleepstudy_intercept", "reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, "reaction"),
      Workload("sleepstudy_slopes", "reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame, "reaction"),
      Workload("pastes", "strength ~ 1 + (1 | batch / cask)", Pastes.frame, "strength"),
      Workload("penicillin", "diameter ~ 1 + (1 | plate) + (1 | sample)", Penicillin.frame, "diameter")
    )
    val original = right(
      Lmm.fit(workloads.head.formula, Sleepstudy.frame, FitOptions.ml.copy(optimizer = LmmOptimizer.AnalyticGradient))
    )
    for index <- 0 until 39 do
      val seed = Simulation.seedFor(20261026L, index.toLong)
      val y = right(Simulation.response(original, seed))
      println(
        s"FLOW_RESPONSE {\"platform\":\"${BenchmarkPlatform.name}\",\"index\":$index,\"seed\":\"$seed\",\"y\":${y.mkString("[", ",", "]")}}"
      )
    workloads.foreach: workload =>
      val plan = right(Lmm.prepare(workload.formula, workload.frame))
      val canonical = right(plan.fit(accurate))
      val responses =
        Vector.tabulate(15)(index => right(Simulation.response(canonical, Simulation.seedFor(834927L, index.toLong))))
      for strategy <- LmmOptimizer.values do
        val options = if strategy == LmmOptimizer.DerivativeFree then FitOptions.ml else accurate
        val session = right(plan.session(options))
        // Unmeasured work warms each fitting path before its samples.
        (0 until 3).foreach: i =>
          session.refit(responses(i), RefitStart.Cold)
          Lmm.fit(workload.formula, replace(workload.frame, workload.response, responses(i)), options)
        responses.zipWithIndex.foreach: (response, index) =>
          val rebuild = measured(workload.name, strategy, index, "rebuild")(
            Lmm.fit(workload.formula, replace(workload.frame, workload.response, response), options)
          )
          val refit = measured(workload.name, strategy, index, "refit")(session.refit(response, RefitStart.Cold))
          (rebuild, refit) match
            case (Right(a), Right(b)) =>
              assertEqualsDouble(a.objective, b.objective, 1e-6)
              a.beta.zip(b.beta).foreach((x, y) => assertEqualsDouble(x, y, 1e-7))
              assertEquals(a.options, b.options)
              assertEquals(a.optimizer.stopReason, b.optimizer.stopReason)
            case (Left(a), Left(b)) => assertEquals(a.code, b.code)
            case other              => fail(s"$workload/$strategy full-rebuild mismatch: $other")
