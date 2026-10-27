package mixeff4s.benchmark

import mixeff4s.data.ModelFrame
import mixeff4s.error.FitResult
import mixeff4s.fixtures.{GlmmFrames, GlmmAnalysisFrames}
import mixeff4s.formula.Formula
import mixeff4s.glmm.*
import mixeff4s.lmm.Simulation
import mixeff4s.model.Family
import scala.concurrent.duration.*

/** Opt-in matched rebuild/cold-refit measurements with fresh simulated responses. */
class GlmmWorkflowBenchmarkSuite extends munit.FunSuite:
  override val munitTimeout = 300.seconds
  override def munitTests(): Seq[munit.Test] = if BenchmarkPlatform.enabled then super.munitTests() else Seq.empty
  private def right[A](value: FitResult[A]): A = value.fold(e => fail(s"${e.code}: ${e.message}"), identity)
  private final case class Workload(
      name: String,
      formula: String,
      frame: ModelFrame,
      response: String,
      family: Family,
      options: GlmmOptions,
      nodes: Vector[Int],
      responses: Int
  )
  private def replace(frame: ModelFrame, name: String, values: Vector[Double]): ModelFrame =
    right(ModelFrame.of(frame.columns.map { case (key, value) =>
      key -> (if key == name then ModelFrame.numeric(values) else value)
    }*))

  private def measured(workload: String, nodes: Int, index: Int, phase: String, budget: Int)(
      operation: => FitResult[GlmmFit]
  ): FitResult[GlmmFit] =
    val allocation0 = BenchmarkPlatform.allocatedBytes
    val start = System.nanoTime()
    val result = operation
    val elapsed = System.nanoTime() - start
    val allocation1 = BenchmarkPlatform.allocatedBytes
    val allocated = if allocation0 < 0 || allocation1 < 0 then "null" else (allocation1 - allocation0).toString
    val assess0 = BenchmarkPlatform.allocatedBytes
    val assessStart = System.nanoTime()
    val assessment = result.toOption.map(GlmmAssessment.assess(_))
    val assessmentNs = System.nanoTime() - assessStart
    val assess1 = BenchmarkPlatform.allocatedBytes
    val assessAllocated = if assess0 < 0 || assess1 < 0 then "null" else (assess1 - assess0).toString
    val objective = result.toOption.map(_.deviance.toString).getOrElse("null")
    val error = result.left.toOption.map(e => "\"" + e.code + "\"").getOrElse("null")
    val evaluations = result.toOption.map(_.optimizer.fevals.toString).getOrElse("null")
    val reasons =
      assessment.map(_.refusalCodes).getOrElse(Vector.empty).map(x => "\"" + x + "\"").mkString("[", ",", "]")
    val risks = assessment.map(_.riskFlags).getOrElse(Vector.empty).map(x => "\"" + x + "\"").mkString("[", ",", "]")
    val gradient = assessment
      .flatMap(_.outer)
      .map(_.gradient)
      .collect { case JointGradientCheck.Available(_, norm, _) =>
        norm.toString
      }
      .getOrElse("null")
    println(
      s"GLMM_BENCH {\"platform\":\"${BenchmarkPlatform.name}\",\"runtime\":\"${BenchmarkPlatform.runtime}\",\"workload\":\"$workload\",\"nAgq\":$nodes,\"max_evaluations\":$budget,\"index\":$index,\"phase\":\"$phase\",\"ns\":$elapsed,\"allocated_bytes\":$allocated,\"fresh_assessment_ns\":$assessmentNs,\"fresh_assessment_allocated_bytes\":$assessAllocated,\"objective\":$objective,\"evaluations\":$evaluations,\"error\":$error,\"qualified\":${assessment.exists(_.eligibleForBootstrap)},\"projected_norm\":$gradient,\"reasons\":$reasons,\"risks\":$risks}"
    )
    result

  test("new responses preserve estimators and matched outcomes while reporting allocation and failures"):
    val count = BenchmarkPlatform.glmmResponseCount
    require(count > 0 && count <= 64, "benchmark response count must be between1 and64")
    val workloads = Vector(
      Workload(
        "binomial",
        "y ~ 1 + x + (1 | g)",
        GlmmAnalysisFrames.binomial,
        "y",
        Family.Binomial,
        GlmmOptions.jointAgq(5).copy(trials = Some(GlmmAnalysisFrames.trials)),
        Vector(1, 5),
        count
      ),
      Workload(
        "poisson_offset",
        "count ~ 1 + x + (1 | group)",
        GlmmFrames.synthetic_poisson_offset,
        "count",
        Family.Poisson,
        GlmmOptions
          .jointAgq(5)
          .copy(offset = Some(GlmmFrames.synthetic_poisson_offset.numeric("exposure").get.map(math.log))),
        Vector(1, 5),
        count
      ),
      Workload(
        "rare_bernoulli",
        "y ~ 1 + exposure + (1 | group)",
        GlmmAnalysisFrames.rare,
        "y",
        Family.Bernoulli,
        GlmmOptions.jointAgq(5),
        Vector(1, 5),
        count
      ),
      Workload(
        "crossed_poisson",
        "count ~ 1 + x + (1 + x | site) + (1 | batch)",
        GlmmFrames.synthetic_poisson_crossed,
        "count",
        Family.Poisson,
        GlmmOptions.jointLaplace
          .copy(offset = Some(GlmmFrames.synthetic_poisson_crossed.numeric("exposure").get.map(math.log))),
        Vector(1),
        math.min(8, count)
      )
    )
    workloads.foreach: workload =>
      val formula = right(Formula.parse(workload.formula).left.map(mixeff4s.error.MixedModelError.Formula.apply))
      val plan = right(Glmm.prepare(formula, workload.frame, workload.family))
      val original = right(plan.fit(workload.options))
      val responses = Vector.tabulate(workload.responses): index =>
        val seed = Simulation.seedFor(794623L, index.toLong)
        val response = right(GlmmSimulation.response(original, seed))
        println(
          s"GLMM_BENCH_RESPONSE {\"platform\":\"${BenchmarkPlatform.name}\",\"workload\":\"${workload.name}\",\"index\":$index,\"seed\":\"$seed\",\"y\":${response.mkString("[", ",", "]")}}"
        )
        response
      workload.nodes.foreach: nodes =>
        val options = workload.options.copy(
          nAgq = nodes,
          maxEvaluations = if workload.name == "crossed_poisson" then 2000 else workload.options.maxEvaluations
        )
        val session = right(plan.session(options))
        (0 until math.min(3, responses.length)).foreach: i =>
          val _ = session.refit(responses(i))
          val _ = Glmm.fit(
            formula,
            replace(workload.frame, workload.response, responses(i)),
            workload.family,
            options = options
          )
        responses.zipWithIndex.foreach: (response, index) =>
          val rebuilt = measured(workload.name, nodes, index, "rebuild", options.maxEvaluations)(
            Glmm.fit(formula, replace(workload.frame, workload.response, response), workload.family, options = options)
          )
          val refit = measured(workload.name, nodes, index, "refit", options.maxEvaluations)(session.refit(response))
          (rebuilt, refit) match
            case (Right(a), Right(b)) =>
              assertEqualsDouble(a.deviance, b.deviance, 1e-6)
              a.beta.zip(b.beta).foreach((x, y) => assertEqualsDouble(x, y, 1e-5))
              a.theta.zip(b.theta).foreach((x, y) => assertEqualsDouble(x, y, 1e-5))
              assertEquals(a.options, b.options)
              assertEquals(a.nAgq, b.nAgq)
              assertEquals(a.approximation, b.approximation)
            case (Left(a), Left(b)) => assertEquals(a.code, b.code)
            case other              => fail(s"matched outcome mismatch for ${workload.name}/$nodes/$index: $other")
