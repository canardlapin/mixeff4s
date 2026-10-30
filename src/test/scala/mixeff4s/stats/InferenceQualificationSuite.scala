package mixeff4s.stats

import mixeff4s.benchmark.BenchmarkPlatform
import mixeff4s.error.FitResult
import mixeff4s.fixtures.{GlmmAnalysisFrames, Sleepstudy}
import mixeff4s.glmm.*
import mixeff4s.lmm.{Lmm, Simulation}
import mixeff4s.model.Family
import mixeff4s.examples.LikelihoodInference
import scala.concurrent.duration.*

/** Opt-in bounded calibration/coverage observations. Rates are evidence, not test pass criteria. */
class InferenceQualificationSuite extends munit.FunSuite:
  override val munitTimeout = 600.seconds
  override def munitTests(): Seq[munit.Test] = if BenchmarkPlatform.enabled then super.munitTests() else Seq.empty
  private def right[A](value: FitResult[A]): A = value.fold(e => fail(s"${e.code}: ${e.message}"), identity)
  private def nullable(value: Option[String]): String = value.fold("null")(s => "\"" + s + "\"")
  test("fitted-null simulation reports unconditional failures, conditional calibration and interval coverage"):
    val controls = GlmmOptions.jointLaplace.copy(trials = Some(GlmmAnalysisFrames.trials))
    val smallPlan = right(Glmm.prepare("y ~ 1 + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial))
    val largePlan = right(Glmm.prepare("y ~ 1 + x + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial))
    val source = right(smallPlan.fit(controls))
    println(
      s"STUDY_GENERATOR family=Binomial nAgq=1 beta=${source.beta} theta=${source.theta} omitted_slope=0 outer_seed=2026102901 inner_pairs=39 nominal_level=.9"
    )
    val smallSession = right(smallPlan.session(controls)); val largeSession = right(largePlan.session(controls))
    var requested = 0
    (0 until 16).foreach: index =>
      val seed = Simulation.seedFor(2026102901L, index.toLong)
      val y = right(GlmmSimulation.response(source, seed))
      println(
        s"INFERENCE_STUDY_RESPONSE {\"platform\":\"${BenchmarkPlatform.name}\",\"index\":$index,\"seed\":\"$seed\",\"y\":${y.mkString("[", ",", "]")}}"
      )
      val before = BenchmarkPlatform.allocatedBytes; val start = System.nanoTime()
      val a = smallSession.refit(y); val b = largeSession.refit(y)
      val covariance = b.flatMap(f =>
        GlmmInference.coefficients(f, GlmmInferenceOptions(Vector(BootstrapTarget.Coefficient("x")), level = .9))
      )
      val bootstrap = for
        small <- a
        large <- b
        result <- LikelihoodBootstrap.compare(
          small,
          large,
          LikelihoodBootstrapOptions(39, Simulation.seedFor(seed, 500L), monteCarloLevel = .95)
        )
      yield result
      val calibration = bootstrap.flatMap(_.calibration)
      val covered = covariance.toOption.map(_.intervals.head).map(r => r.lower <= 0 && r.upper >= 0)
      val bytes = BenchmarkPlatform.allocatedBytes
      val allocation = if before < 0 || bytes < 0 then "null" else (bytes - before).toString
      val code = nullable(calibration.left.toOption.map(_.code))
      val p = calibration.toOption.map(_.pvalue.toString).getOrElse("null")
      val coverage = covered.map(_.toString).getOrElse("null")
      println(
        s"INFERENCE_STUDY {\"platform\":\"${BenchmarkPlatform.name}\",\"runtime\":\"${BenchmarkPlatform.runtime}\",\"index\":$index,\"ns\":${System.nanoTime() - start},\"allocated_bytes\":$allocation,\"calibration_error\":$code,\"pvalue\":$p,\"coefficient_covered\":$coverage,\"coefficient_error\":${nullable(covariance.left.toOption.map(_.code))},\"pair_attempts\":${bootstrap.toOption.map(_.attempts.length).getOrElse(0)},\"pair_failures\":${bootstrap.toOption.map(_.failedOrUnqualifiedIndices.length).getOrElse(0)}}"
      )
      bootstrap.toOption.foreach: run =>
        run.attempts.foreach: pair =>
          println(
            s"INFERENCE_STUDY_PAIR platform=${BenchmarkPlatform.name} outer=$index index=${pair.index} seed=${pair.seed} statistic=${pair.statistic.map(_.likelihoodRatio)}"
          )
      requested += 1
    assertEquals(requested, 16)
  test("Gaussian coefficient profile coverage retains unavailable endpoints"):
    val source = right(Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, LikelihoodInference.tightML))
    println(
      s"STUDY_GENERATOR family=Gaussian beta=${source.beta} theta=${source.theta} sigma=${source.sigma} outer_seed=2026102902 nominal_level=.9"
    )
    val plan = source.validatedContext.toOption.get.compiled
    val session = right(plan.session(LikelihoodInference.tightML))
    (0 until 16).foreach: index =>
      val seed = Simulation.seedFor(2026102902L, index.toLong)
      val y = right(Simulation.response(source, seed))
      val start = System.nanoTime(); val before = BenchmarkPlatform.allocatedBytes
      val result = session
        .refit(y)
        .flatMap(f =>
          LmmProfile.confint(
            f,
            LmmProfileTarget.Coefficient("days"),
            LmmProfileOptions(
              level = .9,
              maxEvaluations = 2000,
              optimizerControls = LikelihoodInference.tightML.optimizerControls.get
            )
          )
        )
      val endpoints = result.toOption.filter(p => p.lower.value.nonEmpty && p.upper.value.nonEmpty)
      val covered = endpoints.map(p => p.lower.value.get <= source.beta(1) && p.upper.value.get >= source.beta(1))
      val error = result.left.toOption
        .map(_.code)
        .orElse(result.toOption.filter(_ => endpoints.isEmpty).map(_ => "profile_endpoint_unqualified"))
      val bytes = BenchmarkPlatform.allocatedBytes
      val allocation = if before < 0 || bytes < 0 then "null" else (bytes - before).toString
      println(
        s"PROFILE_STUDY {\"platform\":\"${BenchmarkPlatform.name}\",\"index\":$index,\"seed\":\"$seed\",\"ns\":${System.nanoTime() - start},\"allocated_bytes\":$allocation,\"covered\":${covered.map(_.toString).getOrElse("null")},\"error\":${nullable(error)}}"
      )
