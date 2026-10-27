package mixeff4s.glmm

import mixeff4s.examples.GlmmAnalysis
import mixeff4s.comparison.JsonValue
import mixeff4s.lmm.Simulation
import scala.concurrent.duration.*

class GlmmAnalysisSuite extends munit.FunSuite:
  override val munitTimeout = 120.seconds
  private lazy val output = GlmmAnalysis.run().fold(e => fail(s"${e.code}: ${e.message}"), identity)

  test("documented compile fit assess predict simulate bootstrap workflow retains 39 qualified AGQ attempts"):
    val result = output.bootstrap
    assert(output.assessment.eligibleForBootstrap, clues(output.assessment))
    assertEquals(output.fit.approximation, Approximation.JointAgq)
    assertEquals(output.fit.nAgq, 5)
    assertEquals(output.prediction.length, 3)
    assertEquals(output.residuals.length, output.fit.n)
    assertEquals(output.simulated.length, output.fit.n)
    assertEquals(result.attempts.length, 39)
    assertEquals(
      result.qualifiedCount,
      39,
      clues(result.attempts.filterNot(_.qualified).map(a => (a.index, a.refit.left.toOption, a.assessment)))
    )
    assert(result.intervals.forall(_.isRight))
    val first = result.intervals.head.toOption.get
    val contrast = result.intervals(1).toOption.get
    assert(first.lower < first.upper)
    assertEqualsDouble(contrast.lower, 2 * first.lower, 1e-12)
    assertEqualsDouble(contrast.upper, 2 * first.upper, 1e-12)
    result.attempts.foreach: attempt =>
      assertEquals(attempt.seed, Simulation.seedFor(result.options.seed, attempt.index.toLong))
      val fit = attempt.refit.toOption.get
      assertEquals(fit.options, output.fit.options)
      assertEquals(fit.approximation, Approximation.JointAgq)
      assertEquals(fit.nAgq, 5)
      assertEquals(fit.y, attempt.response.toOption.get)
      println(
        s"GLMM_BOOT_RESPONSE {\"index\":${attempt.index},\"seed\":\"${attempt.seed}\",\"y\":${attempt.response.toOption.get.mkString("[", ",", "]")}}"
      )
      println(
        s"GLMM_BOOT_FIT {\"index\":${attempt.index},\"beta\":${fit.beta.mkString("[", ",", "]")},\"theta\":${fit.theta.mkString("[", ",", "]")},\"objective\":${fit.deviance}}"
      )
    println(s"GLMM_BOOT_INTERVAL qualified=${result.qualifiedCount} lower=${first.lower} upper=${first.upper}")

  test("independent lme4 refits of the actual39 responses reproduce coefficients and percentile intervals"):
    val frozen = JsonValue.parse(GlmmBootstrapFrozen.json)
    val base = frozen.req("base_fit")
    output.fit.beta.zip(base.req("beta").asArray.map(_.asDouble)).foreach((a, b) => assertEqualsDouble(a, b, 1e-5))
    output.fit.theta.zip(base.req("theta").asArray.map(_.asDouble)).foreach((a, b) => assertEqualsDouble(a, b, 1e-5))
    output.prediction
      .zip(base.req("population_mean_x").req("mean").asArray.map(_.asDouble))
      .foreach((a, b) => assertEqualsDouble(a, b, 1e-6))
    output.fit.conditionalMean
      .zip(base.req("conditional_mean").asArray.map(_.asDouble))
      .foreach((a, b) => assertEqualsDouble(a, b, 1e-5))
    val reference = frozen.req("bootstrap")
    val rows = reference.req("replicates").asArray
    assertEquals(rows.length, output.bootstrap.attempts.length)
    rows
      .zip(output.bootstrap.attempts)
      .foreach: (row, attempt) =>
        assertEquals(row.req("index").asInt, attempt.index)
        assertEquals(row.req("seed").asString, attempt.seed.toString)
        assertEquals(row.req("status").asString, "ok")
        val fitted = attempt.refit.toOption.get
        fitted.beta.zip(row.req("beta").asArray.map(_.asDouble)).foreach((a, b) => assertEqualsDouble(a, b, 1e-5))
        fitted.theta.zip(row.req("theta").asArray.map(_.asDouble)).foreach((a, b) => assertEqualsDouble(a, b, 1e-5))
        assertEqualsDouble(fitted.deviance, row.req("minus2log_full_density_normalized").asDouble, 1e-4)
    val interval = output.interval.toOption.get
    val limits = reference.req("percentile_95").req("x").asArray.map(_.asDouble)
    assertEqualsDouble(interval.lower, limits(0), 1e-6)
    assertEqualsDouble(interval.upper, limits(1), 1e-6)
