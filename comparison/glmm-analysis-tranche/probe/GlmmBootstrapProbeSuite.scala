package mixeff4s.glmm

import mixeff4s.fixtures.GlmmAnalysisFrames
import mixeff4s.formula.Formula
import mixeff4s.stats.{BootstrapTarget, GlmmBootstrap, GlmmBootstrapOptions}
import scala.concurrent.duration.*

/** Initial end-to-end qualification probe; all outcomes are printed and retained. */
class GlmmBootstrapProbeSuite extends munit.FunSuite:
  override val munitTimeout = 120.seconds
  test("grouped-binomial chosen-estimator bootstrap records all attempts"):
    val options = GlmmOptions.jointAgq(5).copy(trials = Some(GlmmAnalysisFrames.trials))
    val fit = Glmm.fit(Formula.parse("y ~ 1 + x + (1 | g)").toOption.get,
      GlmmAnalysisFrames.binomial, mixeff4s.model.Family.Binomial, options = options).fold(e => fail(e.message), identity)
    println(s"GLMM_BOOT_ORIGINAL beta=${fit.beta} theta=${fit.theta} objective=${fit.deviance} assessment=${GlmmAssessment.assess(fit)}")
    val result = GlmmBootstrap.parametric(fit, GlmmBootstrapOptions(39, 20261027L,
      Vector(BootstrapTarget.Coefficient("x")))).fold(e => fail(e.message), identity)
    result.attempts.foreach: attempt =>
      println(s"GLMM_BOOT_ATTEMPT index=${attempt.index} seed=${attempt.seed} qualified=${attempt.qualified} error=${attempt.refit.left.toOption.map(_.code)} assessment=${attempt.assessment.map(_.refusalCodes)}")
      println(s"GLMM_BOOT_RESPONSE {\"index\":${attempt.index},\"seed\":\"${attempt.seed}\",\"y\":${attempt.response.toOption.get.mkString("[",",","]")}}")
      attempt.refit.foreach(f => println(s"GLMM_BOOT_FIT {\"index\":${attempt.index},\"beta\":${f.beta.mkString("[",",","]")},\"theta\":${f.theta.mkString("[",",","]")},\"objective\":${f.deviance}}"))
    println(s"GLMM_BOOT_INTERVAL qualified=${result.qualifiedCount} interval=${result.intervals}")
    assertEquals(result.attempts.length, 39)
