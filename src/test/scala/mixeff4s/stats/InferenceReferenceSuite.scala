package mixeff4s.stats

import mixeff4s.comparison.JsonValue
import mixeff4s.examples.LikelihoodInference
import mixeff4s.fixtures.{GlmmAnalysisFrames, Sleepstudy}
import mixeff4s.glmm.*
import mixeff4s.lmm.Lmm
import mixeff4s.model.Family

class InferenceReferenceSuite extends munit.FunSuite:
  private val reference = JsonValue.parse(InferenceFrozen.json)
  private def array(x: JsonValue) = x.asArray.map(_.asDouble)
  private lazy val gaussian =
    Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, LikelihoodInference.tightML).toOption.get
  private val profileOptions =
    LmmProfileOptions(maxEvaluations = 2000, optimizerControls = LikelihoodInference.tightML.optimizerControls.get)
  Vector(1, 5).foreach: nodes =>
    test(s"independent R selected-objective curvature and nested likelihood ratio nAgq=$nodes"):
      val options = GlmmOptions.jointAgq(nodes).copy(trials = Some(GlmmAnalysisFrames.trials))
      val large =
        Glmm.fit("y ~ 1 + x + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options).toOption.get
      val small = Glmm.fit("y ~ 1 + (1 | g)", GlmmAnalysisFrames.binomial, Family.Binomial, None, options).toOption.get
      val result = GlmmInference.coefficients(large).fold(e => fail(e.message), identity)
      val root = reference.req("binomial_workflow")
      val curvature = root.req(if nodes == 1 then "central_difference" else "selected_nAGQ5_central_difference")
      val matrix = curvature.req("covariance_2_hessian_inverse_beta_block").asArray.map(array)
      result.covariance.zip(matrix).foreach((a, b) => a.zip(b).foreach((x, y) => assertEqualsDouble(x, y, 1e-7)))
      val expected = root.req("likelihood_ratio_small_null").req(s"nAGQ$nodes").req("statistic").asDouble
      assertEqualsDouble(GlmmLrt.statistic(small, large).toOption.get.likelihoodRatio, expected, 1e-5)
      val full = root.req("full_objective_by_estimator").req(s"nAGQ$nodes").req("full_normalized_minus2loglik").asDouble
      assertEqualsDouble(large.deviance, full, 1e-5)
      println(s"R_CURVATURE nodes=$nodes covariance=${result.covariance} reference=$matrix LR=$expected")
  test("ML constrained profile objectives match an independent dense Gaussian likelihood"):
    val expected = reference.req("sleepstudy_ml_profile")
    assertEqualsDouble(gaussian.objective, expected.req("fit").req("exact_marginal_normal_minus2loglik").asDouble, 1e-7)
    expected
      .req("constrained_ml_trace")
      .asArray
      .foreach: row =>
        val name = row.req("parameter").asString
        val target =
          if name == "days" then LmmProfileTarget.Coefficient("days")
          else if name == "random_intercept_sd" then LmmProfileTarget.RandomEffectStdDev("subj")
          else LmmProfileTarget.ResidualStdDev
        val value = row.req("constrained_value").asDouble
        val point = LmmProfile.evaluate(gaussian, target, value, profileOptions).toOption.get
        val actual = point.objective.fold(e => fail(s"${e.code}: ${e.message}"), identity)
        assertEqualsDouble(actual, row.req("minus2loglik").asDouble, 1e-5)
        println(
          s"R_PROFILE_POINT target=$target value=$value objective=$actual reference=${row.req("minus2loglik").asDouble}"
        )
  test("ML profile intervals agree with frozen lme4 and variance is the squared SD axis"):
    val rows = reference.req("sleepstudy_ml_profile").req("profile").req("intervals").asArray.map(array)
    val targets = Vector(
      LmmProfileTarget.RandomEffectStdDev("subj"),
      LmmProfileTarget.ResidualStdDev,
      LmmProfileTarget.Coefficient("days")
    )
    val profiles = targets.map(target => LmmProfile.confint(gaussian, target, profileOptions).toOption.get)
    profiles
      .zip(rows)
      .foreach: (profile, expected) =>
        assertEqualsDouble(profile.lower.value.get, expected(0), 1e-3)
        assertEqualsDouble(profile.upper.value.get, expected(1), 1e-3)
        Vector(profile.lower, profile.upper).foreach: endpoint =>
          assert(endpoint.bracket.nonEmpty)
          assertEqualsDouble(endpoint.objectiveDifference.get, profile.cutoff, 2e-4)
    val variance =
      LmmProfile.confint(gaussian, LmmProfileTarget.RandomEffectVariance("subj"), profileOptions).toOption.get
    assertEqualsDouble(math.sqrt(variance.lower.value.get), profiles.head.lower.value.get, 1e-5)
    assertEqualsDouble(math.sqrt(variance.upper.value.get), profiles.head.upper.value.get, 1e-5)
