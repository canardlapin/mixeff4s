package mixeff4s.glmm

import mixeff4s.comparison.{FitCompare, JsonValue}
import mixeff4s.fixtures.GlmmFrames
import mixeff4s.formula.Formula
import mixeff4s.model.Family

class GlmmParitySuite extends munit.FunSuite:
  private val references = JsonValue.parse(GlmmFrozen.json).req("results").asArray
  private def numbers(row: JsonValue, name: String) = row.req(name).asArray.map(_.asDouble)
  private def strings(value: JsonValue): Vector[String] = value match
    case JsonValue.Str(text) => Vector(text)
    case other               => other.asArray.map(_.asString)

  references.foreach: reference =>
    val name = reference.req("dataset").asString
    val joint = reference.req("n_agq").asInt == 1
    test(s"$name ${if joint then "joint-Laplace" else "fast-PIRLS"} against frozen R fit"):
      val original = name match
        case "cbpp"                      => GlmmFrames.cbpp
        case "culcitalogreg"             => GlmmFrames.culcitalogreg
        case "synthetic_poisson_offset"  => GlmmFrames.synthetic_poisson_offset
        case "synthetic_poisson_crossed" => GlmmFrames.synthetic_poisson_crossed
        case other                       => fail(s"unknown frame $other")
      val (frame, source, family, trials, offset) = name match
        case "cbpp" =>
          val sizes = original.numeric("size").get
          val ys = original.numeric("incidence").get.zip(sizes).map(_ / _)
          (
            original.addNumeric("proportion", ys).toOption.get,
            "proportion ~ 1 + period + (1 | herd)",
            Family.Binomial,
            Some(sizes.map(_.toInt)),
            None
          )
        case "culcitalogreg"            => (original, "predation ~ 1 + ttt + (1 | block)", Family.Bernoulli, None, None)
        case "synthetic_poisson_offset" =>
          (
            original,
            "count ~ 1 + x + (1 | group)",
            Family.Poisson,
            None,
            Some(original.numeric("exposure").get.map(math.log))
          )
        case _ =>
          (
            original,
            "count ~ 1 + x + (1 + x | site) + (1 | batch)",
            Family.Poisson,
            None,
            Some(original.numeric("exposure").get.map(math.log))
          )
      val options = (if joint then GlmmOptions.jointLaplace else GlmmOptions.fastLaplace).copy(
        trials = trials,
        offset = offset,
        maxEvaluations = 20000
      )
      val formula = Formula.parse(source).toOption.get
      val fit = Glmm.fit(formula, frame, family, options = options).fold(e => fail(e.message), identity)
      val beta = numbers(reference, "beta")
      assertEquals(fit.feNames.map(FitCompare.contrastName), reference.req("coef_names").asArray.map(_.asString))
      fit.beta.zip(beta).foreach((a, b) => assertEqualsDouble(a, b, 1e-3))
      val thetaTolerance = if !joint && name == "synthetic_poisson_offset" then 0.02 else 1e-3
      assertEquals(fit.theta.length, numbers(reference, "theta").length)
      fit.theta.zip(numbers(reference, "theta")).foreach((a, b) => assertEqualsDouble(a, b, thetaTolerance))
      val objective =
        if joint then fit.deviance
        else
          val data = GlmmResponse
            .prepare(
              family,
              frame
                .numeric(
                  if name == "cbpp" then "proportion" else if name == "culcitalogreg" then "predation" else "count"
                )
                .get,
              options
            )
            .toOption
            .get
          val raw = data.y.indices
            .map(i => data.workingWeights(i) * Pirls.devianceComponent(family, data.y(i), fit.conditionalMean(i)))
            .sum
          fit.deviance - raw + data.minusTwoLogDensity(fit.conditionalEta.toArray)
      println(
        s"GLMM_RESULT name=$name joint=$joint objective=$objective beta=${fit.beta} theta=${fit.theta} stop=${fit.optimizer.stopReason} evals=${fit.optimizer.fevals} assessment=${fit.jointAssessment}"
      )
      val objectiveTolerance = if !joint && name == "synthetic_poisson_offset" then 0.01 else 1e-4
      assertEqualsDouble(objective, reference.req("objective_minus2loglik_with_constants").asDouble, objectiveTolerance)
      val etaTolerance = if !joint && name == "synthetic_poisson_offset" then 0.02 else 2e-3
      fit.conditionalEta
        .zip(numbers(reference, "conditional_eta"))
        .foreach((a, b) => assertEqualsDouble(a, b, etaTolerance))
      assertEquals(fit.conditionalEta.length, reference.req("n").asInt)
      assertEquals(fit.conditionalMean.length, reference.req("n").asInt)
      fit.conditionalMean
        .zip(numbers(reference, "conditional_mu"))
        .foreach: (a, b) =>
          val tolerance =
            if family == Family.Poisson then math.abs(b) * math.expm1(etaTolerance) else etaTolerance / 4.0
          assertEqualsDouble(a, b, tolerance)
      val etaDelta = fit.conditionalEta.zip(numbers(reference, "conditional_eta")).map((a, b) => math.abs(a - b)).max
      val muDelta = fit.conditionalMean.zip(numbers(reference, "conditional_mu")).map((a, b) => math.abs(a - b)).max
      println(s"GLMM_DELTA name=$name joint=$joint eta=$etaDelta mu=$muDelta")
      fit.randomEffects.foreach: term =>
        val modes = reference.req("modes").req(term.groupingName)
        val referenceLevels = strings(modes.req("levels"))
        val referenceNames = strings(modes.req("names"))
        val referenceValues = modes.req("values").asArray.map(_.asArray.map(_.asDouble))
        term.levels
          .zip(term.values)
          .foreach: (level, values) =>
            values
              .zip(term.coefficientNames)
              .foreach: (value, coefficient) =>
                val expected = referenceValues(referenceLevels.indexOf(level))(
                  referenceNames.indexOf(FitCompare.contrastName(coefficient))
                )
                assertEqualsDouble(value, expected, etaTolerance)
      if joint then
        assert(fit.jointAssessment.exists(_.firstOrderPassed), clues(fit.jointAssessment))
        fit.inner match
          case GlmmInnerResult.ConditionalMode(result) => assert(result.scoreMaxAbs <= 1e-8)
          case other                                   => fail(s"unexpected $other")
