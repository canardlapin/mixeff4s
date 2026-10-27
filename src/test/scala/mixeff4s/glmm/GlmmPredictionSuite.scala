package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.fixtures.GlmmFrames
import mixeff4s.formula.Formula
import mixeff4s.lmm.{NewGroupPolicy, PredictionMode}
import mixeff4s.model.Family

class GlmmPredictionSuite extends munit.FunSuite:
  private def right[A](value: Either[?, A]): A = value.fold(error => fail(error.toString), identity)
  private lazy val fit =
    right(
      Glmm.fit(
        Formula.parse("predation ~ 1 + ttt + (1 | block)").toOption.get,
        GlmmFrames.culcitalogreg,
        Family.Bernoulli,
        options = GlmmOptions.fastLaplace.copy(maxEvaluations = 10000)
      )
    )

  test("fitted link and response predictions preserve the fitted conditional snapshot"):
    val link = right(GlmmPrediction.fitted(fit, scale = GlmmPredictionScale.Link))
    val response = right(GlmmPrediction.fitted(fit))
    assertEquals(link.length, fit.n)
    link.zip(fit.conditionalEta).foreach((actual, expected) => assertEqualsDouble(actual, expected, 1e-8))
    response.zip(fit.conditionalMean).foreach((actual, expected) => assertEqualsDouble(actual, expected, 1e-8))

  test("population prediction removes modes and frozen encoding refuses unknown groups unless zero-effect is explicit"):
    val data = GlmmFrames.culcitalogreg
    val conditional = right(GlmmPrediction.predict(fit, data))
    val population = right(GlmmPrediction.predict(fit, data, PredictionMode.Population))
    assertEquals(conditional.length, population.length)
    assert(conditional.zip(population).exists((a, b) => math.abs(a - b) > 1e-12))
    val ttt = data.factor("ttt").get
    val unknown = right(
      ModelFrame.of(
        "ttt" -> ModelFrame.factor(ttt.values, ttt.levels).toOption.get,
        "block" -> ModelFrame.factor(Vector.fill(data.nRows)("new-block"))
      )
    )
    assertEquals(GlmmPrediction.predict(fit, unknown).left.toOption.get.code, "prediction_new_group")
    assert(GlmmPrediction.predict(fit, unknown, unknownGroups = NewGroupPolicy.ZeroEffect).isRight)

  test("raw and Pearson residuals are finite"):
    assert(right(GlmmPrediction.residuals(fit)).forall(_.isFinite))
    assert(right(GlmmPrediction.residuals(fit, GlmmResidualScale.Pearson)).forall(_.isFinite))

  test("new-data factor level order is aligned to training and incompatible columns refuse"):
    val raw = GlmmFrames.culcitalogreg
    val ttt = raw.factor("ttt").get
    val reordered = right(
      ModelFrame.of(
        "ttt" -> ModelFrame.factor(ttt.values, ttt.levels.reverse).toOption.get,
        "block" -> raw.column("block").get
      )
    )
    val expected = right(GlmmPrediction.predict(fit, raw))
    right(GlmmPrediction.predict(fit, reordered)).zip(expected).foreach((a, b) => assertEqualsDouble(a, b, 1e-12))
    val wrongType =
      right(ModelFrame.of("ttt" -> ModelFrame.numeric(Vector.fill(raw.nRows)(1.0)), "block" -> raw.column("block").get))
    assertEquals(GlmmPrediction.predict(fit, wrongType).left.toOption.get.code, "prediction_column_type")
    val unknownFixed = right(
      ModelFrame.of("ttt" -> ModelFrame.factor(Vector.fill(raw.nRows)("unknown")), "block" -> raw.column("block").get)
    )
    assertEquals(GlmmPrediction.predict(fit, unknownFixed).left.toOption.get.code, "prediction_unknown_fixed_level")

  test("offset fits require explicit compatible offsets for new data and b-zero differs from marginal mean"):
    val data = GlmmFrames.synthetic_poisson_offset
    val offsets = data.numeric("exposure").get.map(math.log)
    val fitted = right(
      Glmm.fit(
        Formula.parse("count ~ 1 + x + (1 | group)").toOption.get,
        data,
        Family.Poisson,
        options = GlmmOptions.jointLaplace.copy(offset = Some(offsets))
      )
    )
    assertEquals(GlmmPrediction.predict(fitted, data).left.toOption.get.code, "prediction_offset_required")
    val predicted = right(GlmmPrediction.predict(fitted, data, offset = Some(offsets)))
    predicted.zip(fitted.conditionalMean).foreach((a, b) => assertEqualsDouble(a, b, 1e-10))
    assert(GlmmPrediction.predict(fitted, data, offset = Some(Vector(0.0))).isLeft)
    assert(GlmmPrediction.predict(fitted, data, offset = Some(Vector.fill(data.nRows)(Double.NaN))).isLeft)
    val population = right(
      GlmmPrediction.predict(
        fitted,
        data,
        PredictionMode.Population,
        GlmmPredictionScale.Response,
        offset = Some(offsets)
      )
    )
    population
      .zip(data.numeric("x").get)
      .zip(offsets)
      .foreach: (pair, offset) =>
        val (value, x) = pair
        val zeroEffect = math.exp(fitted.beta(0) + fitted.beta(1) * x + offset)
        assertEqualsDouble(value, zeroEffect, 1e-12)
        val marginal = zeroEffect * math.exp(0.5 * fitted.theta.head * fitted.theta.head)
        assert(marginal > value)
