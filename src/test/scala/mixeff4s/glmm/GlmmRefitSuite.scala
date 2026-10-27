package mixeff4s.glmm

import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.fixtures.GlmmFrames
import mixeff4s.formula.Formula
import mixeff4s.lmm.SimulationMode
import mixeff4s.model.Family

class GlmmRefitSuite extends munit.FunSuite:
  private val original = GlmmFrames.synthetic_poisson_offset
  private val formula = Formula.parse("count ~ 1 + x + (1 | group)").toOption.get
  private val a = original.numeric("count").get
  private val b = a.zipWithIndex.map((y, i) => y + (i % 3).toDouble)
  private val offsets = original.numeric("exposure").get.map(math.log)
  private def frame(y: Vector[Double]): ModelFrame = ModelFrame
    .of(original.columns.map:
      case ("count", _) => "count" -> Column.Numeric(y)
      case entry        => entry
    *)
    .toOption
    .get
  private def equalFit(a: GlmmFit, b: GlmmFit): Unit =
    assertEqualsDouble(a.deviance, b.deviance, 1e-8)
    a.beta.zip(b.beta).foreach((x, y) => assertEqualsDouble(x, y, 1e-6))
    a.theta.zip(b.theta).foreach((x, y) => assertEqualsDouble(x, y, 1e-6))
    assertEquals(a.options, b.options)
    assertEquals(a.y, b.y)
    assertEquals(a.approximation, b.approximation)
    assertEquals(a.fixedEffects, b.fixedEffects)

  Vector(GlmmOptions.fastLaplace, GlmmOptions.jointLaplace, GlmmOptions.jointAgq(5)).foreach: chosen =>
    test(s"${chosen.fast}/${chosen.nAgq}: A B invalid A cold refits equal independent fresh fits"):
      val options = chosen.copy(offset = Some(offsets))
      val plan = Glmm.prepare(formula, original, Family.Poisson).toOption.get
      val session = plan.session(options).toOption.get
      val first = session.refit(a).fold(e => fail(e.message), identity)
      val oldSnapshot = (
        first.beta,
        first.theta,
        first.deviance,
        first.conditionalMean,
        first.randomEffects,
        GlmmPrediction.fitted(first),
        GlmmSimulation.response(first, 42L, SimulationMode.ConditionOnModes)
      )
      val second = session.refit(b).fold(e => fail(e.message), identity)
      assert(first.y != second.y)
      equalFit(
        second,
        Glmm.fit(formula, frame(b), Family.Poisson, options = options).fold(e => fail(e.message), identity)
      )
      assert(session.refit(Vector(1.0)).left.toOption.exists(_.code == "dimension_mismatch"))
      assert(session.refit(a.updated(0, Double.NaN)).isLeft)
      val again = session.refit(a).fold(e => fail(e.message), identity)
      equalFit(
        again,
        Glmm.fit(formula, frame(a), Family.Poisson, options = options).fold(e => fail(e.message), identity)
      )
      assertEquals(
        oldSnapshot,
        (
          first.beta,
          first.theta,
          first.deviance,
          first.conditionalMean,
          first.randomEffects,
          GlmmPrediction.fitted(first),
          GlmmSimulation.response(first, 42L, SimulationMode.ConditionOnModes)
        )
      )
      val batch = session.refitAll(Vector(a, Vector(1.0), b))
      assertEquals(batch.map(_.isRight), Vector(true, false, true))
      equalFit(batch(2).toOption.get, second)

  test("separate sessions cannot contaminate each other and controls remain authoritative"):
    val plan = Glmm.prepare(formula, original, Family.Poisson).toOption.get
    val options = GlmmOptions.jointAgq(5).copy(offset = Some(offsets))
    val left = plan.session(options).toOption.get
    val right = plan.session(options).toOption.get
    val first = left.refit(a).toOption.get
    assert(right.refit(b).isRight)
    equalFit(left.refit(a).toOption.get, first)
    val capped = plan.session(options.copy(maxEvaluations = 1)).toOption.get
    assert(capped.refit(a).left.toOption.exists(_.code == "optimizer_budget_exhausted"))
    assertEquals(options, left.options)

  test("binomial refits preserve trials and replace cached response constants"):
    val raw = GlmmFrames.cbpp
    val sizes = raw.numeric("size").get.map(_.toInt)
    val ys = raw.numeric("incidence").get.zip(sizes).map((k, n) => k / n.toDouble)
    val changed = raw.numeric("incidence").get.zip(sizes).map((k, n) => math.min(k + 1.0, n.toDouble) / n.toDouble)
    val data = raw.addNumeric("proportion", ys).toOption.get
    val f = Formula.parse("proportion ~ 1 + period + (1 | herd)").toOption.get
    val options = GlmmOptions.jointAgq(5).copy(trials = Some(sizes))
    val plan = Glmm.prepare(f, data, Family.Binomial).toOption.get
    val session = plan.session(options).toOption.get
    assert(session.refit(ys).isRight)
    val result = session.refit(changed).fold(e => fail(e.message), identity)
    equalFit(result, plan.fit(changed, options).fold(e => fail(e.message), identity))
    assertEquals(result.options.trials, Some(sizes))
    assertEquals(result.nAgq, 5)
    assert(session.refit(ys.updated(0, 0.12345)).isLeft)
    equalFit(session.refit(changed).toOption.get, result)

  test("large-trial simulated proportions roundtrip through response preparation and refit validation"):
    val sizes = Vector(1000000000, 1000000007, Int.MaxValue, Int.MaxValue)
    val counts = Vector(130850654, 257111068, 1000000000, 1200000000)
    val proportions = counts.zip(sizes).map((k, n) => k.toDouble / n.toDouble)
    val options = GlmmOptions.jointLaplace.copy(trials = Some(sizes), maxEvaluations = 1)
    val prepared = GlmmResponse.prepare(Family.Binomial, proportions, options).fold(e => fail(e.message), identity)
    assertEquals(prepared.y, proportions)
    val raw = ModelFrame
      .of("y" -> ModelFrame.numeric(proportions), "g" -> ModelFrame.factor(Vector("a", "a", "b", "b")))
      .toOption
      .get
    val plan = Glmm.prepare("y ~ 1 + (1 | g)", raw, Family.Binomial).toOption.get
    val session = plan.session(options).fold(e => fail(e.message), identity)
    val attempt = session.refit(proportions)
    assert(!attempt.left.toOption.exists(_.code == "invalid_argument"), clues(attempt.left.toOption))
    assert(
      GlmmResponse.prepare(Family.Binomial, proportions.updated(0, proportions(0) + 0.25 / sizes.head), options).isLeft
    )

  test("reusable and direct front doors share the same explicit fast-PIRLS default"):
    val data = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(0.0, 1.0, 1.0, 0.0, 0.0, 1.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c"))
      )
      .toOption
      .get
    val plan = Glmm.prepare("y ~ 1 + (1 | g)", data, Family.Bernoulli).toOption.get
    val reusable = plan.fit().toOption.get
    val direct = Glmm.fit("y ~ 1 + (1 | g)", data, Family.Bernoulli).toOption.get
    assertEquals(reusable.approximation, Approximation.FastPirls)
    equalFit(reusable, direct)
