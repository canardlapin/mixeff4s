package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.formula.Formula
import mixeff4s.model.{Family, Link}

class GlmmInputSuite extends munit.FunSuite:
  private val formula = Formula.parse("y ~ 1 + (1 | g)").toOption.get
  private def frame(y: Vector[Double]) = ModelFrame
    .of("y" -> ModelFrame.numeric(y), "g" -> ModelFrame.factor(y.indices.map(i => s"g${i / 2}").toVector))
    .toOption
    .get

  test("binomial counts and ambiguous fractional responses refuse in both estimators"):
    Vector(GlmmOptions.fastLaplace, GlmmOptions.jointLaplace).foreach: options =>
      Vector(Vector(0.0, 2.0, 1.0, 1.0), Vector(0.0, 0.5, 1.0, 0.5)).foreach: ys =>
        assert(Glmm.fit(formula, frame(ys), Family.Binomial, options = options).isLeft)
      val ys = Vector(0.0, 0.5, 1.0, 0.5)
      assert(
        Glmm.fit(formula, frame(ys), Family.Binomial, options = options.copy(trials = Some(Vector.fill(4)(2)))).isRight
      )

  test("invalid trials, cases, offsets, controls and fractional Poisson counts refuse"):
    val ys = Vector(0.0, 1.0, 0.0, 1.0)
    Vector(Some(Vector(1, 0, 2, 2)), Some(Vector(1, 2))).foreach: trials =>
      assert(Glmm.fit(formula, frame(ys), Family.Binomial, options = GlmmOptions(trials = trials)).isLeft)
    assert(
      Glmm
        .fit(
          formula,
          frame(Vector(0.0, 0.2, 1.0, 0.5)),
          Family.Binomial,
          options = GlmmOptions(trials = Some(Vector.fill(4)(2)))
        )
        .isLeft
    )
    assert(
      Glmm.fit(formula, frame(ys), Family.Bernoulli, options = GlmmOptions(trials = Some(Vector.fill(4)(1)))).isLeft
    )
    assert(Glmm.fit(formula, frame(Vector(0.0, 1.5, 2.0, 1.0)), Family.Poisson).isLeft)
    Vector(
      GlmmOptions(weights = Some(Vector(1.0, 0.0, 1.0, 1.0))),
      GlmmOptions(offset = Some(Vector.fill(4)(Double.NaN))),
      GlmmOptions(maxEvaluations = 0),
      GlmmOptions(conditionalScoreTolerance = Double.NaN)
    ).foreach: options =>
      assert(Glmm.fit(formula, frame(ys), Family.Bernoulli, options = options).isLeft)

  test("estimated-dispersion families refuse until dispersion is explicitly implemented"):
    Vector(Family.Gamma, Family.InverseGaussian, Family.Normal).foreach: family =>
      Vector(GlmmOptions.fastLaplace, GlmmOptions.jointLaplace).foreach: options =>
        val result = Glmm.fit(formula, frame(Vector(0.5, 1.0, 1.5, 2.0)), family, Some(Link.Log), options)
        assert(result.left.toOption.exists(_.code == "unsupported"))

  test("logit inverse and derivative stay finite at large eta"):
    Vector(-1000.0, -50.0, 0.0, 50.0, 1000.0).foreach: eta =>
      assert(Link.Logit.linkinv(eta).isFinite)
      assert(Link.Logit.muEta(eta).isFinite)
    assertEqualsDouble(Link.Logit.linkinv(1000.0), 1.0, 0.0)
    assertEqualsDouble(Link.Logit.linkinv(-1000.0), 0.0, 0.0)
