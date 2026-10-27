package mixeff4s.glmm

import mixeff4s.comparison.{FitCompare, JsonValue}
import mixeff4s.data.ModelFrame
import mixeff4s.fixtures.GlmmFrames
import mixeff4s.formula.Formula
import mixeff4s.model.Family

class AdaptiveQuadratureSuite extends munit.FunSuite:
  private val source = Formula.parse("y ~ 1 + x + (1 | g)").toOption.get
  private val frame = ModelFrame
    .of(
      "y" -> ModelFrame.numeric(Vector(0.0, 1.0, 1.0, 0.0, 1.0, 0.0)),
      "x" -> ModelFrame.numeric(Vector(-1.0, 0.0, 1.0, -1.0, 0.0, 1.0)),
      "g" -> ModelFrame.factor(Vector("a", "a", "a", "b", "b", "b"))
    )
    .toOption
    .get

  test("every admitted normal quadrature order has independent Gaussian polynomial moments"):
    (1 to AdaptiveQuadrature.MaxNodes).foreach: count =>
      val rule = AdaptiveQuadrature.rule(count).fold(e => fail(e.message), identity)
      val weights = rule.logWeights.map(math.exp)
      assertEqualsDouble(weights.sum, 1.0, 2e-14)
      assertEqualsDouble(rule.nodes.zip(weights).map(_ * _).sum, 0.0, 1e-13)
      (1 to math.min(count - 1, 6)).foreach: power =>
        val expected = (1 to power).map(k => (2 * k - 1).toDouble).product
        val actual = rule.nodes.zip(weights).map((x, w) => w * math.pow(x, 2 * power)).sum
        assertEqualsDouble(actual, expected, 1e-9 * expected)
    val two = AdaptiveQuadrature.rule(2).toOption.get
    two.nodes.zip(Vector(-1.0, 1.0)).foreach((a, b) => assertEqualsDouble(a, b, 1e-14))
    two.logWeights.foreach(w => assertEqualsDouble(math.exp(w), 0.5, 1e-14))

  test("one node is exactly Laplace; zero covariance is exactly the GLM density"):
    val beta = Vector(-0.4, 0.3)
    val laplace = Glmm.jointDeviance(source, frame, Family.Bernoulli, beta, Vector(0.8)).toOption.get
    val one = Glmm
      .jointDeviance(source, frame, Family.Bernoulli, beta, Vector(0.8), options = GlmmOptions.jointAgq(1))
      .toOption
      .get
    assertEquals(one, laplace)
    Vector(2, 5, 9, 15, 31).foreach: nodes =>
      val actual = Glmm
        .jointDeviance(source, frame, Family.Bernoulli, beta, Vector(0.0), options = GlmmOptions.jointAgq(nodes))
        .toOption
        .get
      val expected = Glmm.jointDeviance(source, frame, Family.Bernoulli, beta, Vector(0.0)).toOption.get
      assertEquals(actual, expected)

  test("AGQ refuses fast mode, vector and multiple terms, and out-of-domain orders"):
    Vector(0, 32, -1).foreach: nodes =>
      assert(
        Glmm
          .jointDeviance(
            source,
            frame,
            Family.Bernoulli,
            Vector(0.0, 0.0),
            Vector(1.0),
            options = GlmmOptions.jointAgq(nodes)
          )
          .isLeft
      )
    assert(Glmm.fit(source, frame, Family.Bernoulli, options = GlmmOptions(nAgq = 5)).isLeft)
    val vector = Formula.parse("y ~ 1 + x + (1 + x | g)").toOption.get
    assert(
      Glmm
        .jointDeviance(
          vector,
          frame,
          Family.Bernoulli,
          Vector(0.0, 0.0),
          Vector(1.0, 0.0, 1.0),
          options = GlmmOptions.jointAgq(5)
        )
        .isLeft
    )
    val multi = Formula.parse("y ~ 1 + x + (1 | g) + (0 + x | g)").toOption.get
    assert(
      Glmm
        .jointDeviance(
          multi,
          frame,
          Family.Bernoulli,
          Vector(0.0, 0.0),
          Vector(1.0, 1.0),
          options = GlmmOptions.jointAgq(5)
        )
        .isLeft
    )

  private def values(j: JsonValue, key: String): Vector[Double] = j.req(key).asArray.map(_.asDouble)
  private val references = JsonValue.parse(GlmmAnalysisFrozen.json).req("references")
  private val integrals = JsonValue.parse(GlmmAnalysisFrozen.json).req("integrals")

  for name <- Vector("cbpp", "culcitalogreg", "synthetic_poisson_offset"); nodes <- Vector(5, 9, 15) do
    test(s"$name AGQ$nodes full-density fit against normalized frozen lme4"):
      val (data, formula, family, options) = if name == "cbpp" then
        val raw = GlmmFrames.cbpp
        val n = raw.numeric("size").get
        (
          raw.addNumeric("proportion", raw.numeric("incidence").get.zip(n).map(_ / _)).toOption.get,
          Formula.parse("proportion ~ 1 + period + (1 | herd)").toOption.get,
          Family.Binomial,
          GlmmOptions.jointAgq(nodes).copy(trials = Some(n.map(_.toInt)))
        )
      else if name == "culcitalogreg" then
        (
          GlmmFrames.culcitalogreg,
          Formula.parse("predation ~ 1 + ttt + (1 | block)").toOption.get,
          Family.Bernoulli,
          GlmmOptions.jointAgq(nodes)
        )
      else
        val raw = GlmmFrames.synthetic_poisson_offset
        (
          raw,
          Formula.parse("count ~ 1 + x + (1 | group)").toOption.get,
          Family.Poisson,
          GlmmOptions.jointAgq(nodes).copy(offset = Some(raw.numeric("exposure").get.map(math.log)))
        )
      val reference = references.req("fits").req(name).req(s"nAGQ$nodes")
      val shift = integrals.req("integrals").req(name).req("response_only_saturated_minus2loglik").asDouble
      val fit = Glmm.fit(formula, data, family, options = options).fold(e => fail(e.message), identity)
      assertEquals(fit.approximation, Approximation.JointAgq)
      assertEquals(fit.nAgq, nodes)
      assertEquals(fit.feNames.map(FitCompare.contrastName), reference.req("beta_names").asArray.map(_.asString))
      val full = reference.req("minus2log_lme4_raw").asDouble + shift
      println(
        s"AGQ_FIT name=$name nodes=$nodes objective=${fit.deviance} reference=$full beta=${fit.beta} theta=${fit.theta} assessment=${fit.jointAssessment}"
      )
      assertEqualsDouble(fit.deviance, full, 1e-4)
      fit.beta.zip(values(reference, "beta")).foreach((a, b) => assertEqualsDouble(a, b, 1e-3))
      fit.theta.zip(values(reference, "theta")).foreach((a, b) => assertEqualsDouble(a, b, 1e-3))
      fit.conditionalEta.zip(values(reference, "conditional_eta")).foreach((a, b) => assertEqualsDouble(a, b, 2e-3))
      assert(fit.jointAssessment.exists(_.firstOrderPassed), clues(fit.jointAssessment))
      val evaluated = Glmm.jointDeviance(formula, data, family, fit.beta, fit.theta, options = options).toOption.get
      assertEqualsDouble(
        evaluated.objective,
        evaluated.minusTwoLogDensity + evaluated.penalty + evaluated.logdet + evaluated.quadratureCorrection,
        1e-10
      )

  test("AGQ at independently integrated fixed parameters preserves full response constants"):
    integrals
      .req("fixed_parameter_oracles")
      .asArray
      .foreach: oracle =>
        val xs = values(oracle, "x")
        val counts = values(oracle, "y")
        val family = if oracle.req("family").asString == "binomial" then Family.Binomial else Family.Poisson
        val ns = oracle.field("trials").collect { case JsonValue.Arr(values) => values.map(_.asInt) }
        val ys = ns.fold(counts)(n => counts.zip(n).map((k, size) => k / size.toDouble))
        val data = ModelFrame
          .of(
            "y" -> ModelFrame.numeric(ys),
            "x" -> ModelFrame.numeric(xs),
            "z" -> ModelFrame.numeric(values(oracle, "z")),
            "g" -> ModelFrame.factor(Vector.fill(xs.length)("all"))
          )
          .toOption
          .get
        val options = GlmmOptions
          .jointAgq(31)
          .copy(
            trials = ns,
            offset =
              oracle.field("exposure").collect { case JsonValue.Arr(values) => values.map(v => math.log(v.asDouble)) }
          )
        val scalarSlope = Formula.parse("y ~ 1 + x + (0 + z | g)").toOption.get
        val actual = Glmm
          .jointDeviance(
            scalarSlope,
            data,
            family,
            values(oracle, "beta"),
            Vector(oracle.req("theta").asDouble),
            options = options
          )
          .fold(e => fail(e.message), identity)
        assertEqualsDouble(actual.objective, oracle.req("minus2log_full_density").asDouble, 1e-7)

  test("maximum-order quadrature agrees with independent full-density integrals at frozen parameters"):
    Vector("cbpp", "culcitalogreg", "synthetic_poisson_offset").foreach: name =>
      val reference = references.req("fits").req(name).req("nAGQ15")
      val (data, formula, family, options) = if name == "cbpp" then
        val raw = GlmmFrames.cbpp
        val sizes = raw.numeric("size").get
        (
          raw.addNumeric("proportion", raw.numeric("incidence").get.zip(sizes).map(_ / _)).toOption.get,
          Formula.parse("proportion ~ 1 + period + (1 | herd)").toOption.get,
          Family.Binomial,
          GlmmOptions.jointAgq(31).copy(trials = Some(sizes.map(_.toInt)))
        )
      else if name == "culcitalogreg" then
        (
          GlmmFrames.culcitalogreg,
          Formula.parse("predation ~ 1 + ttt + (1 | block)").toOption.get,
          Family.Bernoulli,
          GlmmOptions.jointAgq(31)
        )
      else
        val raw = GlmmFrames.synthetic_poisson_offset
        (
          raw,
          Formula.parse("count ~ 1 + x + (1 | group)").toOption.get,
          Family.Poisson,
          GlmmOptions.jointAgq(31).copy(offset = Some(raw.numeric("exposure").get.map(math.log)))
        )
      val actual = Glmm
        .jointDeviance(formula, data, family, values(reference, "beta"), values(reference, "theta"), options = options)
        .fold(e => fail(e.message), identity)
      val expected = integrals.req("integrals").req(name).req("minus2log_full_density").asDouble
      println(
        s"AGQ_INTEGRAL name=$name nodes=31 objective=${actual.objective} independent=$expected error=${math.abs(actual.objective - expected)}"
      )
      assertEqualsDouble(actual.objective, expected, 1e-4)
