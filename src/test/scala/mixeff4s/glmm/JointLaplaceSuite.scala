package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.error.MixedModelError
import mixeff4s.formula.Formula
import mixeff4s.model.{Family, Link}

class JointLaplaceSuite extends munit.FunSuite:
  private val y = Vector(0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0)
  private val x = Vector(-1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 2.0)
  private val g = Vector.tabulate(12)(i => s"g${i / 4}")
  private val formula = Formula.parse("y ~ 1 + x + (1 | g)").toOption.get
  private val frame =
    ModelFrame.of("y" -> ModelFrame.numeric(y), "x" -> ModelFrame.numeric(x), "g" -> ModelFrame.factor(g)).toOption.get
  private val offset = Vector.tabulate(12)(i => (i % 3 - 1) * 0.15)

  private def evaluate(
      beta: Vector[Double],
      theta: Double,
      family: Family = Family.Bernoulli,
      data: ModelFrame = frame,
      options: GlmmOptions
  ): JointEvaluation =
    Glmm
      .jointDeviance(formula, data, family, beta, Vector(theta), options = options)
      .fold(e => fail(e.message), identity)

  // Independent scalar group oracle: raw row score bisection, no engine design/factors/PIRLS.
  private def oracle(
      beta: Vector[Double],
      theta: Double,
      ys: Vector[Double],
      trials: Vector[Int],
      weights: Vector[Double],
      offsets: Vector[Double],
      poisson: Boolean = false
  ): (Double, Vector[Double], Double) =
    var total = 0.0
    var logdet = 0.0
    val modes = Vector.tabulate(3): group =>
      val indices = (group * 4 until (group + 1) * 4)
      def mean(i: Int, u: Double): Double =
        val eta = beta(0) + beta(1) * x(i) + offsets(i) + theta * u
        if poisson then math.exp(eta) else 1.0 / (1.0 + math.exp(-eta))
      def score(u: Double): Double = u + theta * indices.map(i => weights(i) * trials(i) * (mean(i, u) - ys(i))).sum
      var lo = -1.0
      var hi = 1.0
      while score(lo) > 0.0 do lo *= 2.0
      while score(hi) < 0.0 do hi *= 2.0
      (0 until 100).foreach: _ =>
        val mid = (lo + hi) / 2.0
        if score(mid) > 0.0 then hi = mid else lo = mid
      val u = (lo + hi) / 2.0
      var curvature = 1.0
      total += u * u
      indices.foreach: i =>
        val m = mean(i, u)
        val n = trials(i)
        val k = math.rint(n * ys(i)).toInt
        val combinatorial = (1 to k).map(j => math.log((n - k + j).toDouble / j)).sum
        val nll =
          if poisson then m - ys(i) * math.log(m) + (1 to ys(i).toInt).map(j => math.log(j.toDouble)).sum
          else -k * math.log(m) - (n - k) * math.log1p(-m) - combinatorial
        total += 2.0 * weights(i) * nll
        curvature += theta * theta * weights(i) * n * (if poisson then m else m * (1.0 - m))
      logdet += math.log(curvature)
      u
    (total + logdet, modes, logdet)

  test("fixed beta modes and final-mean determinant agree with independent scalar oracle"):
    val coefficients = Vector(-0.7, 0.4)
    Vector(0.0, 0.2, 1.3, 3.0).foreach: theta =>
      val actual = evaluate(coefficients, theta, options = GlmmOptions.jointLaplace.copy(offset = Some(offset)))
      val expected = oracle(coefficients, theta, y, Vector.fill(12)(1), Vector.fill(12)(1.0), offset)
      assertEqualsDouble(actual.objective, expected._1, 1e-8)
      assertEqualsDouble(actual.logdet, expected._3, 1e-9)
      actual.sphericalModes.head.zip(expected._2).foreach((a, b) => assertEqualsDouble(a, b, 1e-8))
      assert(actual.conditional.scoreMaxAbs <= 1e-8)
      assertEqualsDouble(actual.objective, actual.minusTwoLogDensity + actual.penalty + actual.logdet, 1e-12)

  test("cold A B A probes are independent of previous modes and weights"):
    val design = Glmm.compile(formula, frame).toOption.get
    val options = GlmmOptions.jointLaplace
    val response = GlmmResponse.prepare(Family.Bernoulli, y, options).toOption.get
    val workspace = JointWorkspace.create(design, Family.Bernoulli, Link.Logit, response, options).toOption.get
    val a = workspace.evaluate(Vector(-0.5, 0.2), Vector(1.2)).toOption.get
    assert(workspace.evaluate(Vector(1.5, -0.4), Vector(0.1)).isRight)
    assertEquals(workspace.evaluate(Vector(-0.5, 0.2), Vector(1.2)).toOption.get, a)

  test("binomial trials and case weights affect constants, score and curvature consistently"):
    val trials = Vector(2, 4, 5, 3, 2, 4, 5, 3, 2, 4, 5, 3)
    val successes = Vector(0, 1, 3, 1, 1, 3, 5, 2, 0, 2, 1, 0)
    val proportions = successes.zip(trials).map((k, n) => k.toDouble / n)
    val data = ModelFrame
      .of("y" -> ModelFrame.numeric(proportions), "x" -> ModelFrame.numeric(x), "g" -> ModelFrame.factor(g))
      .toOption
      .get
    val weights = Vector.tabulate(12)(i => 0.5 + (i % 4) * 0.4)
    val options = GlmmOptions.jointLaplace.copy(trials = Some(trials), weights = Some(weights), offset = Some(offset))
    val beta = Vector(-0.4, 0.3)
    val actual = evaluate(beta, 0.9, Family.Binomial, data, options)
    val expected = oracle(beta, 0.9, proportions, trials, weights, offset)
    assertEqualsDouble(actual.objective, expected._1, 1e-8)
    assertEqualsDouble(actual.logdet, expected._3, 1e-9)
    actual.sphericalModes.head.zip(expected._2).foreach((a, b) => assertEqualsDouble(a, b, 1e-8))

  test("Poisson constants, zero counts, offsets and power weights agree with scalar oracle"):
    val counts = Vector(0.0, 1.0, 0.0, 2.0, 3.0, 2.0, 7.0, 4.0, 0.0, 0.0, 1.0, 2.0)
    val data = ModelFrame
      .of("y" -> ModelFrame.numeric(counts), "x" -> ModelFrame.numeric(x), "g" -> ModelFrame.factor(g))
      .toOption
      .get
    val weights = Vector.tabulate(12)(i => 0.7 + (i % 3) * 0.3)
    val beta = Vector(0.2, 0.3)
    val options = GlmmOptions.jointLaplace.copy(offset = Some(offset), weights = Some(weights))
    val actual = evaluate(beta, 0.8, Family.Poisson, data, options)
    val expected = oracle(beta, 0.8, counts, Vector.fill(12)(1), weights, offset, poisson = true)
    assertEqualsDouble(actual.objective, expected._1, 1e-8)
    assertEqualsDouble(actual.logdet, expected._3, 1e-9)

  test("joint API refuses wrong dimensions and invalid theta"):
    assert(
      Glmm
        .jointDeviance(formula, frame, Family.Bernoulli, Vector(0.0), Vector(1.0))
        .left
        .toOption
        .exists(_.isInstanceOf[MixedModelError.DimensionMismatch])
    )
    Vector(Double.NaN, -0.1).foreach: theta =>
      assert(Glmm.jointDeviance(formula, frame, Family.Bernoulli, Vector(0.0, 0.0), Vector(theta)).isLeft)

  test("conditional and outer exhaustion are typed failures"):
    val options = GlmmOptions.jointLaplace.copy(jointMaxIterations = 1)
    Glmm.jointDeviance(formula, frame, Family.Bernoulli, Vector(-0.4, 0.3), Vector(2.0), options = options) match
      case Left(_: MixedModelError.PirlsNotConverged) => ()
      case other                                      => fail(s"expected conditional exhaustion, got $other")
    Glmm.fit(formula, frame, Family.Bernoulli, options = GlmmOptions.jointLaplace.copy(maxEvaluations = 1)) match
      case Left(_: MixedModelError.OptimizerBudgetExhausted) => ()
      case other                                             => fail(s"expected outer exhaustion, got $other")

  test("joint fit is separately labelled and records conditional and outer evidence"):
    val fit = Glmm
      .fit(formula, frame, Family.Bernoulli, options = GlmmOptions.jointLaplace)
      .fold(e => fail(e.message), identity)
    assertEquals(fit.approximation, Approximation.JointLaplace)
    assertEquals(fit.objectiveDefinition, "minus_two_loglikelihood_with_response_constants")
    assertEquals(fit.optimizer.x, fit.beta ++ fit.theta)
    fit.inner match
      case GlmmInnerResult.ConditionalMode(result) => assert(result.scoreMaxAbs <= result.scoreTolerance)
      case other                                   => fail(s"wrong inner evidence $other")
    assert(fit.jointAssessment.nonEmpty)
    assert(fit.conditionalMean.forall(p => p >= 0.0 && p <= 1.0))

  test("near-integer accepted input is canonicalized once even under large case weights"):
    val options = GlmmOptions.jointLaplace.copy(weights = Some(Vector.fill(12)(1e8)), trials = Some(Vector.fill(12)(1)))
    val canonical = y
    val near = y.map(v => if v == 0.0 then 1e-9 else v)
    def data(ys: Vector[Double]) = ModelFrame
      .of("y" -> ModelFrame.numeric(ys), "x" -> ModelFrame.numeric(x), "g" -> ModelFrame.factor(g))
      .toOption
      .get
    val a = evaluate(Vector(-0.4, 0.3), 0.0, Family.Binomial, data(canonical), options)
    val b = evaluate(Vector(-0.4, 0.3), 0.0, Family.Binomial, data(near), options)
    assertEquals(a, b)

  test("zero theta at extreme finite logit eta retains finite likelihood and identity Hessian"):
    Vector(-40.0, 40.0, -1000.0, 1000.0).foreach: eta =>
      val actual = evaluate(Vector(eta, 0.0), 0.0, options = GlmmOptions.jointLaplace)
      val logPartition = math.max(eta, 0.0) + math.log1p(math.exp(-math.abs(eta)))
      val expected = y.map(yi => 2.0 * (logPartition - yi * eta)).sum
      assertEqualsDouble(actual.objective, expected, 1e-10)
      assertEqualsDouble(actual.logdet, 0.0, 0.0)
      assertEqualsDouble(actual.conditional.scoreMaxAbs, 0.0, 0.0)

  test("stable logistic score retains small residuals under large case weights"):
    val response = GlmmResponse
      .prepare(Family.Bernoulli, Vector(1.0), GlmmOptions.jointLaplace.copy(weights = Some(Vector(1e16))))
      .toOption
      .get
    val expected = 1e16 * math.exp(-40.0) / (1.0 + math.exp(-40.0))
    assertEqualsDouble(response.scoreResidual(0, 40.0, 1.0), expected, 1e-16)
    assert(expected > 0.04)

  test("grouped and expanded binomial have equivalent modes and assessment despite response constants"):
    val trials = Vector.fill(12)(2)
    val successes = Vector(0, 1, 1, 0, 1, 2, 2, 2, 0, 1, 1, 0)
    val proportions = successes.map(_ / 2.0)
    val grouped = ModelFrame
      .of("y" -> ModelFrame.numeric(proportions), "x" -> ModelFrame.numeric(x), "g" -> ModelFrame.factor(g))
      .toOption
      .get
    val expandedY = successes.flatMap(k => Vector.tabulate(2)(j => if j < k then 1.0 else 0.0))
    val expanded = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(expandedY),
        "x" -> ModelFrame.numeric(x.flatMap(v => Vector(v, v))),
        "g" -> ModelFrame.factor(g.flatMap(v => Vector(v, v)))
      )
      .toOption
      .get
    val options = GlmmOptions.jointLaplace.copy(trials = Some(trials))
    def inspect(data: ModelFrame, family: Family, config: GlmmOptions) =
      val design = Glmm.compile(formula, data).toOption.get
      val response = GlmmResponse.prepare(family, data.numeric("y").get, config).toOption.get
      val workspace = JointWorkspace.create(design, family, Link.Logit, response, config).toOption.get
      val parameters = Vector(-0.4, 0.3, 0.8)
      val result = workspace.evaluateParameters(parameters).toOption.get
      (
        result,
        JointLaplace.assess(
          workspace,
          parameters,
          result.objective,
          Vector(Double.NegativeInfinity, Double.NegativeInfinity, 0.0),
          config.stationarityTolerance
        )
      )
    val (a, assessmentA) = inspect(grouped, Family.Binomial, options)
    val (b, assessmentB) = inspect(expanded, Family.Bernoulli, GlmmOptions.jointLaplace)
    assertEqualsDouble(b.objective - a.objective, 2.0 * successes.count(_ == 1) * math.log(2.0), 1e-9)
    a.sphericalModes.flatten.zip(b.sphericalModes.flatten).foreach((u, v) => assertEqualsDouble(u, v, 1e-8))
    assertEquals(assessmentA.gradient.passed, assessmentB.gradient.passed)
    assertEquals(assessmentA.boundary, assessmentB.boundary)
    (assessmentA.gradient, assessmentB.gradient) match
      case (JointGradientCheck.Available(ga, _, ta), JointGradientCheck.Available(gb, _, tb)) =>
        assertEquals(ta, tb)
        ga.zip(gb).foreach((u, v) => assertEqualsDouble(u, v, 1e-7))
      case other => fail(s"missing gradient $other")

  test("a zero-theta saddle is not accepted on projected gradient alone"):
    val interceptFormula = Formula.parse("y ~ 1 + (1 | g)").toOption.get
    val design = Glmm.compile(interceptFormula, frame).toOption.get
    val options = GlmmOptions.jointLaplace
    val response = GlmmResponse.prepare(Family.Bernoulli, y, options).toOption.get
    val workspace = JointWorkspace.create(design, Family.Bernoulli, Link.Logit, response, options).toOption.get
    val parameters = Vector(0.0, 0.0)
    val result = workspace.evaluateParameters(parameters).toOption.get
    val assessment = JointLaplace.assess(
      workspace,
      parameters,
      result.objective,
      Vector(Double.NegativeInfinity, 0.0),
      options.stationarityTolerance
    )
    assert(assessment.gradient.passed, clues(assessment))
    assert(!assessment.firstOrderPassed)
    assessment.boundary match
      case JointBoundaryCheck.Descent(_, _, objective) => assert(objective < result.objective)
      case other                                       => fail(s"expected feasible boundary descent, got $other")

  test("degenerate response constants are exactly zero rather than amplified under large case weights"):
    val options =
      GlmmOptions.jointLaplace.copy(weights = Some(Vector.fill(12)(1e16)), trials = Some(Vector.fill(12)(7)))
    val ones = ModelFrame
      .of("y" -> ModelFrame.numeric(Vector.fill(12)(1.0)), "x" -> ModelFrame.numeric(x), "g" -> ModelFrame.factor(g))
      .toOption
      .get
    val binary = evaluate(Vector(1000.0, 0.0), 0.0, Family.Binomial, ones, options)
    assertEqualsDouble(binary.objective, 0.0, 0.0)
    val zeros = ModelFrame
      .of("y" -> ModelFrame.numeric(Vector.fill(12)(0.0)), "x" -> ModelFrame.numeric(x), "g" -> ModelFrame.factor(g))
      .toOption
      .get
    val count = evaluate(Vector(-1000.0, 0.0), 0.0, Family.Poisson, zeros, options.copy(trials = None))
    assertEqualsDouble(count.objective, 0.0, 0.0)
