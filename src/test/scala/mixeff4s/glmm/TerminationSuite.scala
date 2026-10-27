package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.error.MixedModelError
import mixeff4s.formula.Formula
import mixeff4s.model.Family

class TerminationSuite extends munit.FunSuite:
  private val formula = Formula.parse("y ~ 1 + (1 | g)").toOption.get
  private val balanced = ModelFrame
    .of(
      "y" -> ModelFrame.numeric(Vector(0.0, 1.0, 1.0, 0.0, 0.0, 1.0)),
      "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c"))
    )
    .toOption
    .get
  private val unequal = ModelFrame
    .of(
      "y" -> ModelFrame.numeric(Vector(0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0)),
      "g" -> ModelFrame.factor(Vector.fill(4)("a") ++ Vector.fill(4)("b") ++ Vector.fill(4)("c"))
    )
    .toOption
    .get

  test("profiled deviance preserves wrong theta dimension as a typed error"):
    Glmm.profiledDeviance(formula, balanced, Family.Bernoulli, Vector.empty) match
      case Left(err: MixedModelError.DimensionMismatch) => assertEquals(err.code, "dimension_mismatch")
      case other                                        => fail(s"expected dimension mismatch, got $other")

  test("nonfinite or negative diagonal theta is refused"):
    Vector(Double.NaN, Double.PositiveInfinity, -0.1).foreach: theta =>
      Glmm.profiledDeviance(formula, balanced, Family.Bernoulli, Vector(theta)) match
        case Left(_: MixedModelError.InvalidArgument) => ()
        case other                                    => fail(s"expected invalid argument for $theta, got $other")

  test("an unfinished PIRLS probe returns a refusal, never a finite penalty"):
    Glmm.profiledDeviance(
      formula,
      unequal,
      Family.Bernoulli,
      Vector(1.0),
      options = GlmmOptions(pirlsMaxIterations = 1)
    ) match
      case Left(err @ MixedModelError.PirlsNotConverged(iterations)) =>
        assertEquals(iterations, 1)
        assertEquals(err.code, "pirls_not_converged")
      case other => fail(s"expected PIRLS exhaustion, got $other")

  test("a public GLMM refuses unfinished PIRLS"):
    Glmm.fit(formula, unequal, Family.Bernoulli, options = GlmmOptions(pirlsMaxIterations = 1)) match
      case Left(_: MixedModelError.PirlsNotConverged) => ()
      case other                                      => fail(s"expected PIRLS exhaustion, got $other")

  test("a public GLMM refuses an exhausted optimizer budget"):
    Glmm.fit(formula, balanced, Family.Bernoulli, options = GlmmOptions(maxEvaluations = 1)) match
      case Left(err @ MixedModelError.OptimizerBudgetExhausted(evaluations)) =>
        assertEquals(evaluations, 1)
        assertEquals(err.code, "optimizer_budget_exhausted")
      case other => fail(s"expected optimizer exhaustion, got $other")

  test("an accepted GLMM records the final inner and outer termination separately"):
    val fit = Glmm.fit(formula, balanced, Family.Bernoulli).fold(err => fail(err.message), identity)
    assertEquals(fit.optimizer.x, fit.theta)
    assert(fit.optimizer.requireConvergence.isRight)
    fit.inner match
      case GlmmInnerResult.FastPirls(result) =>
        assert(result.iterations >= 1 && result.iterations <= Pirls.MaxIter)
        assert(result.objectiveChange < result.tolerance)
      case other => fail(s"wrong inner criterion $other")
    // Equal successes/failures in every group: beta=0, theta=0, deviance=2n log(2).
    assertEqualsDouble(fit.beta.head, 0.0, 1e-10)
    assertEqualsDouble(fit.theta.head, 0.0, 1e-5)
    assertEqualsDouble(fit.deviance, 12.0 * math.log(2.0), 1e-8)
