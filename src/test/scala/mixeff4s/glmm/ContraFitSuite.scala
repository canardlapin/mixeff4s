package mixeff4s.glmm

import mixeff4s.fixtures.Contra
import mixeff4s.formula.Formula
import mixeff4s.model.Family

class ContraFitSuite extends munit.FunSuite:
  private val formula =
    Formula.parse("use_num ~ 1 + age + age2 + urban + livch + (1 | urban_dist)").getOrElse(fail("parse"))

  test("contra Laplace deviance at the Julia theta matches MixedModels.jl"):
    val juliaTheta = Vector(0.5720746212924732)
    val dev =
      Glmm.profiledDeviance(formula, Contra.frame, Family.Bernoulli, juliaTheta).getOrElse(fail("deviance"))
    assertEqualsDouble(dev, 2361.657202855648, 1.0)

  test("contra fast-PIRLS fit matches MixedModels.jl"):
    val fit = Glmm.fit(formula, Contra.frame, Family.Bernoulli).getOrElse(fail("fit"))
    assertEquals(fit.approximation, Approximation.FastPirls)
    assertEquals(fit.algorithmLabel, "fast-PIRLS")
    assertEquals(fit.theta.length, 1)
    assertEqualsDouble(fit.theta(0), 0.5720746212924732, 0.01)
    assertEqualsDouble(fit.deviance, 2361.657202855648, 1.0)
