package mixeff4s.glmm

import mixeff4s.fixtures.GlmmFrames
import mixeff4s.formula.Formula
import mixeff4s.lmm.SimulationMode
import mixeff4s.model.Family

class GlmmSimulationSuite extends munit.FunSuite:
  private def right[A](value: Either[?, A]): A = value.fold(error => fail(error.toString), identity)
  private lazy val fit = right(
    Glmm.fit(
      Formula.parse("count ~ 1 + x + (1 | group)").toOption.get,
      GlmmFrames.synthetic_poisson_offset,
      Family.Poisson,
      options = GlmmOptions.fastLaplace.copy(
        offset = Some(GlmmFrames.synthetic_poisson_offset.numeric("exposure").get.map(math.log)),
        maxEvaluations = 10000
      )
    )
  )

  test("GLMM simulations are seeded, integral for Poisson, and condition modes without residual sigma"):
    val a = right(GlmmSimulation.response(fit, 991L))
    val b = right(GlmmSimulation.response(fit, 991L))
    val conditional = right(GlmmSimulation.response(fit, 991L, SimulationMode.ConditionOnModes))
    assertEquals(a, b)
    assert(a.forall(value => value >= 0.0 && value == math.rint(value)))
    assertEquals(conditional.length, fit.n)

  test("power-likelihood fits refuse response simulation"):
    val weighted = right(
      Glmm.fit(
        Formula.parse("count ~ 1 + x + (1 | group)").toOption.get,
        GlmmFrames.synthetic_poisson_offset,
        Family.Poisson,
        options = GlmmOptions.fastLaplace.copy(
          weights = Some(Vector.fill(GlmmFrames.synthetic_poisson_offset.nRows)(2.0)),
          offset = Some(GlmmFrames.synthetic_poisson_offset.numeric("exposure").get.map(math.log)),
          maxEvaluations = 10000
        )
      )
    )
    assertEquals(GlmmSimulation.response(weighted, 2L).left.toOption.get.code, "simulation_case_weights")
