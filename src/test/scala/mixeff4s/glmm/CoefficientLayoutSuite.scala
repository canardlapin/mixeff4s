package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.model.Family

class CoefficientLayoutSuite extends munit.FunSuite:
  test("GLMM discloses dropped columns using the same retained coefficient layout"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector.fill(3)(Vector(0.0, 1.0, 0.0, 1.0)).flatten),
        "x" -> ModelFrame.numeric(Vector.fill(3)(Vector(0.0, 0.0, 1.0, 1.0)).flatten),
        "zero" -> ModelFrame.numeric(Vector.fill(12)(0.0)),
        "g" -> ModelFrame.factor(Vector.fill(4)("a") ++ Vector.fill(4)("b") ++ Vector.fill(4)("c"))
      )
      .toOption
      .get
    val fit = Glmm.fit("y ~ 1 + zero + x + (1 | g)", frame, Family.Bernoulli).fold(err => fail(err.message), identity)
    val reference = Glmm.fit("y ~ 1 + x + (1 | g)", frame, Family.Bernoulli).fold(err => fail(err.message), identity)
    assertEquals(fit.fixedEffects.requestedNames, Vector("(Intercept)", "zero", "x"))
    assertEquals(fit.fixedEffects.retainedIndices, Vector(0, 2))
    assertEquals(fit.fixedEffects.droppedNames, Vector("zero"))
    assertEquals(fit.feNames, reference.feNames)
    assertEquals(fit.beta.length, 2)
    fit.beta.indices.foreach(i => assertEqualsDouble(fit.beta(i), reference.beta(i), 1e-10))
    assertEqualsDouble(fit.deviance, reference.deviance, 1e-10)
    assert(compileErrors("def change(fit: mixeff4s.glmm.GlmmFit) = fit.copy(beta = Vector.empty)").nonEmpty)
