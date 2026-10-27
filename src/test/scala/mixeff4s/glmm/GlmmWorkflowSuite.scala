package mixeff4s.glmm

import mixeff4s.examples.GlmmWorkflow

class GlmmWorkflowSuite extends munit.FunSuite:
  test("documented grouped-binomial workflow runs with qualified inner and outer evidence"):
    val fit = GlmmWorkflow.run().fold(e => fail(e.message), identity)
    assertEquals(fit.approximation, Approximation.JointLaplace)
    assert(fit.jointAssessment.exists(_.firstOrderPassed), clues(fit.jointAssessment))
    assertEquals(fit.randomEffects.head.levels.length, 15)
    assertEquals(fit.conditionalMean.length, 56)
    assertEqualsDouble(fit.deviance, 184.05256374299438, 1e-4)
