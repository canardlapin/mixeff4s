package mixeff4s.lmm

import mixeff4s.data.ModelFrame
import mixeff4s.error.MixedModelError
import mixeff4s.formula.Formula

class LmmSuite extends munit.FunSuite:
  private val frame = ModelFrame
    .of(
      "y" -> ModelFrame.numeric(Vector(1.0, 2.1, 3.0, 4.2, 5.1, 6.0)),
      "x" -> ModelFrame.numeric(Vector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)),
      "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c"))
    )
    .getOrElse(fail("frame"))

  test("refuses a formula with no random effects"):
    val formula = Formula.parse("y ~ 1 + x").getOrElse(fail("parse"))
    Lmm.fit(formula, frame, FitOptions.reml) match
      case Left(MixedModelError.NoRandomEffects) => ()
      case other                                 => fail(s"expected NoRandomEffects, got $other")

  test("parses then refuses structured covariance"):
    Lmm.fit("y ~ x + ar1(1 | g)", frame, FitOptions.reml) match
      case Left(err: MixedModelError.Unsupported) =>
        assert(err.details.contains("ar1"), clues(err.details))
        assertEquals(err.code, "unsupported")
      case other =>
        fail(s"expected Unsupported, got $other")

  test("single-term random intercept fits"):
    Lmm.fit("y ~ 1 + x + (1 | g)", frame, FitOptions.reml) match
      case Right(fit) =>
        assertEquals(fit.beta.length, 2)
        assertEquals(fit.theta.length, 1)
        assert(fit.sigma > 0.0, clues(fit.sigma))
        assert(fit.objective.isFinite, clues(fit.objective))
      case other =>
        fail(s"expected a certified fit, got $other")

  test("zerocorr fits two diagonal theta slots and no correlation"):
    val design = Lmm.compile("y ~ 1 + x + (1 + x || g)", frame).getOrElse(fail("compile"))
    assertEquals(design.nTheta, 2)
    assertEquals(design.parmap, Vector((0, 0, 0), (0, 1, 1)))
    Lmm.fit("y ~ 1 + x + (1 + x || g)", frame, FitOptions.reml) match
      case Right(fit) =>
        assertEquals(fit.theta.length, 2)
        assert(fit.theta.forall(_ >= 0.0), clues(fit.theta))
        assertEquals(fit.varcorr.components.head.correlations, Vector(0.0))
        assert(fit.objective.isFinite, clues(fit.objective))
      case other =>
        fail(s"expected a certified zerocorr fit, got $other")

  test("parses then refuses compound symmetry"):
    Lmm.fit("y ~ x + cs(1 + x | g)", frame, FitOptions.reml) match
      case Left(err: MixedModelError.Unsupported) =>
        assert(err.details.contains("cs"), clues(err.details))
      case other =>
        fail(s"expected Unsupported, got $other")
