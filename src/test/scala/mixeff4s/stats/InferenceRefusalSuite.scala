package mixeff4s.stats

import mixeff4s.error.MixedModelError
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.glmm.Glmm
import mixeff4s.lmm.{FitOptions, Lmm}
import mixeff4s.model.Family

class InferenceRefusalSuite extends munit.FunSuite:
  test("profile intervals are a typed refusal"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml).getOrElse(fail("fit"))
    Profile.confint(fit) match
      case Left(err: MixedModelError.InferenceUnavailable) =>
        assertEquals(err.code, "profile_unavailable")
      case other =>
        fail(s"expected profile refusal, got $other")

  test("parametric bootstrap is a typed refusal"):
    val fit = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml).getOrElse(fail("fit"))
    Bootstrap.parametric(fit, 100) match
      case Left(err: MixedModelError.InferenceUnavailable) =>
        assertEquals(err.code, "bootstrap_unavailable")
      case other =>
        fail(s"expected bootstrap refusal, got $other")

  test("GLMM profile intervals stay refused"):
    val frame = mixeff4s.data.ModelFrame
      .of(
        "y" -> mixeff4s.data.ModelFrame.numeric(Vector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)),
        "g" -> mixeff4s.data.ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c"))
      )
      .getOrElse(fail("frame"))
    val fit = Glmm.fit("y ~ 1 + (1 | g)", frame, Family.Bernoulli).getOrElse(fail("glmm"))
    Profile.confint(fit) match
      case Left(err: MixedModelError.InferenceUnavailable) =>
        assertEquals(err.code, "profile_unavailable")
        assert(err.details.contains("GLMM"), clues(err.details))
      case other =>
        fail(s"expected GLMM profile refusal, got $other")
