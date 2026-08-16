package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.error.MixedModelError
import mixeff4s.formula.Formula
import mixeff4s.model.{Family, Link}

class GlmmSuite extends munit.FunSuite:
  private val frame = ModelFrame
    .of(
      "y" -> ModelFrame.numeric(Vector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)),
      "x" -> ModelFrame.numeric(Vector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)),
      "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c"))
    )
    .getOrElse(fail("frame"))

  test("Normal+Identity is refused as a GLMM"):
    Glmm.fit("y ~ 1 + x + (1 | g)", frame, Family.Normal, Some(Link.Identity), GlmmOptions.fastLaplace) match
      case Left(err: MixedModelError.UnsupportedFamilyLink) =>
        assertEquals(err.family, "normal")
        assertEquals(err.link, "identity")
        assertEquals(err.code, "unsupported_family_link")
      case other =>
        fail(s"expected UnsupportedFamilyLink, got $other")

  test("fast=false is refused rather than sold as glmer"):
    Glmm.fit(
      Formula.parse("y ~ 1 + x + (1 | g)").getOrElse(fail("parse")),
      frame,
      Family.Bernoulli,
      None,
      GlmmOptions(fast = false)
    ) match
      case Left(err: MixedModelError.Unsupported) =>
        assert(err.details.contains("fast-PIRLS"), clues(err.details))
        assert(!err.details.contains("glmer"), clues(err.details))
      case other =>
        fail(s"expected Unsupported, got $other")

  test("nAGQ > 1 is refused in this slice"):
    Glmm.fit(
      Formula.parse("y ~ 1 + x + (1 | g)").getOrElse(fail("parse")),
      frame,
      Family.Bernoulli,
      None,
      GlmmOptions(nAgq = 7)
    ) match
      case Left(err: MixedModelError.Unsupported) =>
        assert(err.details.contains("fast-PIRLS"), clues(err.details))
      case other =>
        fail(s"expected Unsupported, got $other")

  test("a small Bernoulli RE intercept is labelled fast-PIRLS"):
    Glmm.fit("y ~ 1 + x + (1 | g)", frame, Family.Bernoulli) match
      case Right(fit) =>
        assertEquals(fit.approximation, Approximation.FastPirls)
        assertEquals(fit.algorithmLabel, "fast-PIRLS")
        assertEquals(fit.family, Family.Bernoulli)
        assertEquals(fit.link, Link.Logit)
        assertEquals(fit.theta.length, 1)
        assertEquals(fit.beta.length, 2)
        assert(fit.deviance.isFinite, clues(fit.deviance))
      case other =>
        fail(s"expected a labelled fast-PIRLS fit, got $other")
