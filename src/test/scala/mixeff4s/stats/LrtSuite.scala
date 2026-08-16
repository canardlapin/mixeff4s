package mixeff4s.stats

import mixeff4s.error.MixedModelError
import mixeff4s.fixtures.{Pastes, Sleepstudy}
import mixeff4s.lmm.{FitOptions, Lmm}

class LrtSuite extends munit.FunSuite:
  test("sleepstudy nested FE LRT matches MixedModels.jl"):
    val smaller =
      Lmm.fit("reaction ~ 1 + (1 + days | subj)", Sleepstudy.frame, FitOptions.ml).getOrElse(fail("fm0"))
    val larger =
      Lmm.fit("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame, FitOptions.ml).getOrElse(fail("fm1"))
    val lrt = Lrt.compare(smaller, larger).getOrElse(fail("lrt"))
    assertEquals(smaller.dof, 5)
    assertEquals(larger.dof, 6)
    assertEquals(lrt.chisqDof, 1)
    assertEqualsDouble(lrt.chisq, 23.5365, 0.05)
    assert(lrt.pvalue < 1e-5, clues(lrt.pvalue))

  test("pastes nested RE LRT matches MixedModels.jl"):
    val smaller =
      Lmm.fit("strength ~ 1 + (1 | batch:cask)", Pastes.frame, FitOptions.ml).getOrElse(fail("m1"))
    val larger =
      Lmm.fit("strength ~ 1 + (1 | batch / cask)", Pastes.frame, FitOptions.ml).getOrElse(fail("m2"))
    val lrt = Lrt.compare(smaller, larger).getOrElse(fail("lrt"))
    assertEqualsDouble(lrt.pvalue, 0.5233767965780878, 0.05)

  test("REML LRT with different fixed effects is refused"):
    val smaller =
      Lmm.fit("reaction ~ 1 + (1 + days | subj)", Sleepstudy.frame, FitOptions.reml).getOrElse(fail("fm0"))
    val larger =
      Lmm.fit("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame, FitOptions.reml).getOrElse(fail("fm1"))
    Lrt.compare(smaller, larger) match
      case Left(err: MixedModelError.InferenceUnavailable) =>
        assertEquals(err.reasonCode, "lrt_unavailable")
        assertEquals(err.code, "lrt_unavailable")
        assert(err.details.contains("REML"), clues(err.details))
      case other =>
        fail(s"expected REML refusal, got $other")

  test("mixing ML and REML is refused"):
    val ml = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml).getOrElse(fail("ml"))
    val reml = Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.reml).getOrElse(fail("reml"))
    Lrt.compare(ml, reml) match
      case Left(err: MixedModelError.InferenceUnavailable) =>
        assertEquals(err.code, "lrt_unavailable")
      case other =>
        fail(s"expected criterion refusal, got $other")
