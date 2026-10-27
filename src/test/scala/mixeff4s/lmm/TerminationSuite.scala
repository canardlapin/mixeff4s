package mixeff4s.lmm

import gale.linalg.Matrix
import mixeff4s.design.FeMat
import mixeff4s.error.MixedModelError
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.optimizer.TrustBqStopReason

class TerminationSuite extends munit.FunSuite:
  test("a public LMM refuses an exhausted main optimizer budget"):
    Lmm.fit("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame, FitOptions.ml.copy(maxEvaluations = 1)) match
      case Left(err @ MixedModelError.OptimizerBudgetExhausted(evaluations)) =>
        assertEquals(evaluations, 1)
        assertEquals(err.code, "optimizer_budget_exhausted")
      case other => fail(s"expected exhausted budget, got $other")

  test("an accepted fit retains evidence for its actual theta and objective"):
    val fit = Lmm
      .fit("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame, FitOptions.ml)
      .fold(err => fail(err.message), identity)
    assertEquals(fit.optimizer.x, fit.theta)
    assertEqualsDouble(fit.optimizer.fmin, fit.objective, 1e-10)
    assert(fit.optimizer.fevals > 1)
    assert(fit.optimizer.requireConvergence.isRight)
    assertNotEquals(fit.optimizer.stopReason, TrustBqStopReason.MaxEvaluations)

  test("a zero fixed-effect pivot cannot become a finite REML objective"):
    val original = Lmm.compile("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame).toOption.get
    // Inject an inconsistent compiled matrix to exercise the factor guard independently of QR truncation.
    val broken = original.copy(xy = FeMat(Matrix.tabulate(original.n, original.p + 1): (i, j) =>
      if j == 1 then 0.0 else original.xy.xy(i, j)))
    val workspace = PlsWorkspace(broken, reml = true).toOption.get
    workspace.objectiveAt(broken.theta) match
      case Left(_: MixedModelError.Optimization) => ()
      case other                                 => fail(s"expected invalid objective diagnostic, got $other")
    assert(!workspace.objective.isFinite)
