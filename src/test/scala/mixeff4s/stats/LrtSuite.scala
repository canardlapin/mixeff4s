package mixeff4s.stats

import mixeff4s.data.ModelFrame
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.fixtures.{Pastes, Sleepstudy}
import mixeff4s.lmm.{FitOptions, Lmm, LmmFit}

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

  test("pastes preserves the LR statistic but refuses unqualified variance-component calibration"):
    val smaller =
      Lmm.fit("strength ~ 1 + (1 | batch:cask)", Pastes.frame, FitOptions.ml).getOrElse(fail("m1"))
    val larger =
      Lmm.fit("strength ~ 1 + (1 | batch / cask)", Pastes.frame, FitOptions.ml).getOrElse(fail("m2"))
    val lrt = Lrt.statistic(smaller, larger).fold(err => fail(err.message), identity)
    // The old reference p=0.5233767965780878 used an ordinary chi-square(1).
    // Preserve that numerical witness without endorsing it as a valid variance-component law.
    assertEqualsDouble(ChiSquare.sf(lrt.likelihoodRatio, 1), 0.5233767965780878, 0.05)
    assertEquals(lrt.parameterDifference, 1)
    refused(Lrt.compare(smaller, larger), "lrt_reference_unavailable")

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

  private def fit(formula: String, data: ModelFrame = Sleepstudy.frame, options: FitOptions = FitOptions.ml): LmmFit =
    Lmm.fit(formula, data, options).fold(err => fail(err.message), identity)

  private def refused[A](result: FitResult[A], code: String = "lrt_unavailable", detail: String = ""): Unit =
    result match
      case Left(err: MixedModelError.InferenceUnavailable) =>
        assertEquals(err.code, code)
        assert(err.details.contains(detail), clues(err.details))
      case other => fail(s"expected $code, got $other")

  private def frame(days: Vector[Double] = Sleepstudy.days, groups: Vector[String] = Sleepstudy.subj): ModelFrame =
    ModelFrame
      .of(
        "reaction" -> ModelFrame.numeric(Sleepstudy.reaction),
        "days" -> ModelFrame.numeric(days),
        "cube" -> ModelFrame.numeric(Sleepstudy.days.map(x => x * x * x)),
        "subj" -> ModelFrame.factor(groups)
      )
      .fold(err => fail(err.message), identity)

  private lazy val intercept = fit("reaction ~ 1 + days + (1 | subj)")
  private lazy val slope = fit("reaction ~ 1 + days + (1 + days | subj)")

  test("same-named changed predictors cannot certify fixed-effect nesting"):
    val small = fit("reaction ~ 1 + days + (1 | subj)", frame(days = Sleepstudy.days.map(x => x * x)))
    val large = fit("reaction ~ 1 + days + cube + (1 | subj)", frame())
    assert(large.loglikelihood > small.loglikelihood)
    refused(Lrt.compare(small, large), detail = "fixed-effect nesting")

  test("random intercept is nested in intercept plus slope, without a fabricated p-value"):
    val result = Lrt.statistic(intercept, slope).fold(err => fail(err.message), identity)
    assertEquals(result.parameterDifference, 2)
    assert(result.likelihoodRatio > 0)
    refused(Lrt.compare(intercept, slope), "lrt_reference_unavailable")

  test("same-named changed grouping partition cannot certify nesting"):
    val changed = fit("reaction ~ 1 + (1 | subj)", frame(groups = Sleepstudy.subj.reverse))
    // Reverse labels alone preserves this balanced partition; interleave observations instead.
    val interleaved = fit("reaction ~ 1 + (1 | subj)", frame(groups = Vector.tabulate(180)(i => s"g${i % 18}")))
    assert(Lrt.statistic(changed, intercept).isRight)
    refused(Lrt.statistic(interleaved, intercept), detail = "random-effect nesting")

  test("renamed predictors and relabeled reordered group levels preserve numerical design"):
    val relabeled = Sleepstudy.subj.map(s => s"renamed$s")
    val data = ModelFrame
      .of(
        "reaction" -> ModelFrame.numeric(Sleepstudy.reaction),
        "x" -> ModelFrame.numeric(Sleepstudy.days),
        "g" -> ModelFrame.factor(relabeled, relabeled.distinct.reverse).toOption.get
      )
      .toOption
      .get
    val small = fit("reaction ~ 1 + (1 + x | g)", data)
    val result = Lrt.compare(small, slope).fold(err => fail(err.message), identity)
    assertEqualsDouble(result.chisq, 23.5365, 0.05)

  test("diagonal covariance is nested in full covariance but calibration remains unavailable"):
    val diagonal = fit("reaction ~ 1 + days + (1 + days || subj)")
    assert(Lrt.statistic(diagonal, slope).isRight)
    refused(Lrt.compare(diagonal, slope), "lrt_reference_unavailable")
    val largerDiagonal = fit("reaction ~ 1 + days + cube + (1 + days || subj)", frame())
    refused(Lrt.statistic(slope, largerDiagonal), detail = "random-effect nesting")

  test("changed random slope values refuse even when names and partitions match"):
    val small = fit("reaction ~ 1 + (0 + days | subj)", frame(days = Sleepstudy.days.map(x => x * x)))
    refused(Lrt.statistic(small, slope), detail = "random-effect nesting")

  test("group-name inclusion is not grouping-partition inclusion"):
    val small = fit("strength ~ 1 + (1 | batch)", Pastes.frame)
    val large = fit("strength ~ 1 + (1 | batch:cask)", Pastes.frame)
    refused(Lrt.statistic(small, large), detail = "random-effect nesting")

  test("REML permits the same fixed basis in a different order for a covariance statistic"):
    val data = Sleepstudy.frame.addNumeric("one", Vector.fill(180)(1.0)).toOption.get
    val small = fit("reaction ~ 0 + one + days + (1 | subj)", data, FitOptions.reml)
    val large = fit("reaction ~ 0 + days + one + (1 + days | subj)", data, FitOptions.reml)
    assert(Lrt.statistic(small, large).isRight)
    refused(Lrt.compare(small, large), "lrt_reference_unavailable")

  test("REML refuses the same named fixed effects with changed values or scaling"):
    val small = fit("reaction ~ 1 + days + (1 | subj)", options = FitOptions.reml)
    for values <- Vector(Sleepstudy.days.map(_ * 2.0), Sleepstudy.days.map(x => x * x)) do
      val large = fit("reaction ~ 1 + days + (1 + days | subj)", frame(days = values), FitOptions.reml)
      refused(Lrt.statistic(small, large), detail = "REML")

  test("rank-dropped zero columns do not change parameter difference or the LR statistic"):
    val data = Sleepstudy.frame.addNumeric("zero", Vector.fill(180)(0.0)).toOption.get
    val small = fit("reaction ~ 1 + zero + (1 + days | subj)", data)
    val large = fit("reaction ~ 1 + zero + days + (1 + days | subj)", data)
    val result = Lrt.compare(small, large).fold(err => fail(err.message), identity)
    assertEquals(result.chisqDof, 1)
    assertEqualsDouble(result.chisq, 23.5365, 0.05)

  test("unbound or modified fits cannot supply comparison evidence"):
    refused(Lrt.statistic(intercept.copy(designEvidence = None), slope), detail = "evidence")
    refused(Lrt.statistic(intercept.copy(y = intercept.y.reverse), slope), detail = "evidence")
    refused(Lrt.statistic(intercept.copy(objective = Double.NaN), slope), detail = "evidence")
    refused(
      Lrt.statistic(intercept.copy(optimizer = intercept.optimizer.copy(fmin = Double.NaN)), slope),
      detail = "evidence"
    )

  test("different responses and non-increasing parameter counts refuse"):
    val data = ModelFrame
      .of(
        "reaction" -> ModelFrame.numeric(Sleepstudy.reaction.map(_ + 1.0)),
        "days" -> ModelFrame.numeric(Sleepstudy.days),
        "subj" -> ModelFrame.factor(Sleepstudy.subj)
      )
      .toOption
      .get
    val changed = fit("reaction ~ 1 + days + (1 + days | subj)", data)
    refused(Lrt.statistic(intercept, changed), detail = "response")
    refused(Lrt.statistic(intercept, intercept), detail = "more parameters")

  test("an actual fitted zero variance permits a statistic but refuses ordinary calibration"):
    val x = Vector.tabulate(80)(i => (i % 8).toDouble)
    val errors = Vector.fill(10)(Vector(1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0)).flatten
    val data = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(x.zip(errors).map((v, e) => 3.0 + 2.0 * v + e)),
        "x" -> ModelFrame.numeric(x),
        "g" -> ModelFrame.factor(Vector.tabulate(80)(i => s"g${i / 8}"))
      )
      .toOption
      .get
    val small = fit("y ~ 1 + (1 | g)", data)
    val large = fit("y ~ 1 + x + (1 | g)", data)
    assert(large.theta.head <= 1e-8)
    assert(Lrt.statistic(small, large).isRight)
    refused(Lrt.compare(small, large), "lrt_reference_unavailable", "boundary")

  test("two smaller random blocks cannot reuse one larger block"):
    val small = fit("reaction ~ 1 + days + (1 | subj) + (1 | subj)")
    assertEquals(small.designEvidence.get.randomBlocks.length, 2)
    refused(Lrt.statistic(small, slope), detail = "random-effect nesting")

  test("random-block matching finds an embedding that a greedy assignment would miss"):
    val small = fit("reaction ~ 1 + days + (1 | subj) + (0 + days | subj)")
    val large = fit("reaction ~ 1 + days + (1 + days | subj) + (1 | subj)")
    assertEquals(small.designEvidence.get.randomBlocks.map(_.columns.length), Vector(1, 1))
    assertEquals(large.designEvidence.get.randomBlocks.map(_.columns.length), Vector(2, 1))
    assert(Lrt.statistic(small, large).isRight)

  test("equivalent scaled fixed bases are conservatively refused rather than declared nonnested"):
    val small = fit("reaction ~ 1 + days + (1 | subj)", frame(days = Sleepstudy.days.map(_ * 2.0)))
    val large = fit("reaction ~ 1 + days + cube + (1 | subj)", frame())
    refused(Lrt.statistic(small, large), detail = "not certified by exact fitted-column inclusion")
