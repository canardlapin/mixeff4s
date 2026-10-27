package mixeff4s.lmm

import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.error.FitResult
import mixeff4s.fixtures.Sleepstudy

class PredictionSuite extends munit.FunSuite:
  private def right[A](value: FitResult[A]): A = value.fold(e => fail(s"${e.code}: ${e.message}"), identity)
  private val n = 48
  private val frame = right(
    ModelFrame.of(
      "y" -> ModelFrame.numeric(Vector.tabulate(n)(i => 2.0 + 0.3 * (i % 6) + 0.4 * (i / 6) + math.sin(i * 1.4))),
      "x" -> ModelFrame.numeric(Vector.tabulate(n)(i => (i % 6 - 2).toDouble)),
      "f" -> right(ModelFrame.factor(Vector.tabulate(n)(i => Vector("a", "b", "c")(i % 3)), Vector("a", "b", "c"))),
      "g" -> ModelFrame.factor(Vector.tabulate(n)(i => s"g${i / 6}"))
    )
  )
  private lazy val fit = right(Lmm.fit("y ~ 1 + f * x + I(x^2) + (1 | g)", frame, FitOptions.ml))

  private def subset(indices: Vector[Int], names: Set[String] = frame.names.toSet): ModelFrame =
    val columns = frame.columns
      .filter(e => names(e._1))
      .map:
        case (name, Column.Numeric(values)) => name -> ModelFrame.numeric(indices.map(values))
        case (name, factor: Column.Factor)  => name -> ModelFrame.factor(indices.map(factor.values))
    right(ModelFrame.of(columns*))

  test("new data reuses training factor levels, transformations and row order"):
    val indices = Vector(20, 1, 4, 0)
    val predicted = right(Prediction.predict(fit, subset(indices)))
    val training = right(Prediction.fitted(fit))
    predicted.zip(indices).foreach((value, row) => assertEqualsDouble(value, training(row), 1e-9))
    val population = right(Prediction.predict(fit, subset(indices, Set("f", "x")), PredictionMode.Population))
    val trainPopulation = right(Prediction.fitted(fit, PredictionMode.Population))
    population.zip(indices).foreach((value, row) => assertEqualsDouble(value, trainPopulation(row), 1e-9))
    val residuals = right(Prediction.residuals(fit))
    residuals.zip(training).zip(fit.y).foreach { case ((e, mu), y) => assertEqualsDouble(e + mu, y, 1e-12) }

  test("unknown fixed levels and changed predictor types refuse"):
    val unknown = right(ModelFrame.of("f" -> ModelFrame.factor(Vector("new")), "x" -> ModelFrame.numeric(Vector(1.0))))
    assertEquals(
      Prediction.predict(fit, unknown, PredictionMode.Population).left.toOption.get.code,
      "prediction_unknown_fixed_level"
    )
    val wrong = right(ModelFrame.of("f" -> ModelFrame.numeric(Vector(1.0)), "x" -> ModelFrame.numeric(Vector(1.0))))
    assertEquals(
      Prediction.predict(fit, wrong, PredictionMode.Population).left.toOption.get.code,
      "prediction_column_type"
    )
    assert(
      Prediction
        .predict(fit, right(ModelFrame.of("f" -> ModelFrame.factor(Vector("a")))), PredictionMode.Population)
        .isLeft
    )

  test("new grouping levels require explicit zero-effect policy"):
    val data = right(
      ModelFrame.of(
        "f" -> ModelFrame.factor(Vector("b")),
        "x" -> ModelFrame.numeric(Vector(2.0)),
        "g" -> ModelFrame.factor(Vector("new"))
      )
    )
    assertEquals(Prediction.predict(fit, data).left.toOption.get.code, "prediction_new_group")
    val conditional = right(Prediction.predict(fit, data, newGroups = NewGroupPolicy.ZeroEffect))
    val population = right(Prediction.predict(fit, data, PredictionMode.Population))
    assertEquals(conditional, population)

  test("rank-dropped columns require estimable new rows"):
    val base = Sleepstudy.frame
    val withAlias = right(ModelFrame.of((base.columns :+ ("twice" -> ModelFrame.numeric(Sleepstudy.days.map(_ * 2))))*))
    val aliased = right(Lmm.fit("reaction ~ 1 + days + twice + (1 | subj)", withAlias, FitOptions.ml))
    val aligned =
      right(ModelFrame.of("days" -> ModelFrame.numeric(Vector(2.0)), "twice" -> ModelFrame.numeric(Vector(4.0))))
    assert(Prediction.predict(aliased, aligned, PredictionMode.Population).isRight)
    val off =
      right(ModelFrame.of("days" -> ModelFrame.numeric(Vector(2.0)), "twice" -> ModelFrame.numeric(Vector(5.0))))
    assertEquals(
      Prediction.predict(aliased, off, PredictionMode.Population).left.toOption.get.code,
      "prediction_nonestimable"
    )

  test("complete snapshot seal rejects changed public fit fields"):
    val altered = Vector(
      fit.copy(beta = fit.beta.updated(0, fit.beta.head + 1)),
      fit.copy(sigma = fit.sigma + 1),
      fit.copy(vcov = fit.vcov.updated(0, fit.vcov.head.updated(0, fit.vcov.head.head + 1))),
      fit.copy(theta = fit.theta.updated(0, fit.theta.head + 0.1)),
      fit.copy(objective = fit.objective + 1),
      fit.copy(options = fit.options.copy(maxEvaluations = 2)),
      fit.copy(fixedEffects =
        fit.fixedEffects.copy(requestedNames = fit.fixedEffects.requestedNames.map(_ + "changed"))
      ),
      fit.copy(optimizer = fit.optimizer.copy(projectedGradientNorm = Some(0.0))),
      fit.copy(y = fit.y.updated(0, fit.y.head + 1)),
      fit.copy(n = fit.n + 1),
      fit.copy(designEvidence = None),
      fit.copy(workflow = None)
    )
    altered.foreach(changed => assertEquals(Prediction.fitted(changed).left.toOption.get.code, "fit_snapshot_mismatch"))
    assert(Prediction.fitted(fit).isRight)

  test("overflow in a dropped interaction cannot masquerade as an estimable row"):
    val training = right(
      ModelFrame.of(
        "y" -> ModelFrame.numeric(Sleepstudy.reaction),
        "x" -> ModelFrame.numeric(Sleepstudy.days),
        "z" -> ModelFrame.numeric(Vector.fill(180)(0.0)),
        "g" -> ModelFrame.factor(Sleepstudy.subj)
      )
    )
    val zero = right(Lmm.fit("y ~ 1 + x:z + (1|g)", training, FitOptions.ml))
    val huge = right(ModelFrame.of("x" -> ModelFrame.numeric(Vector(1e200)), "z" -> ModelFrame.numeric(Vector(1e200))))
    assertEquals(
      Prediction.predict(zero, huge, PredictionMode.Population).left.toOption.get.code,
      "prediction_nonfinite_design"
    )

  test("alias estimability is invariant to the units of a dropped zero column"):
    val training = right(
      ModelFrame.of(
        "reaction" -> ModelFrame.numeric(Sleepstudy.reaction),
        "days" -> ModelFrame.numeric(Sleepstudy.days),
        "zero" -> ModelFrame.numeric(Vector.fill(180)(0.0)),
        "subj" -> ModelFrame.factor(Sleepstudy.subj)
      )
    )
    val fit = right(Lmm.fit("reaction ~ 1 + days + zero + (1|subj)", training, FitOptions.ml))
    val data =
      right(ModelFrame.of("days" -> ModelFrame.numeric(Vector(2.0)), "zero" -> ModelFrame.numeric(Vector(1e-10))))
    assertEquals(
      Prediction.predict(fit, data, PredictionMode.Population).left.toOption.get.code,
      "prediction_nonestimable"
    )

  test("later response refits leave earlier modes and predictions unchanged"):
    val plan = right(Lmm.prepare("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame))
    val session = right(plan.session(FitOptions.ml))
    val first = right(session.refit(Sleepstudy.reaction))
    val before = right(Prediction.fitted(first))
    right(session.refit(Sleepstudy.reaction.zip(Sleepstudy.days).map((y, x) => y + 0.1 * x * x)))
    assertEquals(right(Prediction.fitted(first)), before)

  test("copied compiled designs cannot substitute another prediction recipe"):
    val ordinary = right(Lmm.compile("reaction ~ 1 + days + (1 | subj)", Sleepstudy.frame))
    val squared = right(Lmm.compile("reaction ~ 1 + I(days^2) + (1 | subj)", Sleepstudy.frame))
    assert(Lmm.fit(ordinary.copy(prediction = squared.prediction), FitOptions.ml).isLeft)

  test("chained derived recipes preserve dependencies and reject contradictory supplied helpers"):
    import mixeff4s.formula.{BinOp, DerivedColumn, Expr, FixedTerm, Formula}
    val base = Formula.parse("reaction ~ 1 + days + (1 | subj)").toOption.get
    val chained = base.copy(
      fixedTerms = Vector(FixedTerm.Intercept, FixedTerm.Column("derived2")),
      derived = Vector(
        DerivedColumn("derived1", Expr.Bin(BinOp.Pow, Expr.Col("days"), Expr.Lit(2.0))),
        DerivedColumn("derived2", Expr.Bin(BinOp.Add, Expr.Col("derived1"), Expr.Lit(1.0)))
      )
    )
    val fit = right(Lmm.fit(chained, Sleepstudy.frame, FitOptions.ml))
    val data = right(ModelFrame.of("days" -> ModelFrame.numeric(Vector(2.0))))
    val predicted = right(Prediction.predict(fit, data, PredictionMode.Population))
    assertEqualsDouble(predicted.head, fit.beta.head + 5.0 * fit.beta(1), 1e-9)
    val contradictory =
      right(ModelFrame.of("days" -> ModelFrame.numeric(Vector(2.0)), "derived1" -> ModelFrame.numeric(Vector(99.0))))
    assert(Prediction.predict(fit, contradictory, PredictionMode.Population).isLeft)
