package mixeff4s.lmm

import mixeff4s.data.ModelFrame
import mixeff4s.design.FixedEffectLayout
import mixeff4s.fixtures.Sleepstudy

class CoefficientLayoutSuite extends munit.FunSuite:
  private val frame = Sleepstudy.frame
    .addNumeric("zero", Vector.fill(180)(0.0))
    .toOption
    .get
    .addNumeric("one", Vector.fill(180)(1.0))
    .toOption
    .get
    .addNumeric("doubleDays", Sleepstudy.days.map(_ * 2.0))
    .toOption
    .get

  private def fitted(columns: String, options: FitOptions): LmmFit =
    Lmm.fit(s"reaction ~ 0 + $columns + (1 | subj)", frame, options).fold(err => fail(err.message), identity)

  private lazy val reducedMl = fitted("one + days", FitOptions.ml)
  private lazy val reducedReml = fitted("one + days", FitOptions.reml)

  for
    (columns, droppedIndex) <- Vector(("zero + one + days", 0), ("one + zero + days", 1), ("one + days + zero", 2))
    options <- Vector(FitOptions.ml, FitOptions.reml)
  do
    test(s"${options.criterion}: a dropped column at index $droppedIndex preserves all coefficient axes"):
      val fit = fitted(columns, options)
      val reference = if options.criterion == Criterion.ML then reducedMl else reducedReml
      assertEquals(fit.fixedEffects.requestedNames, columns.split(" \\+ ").toVector)
      assertEquals(fit.fixedEffects.retainedIndices, (0 until 3).filterNot(_ == droppedIndex).toVector)
      assertEquals(fit.fixedEffects.droppedIndices, Vector(droppedIndex))
      assertEquals(fit.fixedEffects.droppedNames, Vector("zero"))
      assertEquals(fit.feNames, Vector("one", "days"))
      assertAligned(fit)
      assertEquivalent(fit, reference)

  test("original middle-column counterexample has a finite days standard error and z"):
    val fit = Lmm
      .fit("reaction ~ 1 + zero + days + (1 | subj)", frame, FitOptions.ml)
      .fold(err => fail(err.message), identity)
    assertEquals(fit.feNames, Vector("(Intercept)", "days"))
    assertEquals(fit.fixedEffects.droppedIndices, Vector(1))
    assertAligned(fit)
    assertEqualsDouble(fit.stderror(1), reducedMl.stderror(1), 1e-8)
    assert(fit.coefTable.zValues(1).isFinite)

  test("balanced random-intercept covariance agrees with the analytic GLS calculation"):
    val fit = fitted("one + zero + days", FitOptions.ml)
    val sigmaSq = fit.sigma * fit.sigma
    val slopeVariance = sigmaSq / (18.0 * 82.5)
    val interceptVariance = sigmaSq * (fit.theta.head * fit.theta.head / 18.0 + 1.0 / 180.0) +
      4.5 * 4.5 * slopeVariance
    val expected = Vector(Vector(interceptVariance, -4.5 * slopeVariance), Vector(-4.5 * slopeVariance, slopeVariance))
    assertMatrix(fit.vcov, expected, 1e-8)

  test("a nonzero scaled alias at theta zero agrees with independent OLS fractions"):
    val x = (1 to 8).map(_.toDouble).toVector
    val residual = Vector(1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0)
    val data = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(x.zip(residual).map((v, e) => 3.0 + 2.0 * v + e)),
        "x" -> ModelFrame.numeric(x),
        "z" -> ModelFrame.numeric(x.map(_ * 10.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c", "d", "d"))
      )
      .toOption
      .get
    val design = Lmm.compile("y ~ 1 + x + z + (1 | g)", data).toOption.get
    assertEquals(design.fe.fixedEffects.retainedIndices, Vector(0, 1))
    assertEquals(design.fe.fixedEffects.droppedNames, Vector("z"))
    val workspace = PlsWorkspace(design, reml = true).toOption.get
    workspace.objectiveAt(Vector(0.0)).fold(err => fail(err.message), identity)
    assertEquals(workspace.feNames, Vector("(Intercept)", "x"))
    assertEqualsDouble(workspace.beta(0), 3.0, 1e-10)
    assertEqualsDouble(workspace.beta(1), 2.0, 1e-10)
    assertEqualsDouble(workspace.sigma * workspace.sigma, 4.0 / 3.0, 1e-10)
    assertMatrix(workspace.vcov, Vector(Vector(17.0 / 21.0, -1.0 / 7.0), Vector(-1.0 / 7.0, 2.0 / 63.0)), 1e-10)
    assertEquals(workspace.coefTable.names, workspace.feNames)
    assertEquals(workspace.coefTable.stdErrors.length, 2)

  test("a fitted scaled alias matches its selected reduced model in ML and REML"):
    Vector(FitOptions.ml, FitOptions.reml).foreach: options =>
      val fit = fitted("one + days + doubleDays", options)
      val reference = fitted("one + days", options)
      assertEquals(fit.fixedEffects.retainedIndices, Vector(0, 1))
      assertEquals(fit.fixedEffects.droppedNames, Vector("doubleDays"))
      assertAligned(fit)
      assertEquivalent(fit, reference)

  test("full-rank column reordering permutes both covariance axes"):
    val fit = fitted("days + one", FitOptions.ml)
    assertEquals(fit.feNames, Vector("days", "one"))
    assertEquals(fit.fixedEffects.retainedIndices, Vector(0, 1))
    assertEquals(fit.fixedEffects.droppedNames, Vector.empty)
    assertEqualsDouble(fit.beta(0), reducedMl.beta(1), 1e-6)
    assertEqualsDouble(fit.beta(1), reducedMl.beta(0), 1e-6)
    assertMatrix(fit.vcov, reducedMl.vcov.reverse.map(_.reverse), 1e-5)

  test("a rank-zero fixed design has empty active outputs and records its dropped column"):
    val design = Lmm.compile("reaction ~ 0 + zero + (1 | subj)", frame).toOption.get
    assertEquals(design.fe.fixedEffects.retainedIndices, Vector.empty)
    assertEquals(design.fe.fixedEffects.droppedNames, Vector("zero"))
    val workspace = PlsWorkspace(design, reml = false).toOption.get
    workspace.objectiveAt(Vector(0.0)).fold(err => fail(err.message), identity)
    assertEquals(workspace.beta, Vector.empty)
    assertEquals(workspace.vcov, Vector.empty)
    assertEquals(workspace.stderror, Vector.empty)
    assertEquals(workspace.coefTable.names, Vector.empty)

  test("layout construction rejects repeated and out-of-range indices, including copy"):
    val layout = FixedEffectLayout(Vector("a", "b", "c"), Vector(2, 0))
    assertEquals(layout.retainedNames, Vector("c", "a"))
    assertEquals(layout.droppedIndices, Vector(1))
    Vector(Vector(0, 0), Vector(-1), Vector(3)).foreach: indices =>
      intercept[IllegalArgumentException](FixedEffectLayout(layout.requestedNames, indices))
      intercept[IllegalArgumentException](layout.copy(retainedIndices = indices))

  test("fit and coefficient table construction cannot silently truncate inconsistent shapes"):
    intercept[IllegalArgumentException](reducedMl.copy(beta = Vector(1.0)))
    intercept[IllegalArgumentException](reducedMl.copy(vcov = Vector(Vector(1.0), Vector(2.0, 3.0))))
    intercept[IllegalArgumentException](reducedMl.copy(vcov = Vector.empty))
    intercept[IllegalArgumentException](CoefTable.wald(Vector("a", "b"), Vector(1.0), Vector(2.0)))
    intercept[IllegalArgumentException](CoefTable.wald(Vector("a"), Vector(1.0), Vector(2.0, 3.0)))
    intercept[IllegalArgumentException](CoefTable(Vector("a"), Vector(1.0), Vector(2.0), Vector.empty))
    intercept[IllegalArgumentException](reducedMl.coefTable.copy(stdErrors = Vector.empty))

  private def assertAligned(fit: LmmFit): Unit =
    val rank = fit.fixedEffects.rank
    assertEquals(fit.beta.length, rank)
    assertEquals(fit.feNames.length, rank)
    assertEquals(fit.vcov.length, rank)
    assert(fit.vcov.forall(_.length == rank))
    assertEquals(fit.stderror.length, rank)
    assertEquals(fit.coefTable.names, fit.feNames)
    assertEquals(fit.coefTable.stdErrors, fit.stderror)
    assertEquals(fit.coefTable.zValues.length, rank)
    assert(fit.stderror.forall(_.isFinite))
    assert(fit.coefTable.zValues.forall(_.isFinite))

  private def assertEquivalent(fit: LmmFit, reference: LmmFit): Unit =
    assertEquals(fit.feNames, reference.feNames)
    assertEqualsDouble(fit.objective, reference.objective, 1e-8)
    assertEqualsDouble(fit.sigma, reference.sigma, 1e-8)
    fit.beta.indices.foreach(i => assertEqualsDouble(fit.beta(i), reference.beta(i), 1e-8))
    assertMatrix(fit.vcov, reference.vcov, 1e-8)
    assertEquals(fit.dof, reference.dof)

  private def assertMatrix(actual: Vector[Vector[Double]], expected: Vector[Vector[Double]], tolerance: Double): Unit =
    assertEquals(actual.length, expected.length)
    actual.indices.foreach: i =>
      assertEquals(actual(i).length, expected(i).length)
      actual(i).indices.foreach: j =>
        assertEqualsDouble(actual(i)(j), expected(i)(j), tolerance, clues(i, j))
