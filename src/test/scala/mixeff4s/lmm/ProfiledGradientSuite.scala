package mixeff4s.lmm

import mixeff4s.data.ModelFrame
import mixeff4s.design.CompiledDesign

class ProfiledGradientSuite extends munit.FunSuite:
  private val n = 12
  private val data = ModelFrame
    .of(
      "y" -> ModelFrame.numeric(Vector.tabulate(n)(i => 2.0 + 0.7 * (i % 4) + math.sin(i * 1.7) + (i / 4))),
      "x" -> ModelFrame.numeric(Vector.tabulate(n)(i => (i % 4 - 1).toDouble)),
      "g" -> ModelFrame.factor(Vector.tabulate(n)(i => s"g${i / 4}")),
      "h" -> ModelFrame.factor(Vector.tabulate(n)(i => s"h${i % 3}")),
      "nested" -> ModelFrame.factor(Vector.tabulate(n)(i => s"n${i / 8}"))
    )
    .toOption
    .get

  for
    formula <- Vector(
      "y ~ 1 + x + (1 | g)",
      "y ~ 1 + x + (1 + x | g)",
      "y ~ 1 + x + (1 + x | g) + (1 | h)",
      "y ~ 1 + x + (1 | g) + (1 | nested)",
      "y ~ 0 + (1 | g)"
    )
    reml <- Vector(false, true)
  do
    test(s"observation-space derivative oracle: $formula, REML=$reml"):
      val d = Lmm.compile(formula, data).toOption.get
      val w = PlsWorkspace(d, reml).toOption.get
      val interior = d.parmap.zipWithIndex.map { case ((_, r, c), i) => if r == c then 0.6 + 0.1 * i else -0.2 }
      val singular = d.parmap.map { case (_, r, c) => if r == c then 0.0 else 0.8 }
      for theta <- Vector(interior, Vector.fill(d.nTheta)(0.0), singular, interior) do
        val actual = w.objectiveAndGradientAt(theta).fold(err => fail(err.message), identity)
        val expected = denseOracle(d, theta, reml)
        actual.gradient
          .zip(expected)
          .foreach: (got, want) =>
            assertEqualsDouble(got, want, 2e-9 * math.max(1.0, math.abs(want)))
        assert(actual.workingBytes <= ProfiledGradientOptions().maxWorkingBytes)
        assertEquals(
          actual.backend,
          if d.nReTerms == 1 then ProfiledGradientBackend.LevelBlocks else ProfiledGradientBackend.Streamed
        )
      val analytic = w.objectiveAndGradientAt(interior).toOption.get.gradient
      interior.indices.foreach: j =>
        val h = 1e-5
        val lo = w.objectiveAt(interior.updated(j, interior(j) - h)).toOption.get
        val hi = w.objectiveAt(interior.updated(j, interior(j) + h)).toOption.get
        val fd = (hi - lo) / (2 * h)
        assertEqualsDouble(analytic(j), fd, 2e-5 * math.max(1.0, math.abs(fd)))

  test("zero Lambda has zero gradient even at a variance maximum: closed-form witness"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(-3.1, -2.9, -1.1, -0.9, 0.9, 1.1, 2.9, 3.1)),
        "g" -> ModelFrame.factor(Vector("a", "a", "b", "b", "c", "c", "d", "d"))
      )
      .toOption
      .get
    val d = Lmm.compile("y ~ 1 + (1 | g)", frame).toOption.get
    for reml <- Vector(false, true) do
      val w = PlsWorkspace(d, reml).toOption.get
      val zero = w.objectiveAndGradientAt(Vector(0.0)).toOption.get
      assertEqualsDouble(zero.gradient.head, 0.0, 1e-12)
      val point = w.objectiveAndGradientAt(Vector(0.1)).toOption.get
      assert(point.objective < zero.objective)
      val t = 0.1
      val m = 1.0 + 2.0 * t * t
      val s = 0.08 + 40.0 / m
      val dlog = (if reml then 3.0 else 4.0) * 4.0 * t / m
      val ds = -160.0 * t / (m * m)
      val expected = dlog + (if reml then 7.0 else 8.0) * ds / s
      assertEqualsDouble(point.gradient.head, expected, 1e-10)

  test("gradient scratch budget is enforced without a dense inverse allocation"):
    val d = Lmm.compile("y ~ 1 + x + (1 + x | g) + (1 | h)", data).toOption.get
    val w = PlsWorkspace(d, true).toOption.get
    assert(w.objectiveAndGradientAt(d.theta, ProfiledGradientOptions(1)).isLeft)
    val evaluated = w.objectiveAndGradientAt(d.theta).toOption.get
    assert(w.objectiveAndGradientAt(d.theta, ProfiledGradientOptions(evaluated.workingBytes)).isRight)
    assert(w.objectiveAndGradientAt(d.theta, ProfiledGradientOptions(evaluated.workingBytes - 1)).isLeft)

  test("gradient scratch sizing rejects impossible arrays before allocation"):
    assert(ProfiledGradient.numericScratchBytes(Vector(Int.MaxValue), 2, 1, 1, true).isLeft)
    assert(ProfiledGradient.numericScratchBytes(Vector(-1), 1, 1, 1, false).isLeft)
    assertEquals(ProfiledGradient.numericScratchBytes(Vector(4, 6), 2, 3, 2, true), Right(1016L))

  // Independent observation-space oracle: Sigma=I+Z Lambda Lambda' Z',
  // P=Sigma^-1 - Sigma^-1 X (X' Sigma^-1 X)^-1 X' Sigma^-1.
  // Uses pivoted Gauss-Jordan inversion, not the production blocked Cholesky.
  private type Mat = Array[Array[Double]]
  private def transpose(a: Mat): Mat =
    if a.isEmpty then Array.empty else Array.tabulate(a(0).length, a.length)((i, j) => a(j)(i))
  private def multiply(a: Mat, b: Mat): Mat =
    val cols = if b.isEmpty then 0 else b(0).length
    Array.tabulate(a.length, cols)((i, j) => a(i).indices.map(k => a(i)(k) * b(k)(j)).sum)
  private def inverse(a: Mat): Mat =
    val size = a.length
    val aug = Array.tabulate(size, 2 * size)((i, j) => if j < size then a(i)(j) else if j - size == i then 1.0 else 0.0)
    var c = 0
    while c < size do
      val pivot = (c until size).maxBy(i => math.abs(aug(i)(c)))
      val row = aug(c); aug(c) = aug(pivot); aug(pivot) = row
      assert(math.abs(aug(c)(c)) > 1e-14)
      val diag = aug(c)(c)
      aug(c).indices.foreach(j => aug(c)(j) /= diag)
      (0 until size)
        .filter(_ != c)
        .foreach: i =>
          val scale = aug(i)(c)
          aug(i).indices.foreach(j => aug(i)(j) -= scale * aug(c)(j))
      c += 1
    aug.map(_.drop(size))

  private def denseOracle(d: CompiledDesign, theta: Vector[Double], reml: Boolean): Vector[Double] =
    var cursor = 0
    val terms = d.reterms.map: term =>
      val fitted = term.withTheta(theta.slice(cursor, cursor + term.nTheta)).toOption.get
      cursor += term.nTheta
      fitted
    val sigma = Array.tabulate(n, n)((i, j) => if i == j then 1.0 else 0.0)
    terms.foreach: term =>
      for lev <- 0 until term.nLevels; col <- 0 until term.vsize do
        val z = Array.tabulate(n): i =>
          if term.refs(i) != lev then 0.0 else (col until term.vsize).map(r => term.z(r, i) * term.lambda(r, col)).sum
        for i <- 0 until n; j <- 0 until n do sigma(i)(j) += z(i) * z(j)
    val inv = inverse(sigma)
    val x = Array.tabulate(n, d.p)((i, j) => d.xy.xy(i, j))
    val pmat = inv.map(_.clone())
    if d.p > 0 then
      val vx = multiply(inv, x)
      val correction = multiply(multiply(vx, inverse(multiply(transpose(x), vx))), transpose(vx))
      for i <- 0 until n; j <- 0 until n do pmat(i)(j) -= correction(i)(j)
    val y = Array.tabulate(n)(i => d.xy.xy(i, d.p))
    val residual = Array.tabulate(n)(i => (0 until n).map(j => pmat(i)(j) * y(j)).sum)
    val ss = y.indices.map(i => y(i) * residual(i)).sum
    d.parmap.map: (t, r, c) =>
      val term = terms(t)
      val derivative = Array.ofDim[Double](n, n)
      for lev <- 0 until term.nLevels do
        val z = Array.tabulate(n)(i => if term.refs(i) == lev then term.z(r, i) else 0.0)
        val w = Array.tabulate(n)(i =>
          if term.refs(i) != lev then 0.0 else (c until term.vsize).map(j => term.z(j, i) * term.lambda(j, c)).sum
        )
        for i <- 0 until n; j <- 0 until n do derivative(i)(j) += z(i) * w(j) + w(i) * z(j)
      val traceMatrix = if reml then pmat else inv
      val trace = (0 until n).map(i => (0 until n).map(j => traceMatrix(i)(j) * derivative(j)(i)).sum).sum
      val quad = (0 until n).map(i => (0 until n).map(j => residual(i) * derivative(i)(j) * residual(j)).sum).sum
      trace - (if reml then n - d.p else n).toDouble * quad / ss
