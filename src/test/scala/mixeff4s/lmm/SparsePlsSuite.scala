package mixeff4s.lmm

import mixeff4s.data.ModelFrame
import mixeff4s.design.CompiledDesign
import mixeff4s.linalg.MatrixBlock

class SparsePlsSuite extends munit.FunSuite:
  private val n = 24
  private val y = Vector.tabulate(n)(i => 3.0 + 0.4 * (i % 4) + math.sin(i * 1.8) + (i / 4) * 0.7)
  private val frame = ModelFrame
    .of(
      "y" -> ModelFrame.numeric(y),
      "x" -> ModelFrame.numeric(Vector.tabulate(n)(i => (i % 4 - 1).toDouble)),
      "g" -> ModelFrame.factor(Vector.tabulate(n)(i => s"g${i / 4}")),
      "h" -> ModelFrame.factor(Vector.tabulate(n)(i => s"h${i % 4}")),
      "nested" -> ModelFrame.factor(Vector.tabulate(n)(i => s"n${i / 8}")),
      "site" -> ModelFrame.factor(Vector.tabulate(n)(i => s"site${i / 12}")),
      "subject" -> ModelFrame.factor(Vector.tabulate(n)(i => s"s${i / 12}-${(i / 3) % 2}")),
      "item" -> ModelFrame.factor(Vector.tabulate(n)(i => s"i${i / 12}-${i % 2}")),
      "partial" -> ModelFrame.factor(Vector.tabulate(n)(i => s"p${if i < 12 then i / 4 else i % 3}"))
    )
    .toOption
    .get

  for
    random <- Vector(
      "(1 | g) + (1 | nested)",
      "(1 | g) + (1 | h)",
      "(1 | subject) + (1 | item) + (1 | site)",
      "(1 | g) + (1 | partial) + (1 | h)",
      "(1+x | g) + (1+x | nested)",
      "(1+x | g) + (1+x | h)"
    )
    reml <- Vector(false, true)
  do
    test(s"blocked factor equals independent dense factor across refreshes: $random REML=$reml"):
      val design = Lmm.compile(s"y ~ 1+x+$random", frame).toOption.get
      val workspace = PlsWorkspace(design, reml).toOption.get
      val theta = design.parmap.map { case (_, r, c) => if r == c then 0.8 else -0.2 }
      val weights = Vector(
        Array.tabulate(n)(i => if i % 3 == 0 then 0.0 else 1.0),
        Array.fill(n)(1.0),
        Array.tabulate(n)(i => 0.5 + 0.1 * (i % 7)),
        Array.fill(n)(1.0)
      )
      val baselinePattern = pattern(workspace.crossProductBlocksSnapshot)
      for sw <- weights; th <- Vector(theta, Vector.fill(theta.length)(0.0), theta) do
        workspace.updateIrlsWeights(sw, y.toArray).fold(e => fail(e.message), identity)
        val objective = workspace.objectiveAt(th).fold(e => fail(e.message), identity)
        val expected = denseFactor(design, th, sw, y)
        compareFactor(workspace, expected)
        assertEquals(pattern(workspace.crossProductBlocksSnapshot), baselinePattern)
        val q = design.nRanef
        val df = if reml then n - design.p else n
        val logdet = (0 until (if reml then q + design.p else q)).map(i => 2 * math.log(expected(i)(i))).sum
        val ss = expected.last.last * expected.last.last
        val want = logdet + df * (1 + math.log(2 * math.Pi * ss / df))
        assertEqualsDouble(objective, want, 2e-8 * math.max(1.0, math.abs(want)))
      // Response refresh keeps topology and agrees with a fully weighted rebuild.
      val changed = y.zipWithIndex.map((value, i) => value + 0.03 * i * i)
      workspace.updateResponse(changed).fold(e => fail(e.message), identity)
      workspace.objectiveAt(theta).fold(e => fail(e.message), identity)
      compareFactor(workspace, denseFactor(design, theta, Array.fill(n)(1.0), changed))

  test("invalid weight and response refreshes leave the previous cross-products intact"):
    val d = Lmm.compile("y ~ 1+x+(1|g)+(1|nested)", frame).toOption.get
    val w = PlsWorkspace(d, false).toOption.get
    val expected = w.objectiveAt(d.theta).toOption.get
    for bad <- Vector(Double.NaN, Double.PositiveInfinity, -1.0) do
      val weights = Array.fill(n)(1.0)
      weights(3) = bad
      assert(w.updateIrlsWeights(weights, y.toArray).isLeft)
      assertEqualsDouble(w.objectiveAt(d.theta).toOption.get, expected, 1e-12)
    val badY = y.toArray
    badY(3) = Double.NaN
    assert(w.updateIrlsWeights(Array.fill(n)(1.0), badY).isLeft)
    assertEqualsDouble(w.objectiveAt(d.theta).toOption.get, expected, 1e-12)

  private def pattern(blocks: Vector[MatrixBlock]): Vector[Vector[(Int, Int)]] =
    blocks.map:
      case MatrixBlock.Sparse(mat) =>
        val result = Vector.newBuilder[(Int, Int)]
        mat.foreachStored((r, c, _) => result += ((r, c)))
        result.result()
      case _ => Vector.empty

  private def compareFactor(workspace: PlsWorkspace, expected: Array[Array[Double]]): Unit =
    val blocks = workspace.factorBlocksSnapshot
    val sizes = workspace.design.reterms.map(_.nRanef) :+ (workspace.p + 1)
    val offsets = sizes.scanLeft(0)(_ + _)
    for i <- sizes.indices; j <- 0 to i do
      val block = blocks(MatrixBlock.blockIndex(i, j))
      for r <- 0 until sizes(i); c <- 0 until sizes(j) if offsets(i) + r >= offsets(j) + c do
        val want = expected(offsets(i) + r)(offsets(j) + c)
        assertEqualsDouble(block(r, c), want, 1e-8 * math.max(1.0, math.abs(want)), clues(i, j, r, c))

  /** Independently assemble observation rows [Z Lambda | X | y], then use a small dense scalar Cholesky. */
  private def denseFactor(
      d: CompiledDesign,
      theta: Vector[Double],
      sw: Array[Double],
      response: Vector[Double]
  ): Array[Array[Double]] =
    val offsets = d.reterms.map(_.nRanef).scanLeft(0)(_ + _)
    val q = offsets.last
    val size = q + d.p + 1
    val gram = Array.tabulate(size, size)((i, j) => if i == j && i < q then 1.0 else 0.0)
    var ti = 0
    val terms = d.reterms.map: term =>
      val result = term.withTheta(theta.slice(ti, ti + term.nTheta)).toOption.get
      ti += term.nTheta
      result
    for obs <- 0 until n do
      val row = new Array[Double](size)
      terms.zipWithIndex.foreach: (term, t) =>
        for col <- 0 until term.vsize do
          row(offsets(t) + term.refs(obs) * term.vsize + col) =
            sw(obs) * (col until term.vsize).map(r => term.z(r, obs) * term.lambda(r, col)).sum
      for col <- 0 until d.p do row(q + col) = sw(obs) * d.xy.xy(obs, col)
      row(size - 1) = sw(obs) * response(obs)
      for i <- 0 until size; j <- 0 to i do gram(i)(j) += row(i) * row(j)
    val lower = Array.ofDim[Double](size, size)
    for i <- 0 until size; j <- 0 to i do
      val remainder = gram(i)(j) - (0 until j).map(k => lower(i)(k) * lower(j)(k)).sum
      lower(i)(j) = if i == j then math.sqrt(remainder) else remainder / lower(j)(j)
      assert(lower(i)(j).isFinite, clues(i, j, remainder))
    lower
