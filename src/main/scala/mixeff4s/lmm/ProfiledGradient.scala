package mixeff4s.lmm

import mixeff4s.design.{CompiledDesign, ReMat}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.linalg.{BlockedCholesky, MatrixBlock}

/** Backend describes exact computation, not a change of estimator. */
enum ProfiledGradientBackend:
  case LevelBlocks, Streamed

final case class ProfiledGradientOptions(maxWorkingBytes: Long = 64L * 1024L * 1024L)

final case class ProfiledEvaluation(
    objective: Double,
    gradient: Vector[Double],
    backend: ProfiledGradientBackend,
    workingBytes: Long
)

/** Exact random-coefficient-space derivative. Crossed determinant traces are streamed one column at a time; a single
  * random term uses independent level solves. Neither backend inverts Lambda or forms n by n covariance. Buffers are
  * reused between evaluations. The budget bounds numeric scratch, excluding existing A/L/design storage and JVM/JS
  * object headers (reported separately by the measurement harness).
  */
private[lmm] final class ProfiledGradient private (
    design: CompiledDesign,
    reml: Boolean,
    a: Array[MatrixBlock],
    l: Array[MatrixBlock],
    val workingBytes: Long
):
  private val terms = design.reterms
  private val k = terms.length
  private val p = design.p
  private val parameters = design.parmap.zipWithIndex
  private def vectors() = terms.map(t => new Array[Double](t.nRanef)).toArray
  private val e = vectors()
  private val u = vectors()
  private val g = vectors()
  private val lu = vectors()
  private val scratch = vectors()
  private val rhs = vectors()
  private val q = terms.map(t => new Array[Double](if reml then t.nRanef * p else 0)).toArray
  private val qc = terms.map(t => new Array[Double](if reml then t.nRanef * p else 0)).toArray
  private val vMinusT = terms.map(t => new Array[Double](if reml then t.nRanef * p else 0)).toArray
  private val fixed = new Array[Double](p)
  private val result = new Array[Double](design.nTheta)
  private val levelScratch = new Array[Double](terms.map(_.vsize).max)
  val backend: ProfiledGradientBackend =
    if k == 1 then ProfiledGradientBackend.LevelBlocks else ProfiledGradientBackend.Streamed

  private def clear(v: Array[Array[Double]]): Unit = v.foreach(java.util.Arrays.fill(_, 0.0))

  private def multiplyA(x: Array[Array[Double]], out: Array[Array[Double]]): Unit =
    clear(out)
    var i = 0
    while i < k do
      var j = 0
      while j <= i do
        val block = a(MatrixBlock.blockIndex(i, j))
        block.addMatVec(x(j), out(i))
        if j != i then block.addMatVec(x(i), out(j), transpose = true)
        j += 1
      i += 1

  /** Load A times a coordinate vector directly. A full matrix-vector multiply would visit every stored entry for each
    * streamed column, although only this column contributes.
    */
  private def crossColumn(term: Int, column: Int, out: Array[Array[Double]]): Unit =
    var i = 0
    while i < k do
      val block = if i >= term then a(MatrixBlock.blockIndex(i, term)) else a(MatrixBlock.blockIndex(term, i))
      var row = 0
      while row < out(i).length do
        out(i)(row) = if i >= term then block(row, column) else block(column, row)
        row += 1
      i += 1

  private def lambda(
      x: Array[Array[Double]],
      out: Array[Array[Double]],
      fitted: Vector[ReMat],
      transpose: Boolean
  ): Unit =
    clear(out)
    terms.indices.foreach: t =>
      val term = fitted(t)
      var lev = 0
      while lev < term.nLevels do
        val start = lev * term.vsize
        var r = 0
        while r < term.vsize do
          var c = 0
          while c <= r do
            if transpose then out(t)(start + c) += term.lambda(r, c) * x(t)(start + r)
            else out(t)(start + r) += term.lambda(r, c) * x(t)(start + c)
            c += 1
          r += 1
        lev += 1

  private def solveRandom(x: Array[Array[Double]]): Unit =
    var i = 0
    while i < k do
      var j = 0
      while j < i do
        l(MatrixBlock.blockIndex(i, j)).addMatVec(x(j), x(i), scale = -1.0)
        j += 1
      BlockedCholesky.solveLowerAgainstRhs(l(MatrixBlock.blockIndex(i, i)), x(i))
      i += 1
    i = k - 1
    while i >= 0 do
      var j = i + 1
      while j < k do
        l(MatrixBlock.blockIndex(j, i)).addMatVec(x(j), x(i), transpose = true, scale = -1.0)
        j += 1
      BlockedCholesky.solveUpperFromLowerTransposeAgainstRhs(l(MatrixBlock.blockIndex(i, i)), x(i))
      i -= 1

  private def validFactors: Boolean =
    def valid(value: Double): Boolean = value.isFinite && value >= BlockedCholesky.SolveZeroTolerance
    terms.indices.forall: t =>
      l(MatrixBlock.blockIndex(t, t)) match
        case MatrixBlock.Diagonal(values)      => values.forall(valid)
        case MatrixBlock.BlockDiagonal(blocks) => blocks.forall(b => (0 until b.rows).forall(i => valid(b(i, i))))
        case MatrixBlock.Dense(mat)            => (0 until mat.rows).forall(i => valid(mat(i, i)))
        case MatrixBlock.Sparse(_)             => false

  /** C^-1 times a fixed-effects vector, using Lxx without constructing an inverse. */
  private def solveFixed(x: Array[Double]): Unit =
    val last = l(MatrixBlock.blockIndex(k, k))
    var i = 0
    while i < p do
      var j = 0
      while j < i do
        x(i) -= last(i, j) * x(j)
        j += 1
      x(i) /= last(i, i)
      i += 1
    i = p - 1
    while i >= 0 do
      var j = i + 1
      while j < p do
        x(i) -= last(j, i) * x(j)
        j += 1
      x(i) /= last(i, i)
      i -= 1

  def evaluate(
      fitted: Vector[ReMat],
      beta: Vector[Double],
      objective: Double,
      pwrss: Double
  ): FitResult[ProfiledEvaluation] =
    val df = if reml then design.n - p else design.n
    val last = l(MatrixBlock.blockIndex(k, k))
    if !validFactors || !objective.isFinite || !pwrss.isFinite || pwrss <= 0.0 || df <= 0 ||
      (0 until p).exists(i => !last(i, i).isFinite || last(i, i) <= 0.0)
    then
      Left(
        MixedModelError.Optimization("analytic gradient requires finite positive factors and residual sum of squares")
      )
    else
      java.util.Arrays.fill(result, 0.0)
      // e = Z'y - Z'X beta from existing cross-products.
      terms.indices.foreach: t =>
        val cross = a(MatrixBlock.blockIndex(k, t))
        var row = 0
        while row < terms(t).nRanef do
          var value = cross(p, row)
          var j = 0
          while j < p do
            value -= cross(j, row) * beta(j)
            j += 1
          e(t)(row) = value
          row += 1
      lambda(e, u, fitted, transpose = true)
      solveRandom(u)
      lambda(u, lu, fitted, transpose = false)
      multiplyA(lu, g)
      terms.indices.foreach(t => e(t).indices.foreach(i => g(t)(i) = e(t)(i) - g(t)(i)))
      if reml then
        var col = 0
        while col < p do
          terms.indices.foreach: t =>
            val cross = a(MatrixBlock.blockIndex(k, t))
            scratch(t).indices.foreach(i => scratch(t)(i) = cross(col, i))
          lambda(scratch, rhs, fitted, transpose = true)
          solveRandom(rhs)
          terms.indices.foreach(t => rhs(t).indices.foreach(i => q(t)(i * p + col) = rhs(t)(i)))
          lambda(rhs, scratch, fitted, transpose = false)
          multiplyA(scratch, rhs)
          terms.indices.foreach: t =>
            val cross = a(MatrixBlock.blockIndex(k, t))
            rhs(t).indices.foreach(i => vMinusT(t)(i * p + col) = rhs(t)(i) - cross(col, i))
          col += 1
        terms.indices.foreach: t =>
          var row = 0
          while row < terms(t).nRanef do
            System.arraycopy(q(t), row * p, fixed, 0, p)
            solveFixed(fixed)
            System.arraycopy(fixed, 0, qc(t), row * p, p)
            row += 1
      // Trace entries D[b,a] where D=M^-1 Lambda' A. Stream exact columns with O(q) scratch.
      if k == 1 then levelTraces(fitted.head)
      else
        terms.indices.foreach: term =>
          var column = 0
          while column < terms(term).nRanef do
            crossColumn(term, column, rhs)
            lambda(rhs, scratch, fitted, transpose = true)
            solveRandom(scratch)
            val localRow = column % terms(term).vsize
            val start = column - localRow
            parameters.foreach: (param, index) =>
              if param._1 == term && param._2 == localRow then result(index) += 2.0 * scratch(term)(start + param._3)
            column += 1
      val scale = df.toDouble / pwrss
      parameters.foreach: (param, index) =>
        val (t, row, col) = param
        var lev = 0
        while lev < terms(t).nLevels do
          val r = lev * terms(t).vsize + row
          val c = lev * terms(t).vsize + col
          var contribution = -scale * g(t)(r) * u(t)(c)
          if reml then
            var j = 0
            while j < p do
              contribution += vMinusT(t)(r * p + j) * qc(t)(c * p + j)
              j += 1
          result(index) += 2.0 * contribution
          lev += 1
      if result.exists(v => !v.isFinite) then Left(MixedModelError.Optimization("analytic gradient is not finite"))
      else Right(ProfiledEvaluation(objective, result.toVector, backend, workingBytes))

  private def levelTraces(term: ReMat): Unit =
    val s = term.vsize
    val cross = a(0)
    val factor = l(0)
    def entry(block: MatrixBlock, lev: Int, row: Int, col: Int): Double = block match
      case MatrixBlock.Diagonal(values)      => if row == col then values(lev) else 0.0
      case MatrixBlock.BlockDiagonal(blocks) => blocks(lev)(row, col)
      case MatrixBlock.Dense(mat)            => mat(lev * s + row, lev * s + col)
      case MatrixBlock.Sparse(_)             => throw IllegalStateException("random diagonal factors cannot be sparse")
    var lev = 0
    while lev < term.nLevels do
      var column = 0
      while column < s do
        var row = 0
        while row < s do
          var value = 0.0
          var j = row
          while j < s do
            value += term.lambda(j, row) * entry(cross, lev, j, column)
            j += 1
          levelScratch(row) = value
          row += 1
        row = 0
        while row < s do
          var j = 0
          while j < row do
            levelScratch(row) -= entry(factor, lev, row, j) * levelScratch(j)
            j += 1
          levelScratch(row) /= entry(factor, lev, row, row)
          row += 1
        row = s - 1
        while row >= 0 do
          var j = row + 1
          while j < s do
            levelScratch(row) -= entry(factor, lev, j, row) * levelScratch(j)
            j += 1
          levelScratch(row) /= entry(factor, lev, row, row)
          row -= 1
        design.parmap.zipWithIndex.foreach: (param, index) =>
          if param._2 == column then result(index) += 2.0 * levelScratch(param._3)
        column += 1
      lev += 1

private[lmm] object ProfiledGradient:
  def apply(
      design: CompiledDesign,
      reml: Boolean,
      a: Array[MatrixBlock],
      l: Array[MatrixBlock],
      options: ProfiledGradientOptions
  ): FitResult[ProfiledGradient] =
    numericScratchBytes(design.reterms.map(_.nRanef), design.p, design.nTheta, design.reterms.map(_.vsize).max, reml)
      .flatMap: bytes =>
        if options.maxWorkingBytes <= 0 || bytes > options.maxWorkingBytes then
          Left(
            MixedModelError.InvalidArgument(
              s"analytic gradient numeric scratch requires $bytes bytes; budget is ${options.maxWorkingBytes}"
            )
          )
        else Right(new ProfiledGradient(design, reml, a, l, bytes))

  private[lmm] def numericScratchBytes(
      dimensions: Vector[Int],
      p: Int,
      nTheta: Int,
      maxLevelSize: Int,
      reml: Boolean
  ): FitResult[Long] =
    if dimensions.exists(_ < 0) || p < 0 || nTheta < 0 || maxLevelSize < 0 ||
      (reml && dimensions.exists(_.toLong * p.toLong > Int.MaxValue))
    then Left(MixedModelError.InvalidArgument("analytic gradient scratch dimensions exceed array limits"))
    else
      try
        val q = dimensions.foldLeft(0L)((sum, size) => Math.addExact(sum, size.toLong))
        val qp = if reml then Math.multiplyExact(3L, Math.multiplyExact(q, p.toLong)) else 0L
        val vectors = Math.addExact(Math.multiplyExact(6L, q), qp)
        val extras = Math.addExact(p.toLong, Math.addExact(nTheta.toLong, maxLevelSize.toLong))
        Right(Math.multiplyExact(8L, Math.addExact(vectors, extras)))
      catch
        case _: ArithmeticException =>
          Left(MixedModelError.InvalidArgument("analytic gradient scratch byte count overflows"))
