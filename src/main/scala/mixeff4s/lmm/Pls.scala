package mixeff4s.lmm

import mixeff4s.design.{CompiledDesign, FeMat, ReMat}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.linalg.{BlockedCholesky, MatrixBlock, SparseWorkMat, WorkMat}
import mixeff4s.optimizer.{ObjectiveGradient, TrustBq, TrustBqOptions, TrustBqResult}

/** Working PLS state for a compiled LMM. */
final class PlsWorkspace private (
    val design: CompiledDesign,
    val reml: Boolean,
    private var reterms: Vector[ReMat],
    private val aBlocks: Array[MatrixBlock],
    private val lBlocks: Array[MatrixBlock]
):
  def n: Int = design.n
  def p: Int = design.p
  def nTheta: Int = design.nTheta
  def theta: Vector[Double] = reterms.flatMap(_.theta)
  def parmap: Vector[(Int, Int, Int)] = design.parmap
  def fittedReterms: Vector[ReMat] = reterms
  private[lmm] def factorBlocksSnapshot: Vector[MatrixBlock] = lBlocks.iterator.map(_.cloneBlock).toVector
  private[lmm] def crossProductBlocksSnapshot: Vector[MatrixBlock] = aBlocks.iterator.map(_.cloneBlock).toVector

  private val sqrtwts = Array.fill(design.n)(1.0)
  private val workingY = Array.tabulate(design.n)(i => design.xy.xy(i, design.p))

  def response: Vector[Double] = workingY.toVector

  /** Replace only response cross-products. Random design, topology and X'X are unchanged. */
  def updateResponse(values: Vector[Double]): FitResult[Unit] =
    if values.length != n then Left(MixedModelError.DimensionMismatch(s"response length ${values.length}, expected $n"))
    else if values.exists(v => !v.isFinite) then Left(MixedModelError.InvalidArgument("response must be finite"))
    else
      values.indices.foreach(i => workingY(i) = values(i))
      val k = reterms.length
      reterms.zipWithIndex.foreach: (term, index) =>
        val cross = aBlocks(MatrixBlock.blockIndex(k, index)).asDense
        var col = 0
        while col < cross.cols do
          cross(p, col) = 0.0
          col += 1
        var obs = 0
        while obs < n do
          val wy = sqrtwts(obs) * sqrtwts(obs) * workingY(obs)
          var j = 0
          while j < term.vsize do
            cross(p, term.refs(obs) * term.vsize + j) += wy * term.z(j, obs)
            j += 1
          obs += 1
      val cross = aBlocks(MatrixBlock.blockIndex(k, k)).asDense
      var j = 0
      while j <= p do
        var sum = 0.0
        var obs = 0
        while obs < n do
          val other = if j == p then workingY(obs) else design.xy.xy(obs, j)
          sum += sqrtwts(obs) * sqrtwts(obs) * workingY(obs) * other
          obs += 1
        cross(p, j) = sum
        cross(j, p) = sum
        j += 1
      Right(())

  def setTheta(values: Vector[Double]): FitResult[Unit] =
    if values.length != nTheta then
      Left(MixedModelError.DimensionMismatch(s"theta length ${values.length}, expected $nTheta"))
    else if values.exists(v => !v.isFinite) then
      Left(MixedModelError.InvalidArgument("theta must contain only finite values"))
    else if parmap.zip(values).exists { case ((_, row, col), value) => row == col && value < 0.0 } then
      Left(MixedModelError.InvalidArgument("theta diagonal entries must be nonnegative"))
    else
      var offset = 0
      val next = Vector.newBuilder[ReMat]
      var error: Option[MixedModelError] = None
      reterms.foreach: rt =>
        if error.isEmpty then
          rt.withTheta(values.slice(offset, offset + rt.nTheta)) match
            case Left(err)      => error = Some(err)
            case Right(updated) =>
              next += updated
              offset += rt.nTheta
      error.toLeft:
        reterms = next.result()
        ()

  def updateL(): FitResult[Unit] = factor(includeFixed = true)

  /** Factor only I + Lambda' Z' W Z Lambda for a fixed-beta conditional solve. The fixed/response Schur block need not
    * be positive definite in a valid GLMM. Call updateL before reading beta, sigma, vcov or an LMM objective again.
    */
  private[mixeff4s] def updateRandomL(): FitResult[Unit] = factor(includeFixed = false)

  private def factor(includeFixed: Boolean): FitResult[Unit] =
    val k = reterms.length
    var j = 0
    while j < k do
      BlockedCholesky.copyScaleInflate(
        lBlocks(MatrixBlock.blockIndex(j, j)),
        aBlocks(MatrixBlock.blockIndex(j, j)),
        lambdaWork(reterms(j))
      )
      j += 1
    var i = 1
    while i < k do
      j = 0
      while j < i do
        BlockedCholesky.copyAndScaleOffdiag(
          lBlocks(MatrixBlock.blockIndex(i, j)),
          aBlocks(MatrixBlock.blockIndex(i, j)),
          lambdaWork(reterms(i)),
          lambdaWork(reterms(j))
        )
        j += 1
      i += 1
    if includeFixed then
      j = 0
      while j < k do
        BlockedCholesky.copyAndRmulLambda(
          lBlocks(MatrixBlock.blockIndex(k, j)),
          aBlocks(MatrixBlock.blockIndex(k, j)),
          lambdaWork(reterms(j))
        )
        j += 1
      MatrixBlock.copyBlock(lBlocks(MatrixBlock.blockIndex(k, k)), aBlocks(MatrixBlock.blockIndex(k, k)))

    val total = if includeFixed then k + 1 else k
    var col = 0
    var error: Option[MixedModelError] = None
    while col < total && error.isEmpty do
      val diagIdx = MatrixBlock.blockIndex(col, col)
      var jj = 0
      while jj < col do
        val off = lBlocks(MatrixBlock.blockIndex(col, jj))
        lBlocks(diagIdx) = BlockedCholesky.prepareRankKTarget(lBlocks(diagIdx), off)
        BlockedCholesky.rankKDowndate(lBlocks(diagIdx), off)
        jj += 1
      BlockedCholesky.cholesky(lBlocks(diagIdx)) match
        case Left(err) => error = Some(err)
        case Right(_)  =>
          var i = col + 1
          while i < total do
            val targetIdx = MatrixBlock.blockIndex(i, col)
            var prev = 0
            while prev < col do
              lBlocks(targetIdx) = BlockedCholesky.prepareSubtractTarget(
                lBlocks(targetIdx),
                lBlocks(MatrixBlock.blockIndex(i, prev)),
                lBlocks(MatrixBlock.blockIndex(col, prev))
              )
              BlockedCholesky.subtractProduct(
                lBlocks(targetIdx),
                lBlocks(MatrixBlock.blockIndex(i, prev)),
                lBlocks(MatrixBlock.blockIndex(col, prev))
              )
              prev += 1
            lBlocks(targetIdx) = BlockedCholesky.prepareRightSolveTarget(lBlocks(targetIdx), lBlocks(diagIdx))
            BlockedCholesky.rdivLowerTranspose(lBlocks(targetIdx), lBlocks(diagIdx))
            i += 1
      col += 1
    error.toLeft(())

  def objectiveAt(values: Vector[Double]): FitResult[Double] =
    for
      _ <- setTheta(values)
      _ <- updateL()
      value = objective
      _ <-
        if value.isFinite then Right(())
        else
          Left(
            MixedModelError.Optimization(
              "profiled LMM objective requires positive finite residual sum of squares and fixed-effect pivots"
            )
          )
    yield value

  private var gradientScratch: Option[ProfiledGradient] = None

  /** One factorization feeds both the value and its exact derivative; cached buffers never replace factor refresh. */
  def objectiveAndGradientAt(
      values: Vector[Double],
      options: ProfiledGradientOptions = ProfiledGradientOptions()
  ): FitResult[ProfiledEvaluation] =
    for
      value <- objectiveAt(values)
      gradient <- gradientScratch match
        case Some(workspace) =>
          if options.maxWorkingBytes > 0 && workspace.workingBytes <= options.maxWorkingBytes then Right(workspace)
          else
            Left(
              MixedModelError.InvalidArgument(
                s"analytic gradient numeric scratch requires ${workspace.workingBytes} bytes; budget is ${options.maxWorkingBytes}"
              )
            )
        case None =>
          ProfiledGradient(design, reml, aBlocks, lBlocks, options).map: workspace =>
            gradientScratch = Some(workspace)
            workspace
      result <- gradient.evaluate(reterms, beta, value, pwrss)
    yield result

  def objective: Double =
    val denomdf = if reml then (n - p).toDouble else n.toDouble
    val (logdet, pwrss) = determinantAndPwrss
    logdet + denomdf * (1.0 + math.log(2.0 * math.Pi * pwrss / denomdf))

  def beta: Vector[Double] =
    val k = reterms.length
    val lLast = lBlocks(MatrixBlock.blockIndex(k, k)).asDense
    val pp1 = lLast.rows
    val rank = pp1 - 1
    if rank == 0 then Vector.empty
    else
      val coef = Array.tabulate(rank)(j => lLast(pp1 - 1, j))
      var i = rank - 1
      while i >= 0 do
        var s = coef(i)
        var j = i + 1
        while j < rank do
          s -= lLast(j, i) * coef(j)
          j += 1
        coef(i) = s / lLast(i, i)
        i -= 1
      coef.toVector

  def sigma: Double =
    val k = reterms.length
    val lLast = lBlocks(MatrixBlock.blockIndex(k, k)).asDense
    val lastDiag = math.abs(lLast(lLast.rows - 1, lLast.cols - 1))
    val denom = if reml then (n - p).toDouble else n.toDouble
    lastDiag / math.sqrt(denom)

  def pwrss: Double =
    val k = reterms.length
    val lLast = lBlocks(MatrixBlock.blockIndex(k, k)).asDense
    val d = lLast(lLast.rows - 1, lLast.cols - 1)
    d * d

  def logdetRe: Double =
    val k = reterms.length
    var ld = 0.0
    var j = 0
    while j < k do
      ld += BlockedCholesky.logdet(lBlocks(MatrixBlock.blockIndex(j, j)))
      j += 1
    ld

  def updateIrlsWeights(sqrtWeights: Array[Double], yWork: Array[Double]): FitResult[Unit] =
    if sqrtWeights.length != n || yWork.length != n then
      Left(
        MixedModelError.DimensionMismatch(
          s"IRLS weights/response length (${sqrtWeights.length}, ${yWork.length}), expected $n"
        )
      )
    else if sqrtWeights.exists(v => !v.isFinite || v < 0.0) || yWork.exists(v => !v.isFinite) then
      Left(
        MixedModelError.InvalidArgument(
          "IRLS square-root weights must be finite and nonnegative, and working responses finite"
        )
      )
    else
      System.arraycopy(sqrtWeights, 0, sqrtwts, 0, n)
      System.arraycopy(yWork, 0, workingY, 0, n)
      recomputeABlocks()
      Right(())

  /** Conditional modes on the spherical scale, one `nRanef` vector per RE term. */
  def ranefU: Vector[Array[Double]] = ranefUGivenBeta(beta)

  /** Weighted conditional solve using caller-supplied retained-basis beta and current RE factors. */
  private[mixeff4s] def ranefUGivenBeta(coef: Vector[Double]): Vector[Array[Double]] =
    require(coef.length == p, "conditional beta must match retained fixed-effect rank")
    val k = reterms.length
    val wr = Array.ofDim[Double](n)
    var obs = 0
    while obs < n do
      var v = sqrtwts(obs) * workingY(obs)
      var q = 0
      while q < p do
        v -= sqrtwts(obs) * design.xy.xy(obs, q) * coef(q)
        q += 1
      wr(obs) = v
      obs += 1

    val cVecs = Array.ofDim[Array[Double]](k)
    var term = 0
    while term < k do
      val re = reterms(term)
      val c = Array.ofDim[Double](re.nRanef)
      obs = 0
      while obs < n do
        val r = re.refs(obs)
        val sw = sqrtwts(obs)
        var s = 0
        while s < re.vsize do
          c(r * re.vsize + s) += sw * re.z(s, obs) * wr(obs)
          s += 1
        obs += 1
      val scaled = Array.ofDim[Double](re.nRanef)
      var lev = 0
      while lev < re.nLevels do
        var i = 0
        while i < re.vsize do
          var acc = 0.0
          var row = i
          while row < re.vsize do
            acc += re.lambda(row, i) * c(lev * re.vsize + row)
            row += 1
          scaled(lev * re.vsize + i) = acc
          i += 1
        lev += 1
      cVecs(term) = scaled
      term += 1

    val vVecs = Array.ofDim[Array[Double]](k)
    var j = 0
    while j < k do
      val rhs = cVecs(j).clone()
      var m = 0
      while m < j do
        lBlocks(MatrixBlock.blockIndex(j, m)).addMatVec(vVecs(m), rhs, scale = -1.0)
        m += 1
      BlockedCholesky.solveLowerAgainstRhs(lBlocks(MatrixBlock.blockIndex(j, j)), rhs)
      vVecs(j) = rhs
      j += 1

    val uVecs = Array.ofDim[Array[Double]](k)
    j = k - 1
    while j >= 0 do
      val rhs = vVecs(j).clone()
      var m = j + 1
      while m < k do
        lBlocks(MatrixBlock.blockIndex(m, j)).addMatVec(uVecs(m), rhs, transpose = true, scale = -1.0)
        m += 1
      BlockedCholesky.solveUpperFromLowerTransposeAgainstRhs(lBlocks(MatrixBlock.blockIndex(j, j)), rhs)
      uVecs(j) = rhs
      j -= 1
    uVecs.toVector

  private def recomputeABlocks(): Unit =
    val k = reterms.length
    var idx = 0
    var i = 0
    while i < k do
      var j = 0
      while j <= i do
        if i == j then fillReCross(aBlocks(idx), reterms(i))
        else fillReCrossOff(aBlocks(idx), reterms(i), reterms(j))
        idx += 1
        j += 1
      i += 1
    var term = 0
    while term < k do
      fillFeReCross(aBlocks(idx), reterms(term))
      idx += 1
      term += 1
    fillFeCross(aBlocks(idx))

  private def zeroBlock(block: MatrixBlock): Unit =
    block match
      case MatrixBlock.Dense(mat) =>
        mat.fill(0.0)
      case MatrixBlock.Sparse(mat)      => mat.clear()
      case MatrixBlock.Diagonal(values) =>
        java.util.Arrays.fill(values, 0.0)
      case MatrixBlock.BlockDiagonal(blocks) =>
        var b = 0
        while b < blocks.length do
          blocks(b).fill(0.0)
          b += 1

  private def fillReCross(block: MatrixBlock, re: ReMat): Unit =
    zeroBlock(block)
    val s = re.vsize
    var obs = 0
    while obs < n do
      val r = re.refs(obs)
      val sw = sqrtwts(obs)
      block match
        case MatrixBlock.Diagonal(diag) =>
          val z = sw * re.z(0, obs)
          diag(r) += z * z
        case MatrixBlock.BlockDiagonal(blocks) =>
          val blk = blocks(r)
          var si = 0
          while si < s do
            val zsi = sw * re.z(si, obs)
            var sj = 0
            while sj < s do
              blk(si, sj) += zsi * (sw * re.z(sj, obs))
              sj += 1
            si += 1
        case MatrixBlock.Sparse(_)  => throw IllegalStateException("diagonal random blocks cannot be sparse")
        case MatrixBlock.Dense(mat) =>
          var si = 0
          while si < s do
            val zsi = sw * re.z(si, obs)
            var sj = 0
            while sj < s do
              mat(r * s + si, r * s + sj) += zsi * (sw * re.z(sj, obs))
              sj += 1
            si += 1
      obs += 1

  private def fillReCrossOff(block: MatrixBlock, a: ReMat, b: ReMat): Unit =
    zeroBlock(block)
    var obs = 0
    while obs < n do
      val sw = sqrtwts(obs)
      val ri = a.refs(obs)
      val rj = b.refs(obs)
      var si = 0
      while si < a.vsize do
        val za = sw * a.z(si, obs)
        var sj = 0
        while sj < b.vsize do
          val row = ri * a.vsize + si
          val col = rj * b.vsize + sj
          val value = za * (sw * b.z(sj, obs))
          block match
            case MatrixBlock.Sparse(mat) => mat.add(row, col, value)
            case MatrixBlock.Dense(mat)  => mat(row, col) += value
            case _ => throw IllegalStateException("off-diagonal cross-product must be sparse or dense")
          sj += 1
        si += 1
      obs += 1

  private def fillFeReCross(block: MatrixBlock, re: ReMat): Unit =
    val mat = block match
      case MatrixBlock.Dense(m) =>
        m.fill(0.0)
        m
      case other =>
        throw IllegalArgumentException(s"FE×RE A block must stay Dense, got $other")
    val pp1 = p + 1
    var obs = 0
    while obs < n do
      val sw = sqrtwts(obs)
      val r = re.refs(obs)
      var col = 0
      while col < pp1 do
        val wx = if col < p then sw * design.xy.xy(obs, col) else sw * workingY(obs)
        var s = 0
        while s < re.vsize do
          mat(col, r * re.vsize + s) += wx * (sw * re.z(s, obs))
          s += 1
        col += 1
      obs += 1

  private def fillFeCross(block: MatrixBlock): Unit =
    val mat = block match
      case MatrixBlock.Dense(m) =>
        m.fill(0.0)
        m
      case other =>
        throw IllegalArgumentException(s"FE×FE A block must stay Dense, got $other")
    val pp1 = p + 1
    var i = 0
    while i < pp1 do
      var j = 0
      while j <= i do
        var sum = 0.0
        var obs = 0
        while obs < n do
          val sw = sqrtwts(obs)
          val vi = if i < p then sw * design.xy.xy(obs, i) else sw * workingY(obs)
          val vj = if j < p then sw * design.xy.xy(obs, j) else sw * workingY(obs)
          sum += vi * vj
          obs += 1
        mat(i, j) = sum
        if i != j then mat(j, i) = sum
        j += 1
      i += 1

  def varcorr: VarCorr = VarCorr.fromReterms(reterms, sigma)

  /** Active-column vcov: σ² (Lxx^{-1})' Lxx^{-1}, in the same order as beta and feNames. */
  def vcov: Vector[Vector[Double]] =
    val k = reterms.length
    val lLast = lBlocks(MatrixBlock.blockIndex(k, k)).asDense
    val rank = lLast.rows - 1
    if rank == 0 then Vector.empty
    else
      val lInv = Array.fill(rank, rank)(0.0)
      var j = 0
      while j < rank do
        lInv(j)(j) = 1.0
        var i = j
        while i < rank do
          var s = lInv(i)(j)
          var k2 = j
          while k2 < i do
            s -= lLast(i, k2) * lInv(k2)(j)
            k2 += 1
          lInv(i)(j) = s / lLast(i, i)
          i += 1
        j += 1
      val sigmaSq = sigma * sigma
      val active = Array.tabulate(rank, rank): (r, c) =>
        var sum = 0.0
        var t = 0
        while t < rank do
          sum += lInv(t)(r) * lInv(t)(c)
          t += 1
        sigmaSq * sum
      active.iterator.map(_.toVector).toVector

  def stderror: Vector[Double] =
    vcov.zipWithIndex.map((row, i) => math.sqrt(row(i)))

  def feNames: Vector[String] = design.fe.fullRankNames

  def coefTable: CoefTable = CoefTable.wald(feNames, beta, stderror)

  def lowerBounds: Vector[Double] =
    parmap.map: (_, row, col) =>
      if row == col then 0.0 else Double.NegativeInfinity

  def upperBounds: Vector[Double] =
    Vector.fill(nTheta)(Double.PositiveInfinity)

  private def determinantAndPwrss: (Double, Double) =
    val k = reterms.length
    var logdet = 0.0
    var j = 0
    while j < k do
      logdet += BlockedCholesky.logdet(lBlocks(MatrixBlock.blockIndex(j, j)))
      j += 1
    val lDense = lBlocks(MatrixBlock.blockIndex(k, k)).asDense
    val pp1 = lDense.rows
    val lastDiag = lDense(pp1 - 1, pp1 - 1)
    val pwrss = lastDiag * lastDiag
    if reml then
      var logdetLxx = 0.0
      var i = 0
      while i < pp1 - 1 do
        val d = lDense(i, i)
        if !d.isFinite || d <= 0.0 then return (Double.PositiveInfinity, pwrss)
        logdetLxx += math.log(d)
        i += 1
      logdet += 2.0 * logdetLxx
    (logdet, pwrss)

  private def lambdaWork(rt: ReMat): WorkMat =
    WorkMat.tabulate(rt.vsize, rt.vsize)((r, c) => rt.lambda(r, c))

object PlsWorkspace:
  def apply(design: CompiledDesign, reml: Boolean): FitResult[PlsWorkspace] =
    if design.n <= 0 || (reml && design.n <= design.p) then
      Left(
        MixedModelError.InvalidArgument("LMM requires observations and positive residual degrees of freedom for REML")
      )
    else if design.reterms.isEmpty then Left(MixedModelError.NoRandomEffects)
    else
      val k = design.reterms.length
      val blocks = Vector.newBuilder[MatrixBlock]
      var i = 0
      while i < k do
        var j = 0
        while j <= i do
          if i == j then blocks += reCross(design.reterms(i))
          else blocks += reCrossOff(design.reterms(i), design.reterms(j))
          j += 1
        i += 1
      var term = 0
      while term < k do
        blocks += feReCross(design.xy, design.reterms(term))
        term += 1
      blocks += feCross(design.xy)
      val a = blocks.result().toArray
      val l = a.map(_.cloneBlock)
      promoteCrossedFillIn(l, design.reterms)
      Right(new PlsWorkspace(design, reml, design.reterms, a, l))

  private def isNested(a: ReMat, b: ReMat): Boolean =
    if a.refs.length != b.refs.length then false
    else
      val bins = Array.fill(a.nLevels)(-1)
      a.refs.indices.forall: obs =>
        val aref = a.refs(obs)
        val bref = b.refs(obs)
        if bins(aref) < 0 then
          bins(aref) = bref
          true
        else bins(aref) == bref

  private def promoteCrossedFillIn(l: Array[MatrixBlock], reterms: Vector[ReMat]): Unit =
    val k = reterms.length
    var i = 1
    while i < k do
      if (0 until i).exists(j => !isNested(reterms(j), reterms(i))) then
        var row = i
        while row < k do
          val idx = MatrixBlock.blockIndex(row, i)
          l(idx) match
            case MatrixBlock.Dense(_) => ()
            case other                => l(idx) = MatrixBlock.Dense(other.asDense)
          row += 1
      i += 1

  private def reCross(re: ReMat): MatrixBlock =
    val s = re.vsize
    val nLevels = re.nLevels
    if s == 1 then
      val diag = new Array[Double](nLevels)
      var obs = 0
      while obs < re.nObs do
        val r = re.refs(obs)
        val z = re.z(0, obs)
        diag(r) += z * z
        obs += 1
      MatrixBlock.Diagonal(diag)
    else
      val blocks = Array.fill(nLevels)(WorkMat.zeros(s, s))
      var obs = 0
      while obs < re.nObs do
        val k = re.refs(obs)
        val blk = blocks(k)
        var si = 0
        while si < s do
          val zsi = re.z(si, obs)
          var sj = 0
          while sj < s do
            blk(si, sj) += zsi * re.z(sj, obs)
            sj += 1
          si += 1
        obs += 1
      MatrixBlock.BlockDiagonal(blocks)

  private def reCrossOff(a: ReMat, b: ReMat): MatrixBlock =
    // Incidence determines topology, including basis entries currently equal to zero.
    val pairs = a.refs.zip(b.refs).distinct
    // Dense crossings pay CSC indexing costs without enough structural zeros to save work.
    // Choose once from incidence, never from current numeric values or weights.
    val cells = a.nRanef.toLong * b.nRanef
    val block =
      if pairs.length.toLong * 2 >= a.nLevels.toLong * b.nLevels && cells <= Int.MaxValue then
        MatrixBlock.Dense(WorkMat.zeros(a.nRanef, b.nRanef))
      else
        val coordinates = pairs.flatMap: (i, j) =>
          (0 until a.vsize).flatMap(r => (0 until b.vsize).map(c => (i * a.vsize + r, j * b.vsize + c)))
        MatrixBlock.Sparse(SparseWorkMat.fromCoordinates(a.nRanef, b.nRanef, coordinates))
    var obs = 0
    while obs < a.nObs do
      val ri = a.refs(obs)
      val rj = b.refs(obs)
      var si = 0
      while si < a.vsize do
        val za = a.z(si, obs)
        var sj = 0
        while sj < b.vsize do
          val row = ri * a.vsize + si
          val col = rj * b.vsize + sj
          val value = za * b.z(sj, obs)
          block match
            case MatrixBlock.Dense(mat)  => mat(row, col) += value
            case MatrixBlock.Sparse(mat) => mat.add(row, col, value)
            case _ => throw IllegalStateException("off-diagonal cross-products require dense or sparse storage")
          sj += 1
        si += 1
      obs += 1
    block

  private def feReCross(xy: FeMat, re: ReMat): MatrixBlock =
    val pp1 = xy.xy.cols
    val result = WorkMat.zeros(pp1, re.nRanef)
    var obs = 0
    while obs < re.nObs do
      val r = re.refs(obs)
      var col = 0
      while col < pp1 do
        var s = 0
        while s < re.vsize do
          result(col, r * re.vsize + s) += xy.xy(obs, col) * re.z(s, obs)
          s += 1
        col += 1
      obs += 1
    MatrixBlock.Dense(result)

  private def feCross(xy: FeMat): MatrixBlock =
    val pp1 = xy.xy.cols
    val n = xy.n
    val result = WorkMat.zeros(pp1, pp1)
    var i = 0
    while i < pp1 do
      var j = 0
      while j <= i do
        var sum = 0.0
        var obs = 0
        while obs < n do
          sum += xy.xy(obs, i) * xy.xy(obs, j)
          obs += 1
        result(i, j) = sum
        if i != j then result(j, i) = sum
        j += 1
      i += 1
    MatrixBlock.Dense(result)

/** Optimizer evidence is a value, separate from the mutable workspace it evaluated. */
final case class PlsFit(workspace: PlsWorkspace, optimizer: TrustBqResult)

object Pls:
  def smallFamilyOptions(nTheta: Int): TrustBqOptions =
    TrustBqOptions(
      initialRadius = 0.75,
      finalRadius = 1e-5,
      maxEvaluations = 1000,
      ftolAbs = 1e-10,
      ftolRel = 1e-11,
      ftolRequiresLocalRadius = true,
      maxCrossTerms = if nTheta <= 3 then Int.MaxValue else 0,
      stallIterations = 4,
      stallRequiresStableX = true
    )

  def fit(
      design: CompiledDesign,
      reml: Boolean,
      maxEvaluations: Int = 1000,
      optimizer: LmmOptimizer = LmmOptimizer.DerivativeFree,
      gradientOptions: ProfiledGradientOptions = ProfiledGradientOptions(),
      optimizerControls: Option[TrustBqOptions] = None
  ): FitResult[PlsFit] =
    PlsWorkspace(design, reml).flatMap(workspace =>
      optimize(workspace, maxEvaluations, optimizer, gradientOptions, optimizerControls = optimizerControls)
    )

  /** Explicit LMM gradient policy. Large covariance models may plateau while theta drifts along a singular ridge;
    * smaller models keep the strict stable-theta stall rule. Both policies guard objective-based stops with the
    * projected gradient; step/radius stops remain numerical stops without a stationarity guarantee.
    */
  private[lmm] def gradientFamilyOptions(nTheta: Int): TrustBqOptions =
    val small = smallFamilyOptions(nTheta).copy(gradientTolerance = 1e-9, finalRadius = 1e-6)
    val numeric = if nTheta <= 3 then small else small.copy(ftolAbs = 1e-8, ftolRel = 1e-10)
    if nTheta < 7 then numeric
    else numeric.copy(stallIterations = 3, stallFtolRel = 1e-6, stallFtolAbs = 1e-8, stallRequiresStableX = false)

  private[lmm] def optimize(
      workspace: PlsWorkspace,
      maxEvaluations: Int,
      optimizer: LmmOptimizer,
      gradientOptions: ProfiledGradientOptions,
      initial: Option[Vector[Double]] = None,
      optimizerControls: Option[TrustBqOptions] = None
  ): FitResult[PlsFit] =
    val start = initial match
      case Some(values) => workspace.setTheta(values).map(_ => values)
      case None         => workspace.setTheta(workspace.design.theta).flatMap(_ => diagonalFirstStart(workspace))
    start.flatMap: values =>
      val options =
        optimizerControls.getOrElse(smallFamilyOptions(workspace.nTheta)).copy(maxEvaluations = maxEvaluations)
      val optimized = optimizer match
        case LmmOptimizer.DerivativeFree =>
          TrustBq.minimize(values, workspace.lowerBounds, workspace.upperBounds, options)(workspace.objectiveAt)
        case LmmOptimizer.AnalyticGradient =>
          val gradientPolicy =
            optimizerControls.getOrElse(gradientFamilyOptions(workspace.nTheta)).copy(maxEvaluations = maxEvaluations)
          TrustBq.minimizeWithGradient(values, workspace.lowerBounds, workspace.upperBounds, gradientPolicy): point =>
            workspace
              .objectiveAndGradientAt(point, gradientOptions)
              .map(e => ObjectiveGradient(e.objective, e.gradient))
      optimized.flatMap: result =>
        for
          _ <- result.requireConvergence
          value <- workspace.objectiveAt(result.x)
          _ <-
            if value.isFinite then Right(())
            else Left(MixedModelError.Optimization("LMM objective is not finite at the TrustBQ minimizer"))
        yield PlsFit(workspace, result)

  /** Coarse zero-correlation warm start when θ has off-diagonal slots.
    *
    * Single-start TrustBQ from the identity λ can stall on a worse interior (orthodont: 433.40 vs 432.58). mixeff4s has
    * no BOBYQA fallback, so the default fit runs the rust diagonal-first ladder: optimize the diagonal with
    * correlations pinned at 0, then hand that point to the full solve.
    */
  private def diagonalFirstStart(workspace: PlsWorkspace): FitResult[Vector[Double]] =
    val start = workspace.theta
    val diagIdx = workspace.parmap.zipWithIndex.collect:
      case ((_, row, col), i) if row == col => i
    if diagIdx.isEmpty || diagIdx.length == start.length then Right(start)
    else
      val reducedStart = diagIdx.map(start)
      val stageOptions = TrustBqOptions(
        initialRadius = 0.75,
        finalRadius = 1e-3,
        maxEvaluations = 40,
        ftolAbs = 1e-4,
        ftolRel = 1e-6,
        maxCrossTerms = if diagIdx.length <= 3 then Int.MaxValue else 0,
        stallIterations = 3,
        stallRequiresStableX = false
      )
      val stage = TrustBq.minimize(
        reducedStart,
        Vector.fill(diagIdx.length)(0.0),
        Vector.fill(diagIdx.length)(Double.PositiveInfinity),
        stageOptions
      ): reduced =>
        workspace.objectiveAt(expandDiagonal(start, workspace.parmap, diagIdx, reduced))
      stage match
        case Right(result) if result.fmin.isFinite =>
          Right(expandDiagonal(start, workspace.parmap, diagIdx, result.x))
        case _ => Right(start)

  private def expandDiagonal(
      start: Vector[Double],
      parmap: Vector[(Int, Int, Int)],
      diagIdx: Vector[Int],
      reduced: Vector[Double]
  ): Vector[Double] =
    val full = start.toArray
    parmap.zipWithIndex.foreach:
      case ((_, row, col), i) if row != col => full(i) = 0.0
      case _                                => ()
    diagIdx.zip(reduced).foreach((i, value) => full(i) = value)
    full.toVector
