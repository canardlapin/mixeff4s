package mixeff4s.lmm

import mixeff4s.design.{CompiledDesign, FeMat, ReMat}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.linalg.{BlockedCholesky, MatrixBlock, WorkMat}
import mixeff4s.optimizer.{TrustBq, TrustBqOptions}

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

  private val sqrtwts = Array.fill(design.n)(1.0)
  private val workingY = Array.tabulate(design.n)(i => design.xy.xy(i, design.p))

  def setTheta(values: Vector[Double]): FitResult[Unit] =
    if values.length != nTheta then
      Left(MixedModelError.DimensionMismatch(s"theta length ${values.length}, expected $nTheta"))
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

  def updateL(): FitResult[Unit] =
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
    j = 0
    while j < k do
      BlockedCholesky.copyAndRmulLambda(
        lBlocks(MatrixBlock.blockIndex(k, j)),
        aBlocks(MatrixBlock.blockIndex(k, j)),
        lambdaWork(reterms(j))
      )
      j += 1
    MatrixBlock.copyBlock(lBlocks(MatrixBlock.blockIndex(k, k)), aBlocks(MatrixBlock.blockIndex(k, k)))

    val total = k + 1
    var col = 0
    var error: Option[MixedModelError] = None
    while col < total && error.isEmpty do
      val diagIdx = MatrixBlock.blockIndex(col, col)
      var jj = 0
      while jj < col do
        val off = lBlocks(MatrixBlock.blockIndex(col, jj)) match
          case MatrixBlock.Dense(mat) => mat
          case other                  => other.asDense
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
              BlockedCholesky.subtractProduct(
                lBlocks(targetIdx),
                lBlocks(MatrixBlock.blockIndex(i, prev)),
                lBlocks(MatrixBlock.blockIndex(col, prev))
              )
              prev += 1
            BlockedCholesky.rdivLowerTranspose(lBlocks(targetIdx), lBlocks(diagIdx))
            i += 1
      col += 1
    error.toLeft(())

  def objectiveAt(values: Vector[Double]): FitResult[Double] =
    for
      _ <- setTheta(values)
      _ <- updateL()
    yield objective

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
    else
      System.arraycopy(sqrtWeights, 0, sqrtwts, 0, n)
      System.arraycopy(yWork, 0, workingY, 0, n)
      recomputeABlocks()
      Right(())

  /** Conditional modes on the spherical scale, one `nRanef` vector per RE term. */
  def ranefU: Vector[Array[Double]] =
    val k = reterms.length
    val coef = beta
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
        val ljm = lBlocks(MatrixBlock.blockIndex(j, m)).asDense
        val vm = vVecs(m)
        var row = 0
        while row < rhs.length do
          var dot = 0.0
          var col = 0
          while col < vm.length do
            dot += ljm(row, col) * vm(col)
            col += 1
          rhs(row) -= dot
          row += 1
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
        val lmj = lBlocks(MatrixBlock.blockIndex(m, j)).asDense
        val um = uVecs(m)
        var row = 0
        while row < rhs.length do
          var dot = 0.0
          var col = 0
          while col < um.length do
            dot += lmj(col, row) * um(col)
            col += 1
          rhs(row) -= dot
          row += 1
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
    val mat = block match
      case MatrixBlock.Dense(m) =>
        m.fill(0.0)
        m
      case other =>
        throw IllegalArgumentException(s"RE×RE off-diagonal A block must stay Dense, got $other")
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
          mat(ri * a.vsize + si, rj * b.vsize + sj) += za * (sw * b.z(sj, obs))
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

  /** Active-column vcov: σ² (Lxx^{-1})' Lxx^{-1} in pivot order. */
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
      unpivotVcov(active)

  private def unpivotVcov(active: Array[Array[Double]]): Vector[Vector[Double]] =
    val piv = design.fe.piv
    val fullP = piv.length
    val p = active.length
    if p == 0 then Vector.empty
    else if p == fullP then
      val result = Array.fill(fullP, fullP)(0.0)
      var i = 0
      while i < fullP do
        var j = 0
        while j < fullP do
          result(piv(i))(piv(j)) = active(i)(j)
          j += 1
        i += 1
      result.iterator.map(_.toVector).toVector
    else
      val result = Array.fill(fullP, fullP)(Double.NaN)
      var i = 0
      while i < p do
        var j = 0
        while j < p do
          result(piv(i))(piv(j)) = active(i)(j)
          j += 1
        i += 1
      result.iterator.map(_.toVector).toVector

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
        if d > 0.0 then logdetLxx += math.log(d)
        i += 1
      logdet += 2.0 * logdetLxx
    (logdet, pwrss)

  private def lambdaWork(rt: ReMat): WorkMat =
    WorkMat.tabulate(rt.vsize, rt.vsize)((r, c) => rt.lambda(r, c))

object PlsWorkspace:
  def apply(design: CompiledDesign, reml: Boolean): FitResult[PlsWorkspace] =
    if design.reterms.isEmpty then Left(MixedModelError.NoRandomEffects)
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
    val result = WorkMat.zeros(a.nRanef, b.nRanef)
    var obs = 0
    while obs < a.nObs do
      val ri = a.refs(obs)
      val rj = b.refs(obs)
      var si = 0
      while si < a.vsize do
        val za = a.z(si, obs)
        var sj = 0
        while sj < b.vsize do
          result(ri * a.vsize + si, rj * b.vsize + sj) += za * b.z(sj, obs)
          sj += 1
        si += 1
      obs += 1
    MatrixBlock.Dense(result)

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

  def fit(design: CompiledDesign, reml: Boolean): FitResult[PlsWorkspace] =
    PlsWorkspace(design, reml).flatMap: workspace =>
      val start = workspace.theta
      TrustBq
        .minimize(start, workspace.lowerBounds, workspace.upperBounds, smallFamilyOptions(workspace.nTheta))(
          workspace.objectiveAt
        )
        .flatMap: result =>
          workspace.objectiveAt(result.x).map(_ => workspace)
