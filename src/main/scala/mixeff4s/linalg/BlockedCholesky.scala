package mixeff4s.linalg

import mixeff4s.error.{LinAlgError, MixedModelError}

/** Blocked Cholesky helpers ported from mixeff-rs `blocks.rs`. */
object BlockedCholesky:
  val SolveZeroTolerance: Double = 1e-30
  val DefaultZeroPadTolerance: Double = java.lang.Math.ulp(1.0)

  def zeroPadAbs(diagonalScale: Double, relativeTolerance: Double): Double =
    if !diagonalScale.isFinite || !relativeTolerance.isFinite then 0.0
    else relativeTolerance.max(0.0) * diagonalScale.max(0.0)

  def diagonalAbsMax(mat: WorkMat): Double =
    val n = mat.rows.min(mat.cols)
    var max = 0.0
    var i = 0
    while i < n do
      val a = math.abs(mat(i, i))
      if a > max then max = a
      i += 1
    max

  def cholesky(block: MatrixBlock, zeroPadTolerance: Double = DefaultZeroPadTolerance): Either[MixedModelError, Unit] =
    block match
      case MatrixBlock.Diagonal(diag) =>
        val tol = zeroPadAbs(diag.foldLeft(0.0)((m, v) => m.max(math.abs(v))), zeroPadTolerance)
        var i = 0
        while i < diag.length do
          if diag(i) <= 0.0 then
            if diag(i) < -tol then return Left(MixedModelError.PosDefException)
            diag(i) = 0.0
          else diag(i) = math.sqrt(diag(i))
          i += 1
        Right(())

      case MatrixBlock.BlockDiagonal(blocks) =>
        var b = 0
        while b < blocks.length do
          choleskyDense(blocks(b), zeroPadTolerance) match
            case Left(err) => return Left(err)
            case Right(_)  => ()
          b += 1
        Right(())

      case MatrixBlock.Dense(mat) =>
        choleskyDense(mat, zeroPadTolerance)

  private def choleskyDense(mat: WorkMat, zeroPadTolerance: Double): Either[MixedModelError, Unit] =
    val n = mat.rows
    if mat.cols != n then
      Left(
        MixedModelError.LinAlg(LinAlgError.DimensionMismatch(s"Cholesky requires square block, got ${n}x${mat.cols}"))
      )
    else
      val tol = zeroPadAbs(diagonalAbsMax(mat), zeroPadTolerance)
      var j = 0
      while j < n do
        var s = mat(j, j)
        var k = 0
        while k < j do
          s -= mat(j, k) * mat(j, k)
          k += 1
        if s <= 0.0 then
          if s < -tol then return Left(MixedModelError.PosDefException)
          var i = j
          while i < n do
            mat(i, j) = 0.0
            i += 1
        else
          val ljj = math.sqrt(s)
          mat(j, j) = ljj
          var i = j + 1
          while i < n do
            var t = mat(i, j)
            k = 0
            while k < j do
              t -= mat(i, k) * mat(j, k)
              k += 1
            mat(i, j) = t / ljj
            i += 1
          i = 0
          while i < j do
            mat(i, j) = 0.0
            i += 1
        j += 1
      Right(())

  /** C -= A * A' for a dense C. */
  def rankKDowndate(c: MatrixBlock, a: WorkMat): Unit =
    c match
      case MatrixBlock.Dense(cMat) =>
        var i = 0
        while i < cMat.rows do
          var j = 0
          while j <= i do
            var sum = 0.0
            var k = 0
            while k < a.cols do
              sum += a(i, k) * a(j, k)
              k += 1
            cMat(i, j) -= sum
            if i != j then cMat(j, i) -= sum
            j += 1
          i += 1
      case MatrixBlock.Diagonal(diag) =>
        var i = 0
        while i < diag.length do
          var sum = 0.0
          var k = 0
          while k < a.cols do
            val v = a(i, k)
            sum += v * v
            k += 1
          diag(i) -= sum
          i += 1
      case MatrixBlock.BlockDiagonal(blocks) =>
        var rowOffset = 0
        var b = 0
        while b < blocks.length do
          val blk = blocks(b)
          val s = blk.rows
          var i = 0
          while i < s do
            var j = 0
            while j <= i do
              var sum = 0.0
              var k = 0
              while k < a.cols do
                sum += a(rowOffset + i, k) * a(rowOffset + j, k)
                k += 1
              blk(i, j) -= sum
              if i != j then blk(j, i) -= sum
              j += 1
            i += 1
          rowOffset += s
          b += 1

  /** A = A * L^{-T} for lower-triangular L. */
  def rdivLowerTranspose(a: MatrixBlock, l: MatrixBlock): Unit =
    (a, l) match
      case (MatrixBlock.Dense(aMat), MatrixBlock.Diagonal(lDiag)) =>
        var j = 0
        while j < lDiag.length do
          val denom = lDiag(j)
          if math.abs(denom) < SolveZeroTolerance then
            var i = 0
            while i < aMat.rows do
              aMat(i, j) = 0.0
              i += 1
          else
            var i = 0
            while i < aMat.rows do
              aMat(i, j) /= denom
              i += 1
          j += 1

      case (MatrixBlock.Dense(aMat), MatrixBlock.BlockDiagonal(lBlocks)) =>
        var colOffset = 0
        var b = 0
        while b < lBlocks.length do
          val lBlk = lBlocks(b)
          val s = lBlk.rows
          if s == 2 then
            val c0 = colOffset
            val c1 = colOffset + 1
            val l00 = lBlk(0, 0)
            val l10 = lBlk(1, 0)
            val l11 = lBlk(1, 1)
            var i = 0
            while i < aMat.rows do
              val x0 = aMat(i, c0)
              if math.abs(l00) < SolveZeroTolerance then aMat(i, c0) = 0.0
              else aMat(i, c0) = x0 / l00
              if math.abs(l11) < SolveZeroTolerance then aMat(i, c1) = 0.0
              else aMat(i, c1) = (aMat(i, c1) - aMat(i, c0) * l10) / l11
              i += 1
          else
            var j = 0
            while j < s do
              val cj = colOffset + j
              if math.abs(lBlk(j, j)) < SolveZeroTolerance then
                var i = 0
                while i < aMat.rows do
                  aMat(i, cj) = 0.0
                  i += 1
              else
                var i = 0
                while i < aMat.rows do
                  var value = aMat(i, cj)
                  var k = 0
                  while k < j do
                    value -= aMat(i, colOffset + k) * lBlk(j, k)
                    k += 1
                  aMat(i, cj) = value / lBlk(j, j)
                  i += 1
              j += 1
          colOffset += s
          b += 1

      case (MatrixBlock.Dense(aMat), MatrixBlock.Dense(lMat)) =>
        val n = lMat.rows
        var j = 0
        while j < n do
          if math.abs(lMat(j, j)) < SolveZeroTolerance then
            var i = 0
            while i < aMat.rows do
              aMat(i, j) = 0.0
              i += 1
          else
            var i = 0
            while i < aMat.rows do
              var s = aMat(i, j)
              var k = 0
              while k < j do
                s -= aMat(i, k) * lMat(j, k)
                k += 1
              aMat(i, j) = s / lMat(j, j)
              i += 1
          j += 1

      case _ =>
        throw IllegalArgumentException("rdivLowerTranspose: unsupported block pair")

  /** L = Λ' A Λ + I for a same-term diagonal block. `lambda` is vsize × vsize. */
  def copyScaleInflate(l: MatrixBlock, a: MatrixBlock, lambda: WorkMat): Unit =
    val s = lambda.rows
    if s == 1 then
      val lamSq = lambda(0, 0) * lambda(0, 0)
      (l, a) match
        case (MatrixBlock.Diagonal(lDiag), MatrixBlock.Diagonal(aDiag)) =>
          var i = 0
          while i < lDiag.length do
            lDiag(i) = lamSq * aDiag(i) + 1.0
            i += 1
        case (MatrixBlock.Dense(lMat), MatrixBlock.Diagonal(aDiag)) =>
          lMat.fill(0.0)
          var i = 0
          while i < aDiag.length do
            lMat(i, i) = lamSq * aDiag(i) + 1.0
            i += 1
        case (MatrixBlock.Dense(lMat), MatrixBlock.Dense(aMat)) =>
          val n = aMat.rows
          var i = 0
          while i < n do
            var j = 0
            while j < n do
              lMat(i, j) = lamSq * aMat(i, j)
              j += 1
            lMat(i, i) += 1.0
            i += 1
        case _ =>
          throw IllegalArgumentException("copyScaleInflate: scalar block kind mismatch")
    else
      (l, a) match
        case (MatrixBlock.BlockDiagonal(lBlocks), MatrixBlock.BlockDiagonal(aBlocks)) =>
          var b = 0
          while b < aBlocks.length do
            scaleInflateSmall(lBlocks(b), aBlocks(b), lambda)
            b += 1
        case (MatrixBlock.Dense(lMat), MatrixBlock.BlockDiagonal(aBlocks)) =>
          lMat.fill(0.0)
          var b = 0
          while b < aBlocks.length do
            val tmp = WorkMat.zeros(s, s)
            scaleInflateSmall(tmp, aBlocks(b), lambda)
            var i = 0
            while i < s do
              var j = 0
              while j < s do
                lMat(b * s + i, b * s + j) = tmp(i, j)
                j += 1
              i += 1
            b += 1
        case _ =>
          throw IllegalArgumentException("copyScaleInflate: vector block kind mismatch")

  private def scaleInflateSmall(dst: WorkMat, src: WorkMat, lambda: WorkMat): Unit =
    val s = lambda.rows
    if s == 2 then
      val l00 = lambda(0, 0)
      val l01 = lambda(0, 1)
      val l10 = lambda(1, 0)
      val l11 = lambda(1, 1)
      val s00 = src(0, 0)
      val s01 = src(0, 1)
      val s10 = src(1, 0)
      val s11 = src(1, 1)
      val t00 = s00 * l00 + s01 * l10
      val t01 = s00 * l01 + s01 * l11
      val t10 = s10 * l00 + s11 * l10
      val t11 = s10 * l01 + s11 * l11
      dst(0, 0) = l00 * t00 + l10 * t10 + 1.0
      dst(0, 1) = l00 * t01 + l10 * t11
      dst(1, 0) = l01 * t00 + l11 * t10
      dst(1, 1) = l01 * t01 + l11 * t11 + 1.0
    else
      var row = 0
      while row < s do
        var col = 0
        while col < s do
          var sum = 0.0
          var ir = 0
          while ir < s do
            var ic = 0
            while ic < s do
              sum += lambda(ir, row) * src(ir, ic) * lambda(ic, col)
              ic += 1
            ir += 1
          dst(row, col) = sum
          col += 1
        dst(row, row) += 1.0
        row += 1

  /** L = A Λ for an FE×RE block. */
  def copyAndRmulLambda(l: MatrixBlock, a: MatrixBlock, lambda: WorkMat): Unit =
    val sj = lambda.rows
    val aMat = a match
      case MatrixBlock.Dense(mat) => mat
      case other                  => other.asDense
    val lMat = l match
      case MatrixBlock.Dense(mat) => mat
      case _                      => throw IllegalArgumentException("copyAndRmulLambda expects a dense L block")
    val nrows = aMat.rows
    val ncols = aMat.cols
    if sj == 1 then
      val lam = lambda(0, 0)
      var i = 0
      while i < nrows do
        var j = 0
        while j < ncols do
          lMat(i, j) = aMat(i, j) * lam
          j += 1
        i += 1
    else if sj == 2 then
      val l00 = lambda(0, 0)
      val l01 = lambda(0, 1)
      val l10 = lambda(1, 0)
      val l11 = lambda(1, 1)
      val nblocks = ncols / 2
      var b = 0
      while b < nblocks do
        val col0 = b * 2
        val col1 = col0 + 1
        var i = 0
        while i < nrows do
          val x0 = aMat(i, col0)
          val x1 = aMat(i, col1)
          lMat(i, col0) = x0 * l00 + x1 * l10
          lMat(i, col1) = x0 * l01 + x1 * l11
          i += 1
        b += 1
    else
      val nblocks = ncols / sj
      var b = 0
      while b < nblocks do
        var i = 0
        while i < nrows do
          var j = 0
          while j < sj do
            var sum = 0.0
            var inner = 0
            while inner < sj do
              sum += aMat(i, b * sj + inner) * lambda(inner, j)
              inner += 1
            lMat(i, b * sj + j) = sum
            j += 1
          i += 1
        b += 1

  /** L = Λ_i' A Λ_j for an off-diagonal RE×RE block. */
  def copyAndScaleOffdiag(l: MatrixBlock, a: MatrixBlock, lambdaI: WorkMat, lambdaJ: WorkMat): Unit =
    val si = lambdaI.rows
    val sj = lambdaJ.rows
    val aMat = a match
      case MatrixBlock.Dense(mat) => mat
      case other                  => other.asDense
    val lMat = l match
      case MatrixBlock.Dense(mat) => mat
      case _                      => throw IllegalArgumentException("copyAndScaleOffdiag expects a dense L block")
    if si == 1 && sj == 1 then
      val scale = lambdaI(0, 0) * lambdaJ(0, 0)
      var i = 0
      while i < aMat.rows do
        var j = 0
        while j < aMat.cols do
          lMat(i, j) = aMat(i, j) * scale
          j += 1
        i += 1
    else
      val nLevelsI = aMat.rows / si
      val nLevelsJ = aMat.cols / sj
      var bi = 0
      while bi < nLevelsI do
        var bj = 0
        while bj < nLevelsJ do
          var row = 0
          while row < si do
            var col = 0
            while col < sj do
              var sum = 0.0
              var ir = 0
              while ir < si do
                var ic = 0
                while ic < sj do
                  sum += lambdaI(ir, row) * aMat(bi * si + ir, bj * sj + ic) * lambdaJ(ic, col)
                  ic += 1
                ir += 1
              lMat(bi * si + row, bj * sj + col) = sum
              col += 1
            row += 1
          bj += 1
        bi += 1

  /** C -= A * B'. */
  def subtractProduct(c: MatrixBlock, a: MatrixBlock, b: MatrixBlock): Unit =
    val aMat = a match
      case MatrixBlock.Dense(mat) => mat
      case other                  => other.asDense
    val bMat = b match
      case MatrixBlock.Dense(mat) => mat
      case other                  => other.asDense
    c match
      case MatrixBlock.Dense(cMat) =>
        var i = 0
        while i < cMat.rows do
          var j = 0
          while j < cMat.cols do
            var sum = 0.0
            var k = 0
            while k < aMat.cols do
              sum += aMat(i, k) * bMat(j, k)
              k += 1
            cMat(i, j) -= sum
            j += 1
          i += 1
      case MatrixBlock.Diagonal(diag) =>
        var i = 0
        while i < diag.length do
          var sum = 0.0
          var k = 0
          while k < aMat.cols do
            sum += aMat(i, k) * bMat(i, k)
            k += 1
          diag(i) -= sum
          i += 1
      case MatrixBlock.BlockDiagonal(_) =>
        val dense = c.asDense
        var i = 0
        while i < dense.rows do
          var j = 0
          while j < dense.cols do
            var sum = 0.0
            var k = 0
            while k < aMat.cols do
              sum += aMat(i, k) * bMat(j, k)
              k += 1
            dense(i, j) -= sum
            j += 1
          i += 1
        throw IllegalArgumentException("subtractProduct cannot write back into BlockDiagonal")

  /** Forward-solve `L x = rhs` in place for a lower-triangular Cholesky block. */
  def solveLowerAgainstRhs(l: MatrixBlock, rhs: Array[Double]): Unit =
    l match
      case MatrixBlock.Diagonal(diag) =>
        var row = 0
        while row < diag.length do
          val denom = diag(row)
          if math.abs(denom) < SolveZeroTolerance then rhs(row) = 0.0
          else rhs(row) /= denom
          row += 1
      case MatrixBlock.BlockDiagonal(blocks) =>
        var rowOffset = 0
        var b = 0
        while b < blocks.length do
          val block = blocks(b)
          val s = block.rows
          var row = 0
          while row < s do
            val diag = block(row, row)
            if math.abs(diag) < SolveZeroTolerance then rhs(rowOffset + row) = 0.0
            else
              var sum = rhs(rowOffset + row)
              var inner = 0
              while inner < row do
                sum -= block(row, inner) * rhs(rowOffset + inner)
                inner += 1
              rhs(rowOffset + row) = sum / diag
            row += 1
          rowOffset += s
          b += 1
      case MatrixBlock.Dense(mat) =>
        var row = 0
        while row < mat.rows do
          val diag = mat(row, row)
          if math.abs(diag) < SolveZeroTolerance then rhs(row) = 0.0
          else
            var sum = rhs(row)
            var inner = 0
            while inner < row do
              sum -= mat(row, inner) * rhs(inner)
              inner += 1
            rhs(row) = sum / diag
          row += 1

  /** Back-solve `L' x = rhs` in place, treating `L` as lower-triangular. */
  def solveUpperFromLowerTransposeAgainstRhs(l: MatrixBlock, rhs: Array[Double]): Unit =
    l match
      case MatrixBlock.Diagonal(diag) =>
        var row = diag.length - 1
        while row >= 0 do
          val denom = diag(row)
          if math.abs(denom) < SolveZeroTolerance then rhs(row) = 0.0
          else rhs(row) /= denom
          row -= 1
      case MatrixBlock.BlockDiagonal(blocks) =>
        var rowOffset = 0
        var b = 0
        while b < blocks.length do
          val block = blocks(b)
          val s = block.rows
          var row = s - 1
          while row >= 0 do
            val diag = block(row, row)
            if math.abs(diag) < SolveZeroTolerance then rhs(rowOffset + row) = 0.0
            else
              var sum = rhs(rowOffset + row)
              var inner = row + 1
              while inner < s do
                sum -= block(inner, row) * rhs(rowOffset + inner)
                inner += 1
              rhs(rowOffset + row) = sum / diag
            row -= 1
          rowOffset += s
          b += 1
      case MatrixBlock.Dense(mat) =>
        var row = mat.rows - 1
        while row >= 0 do
          val diag = mat(row, row)
          if math.abs(diag) < SolveZeroTolerance then rhs(row) = 0.0
          else
            var sum = rhs(row)
            var inner = row + 1
            while inner < mat.rows do
              sum -= mat(inner, row) * rhs(inner)
              inner += 1
            rhs(row) = sum / diag
          row -= 1

  /** logdet(A) = 2 Σ log(diag(L)) for a Cholesky block. */
  def logdet(block: MatrixBlock): Double =
    val sumLogs =
      block match
        case MatrixBlock.Diagonal(diag) =>
          diag.iterator.filter(_ > 0.0).map(math.log).sum
        case MatrixBlock.BlockDiagonal(blocks) =>
          var ld = 0.0
          var b = 0
          while b < blocks.length do
            val blk = blocks(b)
            var i = 0
            while i < blk.rows.min(blk.cols) do
              val d = blk(i, i)
              if d > 0.0 then ld += math.log(d)
              i += 1
            b += 1
          ld
        case MatrixBlock.Dense(mat) =>
          var ld = 0.0
          val n = mat.rows.min(mat.cols)
          var i = 0
          while i < n do
            val d = mat(i, i)
            if d > 0.0 then ld += math.log(d)
            i += 1
          ld
    sumLogs * 2.0
