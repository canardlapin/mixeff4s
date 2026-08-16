package mixeff4s.linalg

/** One block of the mixed-model A or L array. Contents are mutated in place. */
enum MatrixBlock:
  case Dense(mat: WorkMat)
  case Diagonal(values: Array[Double])
  case BlockDiagonal(blocks: Array[WorkMat])

  def cloneBlock: MatrixBlock =
    this match
      case MatrixBlock.Dense(mat) =>
        MatrixBlock.Dense(mat.cloneMat)
      case MatrixBlock.Diagonal(values) =>
        MatrixBlock.Diagonal(values.clone())
      case MatrixBlock.BlockDiagonal(blocks) =>
        MatrixBlock.BlockDiagonal(blocks.map(_.cloneMat))

  def asDense: WorkMat =
    this match
      case MatrixBlock.Dense(mat) =>
        mat
      case MatrixBlock.Diagonal(values) =>
        val n = values.length
        val mat = WorkMat.zeros(n, n)
        var i = 0
        while i < n do
          mat(i, i) = values(i)
          i += 1
        mat
      case MatrixBlock.BlockDiagonal(blocks) =>
        val n = blocks.map(_.rows).sum
        val mat = WorkMat.zeros(n, n)
        var offset = 0
        var b = 0
        while b < blocks.length do
          val blk = blocks(b)
          var i = 0
          while i < blk.rows do
            var j = 0
            while j < blk.cols do
              mat(offset + i, offset + j) = blk(i, j)
              j += 1
            i += 1
          offset += blk.rows
          b += 1
        mat

object MatrixBlock:
  def copyBlock(dst: MatrixBlock, src: MatrixBlock): Unit =
    (dst, src) match
      case (MatrixBlock.Dense(d), MatrixBlock.Dense(s)) if d.rows == s.rows && d.cols == s.cols =>
        d.copyFrom(s)
      case (MatrixBlock.Diagonal(d), MatrixBlock.Diagonal(s)) if d.length == s.length =>
        System.arraycopy(s, 0, d, 0, d.length)
      case (MatrixBlock.BlockDiagonal(d), MatrixBlock.BlockDiagonal(s)) if d.length == s.length =>
        var i = 0
        while i < d.length do
          d(i).copyFrom(s(i))
          i += 1
      case _ =>
        throw IllegalArgumentException("MatrixBlock.copyBlock: incompatible block kinds")

  def blockIndex(i: Int, j: Int): Int =
    i * (i + 1) / 2 + j
