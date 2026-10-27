package mixeff4s.linalg

/** One block of the mixed-model A or L array. Contents are mutated in place. */
enum MatrixBlock:
  case Dense(mat: WorkMat)
  case Diagonal(values: Array[Double])
  case BlockDiagonal(blocks: Array[WorkMat])

  /** Sparse storage is reserved for off-diagonal blocks. */
  case Sparse(mat: SparseWorkMat)

  lazy val rows: Int = this match
    case MatrixBlock.Dense(mat)        => mat.rows
    case MatrixBlock.Diagonal(values)  => values.length
    case MatrixBlock.BlockDiagonal(bs) => bs.map(_.rows).sum
    case MatrixBlock.Sparse(mat)       => mat.rows

  lazy val cols: Int = this match
    case MatrixBlock.Dense(mat)        => mat.cols
    case MatrixBlock.Diagonal(values)  => values.length
    case MatrixBlock.BlockDiagonal(bs) => bs.map(_.cols).sum
    case MatrixBlock.Sparse(mat)       => mat.cols

  private lazy val uniformBlockSize: Int = this match
    case MatrixBlock.BlockDiagonal(blocks) if blocks.nonEmpty =>
      val size = blocks(0).rows
      if size > 0 && blocks.forall(b => b.rows == size && b.cols == size) then size else -1
    case _ => -1

  def apply(row: Int, col: Int): Double = this match
    case MatrixBlock.Dense(mat)            => mat(row, col)
    case MatrixBlock.Diagonal(values)      => if row == col then values(row) else 0.0
    case MatrixBlock.BlockDiagonal(blocks) =>
      val size = uniformBlockSize
      if size > 0 then
        val index = row / size
        if row >= 0 && col >= 0 && index < blocks.length && col / size == index then
          blocks(index)(row % size, col % size)
        else 0.0
      else MatrixBlock.blockDiagonalEntry(blocks, row, col)
    case MatrixBlock.Sparse(mat) => mat(row, col)

  /** dst += scale * this * x, or dst += scale * this' * x. */
  def addMatVec(x: Array[Double], dst: Array[Double], transpose: Boolean = false, scale: Double = 1.0): Unit =
    val in = if transpose then rows else cols
    val out = if transpose then cols else rows
    require(x.length == in && dst.length == out, s"MatrixBlock.addMatVec dimension mismatch")
    this match
      case MatrixBlock.Sparse(mat) =>
        var col = 0
        while col < mat.cols do
          mat.foreachColumnEntry(col) { (row, value) =>
            if transpose then dst(col) += scale * value * x(row)
            else dst(row) += scale * value * x(col)
          }
          col += 1
      case MatrixBlock.Dense(mat) =>
        var row = 0
        while row < mat.rows do
          val base = row * mat.cols
          var col = 0
          while col < mat.cols do
            val value = mat.data(base + col)
            if transpose then dst(col) += scale * value * x(row)
            else dst(row) += scale * value * x(col)
            col += 1
          row += 1
      case MatrixBlock.Diagonal(values) =>
        var i = 0
        while i < values.length do
          dst(i) += scale * values(i) * x(i)
          i += 1
      case MatrixBlock.BlockDiagonal(blocks) =>
        var offset = 0
        var b = 0
        while b < blocks.length do
          val block = blocks(b)
          var row = 0
          while row < block.rows do
            val base = row * block.cols
            var col = 0
            while col < block.cols do
              val value = block.data(base + col)
              if transpose then dst(offset + col) += scale * value * x(offset + row)
              else dst(offset + row) += scale * value * x(offset + col)
              col += 1
            row += 1
          offset += block.rows
          b += 1

  def cloneBlock: MatrixBlock =
    this match
      case MatrixBlock.Dense(mat) =>
        MatrixBlock.Dense(mat.cloneMat)
      case MatrixBlock.Diagonal(values) =>
        MatrixBlock.Diagonal(values.clone())
      case MatrixBlock.BlockDiagonal(blocks) =>
        MatrixBlock.BlockDiagonal(blocks.map(_.cloneMat))
      case MatrixBlock.Sparse(mat) =>
        MatrixBlock.Sparse(mat.cloneMat)

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
      case MatrixBlock.Sparse(mat) => mat.asDense

object MatrixBlock:
  private[linalg] def blockDiagonalEntry(blocks: Array[WorkMat], row: Int, col: Int): Double =
    if blocks.isEmpty then 0.0
    else
      val first = blocks(0)
      var uniform = true
      var b = 1
      while b < blocks.length && uniform do
        uniform = blocks(b).rows == first.rows && blocks(b).cols == first.cols
        b += 1
      if uniform && first.rows == first.cols then
        val block = row / first.rows
        if block >= 0 && block < blocks.length && col / first.cols == block then
          blocks(block)(row % first.rows, col % first.cols)
        else 0.0
      else
        var offset = 0
        b = 0
        while b < blocks.length do
          val block = blocks(b)
          if row >= offset && row < offset + block.rows && col >= offset && col < offset + block.cols then
            return block(row - offset, col - offset)
          offset += block.rows
          b += 1
        0.0
  def copyBlock(dst: MatrixBlock, src: MatrixBlock): Unit =
    (dst, src) match
      case (MatrixBlock.Dense(d), MatrixBlock.Dense(s)) if d.rows == s.rows && d.cols == s.cols => d.copyFrom(s)
      case (MatrixBlock.Dense(d), src) if d.rows == src.rows && d.cols == src.cols              =>
        d.fill(0.0)
        src match
          case MatrixBlock.Sparse(s) => s.foreachStored((row, col, value) => d(row, col) = value)
          case _                     =>
            var row = 0
            while row < d.rows do
              var col = 0
              while col < d.cols do
                d(row, col) = src(row, col)
                col += 1
              row += 1
      case (MatrixBlock.Diagonal(d), MatrixBlock.Diagonal(s)) if d.length == s.length =>
        System.arraycopy(s, 0, d, 0, d.length)
      case (MatrixBlock.BlockDiagonal(d), MatrixBlock.BlockDiagonal(s)) if d.length == s.length =>
        var i = 0
        while i < d.length do
          d(i).copyFrom(s(i))
          i += 1
      case (MatrixBlock.Sparse(d), MatrixBlock.Sparse(s)) if d.rows == s.rows && d.cols == s.cols =>
        val compatible = d.colOffsets.sameElements(s.colOffsets) && d.rowIndices.sameElements(s.rowIndices)
        if !compatible then throw IllegalArgumentException("MatrixBlock.copyBlock: incompatible sparse topology")
        System.arraycopy(s.values, 0, d.values, 0, d.values.length)
      case _ =>
        throw IllegalArgumentException("MatrixBlock.copyBlock: incompatible block kinds")

  def blockIndex(i: Int, j: Int): Int =
    i * (i + 1) / 2 + j
