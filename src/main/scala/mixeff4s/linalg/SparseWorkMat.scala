package mixeff4s.linalg

/** Mutable values over an immutable, sorted-unique CSC sparsity pattern.
  *
  * Coordinates supplied at construction are structural: a stored zero remains a slot and can subsequently receive a
  * non-zero value.
  */
final class SparseWorkMat private[linalg] (
    val rows: Int,
    val cols: Int,
    private[linalg] val colOffsets: Array[Int],
    private[linalg] val rowIndices: Array[Int],
    private[linalg] val values: Array[Double]
):
  private def offset(row: Int, col: Int): Int =
    if row < 0 || row >= rows || col < 0 || col >= cols then -1
    else
      var lo = colOffsets(col)
      var hi = colOffsets(col + 1) - 1
      while lo <= hi do
        val mid = (lo + hi) >>> 1
        val at = rowIndices(mid)
        if at < row then lo = mid + 1
        else if at > row then hi = mid - 1
        else return mid
      -1

  def apply(row: Int, col: Int): Double =
    val i = offset(row, col)
    if i < 0 then 0.0 else values(i)

  /** Updating a non-structural coordinate is a modelling error, even for zero. */
  def update(row: Int, col: Int, value: Double): Unit =
    val i = offset(row, col)
    if i < 0 then throw IllegalArgumentException(s"SparseWorkMat.update: missing structural slot ($row,$col)")
    values(i) = value

  def add(row: Int, col: Int, value: Double): Unit =
    val i = offset(row, col)
    if i < 0 then throw IllegalArgumentException(s"SparseWorkMat.add: missing structural slot ($row,$col)")
    values(i) += value

  def clear(): Unit = java.util.Arrays.fill(values, 0.0)

  def hasEntry(row: Int, col: Int): Boolean = offset(row, col) >= 0

  def foreachColumnEntry(col: Int)(f: (Int, Double) => Unit): Unit =
    require(col >= 0 && col < cols, s"SparseWorkMat column $col outside 0..${cols - 1}")
    var i = colOffsets(col)
    while i < colOffsets(col + 1) do
      f(rowIndices(i), values(i))
      i += 1

  def foreachStored(f: (Int, Int, Double) => Unit): Unit =
    var col = 0
    while col < cols do
      var i = colOffsets(col)
      while i < colOffsets(col + 1) do
        f(rowIndices(i), col, values(i))
        i += 1
      col += 1

  def cloneMat: SparseWorkMat =
    new SparseWorkMat(rows, cols, colOffsets, rowIndices, values.clone())

  def asDense: WorkMat =
    val result = WorkMat.zeros(rows, cols)
    foreachStored((row, col, value) => result(row, col) = value)
    result

object SparseWorkMat:
  def fromCoordinates(rows: Int, cols: Int, coordinates: Iterable[(Int, Int)]): SparseWorkMat =
    require(rows >= 0 && cols >= 0, s"SparseWorkMat dimensions must be non-negative, got ${rows}x$cols")
    val sorted = coordinates.iterator.toArray.distinct.sortBy { case (row, col) => (col, row) }
    sorted.foreach { case (row, col) =>
      require(
        row >= 0 && row < rows && col >= 0 && col < cols,
        s"SparseWorkMat coordinate ($row,$col) outside ${rows}x$cols"
      )
    }
    val offsets = new Array[Int](cols + 1)
    var col = 0
    var index = 0
    while col < cols do
      offsets(col) = index
      while index < sorted.length && sorted(index)._2 == col do index += 1
      col += 1
    offsets(cols) = index
    val rowIndices = new Array[Int](sorted.length)
    var i = 0
    while i < sorted.length do
      rowIndices(i) = sorted(i)._1
      i += 1
    new SparseWorkMat(rows, cols, offsets, rowIndices, new Array[Double](sorted.length))

  def zeros(rows: Int, cols: Int, coordinates: Iterable[(Int, Int)]): SparseWorkMat =
    fromCoordinates(rows, cols, coordinates)
