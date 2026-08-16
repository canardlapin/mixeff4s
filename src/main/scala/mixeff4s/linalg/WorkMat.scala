package mixeff4s.linalg

/** Mutable row-major work matrix for the blocked PLS hot path. */
final class WorkMat(val rows: Int, val cols: Int, val data: Array[Double]):
  def apply(row: Int, col: Int): Double =
    data(row * cols + col)

  def update(row: Int, col: Int, value: Double): Unit =
    data(row * cols + col) = value

  def fill(value: Double): Unit =
    java.util.Arrays.fill(data, value)

  def copyFrom(src: WorkMat): Unit =
    require(src.rows == rows && src.cols == cols, "WorkMat.copyFrom shape mismatch")
    System.arraycopy(src.data, 0, data, 0, data.length)

  def cloneMat: WorkMat =
    WorkMat(rows, cols, data.clone())

object WorkMat:
  def zeros(rows: Int, cols: Int): WorkMat =
    WorkMat(rows, cols, new Array[Double](rows * cols))

  def tabulate(rows: Int, cols: Int)(f: (Int, Int) => Double): WorkMat =
    val mat = zeros(rows, cols)
    var i = 0
    while i < rows do
      var j = 0
      while j < cols do
        mat(i, j) = f(i, j)
        j += 1
      i += 1
    mat
