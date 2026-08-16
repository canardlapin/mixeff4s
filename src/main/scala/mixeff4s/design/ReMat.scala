package mixeff4s.design

import gale.linalg.*
import gale.sparse.{CSC, Sparse}
import mixeff4s.error.MixedModelError

/** One random-effects term. `z` is stored transposed: `vsize × n`. */
final case class ReMat(
    groupingName: String,
    refs: Vector[Int],
    levels: Vector[String],
    cnames: Vector[String],
    z: DMat,
    lambda: DMat,
    inds: Vector[Int],
    adjA: CSC
):
  def vsize: Int = cnames.length
  def nLevels: Int = levels.length
  def nObs: Int = refs.length
  def nRanef: Int = vsize * nLevels
  def nTheta: Int = inds.length

  def theta: Vector[Double] =
    inds.map: idx =>
      val (row, col) = ReMat.linearToSubscript(idx, vsize)
      lambda(row, col)

  def withTheta(values: Vector[Double]): Either[MixedModelError, ReMat] =
    if values.length != nTheta then
      Left(MixedModelError.DimensionMismatch(s"ReMat.withTheta expected $nTheta values, got ${values.length}"))
    else
      val next = Matrix.tabulate(vsize, vsize)((r, c) => lambda(r, c))
      val builder = Matrix.builderFrom(next)
      values.zip(inds).foreach: (value, idx) =>
        val (row, col) = ReMat.linearToSubscript(idx, vsize)
        builder(row, col) = value
      Right(copy(lambda = builder.result()))

  def zerocorr: ReMat =
    val cleared = Matrix.tabulate(vsize, vsize)((r, c) => if r == c then lambda(r, c) else 0.0)
    copy(lambda = cleared, inds = ReMat.diagonalIndices(vsize))

object ReMat:
  def apply(
      groupingName: String,
      refs: Vector[Int],
      levels: Vector[String],
      cnames: Vector[String],
      z: DMat
  ): Either[MixedModelError, ReMat] =
    val vsize = cnames.length
    val nObs = refs.length
    if z.rows != vsize then
      Left(MixedModelError.DimensionMismatch(s"ReMat: z.rows=${z.rows} must equal vsize=$vsize"))
    else if z.cols != nObs then
      Left(MixedModelError.DimensionMismatch(s"ReMat: z.cols=${z.cols} must equal n=$nObs"))
    else if refs.exists(r => r < 0 || r >= levels.length) then
      Left(MixedModelError.InvalidArgument("ReMat: grouping ref is out of bounds"))
    else
      sparseAdjoint(refs, z, vsize, levels.length, nObs).map: adj =>
        new ReMat(
          groupingName,
          refs,
          levels,
          cnames,
          z,
          Matrix.eye(vsize),
          lowerTriangularIndices(vsize),
          adj
        )

  def buildParmap(reterms: Vector[ReMat]): Vector[(Int, Int, Int)] =
    reterms.zipWithIndex.flatMap: (rt, block) =>
      rt.inds.map: ind =>
        val col = ind / rt.vsize
        val row = ind % rt.vsize
        (block, row, col)

  private def lowerTriangularIndices(s: Int): Vector[Int] =
    (0 until s).flatMap(col => (col until s).map(row => col * s + row)).toVector

  private def diagonalIndices(s: Int): Vector[Int] =
    Vector.tabulate(s)(k => k * s + k)

  private def linearToSubscript(idx: Int, nrows: Int): (Int, Int) =
    (idx % nrows, idx / nrows)

  private def sparseAdjoint(
      refs: Vector[Int],
      z: DMat,
      vsize: Int,
      nLevels: Int,
      nObs: Int
  ): Either[MixedModelError, CSC] =
    val builder = Sparse.coo(vsize * nLevels, nObs)
    refs.zipWithIndex.foreach: (level, obs) =>
      val rowStart = level * vsize
      (0 until vsize).foreach: r =>
        val value = z(r, obs)
        if value != 0.0 then
          val _ = builder.add(rowStart + r, obs, value)
    Right(builder.toCOO().toCSC)
