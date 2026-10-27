package mixeff4s.design

import gale.backend.PureBackend
import gale.linalg.*
import mixeff4s.error.MixedModelError

/** Fixed-effects design with pivoted-QR rank. Columns are stored in pivot order. */
final case class FeTerm(
    x: DMat,
    piv: Vector[Int],
    rank: Int,
    cnames: Vector[String]
):
  def n: Int = x.rows
  def p: Int = x.cols

  def fullRankX: DMat =
    if rank == p then x
    else Matrix.tabulate(n, rank)((i, j) => x(i, j))

  def fullRankNames: Vector[String] =
    cnames.take(rank)

  /** Recover requested column order from the working-column permutation. */
  def fixedEffects: FixedEffectLayout =
    require(cnames.length == p && piv.sorted == Vector.range(0, p), "invalid fixed-effect column permutation")
    require(rank >= 0 && rank <= p, "fixed-effect rank out of bounds")
    val requestedNames = cnames.zip(piv).sortBy(_._2).map(_._1)
    FixedEffectLayout(requestedNames, piv.take(rank))

object FeTerm:
  private val RankTolerance = 1e-8

  def apply(x: DMat, cnames: Vector[String]): Either[MixedModelError, FeTerm] =
    if cnames.length != x.cols then
      Left(
        MixedModelError.DimensionMismatch(
          s"FeTerm: ${cnames.length} names for ${x.cols} columns"
        )
      )
    else if x.cols == 0 then Right(FeTerm(x, Vector.empty, 0, cnames))
    else
      val qr = x.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(RankTolerance)))(using
        PureBackend
      )
      val rawPiv = Vector.tabulate(x.cols)(qr.columnPermutation.apply)
      val rank = qr.diagnostics.rank.getOrElse(x.cols).min(x.cols)
      val kept = rawPiv.take(rank).sorted
      val dropped = rawPiv.drop(rank)
      // Preserve earlier estimable columns, as R's limited-pivot QR does, only
      // when that selection agrees with the rank from the pivoted QR above.
      val piv = if rank < x.cols then inOrderPivot(x, rank).getOrElse(kept ++ dropped) else kept ++ dropped
      val pivoted = Matrix.tabulate(x.rows, x.cols)((i, j) => x(i, piv(j)))
      val names = piv.map(cnames)
      Right(FeTerm(pivoted, piv, rank, names))

  private def inOrderPivot(x: DMat, rank: Int): Option[Vector[Int]] =
    val basis = scala.collection.mutable.ArrayBuffer.empty[Array[Double]]
    val kept = Vector.newBuilder[Int]
    val dropped = Vector.newBuilder[Int]
    var j = 0
    while j < x.cols do
      val residual = Array.tabulate(x.rows)(i => x(i, j))
      val originalNorm = residual.foldLeft(0.0)(math.hypot)
      if !originalNorm.isFinite then return None
      var pass = 0
      while pass < 2 do
        basis.foreach: q =>
          var projection = 0.0
          var i = 0
          while i < x.rows do
            projection += q(i) * residual(i)
            i += 1
          i = 0
          while i < x.rows do
            residual(i) -= projection * q(i)
            i += 1
        pass += 1
      val norm = residual.foldLeft(0.0)(math.hypot)
      if !norm.isFinite then return None
      if originalNorm > 0.0 && norm > 1e-7 * originalNorm && basis.length < x.rows then
        var i = 0
        while i < x.rows do
          residual(i) /= norm
          i += 1
        basis += residual
        kept += j
      else dropped += j
      j += 1
    if basis.length == rank then Some(kept.result() ++ dropped.result()) else None

/** Rank-truncated `[X | y]`. */
final case class FeMat(xy: DMat):
  def n: Int = xy.rows
  def rank: Int = xy.cols - 1

object FeMat:
  def apply(fe: FeTerm, y: Vector[Double]): Either[MixedModelError, FeMat] =
    if y.length != fe.n then Left(MixedModelError.DimensionMismatch(s"FeMat: y has ${y.length} rows, expected ${fe.n}"))
    else
      val x = fe.fullRankX
      Right(FeMat(Matrix.tabulate(fe.n, fe.rank + 1)((i, j) => if j < fe.rank then x(i, j) else y(i))))
