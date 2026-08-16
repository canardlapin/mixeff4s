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
      val qr = x.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(RankTolerance)))(
        using PureBackend
      )
      val rawPiv = Vector.tabulate(x.cols)(qr.columnPermutation.apply)
      val rank = qr.diagnostics.rank.getOrElse(x.cols).min(x.cols)
      val kept = rawPiv.take(rank).sorted
      val dropped = rawPiv.drop(rank)
      val piv = kept ++ dropped
      val pivoted = Matrix.tabulate(x.rows, x.cols)((i, j) => x(i, piv(j)))
      val names = piv.map(cnames)
      Right(FeTerm(pivoted, piv, rank, names))

/** Rank-truncated `[X | y]`. */
final case class FeMat(xy: DMat):
  def n: Int = xy.rows
  def rank: Int = xy.cols - 1

object FeMat:
  def apply(fe: FeTerm, y: Vector[Double]): Either[MixedModelError, FeMat] =
    if y.length != fe.n then
      Left(MixedModelError.DimensionMismatch(s"FeMat: y has ${y.length} rows, expected ${fe.n}"))
    else
      val x = fe.fullRankX
      Right(FeMat(Matrix.tabulate(fe.n, fe.rank + 1)((i, j) => if j < fe.rank then x(i, j) else y(i))))
