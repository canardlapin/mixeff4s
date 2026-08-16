package mixeff4s.data

import mixeff4s.error.MixedModelError

/** Column-oriented table used as the fitting substrate. */
enum Column:
  case Numeric(values: Vector[Double])
  case Factor(levels: Vector[String], refs: Vector[Int], values: Vector[String])

  def length: Int =
    this match
      case Column.Numeric(values)      => values.length
      case Column.Factor(_, _, values) => values.length

enum CategoricalCoding:
  case Treatment, CellMeans

final case class EncodedColumn(name: String, values: Vector[Double])

extension (factor: Column.Factor)
  /** Dummy columns for a factor. Treatment drops the first (reference) level. */
  def encodedColumns(variable: String, coding: CategoricalCoding): Vector[EncodedColumn] =
    val start = if coding == CategoricalCoding.Treatment then 1 else 0
    factor.levels.zipWithIndex
      .drop(start)
      .map: (level, idx) =>
        EncodedColumn(
          s"$variable: $level",
          factor.refs.map(ref => if ref == idx then 1.0 else 0.0)
        )

final class ModelFrame private (
    val columns: Vector[(String, Column)],
    val nRows: Int
):
  def names: Vector[String] = columns.map(_._1)

  def column(name: String): Option[Column] =
    columns.collectFirst { case (n, col) if n == name => col }

  def numeric(name: String): Option[Vector[Double]] =
    column(name).collect { case Column.Numeric(values) => values }

  def factor(name: String): Option[Column.Factor] =
    column(name).collect { case f: Column.Factor => f }

  def contains(name: String): Boolean = column(name).isDefined

  def addNumeric(name: String, values: Vector[Double]): Either[MixedModelError, ModelFrame] =
    if values.length != nRows then
      Left(
        MixedModelError.DimensionMismatch(
          s"column `$name` has ${values.length} rows, expected $nRows"
        )
      )
    else if contains(name) then
      Left(MixedModelError.InvalidArgument(s"duplicate column `$name`"))
    else Right(new ModelFrame(columns :+ (name -> Column.Numeric(values)), nRows))

object ModelFrame:
  def numeric(values: Iterable[Double]): Column.Numeric =
    Column.Numeric(values.toVector)

  def factor(values: Iterable[String]): Column.Factor =
    encodeFactor(values.toVector)

  def of(cols: (String, Column)*): Either[MixedModelError, ModelFrame] =
    if cols.isEmpty then Left(MixedModelError.InvalidArgument("ModelFrame requires at least one column"))
    else
      val n = cols.head._2.length
      val seen = scala.collection.mutable.Set.empty[String]
      var error: Option[MixedModelError] = None
      cols.foreach: (name, col) =>
        if error.isEmpty then
          if name.isBlank then error = Some(MixedModelError.InvalidArgument("column name must not be blank"))
          else if !seen.add(name) then error = Some(MixedModelError.InvalidArgument(s"duplicate column `$name`"))
          else if col.length != n then
            error = Some(
              MixedModelError.DimensionMismatch(
                s"column `$name` has ${col.length} rows, expected $n"
              )
            )
      error.toLeft(new ModelFrame(cols.toVector, n))

  private def encodeFactor(values: Vector[String]): Column.Factor =
    val levels = Vector.newBuilder[String]
    val index = scala.collection.mutable.LinkedHashMap.empty[String, Int]
    val refs = values.map: value =>
      index.getOrElseUpdate(
        value, {
          val id = index.size
          levels += value
          id
        }
      )
    Column.Factor(levels.result(), refs, values)
