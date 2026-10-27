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
    else if values.exists(value => !value.isFinite) then
      Left(MixedModelError.InvalidArgument(s"numeric column `$name` contains a non-finite value"))
    else if contains(name) then Left(MixedModelError.InvalidArgument(s"duplicate column `$name`"))
    else Right(new ModelFrame(columns :+ (name -> Column.Numeric(values)), nRows))

object ModelFrame:
  def numeric(values: Iterable[Double]): Column.Numeric =
    Column.Numeric(values.toVector)

  def factor(values: Iterable[String]): Column.Factor =
    encodeFactor(values.toVector, None).fold(
      err => throw IllegalStateException(err.message),
      identity
    )

  /** Encode a factor with an explicit level order so dummy coding matches a manifest. */
  def factor(values: Iterable[String], levels: Vector[String]): Either[MixedModelError, Column.Factor] =
    encodeFactor(values.toVector, Some(levels))

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
          else validateColumn(name, col).left.foreach(err => error = Some(err))
      error.toLeft(new ModelFrame(cols.toVector, n))

  private def validateColumn(name: String, col: Column): Either[MixedModelError, Unit] =
    col match
      case Column.Numeric(values) =>
        if values.exists(value => !value.isFinite) then
          Left(MixedModelError.InvalidArgument(s"numeric column `$name` contains a non-finite value"))
        else Right(())
      case Column.Factor(levels, refs, values) =>
        if refs.length != values.length then
          Left(
            MixedModelError.DimensionMismatch(
              s"factor column `$name` has ${refs.length} refs for ${values.length} values"
            )
          )
        else if levels.distinct.length != levels.length then
          Left(MixedModelError.InvalidArgument(s"factor column `$name` has duplicate levels"))
        else
          refs
            .zip(values)
            .zipWithIndex
            .collectFirst:
              case ((ref, _), row) if ref < 0 || ref >= levels.length =>
                MixedModelError.InvalidArgument(s"factor column `$name` has out-of-range ref $ref at row $row")
              case ((ref, value), row) if levels(ref) != value =>
                MixedModelError.InvalidArgument(
                  s"factor column `$name` value `$value` at row $row does not match level `${levels(ref)}`"
                )
          match
            case Some(error) => Left(error)
            case None        => Right(())

  private def encodeFactor(
      values: Vector[String],
      levelsOpt: Option[Vector[String]]
  ): Either[MixedModelError, Column.Factor] =
    levelsOpt match
      case None =>
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
        Right(Column.Factor(levels.result(), refs, values))
      case Some(levels) =>
        if levels.isEmpty then Left(MixedModelError.InvalidArgument("factor levels must not be empty"))
        else if levels.distinct.length != levels.length then
          Left(MixedModelError.InvalidArgument("factor levels must be unique"))
        else
          val index = levels.zipWithIndex.toMap
          val refs = Vector.newBuilder[Int]
          var error: Option[MixedModelError] = None
          values.foreach: value =>
            if error.isEmpty then
              index.get(value) match
                case Some(id) => refs += id
                case None     =>
                  error = Some(MixedModelError.InvalidArgument(s"value `$value` is not in the declared factor levels"))
          error.toLeft(Column.Factor(levels, refs.result(), values))
