package mixeff4s.formula

import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.error.MixedModelError

/** A parsed mixed-model formula: response, fixed terms, and random terms. */
final case class Formula(
    response: String,
    fixedTerms: Vector[FixedTerm],
    randomTerms: Vector[RandomTerm],
    derived: Vector[DerivedColumn] = Vector.empty
):
  def hasIntercept: Boolean =
    fixedTerms.exists(_ == FixedTerm.Intercept)

  def columnNames: Vector[String] =
    val fromFixed = fixedTerms.flatMap(_.columnNames)
    val fromRandom = randomTerms.flatMap(_.terms.flatMap(_.columnNames))
    val fromGrouping = randomTerms.flatMap(_.grouping.factorNames)
    (Vector(response) ++ fromFixed ++ fromRandom ++ fromGrouping).distinct

  /** Lower stateless transforms into numeric columns at the data boundary. */
  def materialize(frame: ModelFrame): Either[MixedModelError, ModelFrame] =
    derived.foldLeft[Either[MixedModelError, ModelFrame]](Right(frame)): (acc, derivedCol) =>
      acc.flatMap(current => materializeOne(current, derivedCol))

  override def toString: String =
    val rhs = (fixedTerms.map(_.toString) ++ randomTerms.map(_.toString)).mkString(" + ")
    s"$response ~ $rhs"

object Formula:
  def parse(input: String): Either[FormulaError, Formula] =
    Parser.parse(input)

enum FixedTerm:
  case Intercept
  case NoIntercept
  case Column(name: String)
  case Interaction(names: Vector[String])

  def columnNames: Vector[String] =
    this match
      case FixedTerm.Column(name)         => Vector(name)
      case FixedTerm.Interaction(names)   => names
      case FixedTerm.Intercept | FixedTerm.NoIntercept =>
        Vector.empty

  override def toString: String =
    this match
      case FixedTerm.Intercept          => "1"
      case FixedTerm.NoIntercept        => "0"
      case FixedTerm.Column(name)       => name
      case FixedTerm.Interaction(names) => names.mkString(":")

enum RandomCovariance:
  case Full
  case Diagonal
  case CompoundSymmetry
  case Ar1

  def label: String =
    this match
      case RandomCovariance.Full             => "full"
      case RandomCovariance.Diagonal         => "diagonal"
      case RandomCovariance.CompoundSymmetry => "compound_symmetry"
      case RandomCovariance.Ar1              => "ar1"

  def wrapper: Option[String] =
    this match
      case RandomCovariance.Full             => None
      case RandomCovariance.Diagonal         => Some("diag")
      case RandomCovariance.CompoundSymmetry => Some("cs")
      case RandomCovariance.Ar1              => Some("ar1")

  def isSupportedForFit: Boolean =
    this == RandomCovariance.Full || this == RandomCovariance.Diagonal

enum GroupingFactor:
  case Single(name: String)
  case Interaction(names: Vector[String])
  case Cell(names: Vector[String])

  def factorNames: Vector[String] =
    this match
      case GroupingFactor.Single(name)       => Vector(name)
      case GroupingFactor.Interaction(names) => names
      case GroupingFactor.Cell(names)        => names

  override def toString: String =
    this match
      case GroupingFactor.Single(name)       => name
      case GroupingFactor.Interaction(names) => names.mkString(" & ")
      case GroupingFactor.Cell(names)        => names.mkString(":")

enum RandomTermExpansion:
  case NestedGrouping
  case CrossedGrouping

final case class RandomTermSource(
    written: String,
    expansion: Option[RandomTermExpansion]
)

final case class RandomTerm(
    terms: Vector[FixedTerm],
    grouping: GroupingFactor,
    zerocorr: Boolean,
    covariance: RandomCovariance,
    source: Option[RandomTermSource] = None
):
  override def toString: String =
    val body = terms.map(_.toString).mkString(" + ")
    if zerocorr then s"($body || $grouping)"
    else
      covariance.wrapper match
        case Some(wrapper) => s"$wrapper($body | $grouping)"
        case None          => s"($body | $grouping)"

private def materializeOne(
    frame: ModelFrame,
    derived: DerivedColumn
): Either[MixedModelError, ModelFrame] =
  Transform.materializeColumn(derived, frame).flatMap: computed =>
    frame.column(derived.label) match
      case Some(Column.Numeric(existing)) =>
        if existing.length != computed.length then
          Left(
            MixedModelError.InvalidArgument(
              s"in-formula transform `${derived.label}`: a column with this name already exists in the data but has ${existing.length} rows, expected ${computed.length} — the engine owns this derived column; rename the raw column to avoid the collision"
            )
          )
        else
          val mismatch = existing.indices.find: row =>
            val supplied = existing(row)
            val engine = computed(row)
            val absDiff = math.abs(supplied - engine)
            val relDiff = if math.abs(engine) > 1e-12 then absDiff / math.abs(engine) else absDiff
            absDiff > 1e-12 && relDiff > 1e-10
          mismatch match
            case Some(row) =>
              val supplied = existing(row)
              val engine = computed(row)
              val relDiff =
                if math.abs(engine) > 1e-12 then math.abs(supplied - engine) / math.abs(engine)
                else math.abs(supplied - engine)
              Left(
                MixedModelError.InvalidArgument(
                  s"in-formula transform `${derived.label}` at row $row: the engine computed $engine but the pre-supplied column contains $supplied (relative diff $relDiff). The engine owns this derived-column recipe; there must be exactly one source of truth. Rename the raw column to avoid the collision, or remove it and let the engine compute it."
                )
              )
            case None =>
              Right(frame)
      case Some(_) =>
        Left(
          MixedModelError.InvalidArgument(
            s"in-formula transform `${derived.label}`: a non-numeric column with this name already exists"
          )
        )
      case None =>
        frame.addNumeric(derived.label, computed)
