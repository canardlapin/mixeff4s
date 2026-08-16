package mixeff4s.lmm

import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.{Formula, RandomCovariance}

enum Criterion:
  case ML, REML

final case class FitOptions(criterion: Criterion)

object FitOptions:
  def ml: FitOptions = FitOptions(Criterion.ML)
  def reml: FitOptions = FitOptions(Criterion.REML)

/** Placeholder for a certified LMM fit. The PLS kernel is not implemented yet. */
final case class LmmFit(formula: Formula, options: FitOptions)

/** Linear mixed-model front door: compile a formula against a frame, then fit. */
object Lmm:
  def fit(
      formula: Formula,
      frame: ModelFrame,
      options: FitOptions = FitOptions.reml
  ): FitResult[LmmFit] =
    for
      materialized <- formula.materialize(frame)
      _ <- validate(formula, materialized)
    yield
      // Unreachable until the kernel lands; validate currently ends in Unsupported.
      LmmFit(formula, options)

  def fit(
      source: String,
      frame: ModelFrame,
      options: FitOptions
  ): FitResult[LmmFit] =
    Formula.parse(source) match
      case Left(err)      => Left(MixedModelError.Formula(err))
      case Right(formula) => fit(formula, frame, options)

  private def validate(formula: Formula, frame: ModelFrame): FitResult[Unit] =
    if formula.randomTerms.isEmpty then Left(MixedModelError.NoRandomEffects)
    else
      missingColumn(formula, frame)
        .orElse(wrongResponseType(formula, frame))
        .orElse(refusedCovariance(formula))
        .getOrElse(
          Left(
            MixedModelError.Unsupported(
              "LMM blocked-Cholesky PLS kernel is not implemented yet"
            )
          )
        )

  private def missingColumn(formula: Formula, frame: ModelFrame): Option[FitResult[Unit]] =
    formula.columnNames.find(name => !frame.contains(name)).map: name =>
      Left(MixedModelError.InvalidArgument(s"column `$name` is not present in the model frame"))

  private def wrongResponseType(formula: Formula, frame: ModelFrame): Option[FitResult[Unit]] =
    frame.column(formula.response) match
      case Some(Column.Numeric(_)) => None
      case Some(_) =>
        Some(
          Left(
            MixedModelError.InvalidArgument(
              s"response `${formula.response}` must be numeric"
            )
          )
        )
      case None =>
        Some(
          Left(
            MixedModelError.InvalidArgument(
              s"column `${formula.response}` is not present in the model frame"
            )
          )
        )

  private def refusedCovariance(formula: Formula): Option[FitResult[Unit]] =
    formula.randomTerms.collectFirst:
      case term if !term.covariance.isSupportedForFit =>
        val label = term.covariance match
          case RandomCovariance.CompoundSymmetry => "cs(...)"
          case RandomCovariance.Ar1              => "ar1(...)"
          case other                             => other.label
        Left(
          MixedModelError.Unsupported(
            s"$label random-effect covariance is parsed and refused for fitting"
          )
        )
