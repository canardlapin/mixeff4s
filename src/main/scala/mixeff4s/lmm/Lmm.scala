package mixeff4s.lmm

import mixeff4s.data.ModelFrame
import mixeff4s.design.{CompiledDesign, Design}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.Formula

enum Criterion:
  case ML, REML

final case class FitOptions(criterion: Criterion)

object FitOptions:
  def ml: FitOptions = FitOptions(Criterion.ML)
  def reml: FitOptions = FitOptions(Criterion.REML)

/** Certified profiled LMM fit. */
final case class LmmFit(
    formula: Formula,
    options: FitOptions,
    theta: Vector[Double],
    beta: Vector[Double],
    sigma: Double,
    objective: Double
)

/** Linear mixed-model front door: compile a formula against a frame, then fit. */
object Lmm:
  def compile(formula: Formula, frame: ModelFrame): FitResult[CompiledDesign] =
    Design.compile(formula, frame)

  def compile(source: String, frame: ModelFrame): FitResult[CompiledDesign] =
    Formula.parse(source) match
      case Left(err)      => Left(MixedModelError.Formula(err))
      case Right(formula) => compile(formula, frame)

  def fit(
      formula: Formula,
      frame: ModelFrame,
      options: FitOptions = FitOptions.reml
  ): FitResult[LmmFit] =
    compile(formula, frame).flatMap: design =>
      Pls.fit(design, reml = options.criterion == Criterion.REML).map: workspace =>
        LmmFit(
          formula,
          options,
          workspace.theta,
          workspace.beta,
          workspace.sigma,
          workspace.objective
        )

  def fit(
      source: String,
      frame: ModelFrame,
      options: FitOptions
  ): FitResult[LmmFit] =
    Formula.parse(source) match
      case Left(err)      => Left(MixedModelError.Formula(err))
      case Right(formula) => fit(formula, frame, options)
