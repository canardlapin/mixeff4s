package mixeff4s.compiler

import mixeff4s.data.ModelFrame
import mixeff4s.design.{CompiledDesign, Design}
import mixeff4s.error.FitResult
import mixeff4s.formula.{Formula, FormulaError}
import mixeff4s.pathology.Pathology

/** Front door for the unstable compiled-design artifact. */
object Compiler:
  def artifact(design: CompiledDesign): CompiledArtifact =
    CompiledArtifact.fromDesign(design)

  def semantic(formula: Formula): SemanticModel =
    SemanticModel.from(formula)

  def semantic(source: String): Either[FormulaError, SemanticModel] =
    Formula.parse(source).map(semantic)

  def compile(formula: Formula, frame: ModelFrame): FitResult[CompiledArtifact] =
    Design.compile(formula, frame).map(artifact)

  def compile(source: String, frame: ModelFrame): FitResult[CompiledArtifact] =
    Formula.parse(source) match
      case Left(err)      => Left(mixeff4s.error.MixedModelError.Formula(err))
      case Right(formula) => compile(formula, frame)

  def assess(artifact: CompiledArtifact, theta: Vector[Double]): CompiledArtifact =
    artifact.copy(pathology = Pathology.assessFit(artifact.pathology, theta, artifact.parmap))
