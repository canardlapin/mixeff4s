package mixeff4s.glmm

import mixeff4s.design.{CompiledDesign, FixedEffectLayout}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.Formula
import mixeff4s.model.{Family, Link}

/** Frozen encoding with family/link identity. Each session has independent mutable scratch. */
final class CompiledGlmm private[glmm] (
    private[glmm] val design: CompiledDesign,
    val family: Family,
    val link: Link
):
  def formula: Formula = design.formula
  def n: Int = design.n
  def fixedEffects: FixedEffectLayout = design.fe.fixedEffects
  def fit(options: GlmmOptions = GlmmOptions.fastLaplace): FitResult[GlmmFit] =
    fit(Glmm.responseValues(design), options)
  def fit(response: Vector[Double], options: GlmmOptions): FitResult[GlmmFit] =
    session(options).flatMap(_.refit(response))
  def session(options: GlmmOptions = GlmmOptions.fastLaplace): FitResult[GlmmRefitSession] =
    for
      _ <- Glmm.validate(family, link, options)
      _ <-
        if design.prediction.exists(recipe => !recipe.matches(design)) then
          Left(MixedModelError.InvalidArgument("Compiled design does not match its frozen prediction recipe."))
        else Right(())
      response <- GlmmResponse.prepare(family, Glmm.responseValues(design), options)
      workspace <-
        if options.fast then
          Glmm.fastWorkspace(design, family, link, response, options).map(GlmmSessionWorkspace.Fast.apply)
        else JointWorkspace.create(design, family, link, response, options).map(GlmmSessionWorkspace.Joint.apply)
    yield new GlmmRefitSession(this, options, workspace)

  /** Every attempt is retained; a failed response cannot suppress later cold refits. */
  def fitResponses(
      responses: Vector[Vector[Double]],
      options: GlmmOptions = GlmmOptions.fastLaplace
  ): FitResult[Vector[FitResult[GlmmFit]]] = session(options).map(_.refitAll(responses))

private[glmm] enum GlmmSessionWorkspace:
  case Fast(workspace: GlmmWorkspace)
  case Joint(workspace: JointWorkspace)

/** Serialized response-only refits with the same starts as independent fresh fits. Family/link, trials, offsets,
  * weights, estimator and controls are fixed for the session.
  */
final class GlmmRefitSession private[glmm] (
    private val compiled: CompiledGlmm,
    val options: GlmmOptions,
    private val workspace: GlmmSessionWorkspace
):
  def refit(values: Vector[Double]): FitResult[GlmmFit] = synchronized:
    for
      _ <-
        if values.length != compiled.n then
          Left(MixedModelError.DimensionMismatch(s"response length ${values.length}, expected ${compiled.n}"))
        else Right(())
      response <- GlmmResponse.prepare(compiled.family, values, options)
      _ <-
        if response.y.nonEmpty && response.y.forall(_ == response.y.head) then Left(MixedModelError.ConstantResponse)
        else Right(())
      fit <- workspace match
        case GlmmSessionWorkspace.Fast(scratch) =>
          Glmm.fastFitWorkspace(compiled.design, compiled.family, compiled.link, response, options, scratch)
        case GlmmSessionWorkspace.Joint(scratch) =>
          JointLaplace.fitWorkspace(compiled.design, compiled.family, compiled.link, response, options, scratch)
    yield fit

  def refitAll(responses: Vector[Vector[Double]]): Vector[FitResult[GlmmFit]] = synchronized:
    responses.map(refit)
