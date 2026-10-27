package mixeff4s.lmm

import mixeff4s.data.ModelFrame
import mixeff4s.design.{CompiledDesign, Design, FixedEffectLayout}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.Formula
import mixeff4s.optimizer.{TrustBqOptions, TrustBqResult}

enum Criterion:
  case ML, REML

enum LmmOptimizer:
  case DerivativeFree, AnalyticGradient

/** maxEvaluations bounds the main optimizer, excluding the advisory diagonal start. */
final case class FitOptions(
    criterion: Criterion,
    maxEvaluations: Int = 1000,
    optimizer: LmmOptimizer = LmmOptimizer.DerivativeFree,
    gradientOptions: ProfiledGradientOptions = ProfiledGradientOptions(),
    optimizerControls: Option[TrustBqOptions] = None
):
  /** Explicit controls replace the main-loop family policy. maxEvaluations remains authoritative; the advisory
    * diagonal-start policy is unchanged. Tighter controls do not guarantee assessment success.
    */
  def validate: FitResult[Unit] =
    if maxEvaluations <= 0 then Left(MixedModelError.InvalidArgument("maxEvaluations must be positive"))
    else
      optimizerControls match
        case None    => Right(())
        case Some(c) =>
          val finite = Vector(
            c.initialRadius,
            c.finalRadius,
            c.ftolAbs,
            c.ftolRel,
            c.stallFtolAbs,
            c.stallFtolRel,
            c.etaAccept,
            c.etaExpand,
            c.shrinkFactor,
            c.expandFactor,
            c.gradientTolerance
          ).forall(_.isFinite)
          if !finite || c.initialRadius <= 0 || c.finalRadius <= 0 || c.finalRadius > c.initialRadius ||
            c.ftolAbs < 0 || c.ftolRel < 0 || c.stallFtolAbs < -1 || c.stallFtolRel < -1 || c.stallIterations <= 0 ||
            c.maxCrossTerms < 0 || c.etaAccept < 0 || c.etaAccept >= 1 || c.etaExpand < c.etaAccept || c.etaExpand > 1 ||
            c.shrinkFactor < 0 || c.shrinkFactor >= 1 || c.expandFactor <= 1 || c.gradientTolerance < 0
          then Left(MixedModelError.InvalidArgument("Invalid explicit TrustBQ controls."))
          else Right(())

object FitOptions:
  def ml: FitOptions = FitOptions(Criterion.ML)
  def reml: FitOptions = FitOptions(Criterion.REML)

/** Profiled LMM fit with recorded numerical termination; stationarity is not certified. */
final case class LmmFit(
    formula: Formula,
    options: FitOptions,
    n: Int,
    y: Vector[Double],
    theta: Vector[Double],
    beta: Vector[Double],
    sigma: Double,
    objective: Double,
    fixedEffects: FixedEffectLayout,
    vcov: Vector[Vector[Double]],
    varcorr: VarCorr,
    logdetRe: Double,
    optimizer: TrustBqResult,
    parmap: Vector[(Int, Int, Int)],
    designEvidence: Option[LmmDesignEvidence] = None,
    private[mixeff4s] val workflow: Option[LmmFitContext] = None
):
  require(beta.length == fixedEffects.rank, "beta length must match the retained fixed-effect layout")
  require(
    vcov.length == beta.length && vcov.forall(_.length == beta.length),
    "vcov must be square in retained fixed-effect order"
  )

  def feNames: Vector[String] = fixedEffects.retainedNames

  def dof: Int = beta.length + theta.length + 1

  def loglikelihood: Double = -0.5 * objective

  def stderror: Vector[Double] =
    vcov.zipWithIndex.map((row, i) => math.sqrt(row(i)))

  def coefTable: CoefTable =
    CoefTable.wald(feNames, beta, stderror)

  private[mixeff4s] def validatedContext: FitResult[LmmFitContext] =
    workflow match
      case Some(context) if context.matches(this) => Right(context)
      case _                                      =>
        Left(
          MixedModelError.InferenceUnavailable(
            "fit_snapshot_mismatch",
            "The fit has no matching immutable fitted context."
          )
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
    compile(formula, frame).flatMap(design => fit(design, options))

  /** Fit an already compiled design without parsing, encoding or rank analysis again. */
  def fit(design: CompiledDesign, options: FitOptions): FitResult[LmmFit] =
    if design.prediction.exists(encoding => !encoding.matches(design)) then
      Left(MixedModelError.InvalidArgument("Compiled design does not match its frozen prediction recipe."))
    else
      options.validate.flatMap { _ =>
        Pls
          .fit(
            design,
            options.criterion == Criterion.REML,
            options.maxEvaluations,
            options.optimizer,
            options.gradientOptions,
            options.optimizerControls
          )
          .map(result => snapshot(result, options))
      }

  def prepare(formula: Formula, frame: ModelFrame): FitResult[CompiledLmm] =
    compile(formula, frame).map(new CompiledLmm(_))

  def prepare(source: String, frame: ModelFrame): FitResult[CompiledLmm] =
    compile(source, frame).map(new CompiledLmm(_))

  private[lmm] def snapshot(result: PlsFit, options: FitOptions): LmmFit =
    val workspace = result.workspace
    val design = workspace.design
    val response = workspace.response
    val fit = LmmFit(
      design.formula,
      options,
      design.n,
      response,
      workspace.theta,
      workspace.beta,
      workspace.sigma,
      workspace.objective,
      design.fe.fixedEffects,
      workspace.vcov,
      workspace.varcorr,
      workspace.logdetRe,
      result.optimizer,
      design.parmap,
      Some(LmmDesignEvidence.capture(design, options.criterion, workspace.theta, workspace.objective, response))
    )
    val modes = workspace.ranefU
      .zip(workspace.fittedReterms)
      .map: (u, term) =>
        Vector.tabulate(term.nRanef): index =>
          val row = index % term.vsize
          val start = index - row
          var value = 0.0
          var col = 0
          while col <= row do
            value += term.lambda(row, col) * u(start + col)
            col += 1
          value
    fit.copy(workflow = Some(new LmmFitContext(design, fit, modes)))

  def fit(
      source: String,
      frame: ModelFrame,
      options: FitOptions
  ): FitResult[LmmFit] =
    Formula.parse(source) match
      case Left(err)      => Left(MixedModelError.Formula(err))
      case Right(formula) => fit(formula, frame, options)
