package mixeff4s.lmm

import mixeff4s.design.{CompiledDesign, FixedEffectLayout}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.Formula

/** A previous successful fit is used only when explicitly requested. */
enum RefitStart:
  case Cold, Previous

enum BatchStart:
  case Independent, WarmSequential

/** Immutable compiled plan. Each session owns its numerical workspace; fits are immutable snapshots. */
final class CompiledLmm private[lmm] (private val design: CompiledDesign):
  def formula: Formula = design.formula
  def n: Int = design.n
  def fixedEffects: FixedEffectLayout = design.fe.fixedEffects
  def fit(options: FitOptions = FitOptions.reml): FitResult[LmmFit] = Lmm.fit(design, options)
  def fit(response: Vector[Double], options: FitOptions): FitResult[LmmFit] =
    session(options).flatMap(_.refit(response))
  def session(options: FitOptions = FitOptions.reml): FitResult[LmmRefitSession] =
    options.validate
      .flatMap(_ => PlsWorkspace(design, options.criterion == Criterion.REML))
      .map(new LmmRefitSession(_, options))
  def fitResponses(
      responses: Vector[Vector[Double]],
      options: FitOptions = FitOptions.reml,
      start: BatchStart = BatchStart.Independent
  ): FitResult[Vector[LmmFit]] = session(options).flatMap(_.refitAll(responses, start))

/** Serialized response-only refits. A failed call cannot replace the previous successful warm-start state. The
  * workspace never escapes; prior fit values and design evidence are independent of later calls.
  */
final class LmmRefitSession private[lmm] (private val workspace: PlsWorkspace, val options: FitOptions):
  private var lastSuccessfulTheta: Option[Vector[Double]] = None

  def refit(response: Vector[Double], start: RefitStart = RefitStart.Cold): FitResult[LmmFit] = synchronized:
    val initial = start match
      case RefitStart.Cold     => Right(None)
      case RefitStart.Previous =>
        lastSuccessfulTheta match
          case Some(theta) => Right(Some(theta))
          case None        =>
            Left(MixedModelError.InvalidArgument("previous refit start requires a successful fit in this session"))
    for
      theta <- initial
      _ <- workspace.updateResponse(response)
      result <- Pls.optimize(
        workspace,
        options.maxEvaluations,
        options.optimizer,
        options.gradientOptions,
        theta,
        options.optimizerControls
      )
    yield
      val fit = Lmm.snapshot(result, options)
      lastSuccessfulTheta = Some(fit.theta)
      fit

  /** WarmSequential is explicitly order-dependent. The whole batch is serialized as one session operation. A failure
    * retains the warm-start state of the successful prefix; batches are not transactional.
    */
  def refitAll(
      responses: Vector[Vector[Double]],
      start: BatchStart = BatchStart.Independent
  ): FitResult[Vector[LmmFit]] = synchronized:
    responses.zipWithIndex.foldLeft[FitResult[Vector[LmmFit]]](Right(Vector.empty)): (acc, item) =>
      acc.flatMap: fits =>
        val mode = if start == BatchStart.WarmSequential && item._2 > 0 then RefitStart.Previous else RefitStart.Cold
        refit(item._1, mode).map(fits :+ _)
