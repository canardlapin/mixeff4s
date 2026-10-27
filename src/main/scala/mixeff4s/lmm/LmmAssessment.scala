package mixeff4s.lmm

import gale.linalg.Matrix
import gale.spectral.{Eigen, EigenSelection, EigenVectors}
import mixeff4s.design.ReMat
import mixeff4s.error.FitResult
import mixeff4s.optimizer.TrustBqStopReason

/** Options for an exact, post-fit LMM assessment. They do not change the fit. */
final case class LmmAssessmentOptions(
    gradientOptions: ProfiledGradientOptions = ProfiledGradientOptions(),
    objectiveRelativeTolerance: Double = 1e-9,
    covarianceRankRelativeTolerance: Double = 1e-12,
    covarianceRankAbsoluteTolerance: Double = 0.0,
    boundaryTolerance: Double = 1e-8,
    boundaryProbe: Double = 1e-4
):
  def validate: FitResult[Unit] =
    if objectiveRelativeTolerance < 0.0 || !objectiveRelativeTolerance.isFinite then
      Left(mixeff4s.error.MixedModelError.InvalidArgument("objectiveRelativeTolerance must be finite and nonnegative"))
    else if covarianceRankRelativeTolerance < 0.0 || !covarianceRankRelativeTolerance.isFinite then
      Left(
        mixeff4s.error.MixedModelError.InvalidArgument("covarianceRankRelativeTolerance must be finite and nonnegative")
      )
    else if covarianceRankAbsoluteTolerance < 0.0 || !covarianceRankAbsoluteTolerance.isFinite then
      Left(
        mixeff4s.error.MixedModelError.InvalidArgument("covarianceRankAbsoluteTolerance must be finite and nonnegative")
      )
    else if boundaryTolerance < 0.0 || !boundaryTolerance.isFinite then
      Left(mixeff4s.error.MixedModelError.InvalidArgument("boundaryTolerance must be finite and nonnegative"))
    else if boundaryProbe <= 0.0 || !boundaryProbe.isFinite then
      Left(mixeff4s.error.MixedModelError.InvalidArgument("boundaryProbe must be finite and positive"))
    else Right(())

  def objectiveTolerance(objective: Double): Double = objectiveRelativeTolerance * math.max(1.0, math.abs(objective))

enum GradientCheck:
  case Pass(projectedNorm: Double, tolerance: Double)
  case Failed(projectedNorm: Double, tolerance: Double)
  case Unavailable(reasonCode: String)

enum BoundaryAdequacy:
  case NotApplicable
  case Failed(parameter: Int, step: Double, objective: Double)
  case Unresolved(reasonCode: String)
  case Unavailable(reasonCode: String)

/** Eigenvalues are for the relative covariance Lambda Lambda', not scaled by sigma squared. */
final case class RandomEffectCovarianceRank(
    term: Int,
    eigenvalues: Vector[Double],
    rank: Option[Int],
    tolerance: Double,
    unavailableReason: Option[String]
)

/** Immutable post-fit evidence. A pass is first-order evidence, not a global-optimality claim. */
final case class LmmAssessment(
    numericalTermination: TrustBqStopReason,
    numericalTerminationAccepted: Boolean,
    gradient: GradientCheck,
    boundaryAdequacy: BoundaryAdequacy,
    covarianceRanks: Vector[RandomEffectCovarianceRank],
    reasonCodes: Vector[String],
    bootstrapIneligibilityReason: Option[String],
    eligibleForBootstrap: Boolean
)

object LmmAssessment:
  def assess(fit: LmmFit, options: LmmAssessmentOptions = LmmAssessmentOptions()): LmmAssessment =
    options.validate match
      case Left(error) => unavailable(fit, s"assessment_options_${error.code}")
      case Right(_)    =>
        fit.validatedContext match
          case Left(error)    => unavailable(fit, error.code)
          case Right(context) =>
            val ranks = covarianceRanks(context.design.reterms, fit.theta, options)
            freshEvaluation(context, fit, options) match
              case Left(code)                   => unavailable(fit, code, ranks)
              case Right((objective, gradient)) =>
                val tolerance = options.objectiveTolerance(fit.objective)
                if math.abs(objective - fit.objective) > tolerance then
                  unavailable(fit, "fit_objective_mismatch", ranks)
                else
                  val projected = projectedNorm(fit.theta, gradient, context.design.parmap)
                  val checked =
                    if projected <= tolerance then GradientCheck.Pass(projected, tolerance)
                    else GradientCheck.Failed(projected, tolerance)
                  val boundary = assessBoundary(fit.theta, context.design.parmap, fit.objective, tolerance, options) {
                    theta =>
                      objectiveAt(context, fit, theta)
                  }
                  val accepted = fit.optimizer.stopReason.isAcceptableConvergence
                  val reasonBuilder: scala.collection.mutable.Builder[String, Vector[String]] =
                    Vector.newBuilder[String]
                  if !accepted then reasonBuilder += "numerical_termination_not_accepted"
                  checked match
                    case GradientCheck.Failed(_, _)      => reasonBuilder += "projected_gradient_failed"
                    case GradientCheck.Unavailable(code) => reasonBuilder += code
                    case _                               => ()
                  boundary match
                    case BoundaryAdequacy.Failed(_, _, _)   => reasonBuilder += "boundary_descent_detected"
                    case BoundaryAdequacy.Unresolved(code)  => reasonBuilder += code
                    case BoundaryAdequacy.Unavailable(code) => reasonBuilder += code
                    case _                                  => ()
                  if ranks.exists(_.rank.isEmpty) then reasonBuilder += "covariance_rank_unavailable"
                  else if ranks.exists(rank => rank.rank.exists(_ < rank.eigenvalues.length)) then
                    reasonBuilder += "covariance_reduced_rank"
                  if isBoundary(context.design.parmap, fit.theta, options.boundaryTolerance) then
                    reasonBuilder += "boundary_fit_not_bootstrap_eligible"
                  val result = reasonBuilder.result().distinct
                  LmmAssessment(
                    fit.optimizer.stopReason,
                    accepted,
                    checked,
                    boundary,
                    ranks,
                    result,
                    result.headOption,
                    accepted && checked.isInstanceOf[GradientCheck.Pass] &&
                      !isBoundary(context.design.parmap, fit.theta, options.boundaryTolerance) &&
                      ranks.forall(rank => rank.rank.contains(rank.eigenvalues.length))
                  )

  private[lmm] def assessBoundary(
      theta: Vector[Double],
      parmap: Vector[(Int, Int, Int)],
      baseline: Double,
      tolerance: Double,
      options: LmmAssessmentOptions
  )(objectiveAt: Vector[Double] => Either[String, Double]): BoundaryAdequacy =
    val boundary = parmap.zipWithIndex.collect:
      case ((_, row, col), index) if row == col && math.abs(theta(index)) <= options.boundaryTolerance => index
    if boundary.isEmpty then BoundaryAdequacy.NotApplicable
    else
      var result: Option[BoundaryAdequacy] = None
      var parameter = 0
      while parameter < boundary.length && result.isEmpty do
        val index = boundary(parameter)
        val h = options.boundaryProbe
        val probes = Vector(h, h / 2.0)
        var allDescend = true
        var probe = 0
        while probe < probes.length && result.isEmpty do
          val step = probes(probe)
          objectiveAt(theta.updated(index, step)) match
            case Left(code)                                   => result = Some(BoundaryAdequacy.Unavailable(code))
            case Right(value) if value < baseline - tolerance => ()
            case Right(_)                                     => allDescend = false
          probe += 1
        if result.isEmpty && allDescend then
          objectiveAt(theta.updated(index, h)) match
            case Right(value) => result = Some(BoundaryAdequacy.Failed(index, h, value))
            case Left(code)   => result = Some(BoundaryAdequacy.Unavailable(code))
        parameter += 1
      result.getOrElse(BoundaryAdequacy.Unresolved("boundary_adequacy_unresolved"))

  private def freshEvaluation(
      context: LmmFitContext,
      fit: LmmFit,
      options: LmmAssessmentOptions
  ): Either[String, (Double, Vector[Double])] =
    (for
      workspace <- PlsWorkspace(context.design, fit.options.criterion == Criterion.REML)
      _ <- workspace.updateResponse(fit.y)
      value <- workspace.objectiveAndGradientAt(fit.theta, options.gradientOptions)
    yield (value.objective, value.gradient)).left.map(_.code)

  private def objectiveAt(context: LmmFitContext, fit: LmmFit, theta: Vector[Double]): Either[String, Double] =
    (for
      workspace <- PlsWorkspace(context.design, fit.options.criterion == Criterion.REML)
      _ <- workspace.updateResponse(fit.y)
      value <- workspace.objectiveAt(theta)
    yield value).left.map(_.code)

  private def unavailable(
      fit: LmmFit,
      code: String,
      ranks: Vector[RandomEffectCovarianceRank] = Vector.empty
  ): LmmAssessment =
    LmmAssessment(
      fit.optimizer.stopReason,
      fit.optimizer.stopReason.isAcceptableConvergence,
      GradientCheck.Unavailable(code),
      BoundaryAdequacy.Unavailable(code),
      ranks,
      (Vector(code) ++ ranks.collect { case rank if rank.rank.isEmpty => "covariance_rank_unavailable" }).distinct,
      Some(code),
      false
    )

  private def projectedNorm(theta: Vector[Double], gradient: Vector[Double], parmap: Vector[(Int, Int, Int)]): Double =
    theta.indices.iterator
      .map: index =>
        val projected =
          if parmap(index)._2 == parmap(index)._3 then math.max(0.0, theta(index) - gradient(index))
          else theta(index) - gradient(index)
        math.abs(theta(index) - projected)
      .maxOption
      .getOrElse(0.0)

  private def isBoundary(parmap: Vector[(Int, Int, Int)], theta: Vector[Double], tolerance: Double): Boolean =
    parmap.zip(theta).exists { case ((_, row, col), value) => row == col && math.abs(value) <= tolerance }

  private[mixeff4s] def covarianceRanks(
      terms: Vector[ReMat],
      theta: Vector[Double],
      options: LmmAssessmentOptions
  ): Vector[RandomEffectCovarianceRank] =
    var offset = 0
    terms.zipWithIndex.map: (term, termIndex) =>
      val fitted = term.withTheta(theta.slice(offset, offset + term.nTheta)).toOption.get
      offset += term.nTheta
      val covariance = Matrix.tabulate(fitted.vsize, fitted.vsize): (row, col) =>
        var total = 0.0
        var k = 0
        while k < fitted.vsize do
          total += fitted.lambda(row, k) * fitted.lambda(col, k)
          k += 1
        total
      val decomposition = Eigen.eigSymmetric(covariance, EigenSelection.All, EigenVectors.ValuesOnly)
      val eigenvalues =
        decomposition.toOption.fold(Vector.empty)(eigen => Vector.tabulate(eigen.size)(eigen.eigenvalues(_)))
      val scale = eigenvalues.map(math.abs).foldLeft(0.0)(math.max)
      val tolerance = options.covarianceRankAbsoluteTolerance + options.covarianceRankRelativeTolerance * scale
      RandomEffectCovarianceRank(
        termIndex,
        eigenvalues,
        decomposition.toOption.map(_ => eigenvalues.count(_ > tolerance)),
        tolerance,
        decomposition.left.toOption.map(_ => "covariance_eigendecomposition_failed")
      )
