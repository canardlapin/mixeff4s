package mixeff4s.glmm

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.lmm.{LmmAssessment, LmmAssessmentOptions}
import mixeff4s.model.Family
import mixeff4s.optimizer.TrustBqStopReason

final case class GlmmAssessmentOptions(
    gradientTolerance: Double = 1e-4,
    conditionalScoreTolerance: Double = 1e-8,
    objectiveRelativeTolerance: Double = 1e-9,
    boundaryTolerance: Double = 1e-8,
    covarianceRankRelativeTolerance: Double = 1e-12,
    covarianceRankAbsoluteTolerance: Double = 0.0
):
  def validate: FitResult[Unit] =
    val nonnegative = Vector(
      objectiveRelativeTolerance,
      boundaryTolerance,
      covarianceRankRelativeTolerance,
      covarianceRankAbsoluteTolerance
    )
    if !gradientTolerance.isFinite || gradientTolerance <= 0.0 || !conditionalScoreTolerance.isFinite ||
      conditionalScoreTolerance <= 0.0 || nonnegative.exists(x => !x.isFinite || x < 0.0)
    then Left(MixedModelError.InvalidArgument("invalid GLMM assessment tolerances"))
    else Right(())

/** Sufficient fixed-effect recession evidence. Excluded does not certify a finite joint MLE. */
enum FixedEffectRecession:
  case Excluded(witness: String, rank: Int)
  case Detected(reasonCode: String, direction: Vector[Double])
  case Unavailable(reasonCode: String)

final case class GlmmCovarianceRank(term: String, dimension: Int, rank: Option[Int], eigenvalues: Vector[Double])

/** Fresh numerical evidence for the actual estimator, with conservative bootstrap eligibility. Pure groups are risk
  * flags, not proof that Gaussian-penalized modes do not exist. Finite joint-MLE existence and global optimality remain
  * unassessed.
  */
final case class GlmmAssessment(
    numericalTermination: TrustBqStopReason,
    numericalTerminationAccepted: Boolean,
    conditional: Either[String, ConditionalModeResult],
    outer: Option[JointAssessment],
    covarianceRanks: Vector[GlmmCovarianceRank],
    fixedEffectRecession: FixedEffectRecession,
    riskFlags: Vector[String],
    refusalCodes: Vector[String]
):
  def eligibleForBootstrap: Boolean = refusalCodes.isEmpty
  def bootstrapIneligibilityReason: Option[String] = refusalCodes.headOption

object GlmmAssessment:
  def assess(fit: GlmmFit, options: GlmmAssessmentOptions = GlmmAssessmentOptions()): GlmmAssessment =
    val context = fit.context
    val rankOptions = LmmAssessmentOptions(
      covarianceRankRelativeTolerance = options.covarianceRankRelativeTolerance,
      covarianceRankAbsoluteTolerance = options.covarianceRankAbsoluteTolerance
    )
    val ranks =
      if options.validate.isLeft then Vector.empty
      else
        LmmAssessment
          .covarianceRanks(context.design.reterms, fit.theta, rankOptions)
          .zip(context.design.reterms)
          .map: (evidence, term) =>
            GlmmCovarianceRank(term.groupingName, term.vsize, evidence.rank, evidence.eigenvalues)
    val fixed = fixedEffectEvidence(context)
    val risks = Vector("finite_joint_mle_not_assessed") ++ pureGroupRisks(context)
    val accepted = fit.optimizer.requireConvergence.isRight
    var refusals = Vector.empty[String]
    if options.validate.isLeft then refusals :+= "glmm_assessment_invalid_options"
    if !accepted then refusals :+= "bootstrap_optimizer_unqualified"
    if context.response.caseWeights.exists(_ != 1.0) then refusals :+= "bootstrap_case_weights"
    if fit.approximation == Approximation.FastPirls then refusals :+= "bootstrap_requires_joint"
    if ranks.isEmpty || ranks.exists(_.rank.isEmpty) then refusals :+= "bootstrap_covariance_rank_unavailable"
    else if ranks.exists(r => r.rank.get < r.dimension) then refusals :+= "bootstrap_reduced_rank"
    if context.design.parmap.zip(fit.theta).exists { case ((_, row, col), value) =>
        row == col && value <= options.boundaryTolerance
      }
    then refusals :+= "bootstrap_covariance_boundary"
    fixed match
      case FixedEffectRecession.Excluded(_, _) => ()
      case FixedEffectRecession.Detected(_, _) => refusals :+= "bootstrap_fixed_effect_recession"
      case FixedEffectRecession.Unavailable(_) => refusals :+= "bootstrap_fixed_effect_evidence_unavailable"
    val fresh =
      if fit.approximation == Approximation.FastPirls then
        Left(
          MixedModelError
            .InferenceUnavailable("bootstrap_requires_joint", "Fast-PIRLS is outside GLMM bootstrap qualification.")
        )
      else
        options.validate.flatMap: _ =>
          val controls = context.options.copy(
            conditionalScoreTolerance = options.conditionalScoreTolerance,
            stationarityTolerance = options.gradientTolerance
          )
          for
            workspace <- JointWorkspace.create(context.design, fit.family, fit.link, context.response, controls)
            evaluation <- workspace.evaluate(fit.beta, fit.theta)
            _ <-
              if math.abs(evaluation.objective - fit.deviance) <= options.objectiveRelativeTolerance * math.max(
                  1.0,
                  math.abs(fit.deviance)
                )
              then Right(())
              else
                Left(
                  MixedModelError.InferenceUnavailable(
                    "glmm_assessment_objective_mismatch",
                    "Fresh objective differs from the fitted objective."
                  )
                )
            parameters = fit.beta ++ fit.theta
            lower = Vector.fill(fit.beta.length)(Double.NegativeInfinity) ++ workspace.pls.lowerBounds
            outer = JointLaplace.assess(workspace, parameters, evaluation.objective, lower, options.gradientTolerance)
          yield (evaluation.conditional, outer)
    val (conditional, outer): (Either[String, ConditionalModeResult], Option[JointAssessment]) = fresh match
      case Left(error) =>
        refusals :+= error.code
        (Left(error.code), None)
      case Right((inner, evidence)) =>
        if !evidence.firstOrderPassed then refusals :+= "bootstrap_outer_unqualified"
        (Right(inner), Some(evidence))
    GlmmAssessment(fit.optimizer.stopReason, accepted, conditional, outer, ranks, fixed, risks, refusals.distinct)

  private[glmm] def fixedEffectEvidence(context: GlmmFitContext): FixedEffectRecession =
    val design = context.design
    val response = context.response
    val p = design.p
    val rows = Vector.tabulate(design.n, p)((i, j) => design.xy.xy(i, j))
    if p == 0 then return FixedEffectRecession.Excluded("no_fixed_effect_coordinates", 0)
    if context.family == Family.Poisson then
      val rank = spanningRank(rows.zip(response.y).collect { case (row, y) if y > 0.0 => row }, p)
      if rank == p then FixedEffectRecession.Excluded("positive_count_rows_span_retained_design", rank)
      else FixedEffectRecession.Unavailable("positive_count_rows_do_not_span")
    else if context.family == Family.Binomial || context.family == Family.Bernoulli then
      val mixed = rows.zip(response.y).collect { case (row, y) if y > 0.0 && y < 1.0 => row }
      val mixedRank = spanningRank(mixed, p)
      if mixedRank == p then return FixedEffectRecession.Excluded("mixed_binomial_rows_span_retained_design", mixedRank)
      val cells = rows
        .zip(response.y)
        .groupBy(_._1)
        .iterator
        .collect:
          case (row, observations) if observations.exists(_._2 > 0.0) && observations.exists(_._2 < 1.0) => row
      val cellRank = spanningRank(cells.toVector, p)
      if cellRank == p then return FixedEffectRecession.Excluded("mixed_outcome_cells_span_retained_design", cellRank)
      val signed = rows
        .zip(response.y)
        .flatMap: (row, y) =>
          (if y > 0.0 then Vector(row) else Vector.empty) ++ (if y < 1.0 then Vector(row.map(-_)) else Vector.empty)
      def verifiedDirection(direction: Vector[Double]): Option[FixedEffectRecession] =
        val products = signed.map(row => row.zip(direction).map(_ * _).sum)
        val scale = products.map(math.abs).maxOption.getOrElse(0.0)
        if products.forall(_ >= 0.0) && scale > 1e-10 then
          Some(FixedEffectRecession.Detected("nonnegative_signed_design_direction", direction))
        else None
      if p == 1 then
        val values = signed.map(_.head)
        val scale = values.map(math.abs).maxOption.getOrElse(0.0)
        if values.min < -1e-10 * scale && values.max > 1e-10 * scale then
          FixedEffectRecession.Excluded("one_coordinate_signed_rows_have_both_signs", 1)
        else
          verifiedDirection(Vector(1.0))
            .orElse(verifiedDirection(Vector(-1.0)))
            .getOrElse(FixedEffectRecession.Unavailable("signed_row_geometry_uncertain"))
      else if p == 2 then
        val intercept = (0 until p).find(j => rows.forall(_(j) == 1.0))
        intercept match
          case None           => FixedEffectRecession.Unavailable("two_coordinate_nonintercept_geometry_unassessed")
          case Some(constant) =>
            val variable = 1 - constant
            val success = rows.zip(response.y).collect { case (row, y) if y > 0.0 => row(variable) }
            val failure = rows.zip(response.y).collect { case (row, y) if y < 1.0 => row(variable) }
            if success.isEmpty || failure.isEmpty then FixedEffectRecession.Unavailable("missing_outcome_class")
            else
              val scale = (success ++ failure).map(math.abs).max.max(1e-100)
              val overlap = math.min(success.max - failure.min, failure.max - success.min)
              if overlap > 1e-10 * scale then FixedEffectRecession.Excluded("strict_outcome_hull_overlap", 2)
              else
                val cut = 0.5 * failure.max + 0.5 * success.min
                val direction = Vector.tabulate(2)(j => if j == constant then -cut else 1.0)
                val reverseCut = 0.5 * success.max + 0.5 * failure.min
                val reverse = Vector.tabulate(2)(j => if j == constant then reverseCut else -1.0)
                verifiedDirection(direction)
                  .orElse(verifiedDirection(reverse))
                  .getOrElse(FixedEffectRecession.Unavailable("outcome_hulls_touch_or_are_numerically_uncertain"))
      else FixedEffectRecession.Unavailable("no_sufficient_fixed_effect_witness")
    else FixedEffectRecession.Unavailable("family_recession_unassessed")

  /** Column-scaled, pivoted, twice-reorthogonalized row-span witness. */
  private def spanningRank(rows: Vector[Vector[Double]], p: Int): Int =
    if rows.isEmpty then 0
    else
      val scale = Vector.tabulate(p)(j => rows.iterator.map(row => math.abs(row(j))).max)
      if scale.exists(_ == 0.0) then return 0
      val residuals = rows.map(row => Array.tabulate(p)(j => row(j) / scale(j))).toArray
      var rank = 0
      while rank < p do
        val selected = residuals.indices.maxBy(i => residuals(i).map(x => x * x).sum)
        val norm = math.sqrt(residuals(selected).map(x => x * x).sum)
        if !norm.isFinite || norm <= 1e-8 * math.sqrt(p.toDouble) then return rank
        val basis = residuals(selected).map(_ / norm)
        residuals.foreach: row =>
          var pass = 0
          while pass < 2 do
            var dot = 0.0
            var j = 0
            while j < p do { dot += row(j) * basis(j); j += 1 }
            j = 0
            while j < p do { row(j) -= dot * basis(j); j += 1 }
            pass += 1
        rank += 1
      rank

  private def pureGroupRisks(context: GlmmFitContext): Vector[String] =
    context.design.reterms.flatMap: term =>
      val counts = Array.fill(term.nLevels)(0)
      val successes = Array.fill(term.nLevels)(false)
      val failures = Array.fill(term.nLevels)(false)
      term.refs.indices.foreach: i =>
        val group = term.refs(i)
        counts(group) += 1
        if context.response.y(i) > 0.0 then successes(group) = true
        if context.response.y(i) < 1.0 || context.family == Family.Poisson then failures(group) = true
      val pure = counts.indices.count(g => counts(g) > 0 && (!successes(g) || !failures(g)))
      if pure > 0 then Vector(s"pure_outcome_groups:${term.groupingName}:$pure") else Vector.empty
