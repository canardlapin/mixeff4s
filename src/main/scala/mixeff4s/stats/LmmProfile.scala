package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.lmm.{Criterion, LmmAssessment, LmmAssessmentOptions, LmmFit, PlsWorkspace}
import mixeff4s.optimizer.{TrustBq, TrustBqOptions, TrustBqResult}

enum LmmProfileTarget:
  case Coefficient(name: String)
  case RandomEffectStdDev(group: String)
  case RandomEffectVariance(group: String)
  case ResidualStdDev

final case class LmmProfileOptions(
    level: Double = 0.95,
    refitWithML: Boolean = false,
    maxEvaluations: Int = 300,
    maxBracketSteps: Int = 16,
    maxBisectionSteps: Int = 32,
    endpointTolerance: Double = 1e-5,
    stationarityTolerance: Double = 1e-4,
    optimizerControls: TrustBqOptions =
      TrustBqOptions(initialRadius = 0.5, finalRadius = 1e-5, ftolRequiresLocalRadius = true),
    assessmentOptions: LmmAssessmentOptions = LmmAssessmentOptions()
)

final case class LmmProfileEvaluation(
    value: Double,
    objective: FitResult[Double],
    optimizer: Option[TrustBqResult],
    nuisanceStationarity: Option[Double],
    qualified: Boolean
)
final case class LmmProfileEndpoint(
    value: Option[Double],
    boundary: Boolean,
    reason: Option[String],
    bracket: Option[(Double, Double)] = None,
    objectiveDifference: Option[Double] = None
)
final case class LmmProfileResult(
    target: LmmProfileTarget,
    estimate: Double,
    level: Double,
    cutoff: Double,
    lower: LmmProfileEndpoint,
    upper: LmmProfileEndpoint,
    originalCriterion: Criterion,
    baseline: LmmFit,
    baselineAssessment: LmmAssessment,
    evaluations: Vector[LmmProfileEvaluation],
    options: LmmProfileOptions
)

object LmmProfile:
  private final case class Scalar(
      estimate: Double,
      lower: Double,
      step: Double,
      evaluate: Double => FitResult[(Double, Option[TrustBqResult], Option[Double])],
      lowerOpen: Boolean = false
  )

  private def prepare(
      fit: LmmFit,
      target: LmmProfileTarget,
      options: LmmProfileOptions
  ): FitResult[(LmmFit, LmmAssessment, Scalar)] =
    for
      _ <- valid(options)
      baseline <-
        if fit.options.criterion == Criterion.ML then fit.validatedContext.map(_ => fit)
        else if options.refitWithML then
          fit.validatedContext
            .flatMap(context => context.compiled.fit(fit.y, fit.options.copy(criterion = Criterion.ML)))
        else
          unavailable(
            "profile_reml_requires_ml_refit",
            "REML profiles require refitWithML = true; no criterion switch is implicit."
          )
      scalar <- compile(baseline, target, options)
      assessment = LmmAssessment.assess(baseline, options.assessmentOptions)
      _ <-
        if assessment.eligibleForBootstrap then Right(())
        else
          unavailable(
            "profile_baseline_unqualified",
            assessment.bootstrapIneligibilityReason.getOrElse(assessment.reasonCodes.mkString(", "))
          )
    yield (baseline, assessment, scalar)

  /** A constrained likelihood evaluation; solver failures are retained inside the returned ledger row. */
  def evaluate(
      fit: LmmFit,
      target: LmmProfileTarget,
      value: Double,
      options: LmmProfileOptions = LmmProfileOptions()
  ): FitResult[LmmProfileEvaluation] =
    prepare(fit, target, options).flatMap: (baseline, _, scalar) =>
      if !value.isFinite || value < scalar.lower then
        Left(MixedModelError.InvalidArgument("profile target is outside its finite domain"))
      else Right(evaluation(baseline, scalar, value, options))

  private def evaluation(
      baseline: LmmFit,
      scalar: Scalar,
      x: Double,
      options: LmmProfileOptions
  ): LmmProfileEvaluation =
    scalar.evaluate(x) match
      case Left(e)               => LmmProfileEvaluation(x, Left(e), None, None, false)
      case Right((v, opt, grad)) =>
        val tolerance = options.assessmentOptions.objectiveTolerance(baseline.objective)
        val stopped = opt.forall(_.requireConvergence.isRight)
        val stationary = grad.forall(g => g.isFinite && g <= options.stationarityTolerance)
        val result: FitResult[Double] =
          if !v.isFinite then unavailable("profile_nonfinite", "Nonfinite constrained objective.")
          else if v < baseline.objective - tolerance then
            unavailable(
              "profile_constrained_below_baseline",
              "Constrained likelihood is better than the qualified ML baseline."
            )
          else if !stopped then opt.get.requireConvergence.flatMap(_ => Right(v))
          else if !stationary then
            unavailable(
              "profile_nuisance_unqualified",
              "Constrained nuisance stationarity failed or reached an unqualified boundary."
            )
          else Right(v)
        LmmProfileEvaluation(x, result, opt, grad, result.isRight)

  def confint(
      fit: LmmFit,
      target: LmmProfileTarget,
      options: LmmProfileOptions = LmmProfileOptions()
  ): FitResult[LmmProfileResult] =
    prepare(fit, target, options).map: (baseline, assessment, scalar) =>
      val threshold = baseline.objective + InferenceQuantiles.chiSquare1(options.level)
      val ledger = scala.collection.mutable.ArrayBuffer.empty[LmmProfileEvaluation]
      def at(x: Double): Either[String, Double] =
        val row = evaluation(baseline, scalar, x, options)
        ledger += row
        row.objective.left.map(_.code)
      val center = at(scalar.estimate)
      val centerError = center.left.toOption.orElse(
        center.toOption
          .filter(v =>
            math.abs(v - baseline.objective) > options.assessmentOptions.objectiveTolerance(baseline.objective)
          )
          .map(_ => "profile_center_mismatch")
      )
      val lower = centerError.fold(
        end(scalar.estimate, scalar.lower, -1, scalar.step, threshold, options, at, scalar.lowerOpen)
      )(code => LmmProfileEndpoint(None, false, Some(code)))
      val upper = centerError.fold(
        end(scalar.estimate, Double.PositiveInfinity, 1, scalar.step, threshold, options, at)
      )(code => LmmProfileEndpoint(None, false, Some(code)))
      val baselineInvalid =
        ledger.exists(_.objective.left.toOption.exists(_.code == "profile_constrained_below_baseline"))
      val invalid = LmmProfileEndpoint(None, false, Some("profile_constrained_below_baseline"))
      LmmProfileResult(
        target,
        scalar.estimate,
        options.level,
        threshold - baseline.objective,
        if baselineInvalid then invalid else lower,
        if baselineInvalid then invalid else upper,
        fit.options.criterion,
        baseline,
        assessment,
        ledger.toVector,
        options
      )

  private def compile(fit: LmmFit, target: LmmProfileTarget, options: LmmProfileOptions): FitResult[Scalar] =
    def ws =
      fit.validatedContext.flatMap(c => PlsWorkspace(c.design, false).flatMap(w => w.updateResponse(fit.y).map(_ => w)))
    def thetaProfile(start: Vector[Double])(f: (PlsWorkspace, Vector[Double]) => FitResult[Double]) = ws.flatMap: w =>
      val controls = options.optimizerControls.copy(maxEvaluations = options.maxEvaluations)
      TrustBq
        .minimize(start, w.lowerBounds, w.upperBounds, controls)(x => f(w, x))
        .flatMap: r =>
          f(w, r.x).flatMap(v =>
            gradient(r.x, w.lowerBounds, w.upperBounds, x => f(w, x)).map(g => (v, Some(r), Some(g)))
          )
    target match
      case LmmProfileTarget.Coefficient(name) =>
        val j = fit.feNames.indexOf(name)
        if j < 0 then unavailable("profile_unknown_or_dropped_coefficient", s"'$name' is not retained")
        else
          val se = math.sqrt(fit.vcov(j)(j))
          if !se.isFinite || se <= 0 then
            unavailable("profile_coefficient_variance_unavailable", s"'$name' has no positive variance")
          else
            Right(
              Scalar(
                fit.beta(j),
                Double.NegativeInfinity,
                se.max(.1),
                t =>
                  thetaProfile(fit.theta): (w, x) =>
                    w.objectiveAt(x)
                      .flatMap: _ =>
                        val q = w.pwrss; val aInv = w.vcov(j)(j) / (q / fit.n)
                        val qt = q + (t - w.beta(j)) * (t - w.beta(j)) / aInv
                        val v = w.logdetRe + fit.n * (1 + math.log(2 * math.Pi * qt / fit.n))
                        if v.isFinite && aInv > 0 then Right(v)
                        else Left(MixedModelError.Optimization("non-finite coefficient profile"))
              )
            )
      case LmmProfileTarget.ResidualStdDev =>
        Right(
          Scalar(
            fit.sigma,
            0,
            (fit.sigma * .2).max(.05),
            s =>
              if s <= 0 || !s.isFinite then Left(MixedModelError.InvalidArgument("residual SD must be positive"))
              else
                thetaProfile(fit.theta): (w, x) =>
                  w.objectiveAt(x).map(_ => w.logdetRe + fit.n * math.log(2 * math.Pi * s * s) + w.pwrss / (s * s))
            ,
            lowerOpen = true
          )
        )
      case LmmProfileTarget.RandomEffectStdDev(group)   => random(fit, group, false, options)
      case LmmProfileTarget.RandomEffectVariance(group) => random(fit, group, true, options)

  private def random(fit: LmmFit, group: String, variance: Boolean, options: LmmProfileOptions): FitResult[Scalar] =
    fit.validatedContext.flatMap: c =>
      val hits = c.design.reterms.zipWithIndex.filter(_._1.groupingName == group)
      if hits.size != 1 || hits.head._1.vsize != 1 || hits.head._1.nTheta != 1 then
        unavailable("profile_random_effect_not_scalar_independent", s"'$group' is not scalar independent")
      else
        val term = hits.head._2; val idx = fit.parmap.indexWhere(_._1 == term)
        val estimate = fit.theta(idx) * fit.sigma
        if estimate <= 0 || !estimate.isFinite then
          unavailable("profile_random_effect_boundary", s"'$group' is on the boundary")
        else
          val shown = if variance then estimate * estimate else estimate
          Right(
            Scalar(
              shown,
              0,
              (shown * .2).max(.05),
              value =>
                val tau = if variance then math.sqrt(value) else value
                if tau < 0 || !tau.isFinite then
                  Left(MixedModelError.InvalidArgument("random-effect SD must be nonnegative"))
                else
                  PlsWorkspace(c.design, false).flatMap(w =>
                    w.updateResponse(fit.y)
                      .flatMap: _ =>
                        val initial = fit.theta.indices.filter(_ != idx).map(fit.theta).toVector :+ math.log(fit.sigma)
                        val lower = fit.theta.indices
                          .filter(_ != idx)
                          .map(i => if fit.parmap(i)._2 == fit.parmap(i)._3 then 0.0 else Double.NegativeInfinity)
                          .toVector :+ math.log(1e-12)
                        val upper = Vector.fill(initial.size)(Double.PositiveInfinity)
                        def unpack(x: Vector[Double]) =
                          var k = 0;
                          fit.theta.indices
                            .map(i => if i == idx then tau / math.exp(x.last) else { val z = x(k); k += 1; z })
                            .toVector
                        def f(x: Vector[Double]) = w
                          .objectiveAt(unpack(x))
                          .map(_ =>
                            w.logdetRe + fit.n * math.log(2 * math.Pi * math.exp(2 * x.last)) + w.pwrss / math
                              .exp(2 * x.last)
                          )
                        TrustBq
                          .minimize(
                            initial,
                            lower,
                            upper,
                            options.optimizerControls.copy(maxEvaluations = options.maxEvaluations)
                          )(f)
                          .flatMap(r =>
                            f(r.x).flatMap(v => gradient(r.x, lower, upper, f).map(g => (v, Some(r), Some(g))))
                          )
                  )
            )
          )

  private def gradient(
      x: Vector[Double],
      lo: Vector[Double],
      hi: Vector[Double],
      f: Vector[Double] => FitResult[Double]
  ): FitResult[Double] =
    if x.indices.exists(i => lo(i).isFinite && x(i) <= lo(i) + 1e-8) then
      return Right(Double.PositiveInfinity) // nuisance-boundary qualification is deliberately unavailable
    x.indices.foldLeft[FitResult[Double]](Right(0.0)) { (acc, i) =>
      acc.flatMap { m =>
        val a = math.max(lo(i), x(i) - 1e-5)
        val b = math.min(hi(i), x(i) + 1e-5)
        if a == b then Right(m)
        else
          for
            u <- f(x.updated(i, a))
            v <- f(x.updated(i, b))
          yield m.max(math.abs((v - u) / (b - a)))
      }
    }

  private def end(
      center: Double,
      bound: Double,
      dir: Int,
      step: Double,
      threshold: Double,
      options: LmmProfileOptions,
      at: Double => Either[String, Double],
      lowerOpen: Boolean = false
  ): LmmProfileEndpoint =
    var inside = center; var d = step; var outside: Option[Double] = None; var why: Option[String] = None; var i = 0
    while i < options.maxBracketSteps && outside.isEmpty && why.isEmpty do
      val x =
        if dir < 0 then
          if lowerOpen && center - d <= bound then bound + (inside - bound) / 2 else math.max(bound, center - d)
        else center + d
      at(x) match
        case Left(code)                 => why = Some(code)
        case Right(v) if v >= threshold => outside = Some(x)
        case Right(_) if x == bound     => return LmmProfileEndpoint(Some(bound), true, None)
        case Right(_)                   => inside = x; d *= 2
      i += 1
    outside match
      case None      => LmmProfileEndpoint(None, false, why.orElse(Some("profile_endpoint_not_bracketed")))
      case Some(out) =>
        var a = if dir < 0 then out else inside; var b = if dir < 0 then inside else out; i = 0
        while i < options.maxBisectionSteps && math.abs(b - a) > options.endpointTolerance do
          val m = (a + b) / 2;
          at(m) match
            case Left(code) => return LmmProfileEndpoint(None, false, Some(code))
            case Right(v)   =>
              if v >= threshold then if dir < 0 then a = m else b = m else if dir < 0 then b = m else a = m
          i += 1
        if math.abs(b - a) > options.endpointTolerance then
          LmmProfileEndpoint(None, false, Some("profile_root_unresolved"))
        else
          val midpoint = (a + b) / 2
          at(midpoint) match
            case Left(code)   => LmmProfileEndpoint(None, false, Some(code), Some((a, b)))
            case Right(value) =>
              LmmProfileEndpoint(
                Some(midpoint),
                false,
                None,
                Some((a, b)),
                Some(value - (threshold - InferenceQuantiles.chiSquare1(options.level)))
              )

  private def valid(o: LmmProfileOptions): FitResult[Unit] =
    if !o.level.isFinite || o.level <= 0 || o.level >= 1 || o.maxEvaluations <= 0 || o.maxBracketSteps <= 0 || o.maxBisectionSteps <= 0 || !o.endpointTolerance.isFinite || o.endpointTolerance <= 0 || !o.stationarityTolerance.isFinite || o.stationarityTolerance < 0
    then Left(MixedModelError.InvalidArgument("invalid profile options"))
    else o.assessmentOptions.validate
  private def unavailable[A](c: String, d: String): FitResult[A] = Left(MixedModelError.InferenceUnavailable(c, d))
