package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.glmm.*
import mixeff4s.lmm.*

/** Local likelihood evidence for LRT simulation, deliberately distinct from coefficient-bootstrap eligibility. Boundary
  * nulls require strict, step-stable one-sided VARIANCE scores and positive tangent curvature. This is local numerical
  * evidence, never a finite-MLE/global-optimum certificate.
  */
private[stats] object LikelihoodQualification:
  def glmmBoundary(fit: GlmmFit): Boolean = fit.context.design.parmap.zip(fit.theta).exists { case ((_, r, c), v) =>
    r == c && v <= 1e-8
  }
  def glmm(fit: GlmmFit): FitResult[Unit] =
    val evidence = GlmmAssessment.assess(fit)
    val allowed = Set("bootstrap_covariance_boundary", "bootstrap_reduced_rank", "bootstrap_outer_unqualified")
    val bad = evidence.refusalCodes.filterNot(allowed)
    if bad.nonEmpty then refuse("likelihood_fit_unqualified", bad.mkString(", "))
    else if !glmmBoundary(fit) then
      if evidence.eligibleForBootstrap then
        for
          evaluator <- JointLikelihood(fit)
          _ <- ObservedInformation.compute(
            fit.beta ++ fit.theta,
            Vector.fill(fit.beta.length)(Double.NegativeInfinity) ++ fit.context.design.parmap.map { case (_, r, c) =>
              if r == c then 0.0 else Double.NegativeInfinity
            },
            CurvatureOptions()
          )(evaluator.evaluate)
        yield ()
      else refuse("likelihood_fit_unqualified", evidence.refusalCodes.mkString(", "))
    else if fit.context.design.reterms.exists(_.vsize != 1) then
      refuse(
        "likelihood_boundary_shape",
        "Boundary LRT qualification currently supports independent scalar RE terms only."
      )
    else
      for
        evaluator <- JointLikelihood(fit)
        _ <- boundary(
          fit.beta ++ fit.theta,
          Vector.fill(fit.beta.length)(Double.NegativeInfinity) ++ fit.context.design.parmap.map { case (_, r, c) =>
            if r == c then 0.0 else Double.NegativeInfinity
          },
          fit.deviance
        )(evaluator.evaluate)
      yield ()

  def lmm(fit: LmmFit): FitResult[Unit] =
    for
      context <- fit.validatedContext
      _ <-
        if fit.options.criterion == Criterion.ML then Right(())
        else refuse("lrt_requires_ml", "Null bootstrap comparisons require ML fits.")
      _ <- fit.optimizer.requireConvergence
      workspace <- PlsWorkspace(context.design, false)
      _ <- workspace.updateResponse(fit.y)
      value <- workspace.objectiveAt(fit.theta)
      _ <-
        if math.abs(value - fit.objective) <= 1e-9 * math.max(1, math.abs(value)) then Right(())
        else refuse("likelihood_objective_mismatch", "Fresh LMM objective mismatch.")
      _ <-
        if fit.sigma.isFinite && fit.sigma > 0 then Right(())
        else refuse("likelihood_scale_invalid", "Positive residual scale required.")
      _ <-
        if fit.theta.indices.exists(i => workspace.lowerBounds(i) == 0 && fit.theta(i) <= 1e-8) then
          if context.design.reterms.exists(_.vsize != 1) then
            refuse("likelihood_boundary_shape", "Boundary LRT qualification requires scalar independent RE terms.")
          else boundary(fit.theta, workspace.lowerBounds, fit.objective)(workspace.objectiveAt)
        else
          val assessment = LmmAssessment.assess(fit)
          if assessment.eligibleForBootstrap then
            ObservedInformation
              .compute(fit.theta, workspace.lowerBounds, CurvatureOptions())(workspace.objectiveAt)
              .map(_ => ())
          else refuse("likelihood_fit_unqualified", assessment.reasonCodes.mkString(", "))
    yield ()

  private def boundary(x: Vector[Double], lower: Vector[Double], baseline: Double)(
      f: Vector[Double] => FitResult[Double]
  ): FitResult[Unit] =
    val active = x.indices.filter(i => lower(i) == 0 && x(i) <= 1e-8).toVector
    val free = x.indices.filterNot(active.contains).toVector
    val zero = active.foldLeft(x)((a, i) => a.updated(i, 0.0))
    for
      f0 <- f(zero)
      _ <-
        if math.abs(f0 - baseline) <= 1e-8 * math.max(1, math.abs(baseline)) then Right(())
        else refuse("likelihood_boundary_snap", "Boundary snap changes objective.")
      _ <- active.foldLeft[FitResult[Unit]](Right(())): (acc, index) =>
        acc.flatMap: _ =>
          val h = 1e-4
          def v(step: Double) = f(zero.updated(index, math.sqrt(step)))
          for
            one <- v(h)
            two <- v(2 * h)
            half <- v(h / 2)
            g1 = (-3 * f0 + 4 * one - two) / (2 * h)
            g2 = (-3 * f0 + 4 * half - one) / h
            error = math.abs(g1 - g2)
            _ <-
              if g2.isFinite && g2 > math.max(1e-5, 4 * error) then Right(())
              else
                refuse(
                  "likelihood_boundary_score_unresolved",
                  s"Variance-coordinate score $g2 with step discrepancy $error is not strictly positive."
                )
          yield ()
      _ <-
        if free.isEmpty then Right(())
        else
          val xf = free.map(zero); val lf = free.map(lower)
          def evaluate(p: Vector[Double]) = f(free.zip(p).foldLeft(zero) { case (a, (i, v)) => a.updated(i, v) })
          for
            _ <- free.foldLeft[FitResult[Unit]](Right(())): (acc, j) =>
              acc.flatMap: _ =>
                val h = 1e-4 * math.max(1, math.abs(zero(j)))
                for
                  plus <- f(zero.updated(j, zero(j) + h))
                  minus <- f(zero.updated(j, zero(j) - h))
                  g = (plus - minus) / (2 * h)
                  _ <-
                    if g.isFinite && math.abs(g) <= 1e-4 then Right(())
                    else refuse("likelihood_tangent_gradient", s"Free tangent gradient $g fails.")
                yield ()
            _ <- ObservedInformation.compute(xf, lf, CurvatureOptions())(evaluate)
          yield ()
      _ <- Vector(1e-4, 1e-2, .1, .5, 1.0).foldLeft[FitResult[Unit]](Right(())): (acc, step) =>
        acc.flatMap: _ =>
          val directions = active.map(i => zero.updated(i, step)) :+ active.foldLeft(zero)((a, i) => a.updated(i, step))
          directions.foldLeft[FitResult[Unit]](Right(())): (prior, probe) =>
            prior.flatMap: _ =>
              f(probe).flatMap: value =>
                if value >= f0 - 1e-8 then Right(())
                else refuse("likelihood_boundary_descent", "Feasible boundary probe finds descent.")
    yield ()

  private def refuse(code: String, message: String): FitResult[Nothing] = Left(
    MixedModelError.InferenceUnavailable(code, message)
  )
