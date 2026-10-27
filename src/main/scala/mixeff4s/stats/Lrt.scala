package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.lmm.{Criterion, LmmDesignEvidence, LmmFit, RandomDesignEvidence}

/** An asymptotic chi-square comparison of nested fixed-effect ML models. */
final case class Lrt(
    nobs: Int,
    formulas: Vector[String],
    dof: Vector[Int],
    loglik: Vector[Double],
    deviance: Vector[Double],
    chisq: Double,
    chisqDof: Int,
    pvalue: Double
)

/** Verified nested likelihood comparison. Parameter-count difference is not a reference distribution. */
final case class LrtStatistic(
    nobs: Int,
    formulas: Vector[String],
    dof: Vector[Int],
    loglik: Vector[Double],
    deviance: Vector[Double],
    likelihoodRatio: Double,
    parameterDifference: Int
)

object Lrt:
  private val LogLikTol = 1.0e-10
  private val BoundaryTolerance = 1.0e-8

  /** Ordinary asymptotic calibration is supported only for fixed-effect ML changes with unchanged RE design. Covariance
    * changes and fitted boundaries require a calibration not implemented here.
    */
  def compare(smaller: LmmFit, larger: LmmFit): FitResult[Lrt] =
    for
      evidence <- validated(smaller, larger)
      (small, large) = evidence
      _ <-
        if !sameRandom(small, large) then
          unavailable("random-effect covariance changes require an explicit reference law", "lrt_reference_unavailable")
        else if atBoundary(smaller) || atBoundary(larger) then
          unavailable(
            "fitted covariance is on or near a boundary; ordinary chi-square calibration is not qualified",
            "lrt_reference_unavailable"
          )
        else Right(())
    yield
      val result = summarize(smaller, larger)
      Lrt(
        result.nobs,
        result.formulas,
        result.dof,
        result.loglik,
        result.deviance,
        result.likelihoodRatio,
        result.parameterDifference,
        ChiSquare.sf(result.likelihoodRatio, result.parameterDifference)
      )

  /** Certify nesting by exact inclusion of fitted basis columns and matching grouping partitions. Equivalent
    * reparameterizations not witnessed by column inclusion are conservatively refused. REML additionally requires the
    * same retained fixed-effect basis, up to permutation.
    */
  def statistic(smaller: LmmFit, larger: LmmFit): FitResult[LrtStatistic] =
    validated(smaller, larger).map(_ => summarize(smaller, larger))

  private def validated(smaller: LmmFit, larger: LmmFit): FitResult[(LmmDesignEvidence, LmmDesignEvidence)] =
    for
      small <- evidence(smaller)
      large <- evidence(larger)
      _ <-
        if smaller.n != larger.n || smaller.y != larger.y then
          unavailable("models were not fitted to the same row-aligned response values")
        else if smaller.options.criterion != larger.options.criterion then
          unavailable("models mix REML and ML fit criteria; refit with a common criterion")
        else if smaller.options.criterion == Criterion.REML &&
          !sameColumns(small.fixedColumns, large.fixedColumns)
        then unavailable("REML likelihood comparisons require identical retained fixed-effect bases; refit with ML")
        else if !includes(small.fixedColumns, large.fixedColumns)(_ == _) then
          unavailable("fixed-effect nesting is not certified by exact fitted-column inclusion")
        else if !includes(small.randomBlocks, large.randomBlocks)(blockCovered) then
          unavailable(
            "random-effect nesting is not certified by grouping partitions, basis columns and covariance constraints"
          )
        else if larger.dof <= smaller.dof then
          unavailable("larger model must have more parameters than the smaller model")
        else if larger.loglikelihood - smaller.loglikelihood < -LogLikTol then
          unavailable("log-likelihood must not be lower in the larger model")
        else Right(())
    yield (small, large)

  private def evidence(fit: LmmFit): FitResult[LmmDesignEvidence] =
    fit.designEvidence match
      case Some(value)
          if value.matches(fit) && fit.objective.isFinite && fit.sigma.isFinite && fit.sigma > 0.0 &&
            fit.optimizer.x == fit.theta && fit.optimizer.fmin == fit.objective &&
            fit.optimizer.requireConvergence.isRight =>
        Right(value)
      case _ =>
        Left(
          MixedModelError.InferenceUnavailable(
            "lrt_unavailable",
            "fit lacks matching design and numerical termination evidence"
          )
        )

  private def summarize(smaller: LmmFit, larger: LmmFit): LrtStatistic =
    LrtStatistic(
      smaller.n,
      Vector(smaller.formula.toString, larger.formula.toString),
      Vector(smaller.dof, larger.dof),
      Vector(smaller.loglikelihood, larger.loglikelihood),
      Vector(smaller.objective, larger.objective),
      math.max(0.0, 2.0 * (larger.loglikelihood - smaller.loglikelihood)),
      larger.dof - smaller.dof
    )

  private def sameColumns(a: Vector[Vector[Double]], b: Vector[Vector[Double]]): Boolean =
    a.length == b.length && includes(a, b)(_ == _)

  private def blockCovered(small: RandomDesignEvidence, large: RandomDesignEvidence): Boolean =
    small.partition == large.partition && (small.diagonal || !large.diagonal) &&
      includes(small.columns, large.columns)(_ == _)

  private def sameRandom(a: LmmDesignEvidence, b: LmmDesignEvidence): Boolean =
    a.randomBlocks.length == b.randomBlocks.length &&
      includes(a.randomBlocks, b.randomBlocks): (small, large) =>
        small.diagonal == large.diagonal && small.columns.length == large.columns.length &&
          blockCovered(small, large)

  /** Bipartite matching preserves multiplicity and avoids greedy matching failures. */
  private def includes[A](small: Vector[A], large: Vector[A])(covers: (A, A) => Boolean): Boolean =
    if small.length > large.length then false
    else
      val edges = small.map(s => large.indices.filter(j => covers(s, large(j))).toVector)
      val owner = Array.fill(large.length)(-1)
      def assign(i: Int, seen: Array[Boolean]): Boolean =
        edges(i).exists: j =>
          if seen(j) then false
          else
            seen(j) = true
            if owner(j) < 0 || assign(owner(j), seen) then
              owner(j) = i
              true
            else false
      small.indices.forall(i => assign(i, Array.fill(large.length)(false)))

  private def atBoundary(fit: LmmFit): Boolean =
    fit.parmap.zip(fit.theta).exists { case ((_, row, col), value) =>
      row == col && value <= BoundaryTolerance
    }

  private def unavailable(details: String, code: String = "lrt_unavailable"): FitResult[Unit] =
    Left(MixedModelError.InferenceUnavailable(code, details))
