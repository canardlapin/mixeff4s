package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.{FixedTerm, GroupingFactor, RandomTerm}
import mixeff4s.lmm.{Criterion, LmmFit}

/** Likelihood-ratio comparison of two nested LMM fits. */
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

object Lrt:
  private val LogLikTol = 1.0e-10

  def compare(smaller: LmmFit, larger: LmmFit): FitResult[Lrt] =
    refuse(smaller, larger).map: _ =>
      val ll0 = smaller.loglikelihood
      val ll1 = larger.loglikelihood
      val diff = ll1 - ll0
      val withinTol = diff <= 0.0
      val chi = if withinTol then 0.0 else 2.0 * diff
      val ddof = larger.dof - smaller.dof
      Lrt(
        smaller.n,
        Vector(smaller.formula.toString, larger.formula.toString),
        Vector(smaller.dof, larger.dof),
        Vector(ll0, ll1),
        Vector(smaller.objective, larger.objective),
        chi,
        ddof,
        ChiSquare.sf(chi, ddof)
      )

  private def refuse(smaller: LmmFit, larger: LmmFit): FitResult[Unit] =
    if smaller.n != larger.n || smaller.y != larger.y then
      unavailable("models were not fitted to the same response values")
    else if smaller.options.criterion != larger.options.criterion then
      unavailable("models mix REML and ML fit criteria; refit with a common criterion")
    else if smaller.options.criterion == Criterion.REML && smaller.feNames != larger.feNames then
      unavailable("REML likelihood-ratio tests require identical fixed effects; refit with ML")
    else if !fixedNested(smaller, larger) then unavailable("fixed-effect column spaces are not nested")
    else if !randomNested(smaller, larger) then unavailable("random-effect term structures are not nested")
    else if larger.dof <= smaller.dof then
      unavailable("larger model must have more degrees of freedom than the smaller model")
    else if larger.loglikelihood - smaller.loglikelihood < -LogLikTol then
      unavailable("Log-likelihood must not be lower in models with more degrees of freedom")
    else Right(())

  private def fixedNested(smaller: LmmFit, larger: LmmFit): Boolean =
    smaller.feNames.toSet.subsetOf(larger.feNames.toSet)

  private def randomNested(smaller: LmmFit, larger: LmmFit): Boolean =
    smaller.formula.randomTerms.forall: small =>
      larger.formula.randomTerms.exists(large => termCovers(large, small))

  private def termCovers(large: RandomTerm, small: RandomTerm): Boolean =
    groupingCompatible(large.grouping, small.grouping) &&
      columnsOf(small).subsetOf(columnsOf(large))

  private def groupingCompatible(large: GroupingFactor, small: GroupingFactor): Boolean =
    large == small || groupingNames(small).toSet.subsetOf(groupingNames(large).toSet)

  private def groupingNames(grouping: GroupingFactor): Vector[String] =
    grouping match
      case GroupingFactor.Single(name)       => Vector(name)
      case GroupingFactor.Interaction(names) => names
      case GroupingFactor.Cell(names)        => names

  private def columnsOf(term: RandomTerm): Set[String] =
    val cols = term.terms.flatMap:
      case FixedTerm.Column(name)                      => Vector(name)
      case FixedTerm.Interaction(names)                => Vector(names.mkString(":"))
      case FixedTerm.Intercept | FixedTerm.NoIntercept =>
        Vector.empty
    if cols.isEmpty then Set("1") else cols.toSet

  private def unavailable(details: String): FitResult[Unit] =
    Left(MixedModelError.InferenceUnavailable("lrt_unavailable", details))
