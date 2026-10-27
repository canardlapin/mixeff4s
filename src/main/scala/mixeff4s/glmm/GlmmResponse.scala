package mixeff4s.glmm

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.model.{Family, Link}

/** Binomial responses are proportions; trials are explicit integers. Case weights multiply both the conditional log
  * density and its score/curvature.
  */
private[mixeff4s] final case class GlmmResponse(
    family: Family,
    y: Vector[Double],
    caseWeights: Vector[Double],
    trials: Vector[Int],
    offset: Vector[Double]
):
  val workingWeights: Vector[Double] = caseWeights.zip(trials).map((w, n) => w * n.toDouble)
  private val counts = y.indices
    .map(i => if family == Family.Bernoulli || family == Family.Binomial then math.rint(y(i) * trials(i)) else y(i))
    .toVector
  def scoreResidual(i: Int, eta: Double, mean: Double): Double =
    val residual = family match
      case Family.Bernoulli | Family.Binomial if eta >= 0.0 =>
        val e = math.exp(-eta)
        (counts(i) - trials(i)) + trials(i) * e / (1.0 + e)
      case _ => counts(i) - trials(i) * mean
    caseWeights(i) * residual
  private val constants = y.indices.map: i =>
    family match
      case Family.Binomial =>
        val n = trials(i).toDouble
        val successes = counts(i)
        if successes == 0.0 || successes == n then 0.0
        else logGamma(n + 1.0) - logGamma(successes + 1.0) - logGamma(n - successes + 1.0)
      case Family.Poisson => -logGamma(y(i) + 1.0)
      case _              => 0.0

  /** Parameter-dependent part only, for response-constant-free quadrature differences. */
  private[glmm] def parameterDeviance(i: Int, eta: Double): Double =
    val loss = family match
      case Family.Bernoulli | Family.Binomial =>
        val n = trials(i).toDouble
        val s = counts(i)
        if eta >= 0.0 then (n - s) * eta + n * math.log1p(math.exp(-eta))
        else -s * eta + n * math.log1p(math.exp(eta))
      case Family.Poisson => math.exp(eta) - y(i) * eta
      case _              => Double.NaN
    2.0 * caseWeights(i) * loss

  def minusTwoLogDensity(eta: Array[Double]): Double =
    var result = 0.0
    var i = 0
    while i < y.length do
      result += parameterDeviance(i, eta(i)) - 2.0 * caseWeights(i) * constants(i)
      i += 1
    result

  // Lanczos, g=7. All arguments here are positive, so no reflection is needed.
  private def logGamma(x: Double): Double =
    if x == 1.0 || x == 2.0 then return 0.0
    val coefficients = Array(676.5203681218851, -1259.1392167224028, 771.3234287776531, -176.6150291621406,
      12.507343278686905, -0.13857109526572012, 9.984369578019572e-6, 1.5056327351493116e-7)
    val z = x - 1.0
    var a = 0.9999999999998099
    var j = 0
    while j < coefficients.length do
      a += coefficients(j) / (z + j + 1.0)
      j += 1
    val t = z + 7.5
    0.9189385332046727 + (z + 0.5) * math.log(t) - t + math.log(a)

private[mixeff4s] object GlmmResponse:
  def prepare(family: Family, y: Vector[Double], options: GlmmOptions): FitResult[GlmmResponse] =
    for
      _ <- Pirls.validateResponse(family, y)
      weights = options.weights.getOrElse(Vector.fill(y.length)(1.0))
      _ <- Pirls.validateCaseWeights(weights, y.length)
      offset = options.offset.getOrElse(Vector.fill(y.length)(0.0))
      _ <- Pirls.validateOffset(offset, y.length)
      trials = options.trials.getOrElse(Vector.fill(y.length)(1))
      _ <-
        if options.trials.nonEmpty && family != Family.Binomial then
          Left(MixedModelError.InvalidArgument("trials are accepted only for Family.Binomial"))
        else if trials.length != y.length || trials.exists(_ <= 0) then
          Left(MixedModelError.InvalidArgument("trials must be positive integers, one per observation"))
        else if family == Family.Binomial && y.indices.exists: i =>
            val successes = y(i) * trials(i).toDouble
            math.abs(successes - math.rint(successes)) > math.max(1e-8, 4.0 * math.ulp(successes))
        then
          Left(
            MixedModelError.InvalidArgument(
              "binomial proportions must correspond to integer successes at the supplied trials"
            )
          )
        else if weights.zip(trials).exists((w, n) => !(w * n.toDouble).isFinite) then
          Left(MixedModelError.InvalidArgument("case weight times trials must be finite"))
        else Right(())
      canonical = y.indices.map: i =>
        family match
          case Family.Bernoulli => math.rint(y(i))
          case Family.Binomial  => math.rint(y(i) * trials(i)) / trials(i).toDouble
          case Family.Poisson   => math.rint(y(i))
          case _                => y(i)
    yield GlmmResponse(family, canonical.toVector, weights, trials, offset)

  def validatePath(family: Family, link: Link, options: GlmmOptions): FitResult[Unit] =
    if options.nAgq < 1 || options.nAgq > AdaptiveQuadrature.MaxNodes then
      Left(MixedModelError.Unsupported(s"nAgq must be between 1 and ${AdaptiveQuadrature.MaxNodes}"))
    else if options.fast && options.nAgq != 1 then Left(MixedModelError.Unsupported("AGQ requires fast=false"))
    else if options.maxEvaluations <= 0 || options.pirlsMaxIterations <= 0 || options.jointMaxIterations <= 0 then
      Left(MixedModelError.InvalidArgument("optimizer and PIRLS iteration budgets must be positive"))
    else if !options.conditionalScoreTolerance.isFinite || options.conditionalScoreTolerance <= 0.0 then
      Left(MixedModelError.InvalidArgument("conditionalScoreTolerance must be finite and positive"))
    else if !options.stationarityTolerance.isFinite || options.stationarityTolerance <= 0.0 then
      Left(MixedModelError.InvalidArgument("stationarityTolerance must be finite and positive"))
    else if family.hasDispersion then
      Left(MixedModelError.Unsupported("estimated-dispersion GLMM families require a qualified dispersion contract"))
    else if !options.fast && !((family == Family.Bernoulli || family == Family.Binomial) && link == Link.Logit ||
        family == Family.Poisson && link == Link.Log)
    then
      Left(
        MixedModelError.Unsupported("joint Laplace currently supports binomial/Bernoulli logit and Poisson log only")
      )
    else Right(())
