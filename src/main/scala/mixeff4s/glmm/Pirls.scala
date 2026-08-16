package mixeff4s.glmm

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.model.{Family, Link}

/** Working-response and exponential-family helpers for labelled fast-PIRLS. */
object Pirls:
  val MaxIter: Int = 10
  val ConvergenceTol: Double = 1.0e-5
  val MaxHalvings: Int = 10
  private val BoundedMeanEps: Double = 1e-15
  private val LogLinkEtaBound: Double = 30.0

  def workingObservation(
      family: Family,
      link: Link,
      y: Double,
      eta: Double,
      mu: Double,
      caseWeight: Double,
      offset: Double = 0.0
  ): (Double, Double) =
    val (workingMu, etaForDerivative) = boundedMeanAndEta(family, link, mu, eta)
    val dmuDeta = link.muEta(etaForDerivative)
    val varMu = family.variance(workingMu)
    val weight =
      if dmuDeta.isFinite && varMu.isFinite && varMu > 0.0 then caseWeight * dmuDeta * dmuDeta / varMu
      else 0.0
    val resid =
      if !dmuDeta.isFinite || math.abs(dmuDeta) < 1e-15 then 0.0
      else (y - workingMu) / dmuDeta
    (math.sqrt(math.max(weight, 0.0)), eta + resid - offset)

  def boundedMeanAndEta(family: Family, link: Link, mu: Double, eta: Double): (Double, Double) =
    family match
      case Family.Bernoulli | Family.Binomial =>
        val boundedMu = mu.max(BoundedMeanEps).min(1.0 - BoundedMeanEps)
        (boundedMu, link.link(boundedMu))
      case Family.Poisson | Family.NegativeBinomial =>
        link match
          case Link.Log =>
            val boundedEta = eta.max(-LogLinkEtaBound).min(LogLinkEtaBound)
            (math.exp(boundedEta), boundedEta)
          case Link.Sqrt =>
            val boundedMu = mu.max(BoundedMeanEps)
            val minEta = math.sqrt(boundedMu)
            val boundedEta =
              if math.abs(eta) < minEta then if eta < 0.0 then -minEta else minEta
              else eta
            (boundedEta * boundedEta, boundedEta)
          case _ => (mu, eta)
      case _ => (mu, eta)

  def converged(obj: Double, accepted: Double, tol: Double = ConvergenceTol): Boolean =
    math.abs(obj - accepted) < tol

  def isBinary(value: Double): Boolean =
    math.abs(value) < 1e-12 || math.abs(value - 1.0) < 1e-12

  def isNonnegativeInteger(value: Double): Boolean =
    value >= 0.0 && math.abs(value - value.round.toDouble) < 1e-12

  def isInterceptColumn(name: String): Boolean =
    name == "1" || name == "(Intercept)" || name == "Intercept" || name == "intercept"

  def refuseFamilyLink(family: Family, link: Link): FitResult[Unit] =
    if family == Family.Normal && link == Link.Identity then
      Left(MixedModelError.UnsupportedFamilyLink(family.label, link.label))
    else if family == Family.NegativeBinomial then
      Left(
        MixedModelError.Unsupported(
          "negative-binomial GLMM requires a fixed theta; not in this fast-PIRLS slice"
        )
      )
    else if !family.allows(link) then Left(MixedModelError.UnsupportedFamilyLink(family.label, link.label))
    else Right(())

  def validateResponse(family: Family, y: Vector[Double]): FitResult[Unit] =
    y.zipWithIndex.foldLeft[FitResult[Unit]](Right(())):
      case (Left(err), _)           => Left(err)
      case (Right(_), (value, idx)) =>
        if !value.isFinite then
          Left(MixedModelError.InvalidArgument(s"response at index $idx must be finite (got $value)"))
        else
          family match
            case Family.Bernoulli if !isBinary(value) =>
              Left(
                MixedModelError.InvalidArgument(
                  s"bernoulli GLMM response must be exactly 0 or 1; index $idx has $value"
                )
              )
            case Family.Binomial if !(value >= 0.0 && value <= 1.0) && !isNonnegativeInteger(value) =>
              Left(
                MixedModelError.InvalidArgument(
                  s"binomial GLMM response must be a proportion in [0, 1] or a non-negative integer; index $idx has $value"
                )
              )
            case Family.Poisson if value < 0.0 =>
              Left(
                MixedModelError.InvalidArgument(
                  s"poisson GLMM response must be non-negative; index $idx has $value"
                )
              )
            case Family.Gamma | Family.InverseGaussian if value <= 0.0 =>
              Left(
                MixedModelError.InvalidArgument(
                  s"${family.label} GLMM response must be strictly positive; index $idx has $value"
                )
              )
            case _ => Right(())

  def validateCaseWeights(weights: Vector[Double], n: Int): FitResult[Unit] =
    if weights.length != n then
      Left(
        MixedModelError.InvalidArgument(
          s"case weights length (${weights.length}) does not match number of observations ($n)"
        )
      )
    else
      weights.zipWithIndex.collectFirst:
        case (w, i) if !w.isFinite || w <= 0.0 =>
          MixedModelError.InvalidArgument(s"case weight at index $i must be finite and positive (got $w)")
      match
        case Some(err) => Left(err)
        case None      => Right(())

  def validateOffset(offset: Vector[Double], n: Int): FitResult[Unit] =
    if offset.length != n then
      Left(
        MixedModelError.InvalidArgument(
          s"offset length (${offset.length}) does not match number of observations ($n)"
        )
      )
    else
      offset.zipWithIndex.collectFirst:
        case (v, i) if !v.isFinite =>
          MixedModelError.InvalidArgument(s"offset at index $i must be finite (got $v)")
      match
        case Some(err) => Left(err)
        case None      => Right(())

  def initialResponseMean(family: Family, y: Vector[Double], weights: Vector[Double]): Option[Double] =
    if y.isEmpty then None
    else
      var weightedSum = 0.0
      var weightSum = 0.0
      var i = 0
      while i < y.length do
        val w = if i < weights.length then weights(i) else 1.0
        weightedSum += w * y(i)
        weightSum += w
        i += 1
      if weightSum <= 0.0 then None
      else Some(initialMeanForLink(family, weightedSum / weightSum))

  def initialMeanForLink(family: Family, mean: Double): Double =
    family match
      case Family.Bernoulli | Family.Binomial => mean.max(1e-6).min(1.0 - 1e-6)
      case Family.Poisson | Family.NegativeBinomial | Family.Gamma | Family.InverseGaussian =>
        mean.max(1e-6)
      case Family.Normal => mean.max(0.0)

  def devianceComponent(family: Family, y: Double, mu: Double): Double =
    family match
      case Family.Bernoulli | Family.Binomial =>
        val m = mu.max(1e-15).min(1.0 - 1e-15)
        if y == 1.0 then -2.0 * math.log(m)
        else if y == 0.0 then -2.0 * math.log(1.0 - m)
        else 2.0 * (y * math.log(y / m) + (1.0 - y) * math.log((1.0 - y) / (1.0 - m)))
      case Family.Poisson =>
        if y == 0.0 then 2.0 * mu
        else 2.0 * (y * math.log(y / mu) - (y - mu))
      case Family.Normal =>
        val r = y - mu
        r * r
      case Family.Gamma =>
        val m = mu.max(1e-6)
        if y == 0.0 then 2.0 * math.log(m)
        else -2.0 * (math.log(y / m) - (y - m) / m)
      case Family.InverseGaussian =>
        val m = mu.max(1e-6)
        val r = y - m
        r * r / (y * m * m)
      case Family.NegativeBinomial =>
        Double.NaN
