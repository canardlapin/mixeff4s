package mixeff4s.model

/** Conditional distribution families for GLMMs. */
enum Family:
  case Normal, Bernoulli, Binomial, Poisson, NegativeBinomial, Gamma, InverseGaussian

  def label: String =
    this match
      case Family.Normal            => "normal"
      case Family.Bernoulli         => "bernoulli"
      case Family.Binomial          => "binomial"
      case Family.Poisson           => "poisson"
      case Family.NegativeBinomial  => "negative_binomial"
      case Family.Gamma             => "gamma"
      case Family.InverseGaussian   => "inverse_gaussian"

  def canonicalLink: Link =
    this match
      case Family.Normal           => Link.Identity
      case Family.Bernoulli        => Link.Logit
      case Family.Binomial         => Link.Logit
      case Family.Poisson          => Link.Log
      case Family.NegativeBinomial => Link.Log
      case Family.Gamma            => Link.Inverse
      case Family.InverseGaussian  => Link.Inverse

  def hasDispersion: Boolean =
    this == Family.Normal || this == Family.Gamma || this == Family.InverseGaussian

  def variance(mu: Double): Double =
    this match
      case Family.Normal           => 1.0
      case Family.Bernoulli        => mu * (1.0 - mu)
      case Family.Binomial         => mu * (1.0 - mu)
      case Family.Poisson          => mu
      case Family.NegativeBinomial => mu + mu * mu
      case Family.Gamma            => mu * mu
      case Family.InverseGaussian  => mu * mu * mu

  def allows(link: Link): Boolean =
    this match
      case Family.Normal =>
        link == Link.Log || link == Link.Inverse || link == Link.Sqrt
      case Family.Bernoulli | Family.Binomial =>
        link == Link.Logit || link == Link.Probit || link == Link.Cloglog
      case Family.Poisson =>
        link == Link.Log || link == Link.Sqrt
      case Family.NegativeBinomial =>
        link == Link.Log
      case Family.Gamma | Family.InverseGaussian =>
        link == Link.Log || link == Link.Inverse

enum Link:
  case Identity, Log, Logit, Probit, Cloglog, Inverse, Sqrt

  def label: String =
    this match
      case Link.Identity => "identity"
      case Link.Log      => "log"
      case Link.Logit    => "logit"
      case Link.Probit   => "probit"
      case Link.Cloglog  => "cloglog"
      case Link.Inverse  => "inverse"
      case Link.Sqrt     => "sqrt"

  def link(mu: Double): Double =
    this match
      case Link.Identity => mu
      case Link.Log      => math.log(mu)
      case Link.Logit    => math.log(mu / (1.0 - mu))
      case Link.Probit   => Dist.inverseCdf(mu)
      case Link.Cloglog  => math.log(-math.log1p(-mu))
      case Link.Inverse  => 1.0 / mu
      case Link.Sqrt     => math.sqrt(mu)

  def linkinv(eta: Double): Double =
    this match
      case Link.Identity => eta
      case Link.Log      => math.exp(eta)
      case Link.Logit =>
        val e = math.exp(eta)
        e / (1.0 + e)
      case Link.Probit => Dist.cdf(eta)
      case Link.Cloglog =>
        -math.expm1(-math.exp(eta))
      case Link.Inverse => 1.0 / eta
      case Link.Sqrt    => eta * eta

  def muEta(eta: Double): Double =
    this match
      case Link.Identity => 1.0
      case Link.Log      => math.exp(eta)
      case Link.Logit =>
        val e = math.exp(eta)
        e / math.pow(1.0 + e, 2)
      case Link.Probit => Dist.pdf(eta)
      case Link.Cloglog =>
        if eta == Double.PositiveInfinity then 0.0
        else math.exp(eta - math.exp(eta))
      case Link.Inverse => -1.0 / (eta * eta)
      case Link.Sqrt    => 2.0 * eta

/** Standard-normal helpers for the probit link. */
private object Dist:
  def pdf(x: Double): Double =
    math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.Pi)

  def cdf(x: Double): Double =
    0.5 * (1.0 + erf(x / math.sqrt(2.0)))

  /** Acklam's inverse normal CDF, public-domain coefficients. */
  def inverseCdf(p: Double): Double =
    if p <= 0.0 then Double.NegativeInfinity
    else if p >= 1.0 then Double.PositiveInfinity
    else if p < 0.02425 then
      val q = math.sqrt(-2.0 * math.log(p))
      (((((-7.784894002430293e-3 * q - 3.223964580411365e-1) * q - 2.400758277161838) * q
        - 2.549732539343734) * q + 4.374664141464968) * q + 2.938163982698783)
        / ((((7.784695709041462e-3 * q + 3.224671290700398e-1) * q + 2.445134137142996) * q
          + 3.754408661907416) * q + 1.0)
    else if p > 1.0 - 0.02425 then
      val q = math.sqrt(-2.0 * math.log(1.0 - p))
      -((((( -7.784894002430293e-3 * q - 3.223964580411365e-1) * q - 2.400758277161838) * q
        - 2.549732539343734) * q + 4.374664141464968) * q + 2.938163982698783)
        / ((((7.784695709041462e-3 * q + 3.224671290700398e-1) * q + 2.445134137142996) * q
          + 3.754408661907416) * q + 1.0)
    else
      val q = p - 0.5
      val r = q * q
      (((((( -3.969683028665376e1 * r + 2.209460984213129e2) * r - 2.759285104469687e2) * r
        + 1.383577509590705e2) * r - 3.066479806614716e1) * r + 2.506628277459239) * q)
        / (((((-5.447609879822406e1 * r + 1.615858368580409e2) * r - 1.556989798598866e2) * r
          + 6.680131188771972e1) * r - 1.328068071618342e1) * r + 1.0)

  /** Abramowitz & Stegun 7.1.26. */
  private def erf(x: Double): Double =
    val sign = if x < 0 then -1.0 else 1.0
    val ax = math.abs(x)
    val t = 1.0 / (1.0 + 0.3275911 * ax)
    val y =
      1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * math
        .exp(-ax * ax)
    sign * y
