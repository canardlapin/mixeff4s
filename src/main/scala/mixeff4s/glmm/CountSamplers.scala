package mixeff4s.glmm

import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.lmm.Simulation

/** Deterministic, stateless count samplers.
  *
  * Every proposal is addressed by an indexed SplitMix stream, so rejection in one observation cannot perturb another
  * observation's draw. The BTPE and PTRS branches are direct Scala adaptations of NumPy's distributions.c
  * (NCSA/3-Clause BSD; attribution recorded with the frozen inputs).
  */
private[glmm] object CountSamplers:
  private val MaxAttempts = 100000
  private val LogSqrtTwoPi = 0.9189385332046727

  def bernoulli(p: Double, seed: Long): FitResult[Double] =
    if !p.isFinite || p < 0.0 || p > 1.0 then invalid("Bernoulli probability must be finite and in [0, 1].")
    else Right(if uniform(seed, 0L) < p then 1.0 else 0.0)

  def binomial(n: Int, p: Double, seed: Long): FitResult[Double] =
    if n < 0 then invalid("binomial trials must be non-negative")
    else if !p.isFinite || p < 0.0 || p > 1.0 then invalid("binomial probability must be finite and in [0, 1]")
    else if n == 0 || p == 0.0 then Right(0.0)
    else if p == 1.0 then Right(n.toDouble)
    else
      val r = math.min(p, 1.0 - p)
      val draw = if n.toDouble * r < 30.0 then inversion(n, r, seed) else btpe(n, p, seed)
      draw.map(value => if p > 0.5 then n.toDouble - value else value)

  def poisson(lambda: Double, seed: Long): FitResult[Double] =
    if !lambda.isFinite || lambda < 0.0 then invalid("Poisson mean must be finite and non-negative")
    else if lambda > 1e6 then
      Left(
        MixedModelError
          .InferenceUnavailable("simulation_lambda_limit", "Poisson means above 1e6 are not admitted for simulation.")
      )
    else if lambda == 0.0 then Right(0.0)
    else if lambda < 10.0 then poissonInversion(lambda, seed)
    else poissonPtrs(lambda, seed)

  private def inversion(n: Int, p: Double, seed: Long): FitResult[Double] =
    val q = 1.0 - p
    var px = math.exp(n.toDouble * math.log1p(-p))
    var u = uniform(seed, 0L)
    var x = 0
    var draws = 1
    // Full finite inversion: there is deliberately no tail restart.
    while u > px do
      if draws >= MaxAttempts || x >= n then return exhausted("binomial inversion")
      x += 1
      u -= px
      px = ((n.toDouble - x.toDouble + 1.0) * p * px) / (x.toDouble * q)
      if !px.isFinite || px < 0.0 then return exhausted("binomial inversion")
      draws += 1
    Right(x.toDouble)

  private def btpe(n: Int, p: Double, seed: Long): FitResult[Double] =
    val r = math.min(p, 1.0 - p)
    val q = 1.0 - r
    val fm = n.toDouble * r + r
    val m = math.floor(fm).toLong
    val p1 = math.floor(2.195 * math.sqrt(n.toDouble * r * q) - 4.6 * q) + 0.5
    val xm = m.toDouble + 0.5
    val xl = xm - p1
    val xr = xm + p1
    val c = 0.134 + 20.5 / (15.3 + m.toDouble)
    val left = (fm - xl) / (fm - xl * r)
    val laml = left * (1.0 + left / 2.0)
    val right = (xr - fm) / (xr * q)
    val lamr = right * (1.0 + right / 2.0)
    val p2 = p1 * (1.0 + 2.0 * c)
    val p3 = p2 + c / laml
    val p4 = p3 + c / lamr
    var attempt = 0
    while attempt < MaxAttempts do
      val u = uniform(seed, attempt.toLong * 2L) * p4
      var v = uniform(seed, attempt.toLong * 2L + 1L)
      var y = 0L
      var accepted = false
      if u <= p1 then
        y = math.floor(xm - p1 * v + u).toLong
        accepted = true
      else if u <= p2 then
        val x = xl + (u - p1) / c
        v = v * c + 1.0 - math.abs(m.toDouble - x + 0.5) / p1
        if v <= 1.0 then y = math.floor(x).toLong else y = -1L
      else if u <= p3 then
        y = math.floor(xl + math.log(v) / laml).toLong
        if y >= 0L && v != 0.0 then v = v * (u - p2) * laml else y = -1L
      else
        y = math.floor(xr - math.log(v) / lamr).toLong
        if y <= n.toLong && v != 0.0 then v = v * (u - p3) * lamr else y = -1L
      if !accepted && y >= 0L && y <= n.toLong then
        val k = math.abs(y - m)
        val nrq = n.toDouble * r * q
        if k <= 20L || k >= nrq / 2.0 - 1.0 then
          val s = r / q
          val a = s * (n.toDouble + 1.0)
          var f = 1.0
          var i = m + 1L
          while i <= y do { f *= a / i.toDouble - s; i += 1L }
          i = y + 1L
          while i <= m do { f /= a / i.toDouble - s; i += 1L }
          accepted = v <= f
        else
          val rho = (k.toDouble / nrq) * ((k.toDouble * (k.toDouble / 3.0 + 0.625) + 1.0 / 6.0) / nrq + 0.5)
          val t = -k.toDouble * k.toDouble / (2.0 * nrq)
          val a = math.log(v)
          if a < t - rho then accepted = true
          else if a <= t + rho then
            val x1 = y.toDouble + 1.0; val f1 = m.toDouble + 1.0
            val z = n.toDouble + 1.0 - m.toDouble; val w = n.toDouble - y.toDouble + 1.0
            def stirling(x: Double): Double =
              val x2 = x * x
              (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x / 166320.0
            val threshold = xm * math.log(f1 / x1) + (n.toDouble - m.toDouble + 0.5) * math.log(z / w) +
              (y.toDouble - m.toDouble) * math.log(w * r / (x1 * q)) + stirling(f1) + stirling(z) - stirling(
                x1
              ) - stirling(w)
            accepted = a <= threshold
      if accepted then return Right(y.toDouble)
      attempt += 1
    exhausted("binomial BTPE")

  private def poissonInversion(lambda: Double, seed: Long): FitResult[Double] =
    val threshold = math.exp(-lambda)
    var product = 1.0
    var x = 0L
    while x < MaxAttempts do
      product *= uniform(seed, x)
      if product <= threshold then return Right(x.toDouble)
      x += 1L
    exhausted("Poisson inversion")

  private def poissonPtrs(lambda: Double, seed: Long): FitResult[Double] =
    val slam = math.sqrt(lambda)
    val b = 0.931 + 2.53 * slam; val a = -0.059 + 0.02483 * b
    val invAlpha = 1.1239 + 1.1328 / (b - 3.4); val vr = 0.9277 - 3.6224 / (b - 2.0)
    var attempt = 0
    while attempt < MaxAttempts do
      val u = uniform(seed, attempt.toLong * 2L) - 0.5
      val v = uniform(seed, attempt.toLong * 2L + 1L)
      val us = 0.5 - math.abs(u)
      val k = math.floor((2.0 * a / us + b) * u + lambda + 0.43).toLong
      if us >= 0.07 && v <= vr && k >= 0L then return Right(k.toDouble)
      if k >= 0L && !(us < 0.013 && v > us) &&
        math.log(v) + math.log(invAlpha) - math.log(a / (us * us) + b) <= poissonLogPmf(k.toDouble, lambda)
      then return Right(k.toDouble)
      attempt += 1
    exhausted("Poisson PTRS")

  private def uniform(seed: Long, index: Long): Double = Simulation.uniform(Simulation.seedFor(seed, index))

  /** Stable Poisson log PMF: `-bd0(k, lambda) - stirlerr(k) - log(sqrt(2*pi*k))`. This avoids the cancellation in
    * `-lambda + k*log(lambda) - logGamma(k+1)` near a large mean.
    */
  private[glmm] def poissonLogPmf(k: Double, lambda: Double): Double =
    if k == 0.0 then -lambda
    else -deviancePart(k, lambda) - stirlingError(k) - 0.5 * math.log(2.0 * math.Pi * k)

  private def deviancePart(x: Double, mean: Double): Double =
    val delta = x - mean
    if math.abs(delta) < 0.1 * (x + mean) then
      val v = delta / (x + mean)
      val v2 = v * v
      var term = 2.0 * x * v
      var sum = delta * v
      var j = 1
      while j < 100 do
        term *= v2
        val next = term / (2 * j + 1).toDouble
        val updated = sum + next
        if updated == sum then return updated
        sum = updated
        j += 1
      sum
    else x * math.log(x / mean) + mean - x

  private def stirlingError(x: Double): Double =
    if x < 15.0 then
      val coefficients = Array(676.5203681218851, -1259.1392167224028, 771.3234287776531, -176.6150291621406,
        12.507343278686905, -0.13857109526572012, 9.984369578019572e-6, 1.5056327351493116e-7)
      val z = x - 1.0
      var a = 0.9999999999998099; var j = 0
      while j < coefficients.length do { a += coefficients(j) / (z + j + 1.0); j += 1 }
      val t = z + 7.5
      (LogSqrtTwoPi + (z + 0.5) * math.log(t) - t + math.log(a) + math.log(x)) - (x + 0.5) * math.log(
        x
      ) + x - LogSqrtTwoPi
    else
      val inv = 1.0 / x
      inv * (1.0 / 12.0 + inv * inv * (-1.0 / 360.0 + inv * inv * (1.0 / 1260.0 - inv * inv / 1680.0)))

  private def invalid[A](message: String): FitResult[A] = Left(MixedModelError.InvalidArgument(message))
  private def exhausted[A](algorithm: String): FitResult[A] =
    Left(
      MixedModelError.InferenceUnavailable(
        "simulation_sampler_budget_exhausted",
        s"$algorithm exhausted its finite rejection budget."
      )
    )
