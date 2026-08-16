package mixeff4s.stats

/** χ² survival function for positive integer degrees of freedom. */
private[stats] object ChiSquare:
  def sf(x: Double, df: Int): Double =
    if df < 1 || !x.isFinite then Double.NaN
    else if x <= 0.0 then 1.0
    else
      val z = x / 2.0
      if df == 1 then erfc(math.sqrt(z))
      else if df % 2 == 0 then
        var term = 1.0
        var sum = 1.0
        var k = 1
        val m = df / 2
        while k < m do
          term *= z / k.toDouble
          sum += term
          k += 1
        math.exp(-z) * sum
      else
        var term = math.sqrt(2.0 * x / math.Pi)
        var sum = term
        var k = 3
        while k < df do
          term *= x / k.toDouble
          sum += term
          k += 2
        erfc(math.sqrt(z)) + math.exp(-z) * sum

  /** Complementary error function, Abramowitz & Stegun 7.1.26. */
  private def erfc(x: Double): Double =
    val ax = math.abs(x)
    val t = 1.0 / (1.0 + 0.3275911 * ax)
    val y =
      (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * math
        .exp(-ax * ax)
    if x >= 0.0 then y else 2.0 - y
