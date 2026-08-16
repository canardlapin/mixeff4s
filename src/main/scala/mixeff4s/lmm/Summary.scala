package mixeff4s.lmm

import mixeff4s.design.ReMat

/** Variance and correlation of random effects, plus residual σ. */
final case class VarCorr(
    components: Vector[VarCorrComponent],
    residualSd: Double
)

final case class VarCorrComponent(
    group: String,
    names: Vector[String],
    stdDev: Vector[Double],
    correlations: Vector[Double]
)

object VarCorr:
  def fromReterms(reterms: Vector[ReMat], sigma: Double): VarCorr =
    val components = reterms.map: rt =>
      val s = rt.vsize
      val stdDev = Vector.tabulate(s): i =>
        var rowNormSq = 0.0
        var j = 0
        while j <= i do
          val lam = rt.lambda(i, j)
          rowNormSq += lam * lam
          j += 1
        sigma * math.sqrt(rowNormSq)
      val correlations =
        if s <= 1 then Vector.empty
        else
          val normalized = Array.fill(s, s)(0.0)
          var i = 0
          while i < s do
            val rowNorm = stdDev(i) / sigma
            if rowNorm > 0.0 then
              var j = 0
              while j <= i do
                normalized(i)(j) = rt.lambda(i, j) / rowNorm
                j += 1
            i += 1
          val corr = Vector.newBuilder[Double]
          i = 1
          while i < s do
            var j = 0
            while j < i do
              var dot = 0.0
              var k = 0
              while k <= j do
                dot += normalized(i)(k) * normalized(j)(k)
                k += 1
              corr += dot
              j += 1
            i += 1
          corr.result()
      VarCorrComponent(rt.groupingName, rt.cnames, stdDev, correlations)
    VarCorr(components, sigma)

/** Fixed-effect coefficient table. Wald z is reported; p-values are refused. */
final case class CoefTable(
    names: Vector[String],
    estimates: Vector[Double],
    stdErrors: Vector[Double],
    zValues: Vector[Double]
):
  def pValueCode: String = CoefTable.PValueUnavailableCode
  def pValueReason: String = CoefTable.PValueUnavailableReason

object CoefTable:
  val PValueUnavailableCode = "p_value_unavailable"
  val PValueUnavailableReason =
    "Asymptotic Wald p-values are not reported. Request a certified inference method."

  def wald(names: Vector[String], estimates: Vector[Double], stdErrors: Vector[Double]): CoefTable =
    val z = estimates.zip(stdErrors).map: (est, se) =>
      if se == 0.0 || !se.isFinite then Double.NaN else est / se
    CoefTable(names, estimates, stdErrors, z)
