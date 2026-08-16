package mixeff4s.pathology

import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.error.{FitResult, MixedModelError}

/** Drawn frame plus the formula that matches the spec. Generation does not fit. */
final case class Generated(frame: ModelFrame, formula: String)

/** Deterministic Gaussian LMM draws from a truth spec. q > 2 and GLMM families are refused. */
object Generate:
  def apply(spec: GeneratorSpec): FitResult[Generated] =
    val q = spec.reDim
    val cov = spec.reCovTruth
    if spec.groupSizes.isEmpty || spec.n == 0 then
      Left(MixedModelError.InvalidArgument(s"spec '${spec.label}' has no observations"))
    else if spec.residualSd < 0.0 then
      Left(MixedModelError.InvalidArgument(s"spec '${spec.label}' has a negative residual sd"))
    else if spec.nReSlopes > spec.nFePredictors then
      Left(
        MixedModelError.InvalidArgument(
          s"spec '${spec.label}' requests ${spec.nReSlopes} random slopes but only ${spec.nFePredictors} fixed-effect predictors exist"
        )
      )
    else if spec.feTruth.nonEmpty && spec.feTruth.length != spec.feRank then
      Left(
        MixedModelError.InvalidArgument(
          s"spec '${spec.label}' fe_truth length ${spec.feTruth.length} does not match fe rank ${spec.feRank}"
        )
      )
    else if cov.length != q || cov.exists(_.length != q) then
      Left(
        MixedModelError.InvalidArgument(
          s"re_cov_truth dim mismatch in spec '${spec.label}': expected ${q}×$q got ${cov.length}×${cov.headOption.map(_.length).getOrElse(0)}"
        )
      )
    else if q > 2 then Left(MixedModelError.InvalidArgument(s"generation for re_dim = $q is not implemented"))
    else
      sqrtPsd(cov).flatMap: sqrtSigma =>
        val rng = SplitMix64(spec.seed)
        val primaryRe = Vector.tabulate(spec.groupSizes.length): _ =>
          val z = Vector.fill(q)(rng.nextGaussian())
          mul(sqrtSigma, z)
        val beta = spec.beta
        val y = Vector.newBuilder[Double]
        val groups = Vector.newBuilder[String]
        val predictors = Vector.fill(spec.nFePredictors)(Vector.newBuilder[Double])
        spec.groupSizes.zipWithIndex.foreach: (groupN, gIdx) =>
          val u = primaryRe(gIdx)
          val label = f"g${gIdx + 1}%03d"
          var i = 0
          while i < groupN do
            val x = Vector.fill(spec.nFePredictors)(rng.nextGaussian())
            var eta = if spec.hasIntercept then beta.headOption.getOrElse(0.0) else 0.0
            x.zipWithIndex.foreach: (xj, j) =>
              eta += beta.lift(j + (if spec.hasIntercept then 1 else 0)).getOrElse(0.0) * xj
            var rePos = 0
            if spec.hasIntercept then
              eta += u(rePos)
              rePos += 1
            x.take(spec.nReSlopes)
              .foreach: xj =>
                eta += u(rePos) * xj
                rePos += 1
            y += eta + spec.residualSd * rng.nextGaussian()
            groups += label
            x.zipWithIndex.foreach: (xj, j) =>
              predictors(j) += xj
            i += 1
        val cols = Vector.newBuilder[(String, Column)]
        cols += spec.responseName -> ModelFrame.numeric(y.result())
        predictors.zipWithIndex.foreach: (col, j) =>
          cols += s"x${j + 1}" -> ModelFrame.numeric(col.result())
        cols += spec.groupName -> ModelFrame.factor(groups.result())
        ModelFrame.of(cols.result()*).map(frame => Generated(frame, formula(spec)))

  def formula(spec: GeneratorSpec): String =
    val xs = (1 to spec.nFePredictors).map(i => s"x$i")
    val fePart =
      if xs.isEmpty then "1"
      else if spec.hasIntercept then s"1 + ${xs.mkString(" + ")}"
      else xs.mkString(" + ")
    val slopes = xs.take(spec.nReSlopes)
    val reInner = (spec.hasIntercept, slopes.isEmpty) match
      case (true, true)   => "1"
      case (true, false)  => s"1 + ${slopes.mkString(" + ")}"
      case (false, false) => s"0 + ${slopes.mkString(" + ")}"
      case (false, true)  => "1"
    s"${spec.responseName} ~ $fePart + ($reInner | ${spec.groupName})"

  private def mul(matrix: Vector[Vector[Double]], z: Vector[Double]): Vector[Double] =
    matrix.map(row => row.zip(z).map(_ * _).sum)

  private def sqrtPsd(cov: Vector[Vector[Double]]): FitResult[Vector[Vector[Double]]] =
    cov.length match
      case 0 => Right(Vector.empty)
      case 1 =>
        val a = cov(0)(0)
        if a < -1e-10 then Left(MixedModelError.InvalidArgument("re_cov_truth is not PSD"))
        else Right(Vector(Vector(math.sqrt(math.max(a, 0.0)))))
      case 2 =>
        val a = cov(0)(0)
        val b = cov(1)(1)
        val c = 0.5 * (cov(0)(1) + cov(1)(0))
        val disc = math.sqrt(math.max(0.0, (a - b) * (a - b) + 4.0 * c * c))
        val l1 = (a + b + disc) / 2.0
        val l2 = (a + b - disc) / 2.0
        if l2 < -1e-10 then Left(MixedModelError.InvalidArgument("re_cov_truth is not PSD"))
        else
          val (v1x, v1y) = evec(a, c, l1)
          val (v2x, v2y) = evec(a, c, l2)
          val s1 = math.sqrt(math.max(l1, 0.0))
          val s2 = math.sqrt(math.max(l2, 0.0))
          Right(
            Vector(
              Vector(s1 * v1x * v1x + s2 * v2x * v2x, s1 * v1x * v1y + s2 * v2x * v2y),
              Vector(s1 * v1y * v1x + s2 * v2y * v2x, s1 * v1y * v1y + s2 * v2y * v2y)
            )
          )
      case q =>
        Left(MixedModelError.InvalidArgument(s"generation for re_dim = $q is not implemented"))

  private def evec(a: Double, c: Double, lam: Double): (Double, Double) =
    if math.abs(c) <= 1e-15 then if math.abs(a - lam) <= 1e-12 then (1.0, 0.0) else (0.0, 1.0)
    else
      val y = -(a - lam) / c
      val n = math.hypot(1.0, y)
      (1.0 / n, y / n)

  /** SplitMix64 plus polar Box–Muller. Portable across JVM and Scala.js. */
  private final class SplitMix64(private var state: Long):
    def nextLong(): Long =
      state += 0x9e3779b97f4a7c15L
      var z = state
      z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L
      z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL
      z ^ (z >>> 31)

    def nextUnit(): Double =
      (nextLong() >>> 11).toDouble * (1.0 / (1L << 53).toDouble)

    def nextGaussian(): Double =
      var u = 0.0
      var v = 0.0
      var s = 0.0
      while
        u = nextUnit() * 2.0 - 1.0
        v = nextUnit() * 2.0 - 1.0
        s = u * u + v * v
        s >= 1.0 || s == 0.0
      do ()
      u * math.sqrt(-2.0 * math.log(s) / s)
