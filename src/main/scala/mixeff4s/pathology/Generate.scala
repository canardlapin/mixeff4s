package mixeff4s.pathology

import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.model.{Family, Link}

/** Drawn frame plus the formula that matches the spec. Generation does not fit. */
final case class Generated(frame: ModelFrame, formula: String)

/** Deterministic draws from a truth spec. Bernoulli+Logit is supported; other GLMM families are refused. */
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
    else
      SymmetricPsd.sqrt(cov).flatMap: sqrtSigma =>
        val rng = SplitMix64(spec.seed)
        val nPrimary = spec.groupSizes.length
        val primaryRe = Vector.tabulate(nPrimary): _ =>
          val z = Vector.fill(q)(rng.nextGaussian())
          mul(sqrtSigma, z)
        val secondaryRe = spec.crossed match
          case None => Vector.empty
          case Some(crossed) =>
            val sd = math.sqrt(math.max(crossed.reVar, 0.0))
            Vector.fill(crossed.nLevels)(sd * rng.nextGaussian())
        val beta = spec.beta
        val y = Vector.newBuilder[Double]
        val groups = Vector.newBuilder[String]
        val secondaryGroups = Vector.newBuilder[String]
        val predictors = Vector.fill(spec.nFePredictors)(Vector.newBuilder[Double])
        var responseError: Option[MixedModelError] = None
        def emit(gIdx: Int, hIdx: Option[Int]): Unit =
          if responseError.isDefined || gIdx < 0 || gIdx >= nPrimary then ()
          else
            val x = Vector.fill(spec.nFePredictors)(rng.nextGaussian())
            var eta = if spec.hasIntercept then beta.headOption.getOrElse(0.0) else 0.0
            x.zipWithIndex.foreach: (xj, j) =>
              eta += beta.lift(j + (if spec.hasIntercept then 1 else 0)).getOrElse(0.0) * xj
            val u = primaryRe(gIdx)
            var rePos = 0
            if spec.hasIntercept then
              eta += u(rePos)
              rePos += 1
            x.take(spec.nReSlopes)
              .foreach: xj =>
                eta += u(rePos) * xj
                rePos += 1
            hIdx.foreach: h =>
              if h >= 0 && h < secondaryRe.length then eta += secondaryRe(h)
            sampleResponse(spec, eta, rng) match
              case Left(err) =>
                responseError = Some(err)
              case Right(yi) =>
                y += yi
                groups += f"g${gIdx + 1}%03d"
                hIdx.foreach(h => secondaryGroups += f"h${h + 1}%03d")
                x.zipWithIndex.foreach: (xj, j) =>
                  predictors(j) += xj
        spec.crossedCells match
          case Some(cells) =>
            cells.foreach: (gIdx, hIdx) =>
              emit(gIdx, Some(hIdx))
          case None =>
            spec.groupSizes.zipWithIndex.foreach: (groupN, gIdx) =>
              var i = 0
              while i < groupN && responseError.isEmpty do
                emit(gIdx, None)
                i += 1
        responseError match
          case Some(err) => Left(err)
          case None =>
            val cols = Vector.newBuilder[(String, Column)]
            cols += spec.responseName -> ModelFrame.numeric(y.result())
            predictors.zipWithIndex.foreach: (col, j) =>
              cols += s"x${j + 1}" -> ModelFrame.numeric(col.result())
            cols += spec.groupName -> ModelFrame.factor(groups.result())
            spec.crossed.foreach: crossed =>
              cols += crossed.name -> ModelFrame.factor(secondaryGroups.result())
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
    val secondary = spec.crossed.fold("")(c => s" + (1 | ${c.name})")
    s"${spec.responseName} ~ $fePart + ($reInner | ${spec.groupName})$secondary"

  private def sampleResponse(spec: GeneratorSpec, eta: Double, rng: SplitMix64): FitResult[Double] =
    (spec.family, spec.link) match
      case (Family.Normal, Link.Identity) =>
        Right(eta + spec.residualSd * rng.nextGaussian())
      case (Family.Bernoulli, Link.Logit) =>
        val z = eta + spec.binaryInterceptShift
        val p = 1.0 / (1.0 + math.exp(-z))
        Right(if rng.nextUnit() < p then 1.0 else 0.0)
      case (family, link) =>
        Left(MixedModelError.UnsupportedFamilyLink(family.label, link.label))

  private def mul(matrix: Vector[Vector[Double]], z: Vector[Double]): Vector[Double] =
    matrix.map(row => row.zip(z).map(_ * _).sum)

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
