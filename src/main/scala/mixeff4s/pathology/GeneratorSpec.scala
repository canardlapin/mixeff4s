package mixeff4s.pathology

import mixeff4s.model.{Family, Link}

/** Truth-level synthetic LMM spec. Certification is linear algebra on this record, not a fit. */
final case class GeneratorSpec(
    label: String,
    groupSizes: Vector[Int],
    nFePredictors: Int,
    hasIntercept: Boolean,
    nReSlopes: Int,
    reCovTruth: Vector[Vector[Double]],
    seed: Long = 1L,
    feTruth: Vector[Double] = Vector.empty,
    residualSd: Double = 1.0,
    responseName: String = "y",
    groupName: String = "g",
    family: Family = Family.Normal,
    link: Link = Link.Identity,
    binaryInterceptShift: Double = 0.0,
    feCorrMatrix: Vector[Vector[Double]] = Vector.empty,
    /** Per-predictor scale. Position 0 is x1. Missing entries default to 1. */
    feScales: Vector[Double] = Vector.empty,
    crossed: Option[CrossedSpec] = None
):
  def n: Int = crossedCells.fold(groupSizes.sum)(_.length)
  def reDim: Int = 1 + nReSlopes
  def nTheta: Int = reDim * (reDim + 1) / 2 + crossed.fold(0)(_ => 1)
  def feRank: Int = (if hasIntercept then 1 else 0) + nFePredictors
  def nParams: Int = feRank + nTheta + 1
  def primarySizes: Vector[Int] =
    crossedCells match
      case None        => groupSizes
      case Some(cells) =>
        val counts = Array.fill(groupSizes.length)(0)
        cells.foreach: (i, _) =>
          if i >= 0 && i < counts.length then counts(i) += 1
        counts.toVector
  def minGroupSize: Int =
    val sizes = primarySizes.filter(_ > 0)
    if sizes.isEmpty then 0 else sizes.min
  def maxGroupSize: Int =
    val sizes = primarySizes.filter(_ > 0)
    if sizes.isEmpty then 0 else sizes.max
  def crossedCells: Option[Vector[(Int, Int)]] =
    crossed.map: c =>
      c.cells.getOrElse:
        (0 until groupSizes.length).flatMap(i => (0 until c.nLevels).map(j => (i, j))).toVector
  def crossedSummary: Option[CrossedSummary] =
    crossed.flatMap: c =>
      crossedCells.map(cells => Crossing.summarise(groupSizes.length, c.nLevels, cells))
  def beta: Vector[Double] =
    if feTruth.nonEmpty then feTruth else Vector.fill(feRank)(1.0)
  def predictorCorr: Vector[Vector[Double]] =
    if feCorrMatrix.length == nFePredictors && feCorrMatrix.forall(_.length == nFePredictors) then feCorrMatrix
    else Vector.tabulate(nFePredictors, nFePredictors)((i, j) => if i == j then 1.0 else 0.0)
  def predictorScales: Vector[Double] =
    Vector.tabulate(nFePredictors)(j => feScales.lift(j).getOrElse(1.0))
  def predictorCov: Vector[Vector[Double]] =
    val scales = predictorScales
    val corr = predictorCorr
    Vector.tabulate(nFePredictors, nFePredictors): (i, j) =>
      scales(i) * corr(i)(j) * scales(j)

object GeneratorSpec:
  def lmm(
      label: String,
      groupSizes: Vector[Int],
      nFePredictors: Int,
      nReSlopes: Int,
      reCovTruth: Vector[Vector[Double]],
      seed: Long = 1L,
      feTruth: Vector[Double] = Vector.empty,
      residualSd: Double = 1.0
  ): GeneratorSpec =
    GeneratorSpec(
      label,
      groupSizes,
      nFePredictors,
      hasIntercept = true,
      nReSlopes,
      reCovTruth,
      seed,
      feTruth,
      residualSd
    )

  def bernoulliLogit(spec: GeneratorSpec): GeneratorSpec =
    spec.copy(family = Family.Bernoulli, link = Link.Logit, residualSd = 0.0)

  def extremePrevalence(spec: GeneratorSpec, interceptShift: Double): GeneratorSpec =
    bernoulliLogit(spec).copy(binaryInterceptShift = interceptShift)

  def fullCross(spec: GeneratorSpec, name: String, nLevels: Int, reVar: Double): GeneratorSpec =
    spec.copy(crossed = Some(Crossing.fullCross(name, nLevels, reVar)))

  def emptyCrossings(
      spec: GeneratorSpec,
      name: String,
      nSecondary: Int,
      reVar: Double,
      density: Double,
      seed: Long
  ): GeneratorSpec =
    spec.copy(crossed = Some(Crossing.emptyCrossings(spec.groupSizes.length, name, nSecondary, reVar, density, seed)))

  def sparseConnectedCrossings(spec: GeneratorSpec, name: String, nLevels: Int, reVar: Double): GeneratorSpec =
    spec.copy(groupSizes = Vector.fill(nLevels)(1), crossed = Some(Crossing.sparsePath(name, nLevels, reVar)))

  def blockDiagonalCrossings(
      spec: GeneratorSpec,
      name: String,
      blockSize: Int,
      nBlocks: Int,
      reVar: Double
  ): GeneratorSpec =
    val (sizes, crossed) = Crossing.blockDiagonal(name, blockSize, nBlocks, reVar)
    spec.copy(groupSizes = sizes, crossed = Some(crossed))

  def scaleMismatch(spec: GeneratorSpec, scales: Vector[Double]): GeneratorSpec =
    spec.copy(feScales = scales)

  def setGroupSizes(spec: GeneratorSpec, sizes: Vector[Int]): GeneratorSpec =
    spec.copy(groupSizes = sizes)

  def singletonsWithSlope(spec: GeneratorSpec, nGroups: Int): GeneratorSpec =
    spec.copy(groupSizes = Vector.fill(nGroups.max(0))(1))

  def nearSingularRe(spec: GeneratorSpec, rho: Double): GeneratorSpec =
    val cov = spec.reCovTruth
    if cov.length < 2 then spec
    else
      val off = rho * math.sqrt(cov(0)(0) * cov(1)(1))
      spec.copy(reCovTruth =
        cov
          .updated(0, cov(0).updated(1, off))
          .updated(1, cov(1).updated(0, off))
      )

  /** Inverse-CDF Pareto Type I (`xm = 1`), clamped to `[1, meanSize · 50]`. Seeds do not match rust `StdRng`. */
  def paretoSizes(seed: Long, nGroups: Int, alpha: Double, meanSize: Double): Vector[Int] =
    if nGroups <= 0 then Vector.empty
    else
      val rng = SizeRng(seed)
      val cap = (meanSize * 50.0).max(2.0)
      val denom = alpha / (alpha - 1.0).max(0.5)
      Vector.fill(nGroups):
        val u = rng.nextUnit()
        val x = 1.0 / math.pow(1.0 - u, 1.0 / alpha)
        val raw = x * meanSize / denom
        math.round(raw.max(1.0).min(cap)).toInt

  def collinearFe(spec: GeneratorSpec, i: Int, j: Int, rho: Double): GeneratorSpec =
    val n = spec.nFePredictors
    if i == j || i < 0 || j < 0 || i >= n || j >= n then spec
    else
      val base = spec.predictorCorr
      spec.copy(
        feCorrMatrix = base
          .updated(i, base(i).updated(j, rho))
          .updated(j, base(j).updated(i, rho))
      )

  private final class SizeRng(private var state: Long):
    def nextLong(): Long =
      state += 0x9e3779b97f4a7c15L
      var z = state
      z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L
      z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL
      z ^ (z >>> 31)

    def nextUnit(): Double =
      (nextLong() >>> 11).toDouble * (1.0 / (1L << 53).toDouble)
