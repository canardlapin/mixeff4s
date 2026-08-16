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
    binaryInterceptShift: Double = 0.0
):
  def n: Int = groupSizes.sum
  def reDim: Int = 1 + nReSlopes
  def nTheta: Int = reDim * (reDim + 1) / 2
  def feRank: Int = (if hasIntercept then 1 else 0) + nFePredictors
  def nParams: Int = feRank + nTheta + 1
  def minGroupSize: Int = if groupSizes.isEmpty then 0 else groupSizes.min
  def maxGroupSize: Int = if groupSizes.isEmpty then 0 else groupSizes.max
  def beta: Vector[Double] =
    if feTruth.nonEmpty then feTruth else Vector.fill(feRank)(1.0)

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

  def extremePrevalence(spec: GeneratorSpec, interceptShift: Double): GeneratorSpec =
    spec.copy(
      family = Family.Bernoulli,
      link = Link.Logit,
      binaryInterceptShift = interceptShift,
      residualSd = 0.0
    )
