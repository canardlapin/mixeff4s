package mixeff4s.pathology

/** Truth-level synthetic LMM spec. Certification is linear algebra on this record, not a fit. */
final case class GeneratorSpec(
    label: String,
    groupSizes: Vector[Int],
    nFePredictors: Int,
    hasIntercept: Boolean,
    nReSlopes: Int,
    reCovTruth: Vector[Vector[Double]]
):
  def n: Int = groupSizes.sum
  def reDim: Int = 1 + nReSlopes
  def nTheta: Int = reDim * (reDim + 1) / 2
  def feRank: Int = (if hasIntercept then 1 else 0) + nFePredictors
  def nParams: Int = feRank + nTheta + 1
  def minGroupSize: Int = if groupSizes.isEmpty then 0 else groupSizes.min
  def maxGroupSize: Int = if groupSizes.isEmpty then 0 else groupSizes.max

object GeneratorSpec:
  def lmm(
      label: String,
      groupSizes: Vector[Int],
      nFePredictors: Int,
      nReSlopes: Int,
      reCovTruth: Vector[Vector[Double]]
  ): GeneratorSpec =
    GeneratorSpec(label, groupSizes, nFePredictors, hasIntercept = true, nReSlopes, reCovTruth)
