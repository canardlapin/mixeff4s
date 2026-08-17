package mixeff4s.pathology

/** Correlation-form expected Fisher spectrum and the dimensionless weak-id index. */
object Fisher:
  val WeakIdThreshold = 10.0

  /** Eigenvalues of the predictor correlation block, plus a unit intercept slot, descending. */
  def correlationEigvals(spec: GeneratorSpec): Vector[Double] =
    val nPred = spec.nFePredictors
    if nPred == 0 then Vector(1.0)
    else
      val slope =
        SymmetricPsd.eigvals(spec.predictorCorr).getOrElse(Vector.fill(nPred)(1.0))
      (slope :+ 1.0).sortBy(v => -v)

  /** `n * λ_min / trace`. Scale-invariant because the spectrum is already in correlation form. */
  def weakIdScore(n: Int, eigvals: Vector[Double]): Double =
    if eigvals.isEmpty then Double.PositiveInfinity
    else
      val trace = eigvals.sum
      if !trace.isFinite || math.abs(trace) < 1e-15 then Double.PositiveInfinity
      else n.toDouble * eigvals.min.max(0.0) / trace
