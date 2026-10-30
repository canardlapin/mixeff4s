package mixeff4s.stats

import mixeff4s.error.{FitResult, MixedModelError}

final case class CurvatureOptions(
    relativeStep: Double = 1e-3,
    stabilityTolerance: Double = 0.02,
    minimumScaledPivot: Double = 1e-8,
    maxParameters: Int = 64
):
  private[stats] def validate: FitResult[Unit] =
    if !relativeStep.isFinite || relativeStep <= 0 || !stabilityTolerance.isFinite || stabilityTolerance <= 0 ||
      !minimumScaledPivot.isFinite || minimumScaledPivot <= 0 || minimumScaledPivot >= 1 || maxParameters < 1
    then Left(MixedModelError.InvalidArgument("invalid observed-information controls"))
    else Right(())

/** Hessian of -2 log likelihood and 2 H^-1, in the same parameter order. */
final case class CurvatureEvidence(
    hessian: Vector[Vector[Double]],
    covariance: Vector[Vector[Double]],
    steps: Vector[Double],
    hessianRelativeDifference: Double,
    covarianceRelativeDifference: Double,
    minimumScaledPivot: Double,
    evaluations: Int
)

private[stats] object ObservedInformation:
  private def refuse(code: String, message: String): FitResult[Nothing] =
    Left(MixedModelError.InferenceUnavailable(code, message))

  def compute(parameters: Vector[Double], lower: Vector[Double], options: CurvatureOptions)(
      evaluate: Vector[Double] => FitResult[Double]
  ): FitResult[CurvatureEvidence] =
    for
      _ <- options.validate
      _ <-
        if parameters.nonEmpty && parameters.length <= options.maxParameters && lower.length == parameters.length &&
          parameters.forall(_.isFinite)
        then Right(())
        else refuse("curvature_dimension", "Curvature dimension is empty, invalid or exceeds its explicit budget.")
      steps = parameters.map(x => options.relativeStep * math.max(1.0, math.abs(x)))
      _ <-
        if parameters.indices.forall(i => parameters(i) - steps(i) > lower(i)) then Right(())
        else refuse("curvature_boundary", "Central observed curvature requires interior nuisance parameters.")
      baseline <- evaluate(parameters)
      coarse <- hessian(parameters, steps, baseline)(evaluate)
      fine <- hessian(parameters, steps.map(_ / 2), baseline)(evaluate)
      first <- inverse(coarse, options.minimumScaledPivot)
      second <- inverse(fine, options.minimumScaledPivot)
      scaledCoarse = coarse.indices
        .map(i => coarse.indices.map(j => coarse(i)(j) / (math.sqrt(fine(i)(i)) * math.sqrt(fine(j)(j)))).toVector)
        .toVector
      scaledFine = fine.indices
        .map(i => fine.indices.map(j => fine(i)(j) / (math.sqrt(fine(i)(i)) * math.sqrt(fine(j)(j)))).toVector)
        .toVector
      hd = relative(scaledCoarse, scaledFine)
      scaledFirst = parameters.indices
        .map(i =>
          parameters.indices
            .map(j => first._1(i)(j) / (math.sqrt(second._1(i)(i)) * math.sqrt(second._1(j)(j))))
            .toVector
        )
        .toVector
      scaledSecond = parameters.indices
        .map(i =>
          parameters.indices
            .map(j => second._1(i)(j) / (math.sqrt(second._1(i)(i)) * math.sqrt(second._1(j)(j))))
            .toVector
        )
        .toVector
      cd = math.max(
        relative(scaledFirst, scaledSecond),
        parameters.indices.map(i => math.abs(math.sqrt(first._1(i)(i) / second._1(i)(i)) - 1)).max
      )
      _ <-
        if hd <= options.stabilityTolerance && cd <= options.stabilityTolerance then Right(())
        else
          refuse(
            "curvature_step_unstable",
            s"Step-halving Hessian/covariance differences $hd/$cd exceed ${options.stabilityTolerance}."
          )
    yield CurvatureEvidence(
      fine,
      second._1.map(_.map(_ * 2)),
      steps.map(_ / 2),
      hd,
      cd,
      second._2,
      1 + 4 * parameters.length * parameters.length
    )

  private def hessian(x: Vector[Double], h: Vector[Double], f0: Double)(
      f: Vector[Double] => FitResult[Double]
  ): FitResult[Vector[Vector[Double]]] =
    val n = x.length
    val matrix = Array.ofDim[Double](n, n)
    var error: Option[MixedModelError] = None
    def value(v: FitResult[Double]): Double = v match
      case Right(a) if a.isFinite => a
      case Right(_)               =>
        error = Some(MixedModelError.InferenceUnavailable("curvature_nonfinite", "Nonfinite curvature probe."));
        Double.NaN
      case Left(e) => error = Some(e); Double.NaN
    var i = 0
    while i < n && error.isEmpty do
      matrix(i)(i) =
        (value(f(x.updated(i, x(i) + h(i)))) - 2 * f0 + value(f(x.updated(i, x(i) - h(i))))) / (h(i) * h(i))
      var j = 0
      while j < i && error.isEmpty do
        def probe(a: Double, b: Double) = f(x.updated(i, x(i) + a * h(i)).updated(j, x(j) + b * h(j)))
        val v =
          (value(probe(1, 1)) - value(probe(1, -1)) - value(probe(-1, 1)) + value(probe(-1, -1))) / (4 * h(i) * h(j))
        matrix(i)(j) = v; matrix(j)(i) = v; j += 1
      i += 1
    error.toLeft(matrix.map(_.toVector).toVector)

  /** Scaled Cholesky; no ridge, pseudo-inverse or indefinite fallback. */
  private def inverse(a: Vector[Vector[Double]], pivotTolerance: Double): FitResult[(Vector[Vector[Double]], Double)] =
    val n = a.length
    if a.indices.exists(i => !a(i)(i).isFinite || a(i)(i) <= 0) then
      return refuse("curvature_not_positive_definite", "Observed information has a nonpositive diagonal.")
    val scale = Vector.tabulate(n)(i => math.sqrt(a(i)(i)))
    val l = Array.ofDim[Double](n, n)
    var minPivot = Double.PositiveInfinity
    var i = 0
    while i < n do
      var j = 0
      while j <= i do
        var v = a(i)(j) / (scale(i) * scale(j)); var k = 0
        while k < j do { v -= l(i)(k) * l(j)(k); k += 1 }
        if i == j then
          if !v.isFinite || v <= pivotTolerance then
            return refuse(
              "curvature_not_positive_definite",
              "Observed information is indefinite or numerically ill-conditioned."
            )
          minPivot = math.min(minPivot, v); l(i)(j) = math.sqrt(v)
        else l(i)(j) = v / l(j)(j)
        j += 1
      i += 1
    val inv = Array.ofDim[Double](n, n)
    var col = 0
    while col < n do
      val y = Array.fill(n)(0.0); val x = Array.fill(n)(0.0); i = 0
      while i < n do
        var v = if i == col then 1.0 else 0.0; var k = 0
        while k < i do { v -= l(i)(k) * y(k); k += 1 }
        y(i) = v / l(i)(i); i += 1
      i = n - 1
      while i >= 0 do
        var v = y(i); var k = i + 1
        while k < n do { v -= l(k)(i) * x(k); k += 1 }
        x(i) = v / l(i)(i); i -= 1
      i = 0
      while i < n do { inv(i)(col) = x(i) / (scale(i) * scale(col)); i += 1 }
      col += 1
    if inv.exists(_.exists(x => !x.isFinite)) then
      refuse("curvature_nonfinite", "Nonfinite inverse observed information.")
    else Right((inv.map(_.toVector).toVector, minPivot))

  private def relative(a: Vector[Vector[Double]], b: Vector[Vector[Double]]): Double =
    val numerator = a.indices.map(i => a(i).indices.map(j => math.pow(a(i)(j) - b(i)(j), 2)).sum).sum
    val denominator = b.iterator.flatMap(_.iterator).map(x => x * x).sum
    math.sqrt(numerator / math.max(denominator, 1e-300))
