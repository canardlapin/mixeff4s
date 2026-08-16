package mixeff4s.optimizer

import mixeff4s.error.{FitResult, MixedModelError}

final case class TrustBqOptions(
    initialRadius: Double = 0.75,
    finalRadius: Double = 1e-6,
    maxEvaluations: Int = 1000,
    ftolAbs: Double = 1e-10,
    ftolRel: Double = 1e-10,
    ftolRequiresLocalRadius: Boolean = false,
    etaAccept: Double = 0.05,
    etaExpand: Double = 0.75,
    shrinkFactor: Double = 0.5,
    expandFactor: Double = 1.8,
    maxCrossTerms: Int = Int.MaxValue,
    stallIterations: Int = 4,
    stallFtolRel: Double = -1.0,
    stallFtolAbs: Double = -1.0,
    stallRequiresStableX: Boolean = true,
    reuseSamples: Boolean = false
)

enum TrustBqStopReason:
  case RadiusBelowTolerance, ObjectiveTolerance, MaxEvaluations, StepBelowTolerance,
    ObjectiveStagnation, CertifiedConvergence

  def isAcceptableConvergence: Boolean =
    this != TrustBqStopReason.MaxEvaluations

final case class TrustBqResult(
    x: Vector[Double],
    fmin: Double,
    fevals: Int,
    iterations: Int,
    finalRadius: Double,
    stopReason: TrustBqStopReason
)

object TrustBq:
  private val MinStep = 1e-12

  def minimize(
      initial: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      options: TrustBqOptions
  )(objective: Vector[Double] => FitResult[Double]): FitResult[TrustBqResult] =
    validate(initial, lower, upper, options).flatMap: _ =>
      val cache = scala.collection.mutable.HashMap.empty[Vector[Long], Double]
      var fevals = 0
      def evaluate(x: Vector[Double]): FitResult[Double] =
        if options.reuseSamples then
          val key = x.map(java.lang.Double.doubleToLongBits)
          cache.get(key) match
            case Some(cached) => Right(cached)
            case None =>
              objective(x).flatMap: value =>
                if !value.isFinite then
                  Left(MixedModelError.Optimization("TrustBQ objective returned a non-finite value"))
                else
                  fevals += 1
                  cache.update(key, value)
                  Right(value)
        else
          objective(x).flatMap: value =>
            if !value.isFinite then
              Left(MixedModelError.Optimization("TrustBQ objective returned a non-finite value"))
            else
              fevals += 1
              Right(value)

      var x = project(initial, lower, upper)
      evaluate(x).flatMap: f0 =>
        val initialObjective = f0
        var f = f0
        var bestX = x
        var bestF = f
        var radius = options.initialRadius.max(options.finalRadius)
        var iterations = 0
        var stalled = 0
        var stallBestF = bestF
        var stallBestX = bestX
        val stallFtolRel = if options.stallFtolRel >= 0.0 then options.stallFtolRel else options.ftolRel
        val stallFtolAbs = if options.stallFtolAbs >= 0.0 then options.stallFtolAbs else options.ftolAbs
        val ftolRadius = (options.initialRadius / 16.0).max(options.finalRadius)

        def done(reason: TrustBqStopReason): TrustBqResult =
          TrustBqResult(bestX, bestF, fevals, iterations, radius, reason)

        var result: Option[TrustBqResult] = None
        var error: Option[MixedModelError] = None
        while result.isEmpty && error.isEmpty do
          if fevals >= options.maxEvaluations then result = Some(done(TrustBqStopReason.MaxEvaluations))
          else if radius <= options.finalRadius then result = Some(done(TrustBqStopReason.RadiusBelowTolerance))
          else
            val stallObjTol = stallFtolAbs + stallFtolRel * stallBestF.abs.max(1.0)
            val improvedF = (stallBestF - bestF) > stallObjTol
            val movedX =
              options.stallRequiresStableX &&
                bestX.zip(stallBestX).map((a, b) => (a - b).abs).foldLeft(0.0)(_ max _) > options.finalRadius
            if improvedF || movedX then
              stalled = 0
              stallBestF = bestF
              stallBestX = bestX
            else stalled += 1
            val hasDescended = bestF < initialObjective
            val stallRadiusIsLocal =
              if options.ftolRequiresLocalRadius then radius <= ftolRadius
              else radius < options.initialRadius
            if hasDescended && stalled >= options.stallIterations && stallRadiusIsLocal then
              result = Some(done(TrustBqStopReason.ObjectiveStagnation))
            else
              iterations += 1
              buildQuadratic(x, f, lower, upper, radius, options, () => fevals, evaluate) match
                case Left(err) => error = Some(err)
                case Right(model) =>
                  if fevals >= options.maxEvaluations then
                    result = Some(done(TrustBqStopReason.MaxEvaluations))
                  else
                    val step = trustRegionStep(model, x, lower, upper, radius)
                    val stepNorm = norm(step)
                    if stepNorm <= options.finalRadius then
                      radius *= options.shrinkFactor
                      if radius <= options.finalRadius then
                        result = Some(done(TrustBqStopReason.StepBelowTolerance))
                    else
                      val predicted = model.predictedReduction(step)
                      if !predicted.isFinite || predicted <= 0.0 then radius *= options.shrinkFactor
                      else
                        val trial = project(x.zip(step).map(_ + _), lower, upper)
                        if distance(x, trial) <= options.finalRadius then radius *= options.shrinkFactor
                        else
                          evaluate(trial) match
                            case Left(err) => error = Some(err)
                            case Right(trialF) =>
                              if trialF < bestF then
                                bestF = trialF
                                bestX = trial
                              val actual = f - trialF
                              val ratio = actual / predicted
                              if ratio >= options.etaAccept && actual > 0.0 then
                                val oldF = f
                                x = trial
                                f = trialF
                                if ratio >= options.etaExpand && stepNorm > 0.8 * radius then
                                  radius *= options.expandFactor
                                val objectiveTol = options.ftolAbs + options.ftolRel * oldF.abs.max(1.0)
                                if actual.abs <= objectiveTol &&
                                  (!options.ftolRequiresLocalRadius || radius <= ftolRadius)
                                then result = Some(done(TrustBqStopReason.ObjectiveTolerance))
                              else radius *= options.shrinkFactor

        error.toLeft(result.get)

  private final case class SideSample(delta: Double, value: Double)
  private final case class QuadraticModel(gradient: Array[Double], hessian: Array[Double], n: Int):
    def predictedReduction(step: Vector[Double]): Double =
      var linear = 0.0
      var i = 0
      while i < n do
        linear += gradient(i) * step(i)
        i += 1
      var quadratic = 0.0
      i = 0
      while i < n do
        var acc = 0.0
        var j = 0
        while j < n do
          acc += hessian(i * n + j) * step(j)
          j += 1
        quadratic += step(i) * acc
        i += 1
      -(linear + 0.5 * quadratic)

  private def buildQuadratic(
      x: Vector[Double],
      f0: Double,
      lower: Vector[Double],
      upper: Vector[Double],
      radius: Double,
      options: TrustBqOptions,
      fevals: () => Int,
      evaluate: Vector[Double] => FitResult[Double]
  ): FitResult[QuadraticModel] =
    val n = x.length
    val gradient = new Array[Double](n)
    val hessian = new Array[Double](n * n)
    val side = Array.fill[Option[SideSample]](n)(None)
    val minStep = MinStep * radius.max(1.0)

    def evalAxis(index: Int, delta: Double): FitResult[Double] =
      evaluate(x.updated(index, x(index) + delta))

    var i = 0
    var error: Option[MixedModelError] = None
    while i < n && error.isEmpty && fevals() < options.maxEvaluations do
      val hPlus = feasibleAxisDelta(x(i), lower(i), upper(i), radius, 1.0)
      val hMinus = feasibleAxisDelta(x(i), lower(i), upper(i), radius, -1.0)
      val plus =
        if hPlus.abs > minStep && fevals() < options.maxEvaluations then
          evalAxis(i, hPlus) match
            case Left(err)    => error = Some(err); None
            case Right(value) => Some((hPlus, value))
        else None
      val minus =
        if error.isEmpty && hMinus.abs > minStep && fevals() < options.maxEvaluations then
          evalAxis(i, hMinus) match
            case Left(err)    => error = Some(err); None
            case Right(value) => Some((hMinus, value))
        else None
      (plus, minus) match
        case (Some((hp, fp)), Some((hm, fm))) =>
          val a = hp.abs
          val b = hm.abs
          val denom = a * b * (a + b)
          gradient(i) = (b * b * (fp - f0) - a * a * (fm - f0)) / denom
          hessian(i * n + i) = 2.0 * (b * (fp - f0) + a * (fm - f0)) / denom
          side(i) = Some(selectCrossSide(hp, fp, hm, fm, f0))
        case (Some((delta, value)), None) =>
          applyOneSided(i, delta, value)
        case (None, Some((delta, value))) =>
          applyOneSided(i, delta, value)
        case _ => ()

      def applyOneSided(i: Int, delta: Double, value: Double): Unit =
        var slope = (value - f0) / delta
        var curvature = 0.0
        val secondDelta = 2.0 * delta
        if axisDeltaIsFeasible(x(i), lower(i), upper(i), secondDelta, radius * 2.0) &&
          fevals() < options.maxEvaluations
        then
          evalAxis(i, secondDelta) match
            case Left(err) => error = Some(err)
            case Right(second) =>
              curvature = (second - 2.0 * value + f0) / (delta * delta)
              slope = (value - f0 - 0.5 * curvature * delta * delta) / delta
        gradient(i) = slope
        hessian(i * n + i) = curvature
        side(i) = Some(SideSample(delta, value))
      i += 1

    error match
      case Some(err) => Left(err)
      case None =>
        var crossTerms = 0
        i = 0
        while i < n && error.isEmpty && fevals() < options.maxEvaluations && crossTerms < options.maxCrossTerms do
          var j = i + 1
          while j < n && error.isEmpty && fevals() < options.maxEvaluations && crossTerms < options.maxCrossTerms do
            (side(i), side(j)) match
              case (Some(left), Some(right)) =>
                val trial = project(x.updated(i, x(i) + left.delta).updated(j, x(j) + right.delta), lower, upper)
                if distance(x, trial) > minStep then
                  evaluate(trial) match
                    case Left(err) => error = Some(err)
                    case Right(fij) =>
                      val cross = (fij - left.value - right.value + f0) / (left.delta * right.delta)
                      if cross.isFinite then
                        hessian(i * n + j) = cross
                        hessian(j * n + i) = cross
                        crossTerms += 1
              case _ => ()
            j += 1
          i += 1
        error.toLeft(QuadraticModel(gradient, hessian, n))

  private def selectCrossSide(hp: Double, fp: Double, hm: Double, fm: Double, f0: Double): SideSample =
    if (f0 - fp) >= (f0 - fm) then SideSample(hp, fp) else SideSample(hm, fm)

  private def trustRegionStep(
      model: QuadraticModel,
      x: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      radius: Double
  ): Vector[Double] =
    val n = x.length
    if norm(model.gradient.toVector) <= MinStep then return Vector.fill(n)(0.0)
    var shift = 0.0
    var attempt = 0
    while attempt < 14 do
      tryCholSolve(model.hessian, n, model.gradient.map(v => -v), shift) match
        case Some(raw) =>
          val step = boundAndScale(raw.toVector, x, lower, upper, radius)
          if norm(step) > MinStep && model.predictedReduction(step) > 0.0 then return step
        case None => ()
      shift = if shift == 0.0 then 1e-8 else shift * 10.0
      attempt += 1
    cauchyStep(model, x, lower, upper, radius)

  private def cauchyStep(
      model: QuadraticModel,
      x: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      radius: Double
  ): Vector[Double] =
    val gnorm = norm(model.gradient.toVector)
    if gnorm <= MinStep then return Vector.fill(x.length)(0.0)
    val direction = model.gradient.iterator.map(v => -v / gnorm).toVector
    var hdDot = 0.0
    var i = 0
    while i < model.n do
      var acc = 0.0
      var j = 0
      while j < model.n do
        acc += model.hessian(i * model.n + j) * direction(j)
        j += 1
      hdDot += direction(i) * acc
      i += 1
    val alphaModel = if hdDot > 0.0 then (gnorm / hdDot).min(radius) else radius
    val alpha = alphaModel.min(maxBoundStep(x, direction, lower, upper, radius)).max(0.0)
    direction.map(_ * alpha)

  private def tryCholSolve(h: Array[Double], n: Int, rhs: Array[Double], shift: Double): Option[Array[Double]] =
    val l = h.clone()
    if shift != 0.0 then
      var i = 0
      while i < n do
        l(i * n + i) += shift
        i += 1
    var j = 0
    while j < n do
      var s = l(j * n + j)
      var k = 0
      while k < j do
        s -= l(j * n + k) * l(j * n + k)
        k += 1
      if s <= 0.0 || !s.isFinite then return None
      val ljj = math.sqrt(s)
      l(j * n + j) = ljj
      var i = j + 1
      while i < n do
        var t = l(i * n + j)
        k = 0
        while k < j do
          t -= l(i * n + k) * l(j * n + k)
          k += 1
        l(i * n + j) = t / ljj
        i += 1
      j += 1
    val y = new Array[Double](n)
    var i = 0
    while i < n do
      var s = rhs(i)
      var k = 0
      while k < i do
        s -= l(i * n + k) * y(k)
        k += 1
      y(i) = s / l(i * n + i)
      i += 1
    val x = new Array[Double](n)
    i = n - 1
    while i >= 0 do
      var s = y(i)
      var k = i + 1
      while k < n do
        s -= l(k * n + i) * x(k)
        k += 1
      x(i) = s / l(i * n + i)
      i -= 1
    Some(x)

  private def feasibleAxisDelta(value: Double, lower: Double, upper: Double, radius: Double, sign: Double): Double =
    var delta = math.signum(sign) * radius
    if delta > 0.0 && upper.isFinite then delta = delta.min(upper - value)
    else if delta < 0.0 && lower.isFinite then delta = delta.max(lower - value)
    delta

  private def axisDeltaIsFeasible(
      value: Double,
      lower: Double,
      upper: Double,
      delta: Double,
      maxAbsDelta: Double
  ): Boolean =
    delta.abs <= maxAbsDelta && value + delta >= lower && value + delta <= upper && delta.abs > MinStep

  private def project(x: Vector[Double], lower: Vector[Double], upper: Vector[Double]): Vector[Double] =
    x.indices
      .map: i =>
        var v = x(i)
        if lower(i).isFinite && v < lower(i) then v = lower(i)
        if upper(i).isFinite && v > upper(i) then v = upper(i)
        v
      .toVector

  private def boundAndScale(
      step: Vector[Double],
      x: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      radius: Double
  ): Vector[Double] =
    val clipped = step.indices
      .map: i =>
        var s = step(i)
        val candidate = x(i) + s
        if lower(i).isFinite && candidate < lower(i) then s = lower(i) - x(i)
        if upper(i).isFinite && candidate > upper(i) then s = upper(i) - x(i)
        s
      .toVector
    val stepNorm = norm(clipped)
    if stepNorm > radius && stepNorm > 0.0 then clipped.map(_ * (radius / stepNorm))
    else clipped

  private def maxBoundStep(
      x: Vector[Double],
      direction: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      radius: Double
  ): Double =
    var alpha = radius
    var i = 0
    while i < x.length do
      if direction(i) > 0.0 && upper(i).isFinite then alpha = alpha.min((upper(i) - x(i)) / direction(i))
      else if direction(i) < 0.0 && lower(i).isFinite then alpha = alpha.min((lower(i) - x(i)) / direction(i))
      i += 1
    alpha.max(0.0)

  private def norm(values: Vector[Double]): Double =
    math.sqrt(values.iterator.map(v => v * v).sum)

  private def distance(left: Vector[Double], right: Vector[Double]): Double =
    norm(left.zip(right).map((a, b) => a - b))

  private def validate(
      initial: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      options: TrustBqOptions
  ): FitResult[Unit] =
    val n = initial.length
    if n == 0 then Left(MixedModelError.Optimization("TrustBQ requires at least one parameter"))
    else if lower.length != n || upper.length != n then
      Left(MixedModelError.DimensionMismatch("TrustBQ bounds length does not match parameter length"))
    else if !options.initialRadius.isFinite || options.initialRadius <= 0.0 ||
      !options.finalRadius.isFinite || options.finalRadius <= 0.0 ||
      options.finalRadius > options.initialRadius
    then Left(MixedModelError.Optimization("TrustBQ requires 0 < final_radius <= initial_radius"))
    else if options.maxEvaluations <= 0 then
      Left(MixedModelError.Optimization("TrustBQ max_evaluations must be positive"))
    else if options.stallIterations <= 0 then
      Left(MixedModelError.Optimization("TrustBQ stall_iterations must be positive"))
    else if !(options.etaAccept >= 0.0 && options.etaAccept < 1.0) ||
      !(options.etaExpand >= options.etaAccept && options.etaExpand <= 1.0) ||
      !(options.shrinkFactor >= 0.0 && options.shrinkFactor < 1.0) ||
      options.expandFactor <= 1.0 || !options.expandFactor.isFinite
    then Left(MixedModelError.Optimization("TrustBQ trust-region constants are invalid"))
    else if initial.exists(v => !v.isFinite) then
      Left(MixedModelError.Optimization("TrustBQ initial point must be finite"))
    else if lower.indices.exists(i => lower(i) > upper(i)) then
      Left(MixedModelError.Optimization("TrustBQ lower bound exceeds upper bound"))
    else Right(())
