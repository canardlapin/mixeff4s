package mixeff4s.optimizer

import mixeff4s.error.{FitResult, MixedModelError}

/** Value and exact first derivative supplied together by a gradient oracle. */
final case class ObjectiveGradient(value: Double, gradient: Vector[Double])

/** Gradient-oracle TrustBQ kernel. It deliberately owns no fitted-model state: callers must still reevaluate their
  * model at the returned point.
  */
private[optimizer] object GradientTrustBq:
  private val MinStep = 1e-12
  private val SeedRelativeStep = 1e-4
  private val Sr1Skip = 1e-8

  def minimize(
      initial: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      options: TrustBqOptions
  )(oracle: Vector[Double] => FitResult[ObjectiveGradient]): FitResult[TrustBqResult] =
    validate(initial, lower, upper, options) match
      case Left(error) => Left(error)
      case Right(_)    => minimizeValidated(initial, lower, upper, options)(oracle)

  private def minimizeValidated(
      initial: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      options: TrustBqOptions
  )(oracle: Vector[Double] => FitResult[ObjectiveGradient]): FitResult[TrustBqResult] =
    val n = initial.length
    var evaluations = 0
    def evaluate(point: Vector[Double]): FitResult[ObjectiveGradient] =
      if evaluations >= options.maxEvaluations then Left(MixedModelError.OptimizerBudgetExhausted(evaluations))
      else
        oracle(point).flatMap: answer =>
          evaluations += 1
          if !answer.value.isFinite then
            Left(MixedModelError.Optimization("TrustBQ objective returned a non-finite value"))
          else if answer.gradient.length != n then
            Left(
              MixedModelError.Optimization(
                s"TrustBQ gradient oracle returned ${answer.gradient.length} components for $n parameters"
              )
            )
          else if answer.gradient.exists(v => !v.isFinite) then
            Left(MixedModelError.Optimization("TrustBQ gradient oracle returned a non-finite component"))
          else Right(answer)

    var x = project(initial, lower, upper)
    evaluate(x) match
      case Left(error)          => Left(error)
      case Right(initialAnswer) => {
        var f = initialAnswer.value
        var g = initialAnswer.gradient.toArray
        var bestX = x
        var bestF = f
        var bestGradient = g.clone()
        val initialObjective = f
        var radius = options.initialRadius.max(options.finalRadius)
        var iterations = 0
        var stalled = 0
        var stallBestF = bestF
        var stallBestX = bestX
        val stallFtolRel = if options.stallFtolRel >= 0.0 then options.stallFtolRel else options.ftolRel
        val stallFtolAbs = if options.stallFtolAbs >= 0.0 then options.stallFtolAbs else options.ftolAbs
        val ftolRadius = (options.initialRadius / 16.0).max(options.finalRadius)
        var consecutiveRejections = 0
        var consecutivePoorAccepts = 0
        var consecutiveAccepted = 0
        var radiusAtSeed = radius
        var h = Array.ofDim[Double](n * n)

        // Gradient probes are evidence too: they consume the same budget and
        // may expose a downhill escape from a zero-gradient boundary point.
        var seedBest: Option[(Vector[Double], ObjectiveGradient)] = None
        var axis = 0
        while axis < n && evaluations < options.maxEvaluations do
          val delta0 = SeedRelativeStep * x(axis).abs.max(1.0)
          val forward = feasibleDelta(x(axis), lower(axis), upper(axis), delta0)
          val delta =
            if forward.abs > MinStep then forward else feasibleDelta(x(axis), lower(axis), upper(axis), -delta0)
          if delta.abs > MinStep then
            val probe = x.updated(axis, x(axis) + delta)
            evaluate(probe) match
              case Left(err)     => return Left(err)
              case Right(answer) =>
                if answer.value < bestF then
                  bestF = answer.value; bestX = probe; bestGradient = answer.gradient.toArray;
                  seedBest = Some((probe, answer))
                var row = 0
                while row < n do
                  h(row * n + axis) = (answer.gradient(row) - g(row)) / delta
                  row += 1
          axis += 1
        symmetrize(h, n)
        // A seed that is already better must become the active centre, not
        // merely an output candidate. Its secant model is still valid enough
        // for a trust-region trial and subsequent secants correct it.
        seedBest.foreach: (point, answer) =>
          x = point; f = answer.value; g = answer.gradient.toArray

        // Promotion changes the model centre, so the original seed Hessian is
        // no longer a local derivative. Rebuild it before taking the first
        // trust step; otherwise a high-dimensional start can crawl on a
        // surrogate centred at the discarded point.
        if seedBest.nonEmpty && evaluations + n < options.maxEvaluations then
          java.util.Arrays.fill(h, 0.0)
          axis = 0
          while axis < n && evaluations < options.maxEvaluations do
            val delta0 = SeedRelativeStep * x(axis).abs.max(1.0)
            val forward = feasibleDelta(x(axis), lower(axis), upper(axis), delta0)
            val delta =
              if forward.abs > MinStep then forward else feasibleDelta(x(axis), lower(axis), upper(axis), -delta0)
            if delta.abs > MinStep then
              val point = x.updated(axis, x(axis) + delta)
              evaluate(point) match
                case Left(error)   => return Left(error)
                case Right(answer) =>
                  if answer.value < bestF then
                    bestF = answer.value; bestX = point; bestGradient = answer.gradient.toArray
                  var row = 0
                  while row < n do
                    h(row * n + axis) = (answer.gradient(row) - g(row)) / delta
                    row += 1
            axis += 1
          symmetrize(h, n)

        def reseedHessian(): FitResult[Unit] =
          val fresh = Array.ofDim[Double](n * n)
          var better: Option[(Vector[Double], ObjectiveGradient)] = None
          var coordinate = 0
          while coordinate < n && evaluations < options.maxEvaluations do
            val step = SeedRelativeStep * x(coordinate).abs.max(1.0)
            val forward = feasibleDelta(x(coordinate), lower(coordinate), upper(coordinate), step)
            val delta =
              if forward.abs > MinStep then forward
              else feasibleDelta(x(coordinate), lower(coordinate), upper(coordinate), -step)
            if delta.abs > MinStep then
              val point = x.updated(coordinate, x(coordinate) + delta)
              evaluate(point) match
                case Left(error)   => return Left(error)
                case Right(answer) =>
                  if answer.value < bestF then
                    bestF = answer.value; bestX = point; bestGradient = answer.gradient.toArray;
                    better = Some((point, answer))
                  var row = 0
                  while row < n do
                    fresh(row * n + coordinate) = (answer.gradient(row) - g(row)) / delta
                    row += 1
            coordinate += 1
          symmetrize(fresh, n)
          h = fresh
          better.foreach { case (point, answer) =>
            x = point
            f = answer.value
            g = answer.gradient.toArray
          }
          Right(())

        def result(reason: TrustBqStopReason): TrustBqResult =
          TrustBqResult(
            bestX,
            bestF,
            evaluations,
            iterations,
            radius,
            reason,
            evaluations,
            Some(projectedNorm(bestX, bestGradient, lower, upper))
          )

        /** A zero projected gradient is not a minimum certificate: at a bound, a variance coordinate can have zero
          * first derivative at a maximum; coordinate probes alone miss mixed negative-curvature directions. Refresh a
          * gradient-difference Hessian at this centre, then test only feasible negative-eigenvector directions.
          */
        def findFeasibleEscape(): FitResult[Option[(Vector[Double], ObjectiveGradient)]] =
          val freshH = Array.ofDim[Double](n * n)
          val deltas = Array.ofDim[Double](n)
          var i = 0
          while i < n do
            val step = SeedRelativeStep * x(i).abs.max(1.0)
            val forward = feasibleDelta(x(i), lower(i), upper(i), step)
            deltas(i) = if forward.abs > MinStep then forward else feasibleDelta(x(i), lower(i), upper(i), -step)
            i += 1
          var escape: Option[(Vector[Double], ObjectiveGradient)] = None
          def probe(point: Vector[Double]): FitResult[Unit] =
            if evaluations >= options.maxEvaluations then Left(MixedModelError.OptimizerBudgetExhausted(evaluations))
            else
              evaluate(point).map: answer =>
                if answer.value < bestF then
                  bestF = answer.value; bestX = point; bestGradient = answer.gradient.toArray
                if answer.value < f && escape.forall(_._2.value >= answer.value) then escape = Some((point, answer))
          i = 0
          while i < n do
            if deltas(i).abs > MinStep then
              val point = x.updated(i, x(i) + deltas(i))
              evaluate(point) match
                case Left(error)   => return Left(error)
                case Right(answer) =>
                  if answer.value < bestF then
                    bestF = answer.value; bestX = point; bestGradient = answer.gradient.toArray
                  if answer.value < f && escape.forall(_._2.value >= answer.value) then escape = Some((point, answer))
                  var row = 0
                  while row < n do
                    freshH(row * n + i) = (answer.gradient(row) - g(row)) / deltas(i)
                    row += 1
            i += 1
          if escape.nonEmpty then Right(escape)
          else
            symmetrize(freshH, n)
            val (eigenvalues, eigenvectors) = jacobiEigen(freshH, n)
            i = 0
            while i < n && escape.isEmpty do
              if eigenvalues(i) < -MinStep then
                var sign = -1.0
                while sign <= 1.0 && escape.isEmpty do
                  val direction = Array.tabulate(n)(j => sign * eigenvectors(j * n + i))
                  var j = 0
                  while j < n do
                    if (x(j) <= lower(j) + MinStep && direction(j) < 0.0) || (x(j) >= upper(j) - MinStep && direction(
                        j
                      ) > 0.0)
                    then direction(j) = 0.0
                    j += 1
                  val length = norm(direction)
                  if length > MinStep then
                    val normalized = direction.map(_ / length).toVector
                    val alpha = (SeedRelativeStep).min(maxBound(x, normalized, lower, upper, SeedRelativeStep))
                    if alpha > MinStep then
                      probe(x.indices.map(j => x(j) + alpha * normalized(j)).toVector) match
                        case Left(error) => return Left(error)
                        case Right(_)    => ()
                  sign += 2.0
              i += 1
          Right(escape)

        while true do
          if evaluations >= options.maxEvaluations then return Right(result(TrustBqStopReason.MaxEvaluations))
          if radius <= options.finalRadius then return Right(result(TrustBqStopReason.RadiusBelowTolerance))
          val projected = projectedNorm(x, g, lower, upper)
          val band = options.gradientTolerance * f.abs.max(1.0)
          val gradientAllowsObjectiveStop = options.gradientTolerance <= 0.0 || projected <= 1e3 * band
          if options.gradientTolerance > 0.0 && x == bestX && projected <= band then
            findFeasibleEscape() match
              case Left(_: MixedModelError.OptimizerBudgetExhausted) =>
                return Right(result(TrustBqStopReason.MaxEvaluations))
              case Left(error)                          => return Left(error)
              case Right(Some((escapeX, escapeAnswer))) =>
                x = escapeX; f = escapeAnswer.value; g = escapeAnswer.gradient.toArray
                if f < bestF then bestF = f; bestX = x; bestGradient = g.clone()
              case Right(None) => return Right(result(TrustBqStopReason.GradientBelowTolerance))

          val stallObjectiveTolerance = stallFtolAbs + stallFtolRel * stallBestF.abs.max(1.0)
          val improved = (stallBestF - bestF) > stallObjectiveTolerance
          val moved = options.stallRequiresStableX && distance(bestX, stallBestX) > options.finalRadius
          if improved || moved then
            stalled = 0; stallBestF = bestF; stallBestX = bestX
          else stalled += 1
          val localRadius =
            if options.ftolRequiresLocalRadius then radius <= ftolRadius else radius < options.initialRadius
          if bestF < initialObjective && stalled >= options.stallIterations && localRadius && gradientAllowsObjectiveStop
          then return Right(result(TrustBqStopReason.ObjectiveStagnation))

          iterations += 1
          val step = trustStep(g, h, x, lower, upper, radius)
          if norm(step) <= options.finalRadius then
            radius *= options.shrinkFactor
            if radius <= options.finalRadius then return Right(result(TrustBqStopReason.StepBelowTolerance))
          else
            val predicted = reduction(g, h, step)
            if !predicted.isFinite || predicted <= 0.0 then radius *= options.shrinkFactor
            else
              val trial = project(x.indices.map(i => x(i) + step(i)).toVector, lower, upper)
              if distance(x, trial) <= options.finalRadius then radius *= options.shrinkFactor
              else
                evaluate(trial) match
                  case Left(_: MixedModelError.OptimizerBudgetExhausted) =>
                    return Right(result(TrustBqStopReason.MaxEvaluations))
                  case Left(err)     => return Left(err)
                  case Right(answer) =>
                    if answer.value < bestF then
                      bestF = answer.value
                      bestX = trial
                      bestGradient = answer.gradient.toArray
                    val s = trial.indices.map(i => trial(i) - x(i)).toArray
                    val y = answer.gradient.indices.map(i => answer.gradient(i) - g(i)).toArray
                    secantUpdate(h, s, y, n)
                    val actual = f - answer.value
                    val ratio = actual / predicted
                    if ratio >= options.etaAccept && actual > 0.0 then
                      val oldF = f
                      x = trial; f = answer.value; g = answer.gradient.toArray
                      consecutiveRejections = 0
                      consecutiveAccepted += 1
                      if ratio >= options.etaExpand && norm(step) > 0.8 * radius then radius *= options.expandFactor
                      if ratio < 0.25 then consecutivePoorAccepts += 1
                      else consecutivePoorAccepts = 0
                      val contracted = radius * 16.0 <= radiusAtSeed
                      if (consecutivePoorAccepts >= 6 || consecutiveAccepted >= 6 || contracted) && evaluations + n < options.maxEvaluations
                      then
                        reseedHessian() match
                          case Left(error) => return Left(error)
                          case Right(_)    =>
                            consecutivePoorAccepts = 0
                            consecutiveAccepted = 0
                            radiusAtSeed = radius
                      val objectiveTolerance = options.ftolAbs + options.ftolRel * oldF.abs.max(1.0)
                      val newProjected = projectedNorm(x, g, lower, upper)
                      val newGradientAllowsStop =
                        options.gradientTolerance <= 0.0 || newProjected <= 1e3 * options.gradientTolerance * f.abs.max(
                          1.0
                        )
                      if actual.abs <= objectiveTolerance && (!options.ftolRequiresLocalRadius || radius <= ftolRadius) && newGradientAllowsStop
                      then return Right(result(TrustBqStopReason.ObjectiveTolerance))
                    else
                      radius *= options.shrinkFactor
                      consecutiveRejections += 1
                      consecutiveAccepted = 0
                      consecutivePoorAccepts = 0
                      val contracted = radius * 16.0 <= radiusAtSeed
                      if (consecutiveRejections >= 4 || contracted) && evaluations + n < options.maxEvaluations then
                        reseedHessian() match
                          case Left(error) => return Left(error)
                          case Right(_)    =>
                            consecutiveRejections = 0
                            radiusAtSeed = radius
        Right(result(TrustBqStopReason.MaxEvaluations)) // unreachable; satisfies the expression type
      }

  private def trustStep(
      g: Array[Double],
      h: Array[Double],
      x: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      radius: Double
  ): Vector[Double] =
    val n = g.length
    if norm(g) <= MinStep then return Vector.fill(n)(0.0)
    var shift = 0.0
    var attempt = 0
    while attempt < 14 do
      cholSolve(h, n, g.map(-_), shift) match
        case Some(raw) =>
          val step = boundScale(raw.toVector, x, lower, upper, radius)
          if norm(step) > MinStep && reduction(g, h, step) > 0.0 then
            if shift > 0.0 then
              val length = norm(step)
              val direction = step.map(_ / length)
              val curvature = quadratic(h, direction)
              val limit = maxBound(x, direction, lower, upper, radius)
              val alpha =
                if curvature > 0.0 then (-dot(g, direction.toArray) / curvature).min(limit).max(0.0) else limit
              val ray = direction.map(_ * alpha)
              if reduction(g, h, ray) > reduction(g, h, step) then return ray
            return step
        case None => ()
      shift = if shift == 0.0 then 1e-8 else shift * 10.0
      attempt += 1
    val ng = norm(g)
    val rawDirection = g.map(-_ / ng)
    var i = 0
    while i < n do
      if (x(i) <= lower(i) + MinStep && rawDirection(i) < 0.0) ||
        (x(i) >= upper(i) - MinStep && rawDirection(i) > 0.0)
      then rawDirection(i) = 0.0
      i += 1
    val tangentNorm = norm(rawDirection)
    if tangentNorm <= MinStep then return Vector.fill(n)(0.0)
    val d = rawDirection.map(_ / tangentNorm).toVector
    val curvature = quadratic(h, d)
    val alpha =
      (if curvature > 0.0 then (-dot(g, d.toArray) / curvature).min(radius) else radius)
        .min(maxBound(x, d, lower, upper, radius))
    d.map(_ * alpha.max(0.0))

  private def secantUpdate(h: Array[Double], s: Array[Double], y: Array[Double], n: Int): Unit =
    val ys = dot(y, s)
    val hs = multiply(h, s, n)
    val shs = dot(s, hs)
    if ys > Sr1Skip * norm(y) * norm(s) && ys > 0.0 then
      // Retain measured positive curvature even when the current model is indefinite.
      // Only the old-curvature subtraction requires a positive denominator.
      if shs > 0.0 then rankUpdate(h, hs, -1.0 / shs, n)
      rankUpdate(h, y, 1.0 / ys, n)
    else
      val r = y.indices.map(i => y(i) - hs(i)).toArray
      val denom = dot(r, s)
      if denom.abs >= Sr1Skip * norm(r) * norm(s) && denom != 0.0 then rankUpdate(h, r, 1.0 / denom, n)

  private def rankUpdate(h: Array[Double], v: Array[Double], scale: Double, n: Int): Unit =
    var i = 0; while i < n do { var j = 0; while j < n do { h(i * n + j) += scale * v(i) * v(j); j += 1 }; i += 1 }
  private def multiply(h: Array[Double], v: Array[Double], n: Int): Array[Double] =
    Array.tabulate(n)(i => (0 until n).iterator.map(j => h(i * n + j) * v(j)).sum)
  private def reduction(g: Array[Double], h: Array[Double], s: Vector[Double]): Double =
    -dot(g, s.toArray) - 0.5 * quadratic(h, s)
  private def quadratic(h: Array[Double], v: Vector[Double]): Double = dot(v.toArray, multiply(h, v.toArray, v.length))
  private def dot(a: Array[Double], b: Array[Double]): Double = a.indices.iterator.map(i => a(i) * b(i)).sum
  private def norm(v: Array[Double]): Double = math.sqrt(v.iterator.map(x => x * x).sum)
  private def norm(v: Vector[Double]): Double = math.sqrt(v.iterator.map(x => x * x).sum)
  private def distance(a: Vector[Double], b: Vector[Double]): Double = norm(a.indices.map(i => a(i) - b(i)).toVector)
  private def projectedNorm(x: Vector[Double], g: Array[Double], lower: Vector[Double], upper: Vector[Double]): Double =
    x.indices.iterator.map(i => (x(i) - projectOne(x(i) - g(i), lower(i), upper(i))).abs).max
  private def project(x: Vector[Double], lower: Vector[Double], upper: Vector[Double]): Vector[Double] =
    x.indices.map(i => projectOne(x(i), lower(i), upper(i))).toVector
  private def projectOne(v: Double, low: Double, high: Double): Double = v.max(low).min(high)
  private def feasibleDelta(x: Double, low: Double, high: Double, delta: Double): Double =
    if delta > 0.0 then delta.min(high - x) else delta.max(low - x)
  private def maxBound(
      x: Vector[Double],
      d: Vector[Double],
      low: Vector[Double],
      high: Vector[Double],
      radius: Double
  ): Double =
    x.indices.foldLeft(radius) { (a, i) =>
      if d(i) > 0.0 then a.min((high(i) - x(i)) / d(i))
      else if d(i) < 0.0 then a.min((low(i) - x(i)) / d(i))
      else a
    }
  private def boundScale(
      s: Vector[Double],
      x: Vector[Double],
      low: Vector[Double],
      high: Vector[Double],
      radius: Double
  ): Vector[Double] =
    val clipped = s.indices.map(i => (x(i) + s(i)).max(low(i)).min(high(i)) - x(i)).toVector
    val length = norm(clipped); if length > radius && length > 0.0 then clipped.map(_ * radius / length) else clipped
  private def symmetrize(h: Array[Double], n: Int): Unit = for i <- 0 until n; j <- 0 until i do {
    val v = 0.5 * (h(i * n + j) + h(j * n + i)); h(i * n + j) = v; h(j * n + i) = v
  }

  /** Symmetric Jacobi eigensystem; columns of the second result are eigenvectors. */
  private def jacobiEigen(input: Array[Double], n: Int): (Array[Double], Array[Double]) =
    val a = input.clone()
    val vectors = Array.tabulate(n * n)(k => if k / n == k % n then 1.0 else 0.0)
    var sweep = 0
    while sweep < 32 * n * n do
      var p = 0; var q = 1; var largest = 0.0
      var i = 0
      while i < n do
        var j = i + 1
        while j < n do
          val value = a(i * n + j).abs
          if value > largest then { largest = value; p = i; q = j }
          j += 1
        i += 1
      if largest <= 1e-12 then sweep = 32 * n * n
      else
        val app = a(p * n + p); val aqq = a(q * n + q); val apq = a(p * n + q)
        val angle = 0.5 * math.atan2(2.0 * apq, aqq - app)
        val c = math.cos(angle); val s = math.sin(angle)
        i = 0
        while i < n do
          if i != p && i != q then
            val aip = a(i * n + p); val aiq = a(i * n + q)
            a(i * n + p) = c * aip - s * aiq; a(p * n + i) = a(i * n + p)
            a(i * n + q) = s * aip + c * aiq; a(q * n + i) = a(i * n + q)
          i += 1
        a(p * n + p) = c * c * app - 2.0 * s * c * apq + s * s * aqq
        a(q * n + q) = s * s * app + 2.0 * s * c * apq + c * c * aqq
        a(p * n + q) = 0.0; a(q * n + p) = 0.0
        i = 0
        while i < n do
          val vip = vectors(i * n + p); val viq = vectors(i * n + q)
          vectors(i * n + p) = c * vip - s * viq
          vectors(i * n + q) = s * vip + c * viq
          i += 1
        sweep += 1
    (Array.tabulate(n)(i => a(i * n + i)), vectors)

  private def cholSolve(h: Array[Double], n: Int, rhs: Array[Double], shift: Double): Option[Array[Double]] =
    val l = h.clone(); var i = 0; while i < n do { l(i * n + i) += shift; i += 1 }
    var j = 0
    while j < n do
      var sum = l(j * n + j); var k = 0; while k < j do { sum -= l(j * n + k) * l(j * n + k); k += 1 }
      if sum <= 0.0 || !sum.isFinite then return None
      l(j * n + j) = math.sqrt(sum); i = j + 1
      while i < n do {
        var t = l(i * n + j); k = 0; while k < j do { t -= l(i * n + k) * l(j * n + k); k += 1 };
        l(i * n + j) = t / l(j * n + j); i += 1
      }
      j += 1
    val y = new Array[Double](n); i = 0;
    while i < n do {
      var t = rhs(i); var k = 0; while k < i do { t -= l(i * n + k) * y(k); k += 1 }; y(i) = t / l(i * n + i); i += 1
    }
    val out = new Array[Double](n); i = n - 1;
    while i >= 0 do {
      var t = y(i); var k = i + 1; while k < n do { t -= l(k * n + i) * out(k); k += 1 }; out(i) = t / l(i * n + i);
      i -= 1
    }; Some(out)

  private def validate(
      initial: Vector[Double],
      lower: Vector[Double],
      upper: Vector[Double],
      o: TrustBqOptions
  ): FitResult[Unit] =
    if initial.isEmpty then Left(MixedModelError.Optimization("TrustBQ requires at least one parameter"))
    else if lower.length != initial.length || upper.length != initial.length then
      Left(MixedModelError.DimensionMismatch("TrustBQ bounds length does not match parameter length"))
    else if initial.exists(v => !v.isFinite) || lower.exists(_.isNaN) || upper.exists(_.isNaN) || lower.indices.exists(
        i => lower(i) > upper(i)
      )
    then Left(MixedModelError.Optimization("TrustBQ requires finite initial values and ordered bounds"))
    else if !o.initialRadius.isFinite || o.initialRadius <= 0.0 || !o.finalRadius.isFinite || o.finalRadius <= 0.0 || o.finalRadius > o.initialRadius
    then Left(MixedModelError.Optimization("TrustBQ requires 0 < final_radius <= initial_radius"))
    else if o.maxEvaluations <= 0 || o.stallIterations <= 0 || !o.gradientTolerance.isFinite || o.gradientTolerance < 0.0
    then Left(MixedModelError.Optimization("TrustBQ options are invalid"))
    else if !(o.etaAccept >= 0.0 && o.etaAccept < 1.0) || !(o.etaExpand >= o.etaAccept && o.etaExpand <= 1.0) || !(o.shrinkFactor >= 0.0 && o.shrinkFactor < 1.0) || o.expandFactor <= 1.0 || !o.expandFactor.isFinite
    then Left(MixedModelError.Optimization("TrustBQ trust-region constants are invalid"))
    else Right(())
