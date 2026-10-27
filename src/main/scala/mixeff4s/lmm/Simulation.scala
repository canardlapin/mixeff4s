package mixeff4s.lmm

import mixeff4s.error.{FitResult, MixedModelError}

/** Whether a simulated response draws new random effects or conditions on the fitted modes. */
enum SimulationMode:
  case DrawRandomEffects, ConditionOnModes

/** Seeded Gaussian response simulation from an immutable fitted LMM snapshot. */
object Simulation:
  /** Generate one response on the fitted, retained design. Random-effect draws are `sigma * Lambda * z`; conditional
    * simulations use the captured fitted modes directly and add only residual noise.
    */
  def response(
      fit: LmmFit,
      seed: Long,
      mode: SimulationMode = SimulationMode.DrawRandomEffects
  ): FitResult[Vector[Double]] =
    for
      context <- fit.validatedContext
      _ <-
        if fit.sigma.isFinite && fit.sigma >= 0.0 && fit.beta.forall(_.isFinite) && fit.theta.forall(_.isFinite)
        then Right(())
        else
          Left(
            MixedModelError.InferenceUnavailable("simulation_unavailable", "fit has non-finite simulation parameters")
          )
      simulated =
        val design = context.design
        val random = mode match
          case SimulationMode.DrawRandomEffects => drawModes(fit, design, seed)
          case SimulationMode.ConditionOnModes  => context.modes
        Vector.tabulate(design.n): observation =>
          var value = 0.0
          var column = 0
          while column < fit.beta.length do
            value += design.xy.xy(observation, column) * fit.beta(column)
            column += 1
          var term = 0
          while term < design.reterms.length do
            val re = design.reterms(term)
            val level = re.refs(observation)
            var row = 0
            while row < re.vsize do
              value += re.z(row, observation) * random(term)(level * re.vsize + row)
              row += 1
            term += 1
          value + fit.sigma * normal(seed, design.nRanef.toLong + observation.toLong)
      _ <-
        if simulated.forall(_.isFinite) then Right(())
        else
          Left(
            MixedModelError
              .InferenceUnavailable("simulation_nonfinite_draw", "simulation produced a non-finite response")
          )
    yield simulated

  /** Deterministically derives an independent indexed stream from a base seed. */
  def seedFor(base: Long, index: Long): Long = mix64(base + 0x9e3779b97f4a7c15L * (index + 1L))

  private def drawModes(
      fit: LmmFit,
      design: mixeff4s.design.CompiledDesign,
      seed: Long
  ): Vector[Vector[Double]] =
    val lambdas =
      Array.tabulate(design.reterms.length)(i => Array.fill(design.reterms(i).vsize, design.reterms(i).vsize)(0.0))
    var term = 0
    while term < lambdas.length do
      var diagonal = 0
      while diagonal < lambdas(term).length do
        lambdas(term)(diagonal)(diagonal) = 1.0
        diagonal += 1
      term += 1
    fit.parmap.zip(fit.theta).foreach { case ((block, row, column), value) =>
      lambdas(block)(row)(column) = value
    }
    var offset = 0L
    design.reterms.zipWithIndex.map { case (re, block) =>
      val values = Array.fill(re.nRanef)(0.0)
      var level = 0
      while level < re.nLevels do
        val z = Array.tabulate(re.vsize)(j => normal(seed, offset + level.toLong * re.vsize + j))
        var row = 0
        while row < re.vsize do
          var value = 0.0
          var column = 0
          while column <= row do
            value += lambdas(block)(row)(column) * z(column)
            column += 1
          values(level * re.vsize + row) = fit.sigma * value
          row += 1
        level += 1
      offset += re.nRanef.toLong
      values.toVector
    }

  /** Stateless SplitMix64 plus Box--Muller: reproducible on JVM and Scala.js. */
  private[mixeff4s] def normal(seed: Long, index: Long): Double =
    val pair = index >>> 1
    val u1 = uniform(seedFor(seed, pair * 2L))
    val u2 = uniform(seedFor(seed, pair * 2L + 1L))
    val radius = math.sqrt(-2.0 * math.log(u1))
    val angle = 2.0 * math.Pi * u2
    if (index & 1L) == 0L then radius * math.cos(angle) else radius * math.sin(angle)

  private[mixeff4s] def uniform(bits: Long): Double =
    // Open interval avoids log(0), using the high 53 bits without platform RNGs.
    (((bits >>> 12).toDouble) + 0.5) * math.pow(2.0, -52.0)

  private def mix64(value: Long): Long =
    var z = value
    z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L
    z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL
    z ^ (z >>> 31)
