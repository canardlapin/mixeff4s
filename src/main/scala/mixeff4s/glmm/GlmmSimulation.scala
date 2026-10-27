package mixeff4s.glmm

import mixeff4s.design.CompiledDesign
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.lmm.{Simulation, SimulationMode}
import mixeff4s.model.{Family, Link}

/** Seeded conditional-response simulation from an immutable GLMM fit.
  *
  * A random-effect draw is Lambda*z (there is no Gaussian residual sigma in a GLMM). Each observation owns a separate
  * count-sampler stream, so a rejected proposal cannot shift any later response or a random-effect draw.
  */
object GlmmSimulation:
  def response(
      fit: GlmmFit,
      seed: Long,
      mode: SimulationMode = SimulationMode.DrawRandomEffects
  ): FitResult[Vector[Double]] =
    fit.validatedContext.flatMap: context =>
      for
        _ <- supported(fit)
        _ <-
          if context.response.caseWeights.forall(_ == 1.0) then Right(())
          else
            Left(
              MixedModelError.InferenceUnavailable(
                "simulation_case_weights",
                "Simulation refuses non-unit case weights: power likelihood is not a sampling law."
              )
            )
        _ <-
          if fit.beta.forall(_.isFinite) && fit.theta.forall(_.isFinite) then Right(())
          else
            Left(
              MixedModelError.InferenceUnavailable("simulation_unavailable", "fit has non-finite simulation parameters")
            )
        random = mode match
          case SimulationMode.DrawRandomEffects => drawModes(fit, context.design, seed)
          case SimulationMode.ConditionOnModes  => fit.randomEffects.map(_.values.flatten)
        eta <- eta(fit, context.design, context.response.offset, random)
        mean = eta.map(fit.link.linkinv)
        _ <-
          if mean.forall(_.isFinite) then Right(())
          else
            Left(
              MixedModelError
                .InferenceUnavailable("simulation_nonfinite_draw", "simulation produced a non-finite conditional mean")
            )
        values <- sample(fit.family, context.response.trials, mean, seed)
      yield values

  private def supported(fit: GlmmFit): FitResult[Unit] =
    (fit.family, fit.link) match
      case (Family.Bernoulli | Family.Binomial, Link.Logit) | (Family.Poisson, Link.Log) => Right(())
      case _ => Left(MixedModelError.UnsupportedFamilyLink(fit.family.label, fit.link.label))

  private def eta(
      fit: GlmmFit,
      design: CompiledDesign,
      offsets: Vector[Double],
      random: Vector[Vector[Double]]
  ): FitResult[Vector[Double]] =
    val values = Array.tabulate(design.n)(i =>
      offsets(i) + (0 until design.p).iterator.map(j => design.xy.xy(i, j) * fit.beta(j)).sum
    )
    design.reterms.zipWithIndex.foreach: (term, block) =>
      var i = 0
      while i < design.n do
        var j = 0
        while j < term.vsize do
          values(i) += term.z(j, i) * random(block)(term.refs(i) * term.vsize + j)
          j += 1
        i += 1
    if values.forall(_.isFinite) then Right(values.toVector)
    else
      Left(
        MixedModelError.InferenceUnavailable(
          "simulation_nonfinite_draw",
          "simulation produced a non-finite linear predictor"
        )
      )

  private def sample(family: Family, trials: Vector[Int], mean: Vector[Double], seed: Long): FitResult[Vector[Double]] =
    mean.indices.foldLeft[FitResult[Vector[Double]]](Right(Vector.empty)): (acc, i) =>
      acc.flatMap: values =>
        val stream = Simulation.seedFor(seed, i.toLong + 1L)
        val draw = family match
          case Family.Bernoulli => CountSamplers.bernoulli(mean(i), stream)
          case Family.Binomial  => CountSamplers.binomial(trials(i), mean(i), stream).map(_ / trials(i).toDouble)
          case Family.Poisson   => CountSamplers.poisson(mean(i), stream)
          case _                => Left(MixedModelError.UnsupportedFamilyLink(family.label, ""))
        draw.map(values :+ _)

  private def drawModes(fit: GlmmFit, design: CompiledDesign, seed: Long): Vector[Vector[Double]] =
    val lambdas =
      Array.tabulate(design.reterms.length)(i => Array.fill(design.reterms(i).vsize, design.reterms(i).vsize)(0.0))
    var term = 0
    while term < lambdas.length do
      var diagonal = 0
      while diagonal < lambdas(term).length do { lambdas(term)(diagonal)(diagonal) = 1.0; diagonal += 1 }
      term += 1
    fit.context.design.parmap.zip(fit.theta).foreach { case ((block, row, column), value) =>
      lambdas(block)(row)(column) = value
    }
    var offset = 0L
    design.reterms.zipWithIndex.map: (re, block) =>
      val values = Array.fill(re.nRanef)(0.0)
      var level = 0
      while level < re.nLevels do
        var row = 0
        while row < re.vsize do
          var value = 0.0
          var column = 0
          while column <= row do
            value += lambdas(block)(row)(column) * Simulation.normal(
              Simulation.seedFor(seed, -1L),
              offset + level.toLong * re.vsize + column
            )
            column += 1
          values(level * re.vsize + row) = value
          row += 1
        level += 1
      offset += re.nRanef.toLong
      values.toVector
