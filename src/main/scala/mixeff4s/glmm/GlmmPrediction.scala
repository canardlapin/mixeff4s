package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.lmm.{NewGroupPolicy, PredictionMode}

/** Population predictions set every random effect to zero; they are not marginal means. */
enum GlmmPredictionScale:
  case Link, Response

enum GlmmResidualScale:
  case Raw, Pearson

object GlmmPrediction:
  def fitted(
      fit: GlmmFit,
      mode: PredictionMode = PredictionMode.Conditional,
      scale: GlmmPredictionScale = GlmmPredictionScale.Response
  ): FitResult[Vector[Double]] =
    fit.validatedContext.flatMap: context =>
      linear(fit, context.design, mode, context.response.offset).flatMap(transform(fit, _, scale))

  def predict(
      fit: GlmmFit,
      data: ModelFrame,
      mode: PredictionMode = PredictionMode.Conditional,
      scale: GlmmPredictionScale = GlmmPredictionScale.Response,
      unknownGroups: NewGroupPolicy = NewGroupPolicy.Refuse,
      offset: Option[Vector[Double]] = None
  ): FitResult[Vector[Double]] =
    fit.validatedContext.flatMap: context =>
      for
        offsets <- predictionOffset(fit, data.nRows, offset)
        encoding <- context.design.prediction.toRight(
          MixedModelError.InferenceUnavailable(
            "prediction_encoding_unavailable",
            "This compiled design has no frozen training encoding."
          )
        )
        _ <-
          if encoding.matches(context.design) then Right(())
          else
            Left(
              MixedModelError.InferenceUnavailable(
                "prediction_recipe_mismatch",
                "Prediction recipe does not match the fitted design."
              )
            )
        rows <- encoding.encode(data, mode == PredictionMode.Conditional, unknownGroups == NewGroupPolicy.ZeroEffect)
        eta <- rowsEta(fit, data.nRows, rows.fixed, rows.random, mode, offsets)
        values <- transform(fit, eta, scale)
      yield values

  def residuals(fit: GlmmFit, kind: GlmmResidualScale = GlmmResidualScale.Raw): FitResult[Vector[Double]] =
    fitted(fit, PredictionMode.Conditional, GlmmPredictionScale.Response).flatMap: mean =>
      val response = fit.context.response
      val values = response.y.indices.map: i =>
        val raw = response.y(i) - mean(i)
        kind match
          case GlmmResidualScale.Raw     => raw
          case GlmmResidualScale.Pearson =>
            raw * math.sqrt(response.caseWeights(i) / (fit.family.variance(mean(i)) / response.trials(i).toDouble))
      finite(values.toVector, "residual_nonfinite", "Residual calculation produced a non-finite value.")

  private def predictionOffset(fit: GlmmFit, n: Int, supplied: Option[Vector[Double]]): FitResult[Vector[Double]] =
    supplied.orElse(if fit.options.offset.nonEmpty then None else Some(Vector.fill(n)(0.0))) match
      case None =>
        Left(
          MixedModelError.InferenceUnavailable(
            "prediction_offset_required",
            "New-data offsets are required because the fit used an offset."
          )
        )
      case Some(values) if values.length != n =>
        Left(MixedModelError.DimensionMismatch("prediction offset length must equal new-data rows"))
      case Some(values) if values.exists(value => !value.isFinite) =>
        Left(MixedModelError.InvalidArgument("prediction offsets must be finite"))
      case Some(values) => Right(values)

  private def linear(
      fit: GlmmFit,
      design: mixeff4s.design.CompiledDesign,
      mode: PredictionMode,
      offsets: Vector[Double]
  ): FitResult[Vector[Double]] =
    val values = Array.tabulate(design.n)(i =>
      offsets(i) + (0 until design.p).iterator.map(j => design.xy.xy(i, j) * fit.beta(j)).sum
    )
    if mode == PredictionMode.Conditional then
      design.reterms.zipWithIndex.foreach: (term, block) =>
        var i = 0
        while i < design.n do
          var j = 0
          while j < term.vsize do
            values(i) += term.z(j, i) * fit.randomEffects(block).values(term.refs(i))(j)
            j += 1
          i += 1
    finite(values.toVector, "prediction_nonfinite", "Prediction produced a non-finite value.")

  private def rowsEta(
      fit: GlmmFit,
      n: Int,
      fixed: gale.linalg.DMat,
      random: Vector[mixeff4s.design.PredictionRandomRows],
      mode: PredictionMode,
      offsets: Vector[Double]
  ): FitResult[Vector[Double]] =
    val values =
      Array.tabulate(n)(i => offsets(i) + (0 until fit.beta.length).iterator.map(j => fixed(i, j) * fit.beta(j)).sum)
    if mode == PredictionMode.Conditional then
      random.zipWithIndex.foreach: (term, block) =>
        var i = 0
        while i < n do
          if term.refs(i) >= 0 then
            var j = 0
            while j < term.basis.rows do
              values(i) += term.basis(j, i) * fit.randomEffects(block).values(term.refs(i))(j)
              j += 1
          i += 1
    finite(values.toVector, "prediction_nonfinite", "Prediction produced a non-finite value.")

  private def transform(fit: GlmmFit, eta: Vector[Double], scale: GlmmPredictionScale): FitResult[Vector[Double]] =
    val values = scale match
      case GlmmPredictionScale.Link     => eta
      case GlmmPredictionScale.Response => eta.map(fit.link.linkinv)
    finite(values, "prediction_nonfinite", "Prediction produced a non-finite value.")

  private def finite(values: Vector[Double], code: String, details: String): FitResult[Vector[Double]] =
    if values.forall(_.isFinite) then Right(values) else Left(MixedModelError.InferenceUnavailable(code, details))
