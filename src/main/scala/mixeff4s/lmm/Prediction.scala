package mixeff4s.lmm

import mixeff4s.data.ModelFrame
import mixeff4s.error.{FitResult, MixedModelError}

enum PredictionMode:
  case Population, Conditional

enum NewGroupPolicy:
  case Refuse, ZeroEffect

/** Point predictions on the Gaussian response scale. New groups have no estimated conditional modes. */
object Prediction:
  def fitted(fit: LmmFit, mode: PredictionMode = PredictionMode.Conditional): FitResult[Vector[Double]] =
    fit.validatedContext.flatMap: context =>
      val design = context.design
      val values =
        Array.tabulate(fit.n)(i => (0 until design.p).iterator.map(j => design.xy.xy(i, j) * fit.beta(j)).sum)
      if mode == PredictionMode.Conditional then
        design.reterms.zipWithIndex.foreach: (term, block) =>
          var i = 0
          while i < fit.n do
            var j = 0
            while j < term.vsize do
              values(i) += term.z(j, i) * context.modes(block)(term.refs(i) * term.vsize + j)
              j += 1
            i += 1
      finite(values.toVector)

  def predict(
      fit: LmmFit,
      data: ModelFrame,
      mode: PredictionMode = PredictionMode.Conditional,
      newGroups: NewGroupPolicy = NewGroupPolicy.Refuse
  ): FitResult[Vector[Double]] =
    fit.validatedContext.flatMap: context =>
      context.design.prediction match
        case None =>
          Left(
            MixedModelError.InferenceUnavailable(
              "prediction_encoding_unavailable",
              "This compiled design has no frozen training encoding."
            )
          )
        case Some(encoding) =>
          if !encoding.matches(context.design) then
            Left(
              MixedModelError.InferenceUnavailable(
                "prediction_recipe_mismatch",
                "Prediction recipe does not match the fitted design."
              )
            )
          else
            encoding
              .encode(data, mode == PredictionMode.Conditional, newGroups == NewGroupPolicy.ZeroEffect)
              .flatMap: rows =>
                val values = Array.tabulate(data.nRows)(i =>
                  (0 until fit.beta.length).iterator.map(j => rows.fixed(i, j) * fit.beta(j)).sum
                )
                rows.random.zipWithIndex.foreach: (term, block) =>
                  var i = 0
                  while i < data.nRows do
                    if term.refs(i) >= 0 then
                      var j = 0
                      while j < term.basis.rows do
                        values(i) += term.basis(j, i) * context.modes(block)(term.refs(i) * term.basis.rows + j)
                        j += 1
                    i += 1
                finite(values.toVector)

  def residuals(fit: LmmFit): FitResult[Vector[Double]] =
    fitted(fit).map(values => fit.y.zip(values).map((y, value) => y - value))

  private def finite(values: Vector[Double]): FitResult[Vector[Double]] =
    if values.exists(value => !value.isFinite) then
      Left(MixedModelError.InvalidArgument("Prediction produced non-finite values."))
    else Right(values)
