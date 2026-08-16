package mixeff4s

/** Common imports for describing a mixed model and reading a refusal. */
object prelude:
  export mixeff4s.data.{Column, ModelFrame}
  export mixeff4s.error.{FitResult, LinAlgError, MixedModelError}
  export mixeff4s.formula.{
    FixedTerm,
    Formula,
    FormulaError,
    GroupingFactor,
    RandomCovariance,
    RandomTerm
  }
  export mixeff4s.formula.dsl.dsl.*
  export mixeff4s.design.{CompiledDesign, FeMat, FeTerm, ReMat}
  export mixeff4s.lmm.{Criterion, FitOptions, Lmm, LmmFit}
  export mixeff4s.model.{Family, Link}

  def numeric(values: Iterable[Double]): Column.Numeric = ModelFrame.numeric(values)
  def factorCol(values: Iterable[String]): Column.Factor = ModelFrame.factor(values)
