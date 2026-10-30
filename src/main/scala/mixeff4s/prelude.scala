package mixeff4s

/** Common imports for describing a mixed model and reading a refusal. */
object prelude:
  export mixeff4s.data.{Column, ModelFrame}
  export mixeff4s.error.{FitResult, LinAlgError, MixedModelError}
  export mixeff4s.formula.{FixedTerm, Formula, FormulaError, GroupingFactor, RandomCovariance, RandomTerm}
  export mixeff4s.formula.dsl.dsl.*
  export mixeff4s.design.{CompiledDesign, FeMat, FeTerm, FixedEffectLayout, ReMat}
  export mixeff4s.lmm.{
    BatchStart,
    CoefTable,
    CompiledLmm,
    Criterion,
    FitOptions,
    BoundaryAdequacy,
    GradientCheck,
    Lmm,
    LmmAssessment,
    LmmAssessmentOptions,
    LmmFit,
    LmmOptimizer,
    LmmRefitSession,
    ProfiledGradientOptions,
    Prediction,
    PredictionMode,
    NewGroupPolicy,
    Simulation,
    SimulationMode,
    RefitStart,
    VarCorr,
    VarCorrComponent
  }
  export mixeff4s.glmm.{
    Approximation,
    Glmm,
    GlmmFit,
    CompiledGlmm,
    GlmmRefitSession,
    GlmmPrediction,
    GlmmPredictionScale,
    GlmmResidualScale,
    GlmmSimulation,
    GlmmAssessment,
    GlmmAssessmentOptions,
    FixedEffectRecession,
    GlmmCovarianceRank,
    GlmmOptions,
    GlmmInnerResult,
    GlmmRandomEffects,
    ConditionalModeResult,
    JointEvaluation,
    JointAssessment,
    JointGradientCheck,
    JointBoundaryCheck
  }
  export mixeff4s.stats.{
    Bootstrap,
    BootstrapOptions,
    BootstrapTarget,
    BootstrapResult,
    BootstrapAttempt,
    BootstrapInterval,
    GlmmBootstrap,
    GlmmBootstrapOptions,
    GlmmBootstrapAttempt,
    GlmmBootstrapResult,
    Lrt,
    LrtStatistic,
    GlmmLrt,
    GlmmInference,
    GlmmInferenceOptions,
    GlmmInferenceResult,
    GlmmCoefficientInterval,
    CurvatureOptions,
    CurvatureEvidence,
    LikelihoodBootstrap,
    LikelihoodBootstrapOptions,
    LikelihoodBootstrapAttempt,
    LikelihoodBootstrapResult,
    MonteCarloTail,
    LmmProfile,
    LmmProfileTarget,
    LmmProfileOptions,
    LmmProfileResult,
    LmmProfileEvaluation,
    LmmProfileEndpoint,
    Profile
  }
  export mixeff4s.optimizer.TrustBqOptions
  export mixeff4s.model.{Family, Link}

  def numeric(values: Iterable[Double]): Column.Numeric = ModelFrame.numeric(values)
  def factorCol(values: Iterable[String]): Column.Factor = ModelFrame.factor(values)
