package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.design.{CompiledDesign, Design, FixedEffectLayout}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.Formula
import mixeff4s.lmm.{FitOptions, Pls, PlsWorkspace}
import mixeff4s.model.{Family, Link}
import mixeff4s.optimizer.{TrustBq, TrustBqOptions, TrustBqResult}

/** Binomial y is a proportion; trials are explicit. Positive case weights multiply the entire log density. Joint v1
  * supports canonical binomial/logit and Poisson/log, with nAgq=1.
  */
final case class GlmmOptions(
    fast: Boolean = true,
    nAgq: Int = 1,
    weights: Option[Vector[Double]] = None,
    offset: Option[Vector[Double]] = None,
    maxEvaluations: Int = 1000,
    pirlsMaxIterations: Int = Pirls.MaxIter,
    trials: Option[Vector[Int]] = None,
    jointMaxIterations: Int = 100,
    conditionalScoreTolerance: Double = 1e-8,
    optimizerControls: Option[TrustBqOptions] = None,
    stationarityTolerance: Double = 1e-4
)

object GlmmOptions:
  def fastLaplace: GlmmOptions = GlmmOptions()
  def jointLaplace: GlmmOptions = GlmmOptions(fast = false, maxEvaluations = 10000)

/** Fast and joint objectives have explicit, different definitions. */
enum Approximation:
  case FastPirls, JointLaplace

  def label: String = this match
    case FastPirls    => "fast-PIRLS"
    case JointLaplace => "joint-Laplace"

  def objectiveDefinition: String = this match
    case FastPirls    => "deviance_residuals_plus_penalty_plus_logdet"
    case JointLaplace => "minus_two_loglikelihood_with_response_constants"

/** Inner-loop criteria are distinct types; neither implies outer stationarity. */
enum GlmmInnerResult:
  case FastPirls(result: PirlsResult)
  case ConditionalMode(result: ConditionalModeResult)

/** Link-scale conditional coefficients, with fitted term, level and column order. */
final case class GlmmRandomEffects(
    groupingName: String,
    levels: Vector[String],
    coefficientNames: Vector[String],
    values: Vector[Vector[Double]]
):
  require(values.length == levels.length && values.forall(_.length == coefficientNames.length))

/** An immutable fitted snapshot; construction and copying cannot detach numerical evidence. */
final class GlmmFit private[glmm] (
    val formula: Formula,
    val family: Family,
    val link: Link,
    val approximation: Approximation,
    val theta: Vector[Double],
    val beta: Vector[Double],
    val deviance: Double,
    val fixedEffects: FixedEffectLayout,
    val nAgq: Int,
    val optimizer: TrustBqResult,
    val inner: GlmmInnerResult,
    val conditionalEta: Vector[Double],
    val conditionalMean: Vector[Double],
    val randomEffects: Vector[GlmmRandomEffects],
    val jointAssessment: Option[JointAssessment] = None
):
  require(beta.length == fixedEffects.rank, "beta length must match retained fixed-effect layout")
  def feNames: Vector[String] = fixedEffects.retainedNames
  def algorithmLabel: String = approximation.label
  def objectiveDefinition: String = approximation.objectiveDefinition

/** Generalized mixed models with separately labelled estimators. */
object Glmm:
  def compile(formula: Formula, frame: ModelFrame): FitResult[CompiledDesign] = Design.compile(formula, frame)
  def compile(source: String, frame: ModelFrame): FitResult[CompiledDesign] =
    Formula.parse(source).left.map(MixedModelError.Formula.apply).flatMap(compile(_, frame))

  def fit(
      formula: Formula,
      frame: ModelFrame,
      family: Family,
      link: Option[Link] = None,
      options: GlmmOptions = GlmmOptions.fastLaplace
  ): FitResult[GlmmFit] =
    val resolved = link.getOrElse(family.canonicalLink)
    for
      _ <- validate(family, resolved, options)
      design <- compile(formula, frame)
      response <- GlmmResponse.prepare(family, responseValues(design), options)
      _ <-
        if response.y.nonEmpty && response.y.forall(_ == response.y.head) then Left(MixedModelError.ConstantResponse)
        else Right(())
      fit <-
        if options.fast then fastFit(design, family, resolved, response, options)
        else JointLaplace.fit(design, family, resolved, response, options)
    yield fit

  private def fastFit(
      design: CompiledDesign,
      family: Family,
      link: Link,
      response: GlmmResponse,
      options: GlmmOptions
  ): FitResult[GlmmFit] =
    for
      workspace <- fastWorkspace(design, family, link, response, options)
      result <- TrustBq.minimize(
        workspace.theta,
        workspace.lowerBounds,
        workspace.upperBounds,
        options.optimizerControls
          .getOrElse(Pls.smallFamilyOptions(workspace.nTheta))
          .copy(maxEvaluations = options.maxEvaluations)
      )(th => workspace.evaluate(th).map(_.deviance))
      _ <- result.requireConvergence
      evaluation <- workspace.evaluate(result.x)
    yield new GlmmFit(
      design.formula,
      family,
      link,
      Approximation.FastPirls,
      workspace.theta,
      workspace.currentBeta,
      evaluation.deviance,
      design.fe.fixedEffects,
      options.nAgq,
      result,
      GlmmInnerResult.FastPirls(evaluation.pirls),
      workspace.currentEta,
      workspace.currentMu,
      workspace.currentModes
    )

  /** Fast-PIRLS objective at theta; never dispatches to a different estimator. */
  def profiledDeviance(
      formula: Formula,
      frame: ModelFrame,
      family: Family,
      theta: Vector[Double],
      link: Option[Link] = None,
      options: GlmmOptions = GlmmOptions.fastLaplace
  ): FitResult[Double] =
    val resolved = link.getOrElse(family.canonicalLink)
    for
      _ <- validate(family, resolved, options)
      _ <-
        if options.fast then Right(())
        else Left(MixedModelError.InvalidArgument("profiledDeviance requires fast=true; use jointDeviance"))
      design <- compile(formula, frame)
      response <- GlmmResponse.prepare(family, responseValues(design), options)
      workspace <- fastWorkspace(design, family, resolved, response, options)
      evaluation <- workspace.evaluate(theta)
    yield evaluation.deviance

  /** Evaluate retained-basis beta and parmap theta with a cold, score-qualified conditional solve. */
  def jointDeviance(
      formula: Formula,
      frame: ModelFrame,
      family: Family,
      beta: Vector[Double],
      theta: Vector[Double],
      link: Option[Link] = None,
      options: GlmmOptions = GlmmOptions.jointLaplace
  ): FitResult[JointEvaluation] =
    val resolved = link.getOrElse(family.canonicalLink)
    for
      _ <- validate(family, resolved, options)
      _ <-
        if !options.fast then Right(()) else Left(MixedModelError.InvalidArgument("jointDeviance requires fast=false"))
      design <- compile(formula, frame)
      response <- GlmmResponse.prepare(family, responseValues(design), options)
      workspace <- JointWorkspace.create(design, family, resolved, response, options)
      result <- workspace.evaluate(beta, theta)
    yield result

  def fit(
      source: String,
      frame: ModelFrame,
      family: Family,
      link: Option[Link],
      options: GlmmOptions
  ): FitResult[GlmmFit] =
    Formula.parse(source).left.map(MixedModelError.Formula.apply).flatMap(fit(_, frame, family, link, options))
  def fit(source: String, frame: ModelFrame, family: Family): FitResult[GlmmFit] =
    fit(source, frame, family, None, GlmmOptions.fastLaplace)

  private def validate(family: Family, link: Link, options: GlmmOptions): FitResult[Unit] =
    for
      _ <- Pirls.refuseFamilyLink(family, link)
      _ <- GlmmResponse.validatePath(family, link, options)
      _ <- FitOptions.ml
        .copy(maxEvaluations = options.maxEvaluations, optimizerControls = options.optimizerControls)
        .validate
    yield ()

  private def responseValues(design: CompiledDesign): Vector[Double] =
    Vector.tabulate(design.n)(i => design.xy.xy(i, design.p))

  private def fastWorkspace(
      design: CompiledDesign,
      family: Family,
      link: Link,
      response: GlmmResponse,
      options: GlmmOptions
  ): FitResult[GlmmWorkspace] =
    PlsWorkspace(design, reml = false).map(pls =>
      GlmmWorkspace(pls, family, link, response.y, response.workingWeights, response.offset, options.pirlsMaxIterations)
    )

private final case class GlmmEvaluation(deviance: Double, pirls: PirlsResult)

private final class GlmmWorkspace(
    val pls: PlsWorkspace,
    val family: Family,
    val link: Link,
    val y: Vector[Double],
    val weights: Vector[Double],
    val offset: Vector[Double],
    val maxIterations: Int
):
  private val n = pls.n
  private val p = pls.p
  private val beta = Array.fill(p)(0.0)
  private val u = pls.fittedReterms.map(rt => Array.fill(rt.nRanef)(0.0)).toArray
  private val eta = Array.fill(n)(0.0)
  private val mu = Array.fill(n)(0.0)
  private val sqrtwts = Array.fill(n)(0.0)
  private val workingY = Array.fill(n)(0.0)

  def nTheta: Int = pls.nTheta
  def theta: Vector[Double] = pls.theta
  def feNames: Vector[String] = pls.feNames
  def lowerBounds: Vector[Double] = pls.lowerBounds
  def upperBounds: Vector[Double] = pls.upperBounds
  def currentBeta: Vector[Double] = beta.toVector
  def currentEta: Vector[Double] = eta.toVector
  def currentMu: Vector[Double] = mu.toVector
  def currentModes: Vector[GlmmRandomEffects] =
    pls.fittedReterms
      .zip(u)
      .map: (rt, modes) =>
        val values = Vector.tabulate(rt.nLevels, rt.vsize): (level, row) =>
          (0 until rt.vsize).map(col => rt.lambda(row, col) * modes(level * rt.vsize + col)).sum
        GlmmRandomEffects(rt.groupingName, rt.levels, rt.cnames, values)

  def initializeBeta(): Unit =
    java.util.Arrays.fill(beta, 0.0)
    val intercept = pls.feNames.indexWhere(Pirls.isInterceptColumn)
    if intercept >= 0 then
      Pirls
        .initialResponseMean(family, y, weights)
        .foreach: mean =>
          val offsetMean = offset.sum / offset.length.toDouble
          beta(intercept) = link.link(mean) - offsetMean
    updateEta()

  def evaluate(values: Vector[Double]): FitResult[GlmmEvaluation] =
    for
      _ <- pls.setTheta(values)
      _ = initializeBeta()
      result <- pirls()
      obj = laplaceObjective
      _ <-
        if obj.isFinite then Right(())
        else Left(MixedModelError.Optimization("fast-PIRLS deviance is not finite"))
    yield GlmmEvaluation(obj, result)

  def laplaceObjective: Double =
    var dev = 0.0
    var i = 0
    while i < n do
      dev += weights(i) * Pirls.devianceComponent(family, y(i), mu(i))
      i += 1
    var penalty = 0.0
    var t = 0
    while t < u.length do
      var j = 0
      while j < u(t).length do
        val v = u(t)(j)
        penalty += v * v
        j += 1
      t += 1
    dev + penalty + pls.logdetRe

  def pirls(): FitResult[PirlsResult] =
    var t = 0
    while t < u.length do
      java.util.Arrays.fill(u(t), 0.0)
      t += 1
    updateBAndEta()
    refreshWorking() match
      case Left(err) => return Left(err)
      case Right(_)  => ()
    val uPrev = u.map(_.clone())
    val betaPrev = beta.clone()
    var obj0 = laplaceObjective
    var halvingBound = obj0 * 1.0001
    var iter = 0
    var halvingLimitHits = 0
    while iter < maxIterations do
      refreshWorking() match
        case Left(err) => return Left(err)
        case Right(_)  => ()
      val newBeta = pls.beta
      var q = 0
      while q < p do
        beta(q) = newBeta(q)
        q += 1
      val newU = pls.ranefU
      t = 0
      while t < u.length do
        System.arraycopy(newU(t), 0, u(t), 0, u(t).length)
        t += 1
      updateBAndEta()
      var obj = laplaceObjective
      var nhalf = 0
      while (!obj.isFinite || obj > halvingBound) && nhalf < Pirls.MaxHalvings do
        nhalf += 1
        t = 0
        while t < u.length do
          var j = 0
          while j < u(t).length do
            u(t)(j) = 0.5 * (u(t)(j) + uPrev(t)(j))
            j += 1
          t += 1
        q = 0
        while q < p do
          beta(q) = 0.5 * (beta(q) + betaPrev(q))
          q += 1
        updateBAndEta()
        obj = laplaceObjective
      if !obj.isFinite then return Left(MixedModelError.PirlsStepFailed(iter + 1, nhalf))
      // Fast-PIRLS retains its objective-change rule even after finite halving
      // exhaustion. The stricter conditional-score rule belongs to joint GLMM.
      if obj > halvingBound then halvingLimitHits += 1
      if Pirls.converged(obj, obj0) then
        // The returned determinant must use the accepted mean, including a halved step.
        return refreshWorking().map(_ =>
          PirlsResult(iter + 1, math.abs(obj - obj0), Pirls.ConvergenceTol, halvingLimitHits)
        )
      else
        t = 0
        while t < u.length do
          System.arraycopy(u(t), 0, uPrev(t), 0, u(t).length)
          t += 1
        System.arraycopy(beta, 0, betaPrev, 0, p)
        obj0 = obj
        halvingBound = obj
        iter += 1
    Left(MixedModelError.PirlsNotConverged(iter))

  private def refreshWorking(): FitResult[Unit] =
    var obs = 0
    while obs < n do
      val (sw, wy) = Pirls.workingObservation(family, link, y(obs), eta(obs), mu(obs), weights(obs), offset(obs))
      sqrtwts(obs) = sw
      workingY(obs) = wy
      obs += 1
    pls.updateIrlsWeights(sqrtwts, workingY).flatMap(_ => pls.updateL())

  private def updateEta(): Unit =
    updateBAndEta()

  private def updateBAndEta(): Unit =
    val reterms = pls.fittedReterms
    var obs = 0
    while obs < n do
      var value = offset(obs)
      var q = 0
      while q < p do
        value += pls.design.xy.xy(obs, q) * beta(q)
        q += 1
      eta(obs) = value
      obs += 1
    var t = 0
    while t < reterms.length do
      val rt = reterms(t)
      val ut = u(t)
      obs = 0
      while obs < n do
        val r = rt.refs(obs)
        var s = 0
        while s < rt.vsize do
          var b = 0.0
          var col = 0
          while col < rt.vsize do
            b += rt.lambda(s, col) * ut(r * rt.vsize + col)
            col += 1
          eta(obs) += rt.z(s, obs) * b
          s += 1
        obs += 1
      t += 1
    obs = 0
    while obs < n do
      mu(obs) = link.linkinv(eta(obs))
      obs += 1
