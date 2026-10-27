package mixeff4s.glmm

import mixeff4s.design.CompiledDesign
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.lmm.PlsWorkspace
import mixeff4s.model.{Family, Link}
import mixeff4s.optimizer.{TrustBq, TrustBqOptions}

/** Score is for half the conditional deviance, including the spherical u penalty. */
final case class ConditionalModeResult(iterations: Int, scoreMaxAbs: Double, scoreTolerance: Double, halvings: Int)

/** Retained-constant objective at one [beta; theta], after a cold conditional solve. */
final case class JointEvaluation(
    objective: Double,
    minusTwoLogDensity: Double,
    penalty: Double,
    logdet: Double,
    sphericalModes: Vector[Vector[Double]],
    randomEffects: Vector[GlmmRandomEffects],
    eta: Vector[Double],
    mean: Vector[Double],
    conditional: ConditionalModeResult,
    quadratureCorrection: Double = 0.0
)

enum JointGradientCheck:
  case Available(values: Vector[Double], projectedMaxAbs: Double, tolerance: Double)
  case Unavailable(reasonCode: String)

  def passed: Boolean = this match
    case Available(_, norm, tolerance) => norm <= tolerance
    case Unavailable(_)                => false

enum JointBoundaryCheck:
  case Interior
  case Descent(parameter: Int, step: Double, objective: Double)
  case Unresolved
  case Unavailable(reasonCode: String)

/** Finite-difference first-order evidence, separately from numerical optimizer stopping. No Hessian, global optimum, or
  * finite-MLE existence certificate is implied.
  */
final case class JointAssessment(gradient: JointGradientCheck, boundary: JointBoundaryCheck):
  def firstOrderPassed: Boolean = gradient.passed && boundary == JointBoundaryCheck.Interior

private[glmm] object JointLaplace:
  def fit(
      design: CompiledDesign,
      family: Family,
      link: Link,
      response: GlmmResponse,
      options: GlmmOptions
  ): FitResult[GlmmFit] =
    JointWorkspace
      .create(design, family, link, response, options)
      .flatMap(workspace => fitWorkspace(design, family, link, response, options, workspace))

  private[glmm] def fitWorkspace(
      design: CompiledDesign,
      family: Family,
      link: Link,
      response: GlmmResponse,
      options: GlmmOptions,
      workspace: JointWorkspace
  ): FitResult[GlmmFit] =
    for
      _ <- workspace.pls.setTheta(design.theta)
      _ = workspace.replaceResponse(response)
      start = workspace.initialBeta ++ workspace.pls.theta
      lower = Vector.fill(design.p)(Double.NegativeInfinity) ++ workspace.pls.lowerBounds
      upper = Vector.fill(start.length)(Double.PositiveInfinity)
      controls = options.optimizerControls.getOrElse(
        TrustBqOptions(
          finalRadius = 1e-8,
          ftolAbs = 0.0,
          ftolRel = 0.0,
          stallFtolAbs = 0.0,
          stallFtolRel = 0.0,
          stallIterations = 20
        )
      )
      optimizer <- TrustBq.minimize(start, lower, upper, controls.copy(maxEvaluations = options.maxEvaluations)):
        parameters => workspace.evaluateParameters(parameters).map(_.objective)
      _ <- optimizer.requireConvergence
      evaluation <- workspace.evaluateParameters(optimizer.x)
      assessment = assess(workspace, optimizer.x, evaluation.objective, lower, options.stationarityTolerance)
      // Post-fit probes cannot leave the workspace at a probe or alter the returned snapshot.
      finalEvaluation <- workspace.evaluateParameters(optimizer.x)
    yield new GlmmFit(
      design.formula,
      family,
      link,
      if options.nAgq == 1 then Approximation.JointLaplace else Approximation.JointAgq,
      optimizer.x.drop(design.p),
      optimizer.x.take(design.p),
      finalEvaluation.objective,
      design.fe.fixedEffects,
      options.nAgq,
      optimizer,
      GlmmInnerResult.ConditionalMode(finalEvaluation.conditional),
      finalEvaluation.eta,
      finalEvaluation.mean,
      finalEvaluation.randomEffects,
      Some(assessment),
      new GlmmFitContext(design, response, options, family, link)
    )

  private[glmm] def assess(
      workspace: JointWorkspace,
      parameters: Vector[Double],
      objective: Double,
      lower: Vector[Double],
      tolerance: Double
  ): JointAssessment =
    val gradient = parameters.indices.foldLeft[FitResult[Vector[Double]]](Right(Vector.empty)): (acc, j) =>
      acc.flatMap: values =>
        val h = 1e-4 * math.max(1.0, math.abs(parameters(j)))
        def f(delta: Double) =
          workspace.evaluateParameters(parameters.updated(j, parameters(j) + delta)).map(_.objective)
        val derivative =
          if parameters(j) - h >= lower(j) then
            for
              plus <- f(h)
              minus <- f(-h)
              plusHalf <- f(h / 2.0)
              minusHalf <- f(-h / 2.0)
            yield (4.0 * (plusHalf - minusHalf) / h - (plus - minus) / (2.0 * h)) / 3.0
          else
            for
              plus <- f(h)
              plusTwo <- f(2.0 * h)
            yield (-3.0 * objective + 4.0 * plus - plusTwo) / (2.0 * h)
        derivative.flatMap: value =>
          if value.isFinite then Right(values :+ value)
          else Left(MixedModelError.Optimization("nonfinite joint finite-difference gradient"))
    val check = gradient match
      case Left(error)   => JointGradientCheck.Unavailable(error.code)
      case Right(values) =>
        val projected = values.indices.map: j =>
          if parameters(j) <= lower(j) + 1e-8 && values(j) >= 0.0 then 0.0 else math.abs(values(j))
        JointGradientCheck.Available(values, projected.maxOption.getOrElse(0.0), tolerance)
    val boundaryIndices = parameters.indices.filter(j => lower(j).isFinite && parameters(j) <= lower(j) + 1e-8)
    var boundary: JointBoundaryCheck =
      if boundaryIndices.isEmpty then JointBoundaryCheck.Interior else JointBoundaryCheck.Unresolved
    boundaryIndices.foreach: j =>
      Vector(1e-4, 1e-2, 1.0).foreach: step =>
        workspace.evaluateParameters(parameters.updated(j, parameters(j) + step)) match
          case Left(error) =>
            boundary match
              case JointBoundaryCheck.Descent(_, _, _) => ()
              case _                                   => boundary = JointBoundaryCheck.Unavailable(error.code)
          case Right(probe) if probe.objective < objective - 1e-8 =>
            boundary = JointBoundaryCheck.Descent(j, step, probe.objective)
          case _ => ()
    JointAssessment(check, boundary)

/** Dedicated joint workspace: only RE factors are read. Never reads profiled beta/FE covariance. */
private[glmm] final class JointWorkspace private (
    val pls: PlsWorkspace,
    family: Family,
    link: Link,
    private var response: GlmmResponse,
    options: GlmmOptions,
    quadrature: Option[AdaptiveQuadrature]
):
  private val n = pls.n
  private val p = pls.p
  private val u = pls.fittedReterms.map(rt => Array.fill(rt.nRanef)(0.0)).toArray
  private val eta = Array.fill(n)(0.0)
  private val mu = Array.fill(n)(0.0)
  private val sqrtWeights = Array.fill(n)(0.0)
  private val workingY = Array.fill(n)(0.0)
  private val scores = u.map(_.clone())
  private var beta = Vector.fill(p)(0.0)

  private[glmm] def replaceResponse(values: GlmmResponse): Unit =
    require(values.y.length == n)
    response = values

  def initialBeta: Vector[Double] =
    val intercept = pls.feNames.indexWhere(Pirls.isInterceptColumn)
    if intercept < 0 then Vector.fill(p)(0.0)
    else
      val mean = Pirls.initialResponseMean(family, response.y, response.workingWeights).getOrElse(0.5)
      Vector.fill(p)(0.0).updated(intercept, link.link(mean) - response.offset.sum / n.toDouble)

  def evaluateParameters(parameters: Vector[Double]): FitResult[JointEvaluation] =
    if parameters.length != p + pls.nTheta then
      Left(
        MixedModelError.DimensionMismatch(s"joint parameter length ${parameters.length}, expected ${p + pls.nTheta}")
      )
    else evaluate(parameters.take(p), parameters.drop(p))

  def evaluate(coefficients: Vector[Double], theta: Vector[Double]): FitResult[JointEvaluation] =
    if coefficients.length != p then
      Left(MixedModelError.DimensionMismatch(s"beta length ${coefficients.length}, expected $p"))
    else if coefficients.exists(v => !v.isFinite) then Left(MixedModelError.InvalidArgument("beta must be finite"))
    else
      pls
        .setTheta(theta)
        .flatMap: _ =>
          beta = coefficients
          u.foreach(java.util.Arrays.fill(_, 0.0))
          updateEta()
          solve().flatMap: result =>
            // Score and observed curvature use actual, unclipped canonical-link means.
            refresh(forSolve = false).flatMap: _ =>
              val score = scoreMaxAbs()
              val density = response.minusTwoLogDensity(eta)
              val penalty = u.iterator.flatMap(_.iterator).map(x => x * x).sum
              val logdet = pls.logdetRe
              val objective = density + penalty + logdet
              if !objective.isFinite || !score.isFinite then
                Left(MixedModelError.Optimization("joint Laplace accepted state is not finite"))
              else if score > options.conditionalScoreTolerance then
                Left(MixedModelError.PirlsNotConverged(result.iterations))
              else
                val modes = pls.fittedReterms
                  .zip(u)
                  .map: (rt, block) =>
                    val values = Vector.tabulate(rt.nLevels, rt.vsize): (level, row) =>
                      (0 until rt.vsize).map(col => rt.lambda(row, col) * block(level * rt.vsize + col)).sum
                    GlmmRandomEffects(rt.groupingName, rt.levels, rt.cnames, values)
                val base = JointEvaluation(
                  objective,
                  density,
                  penalty,
                  logdet,
                  u.map(_.toVector).toVector,
                  modes,
                  eta.toVector,
                  mu.toVector,
                  result.copy(scoreMaxAbs = score)
                )
                quadrature match
                  case None       => Right(base)
                  case Some(work) => work.correct(base, pls.theta.head, response)

  private def solve(): FitResult[ConditionalModeResult] =
    var iter = 0
    var halvings = 0
    var accepted = conditionalObjective
    if !accepted.isFinite then
      return Left(MixedModelError.Optimization("joint conditional starting objective is not finite"))
    while iter < options.jointMaxIterations do
      val oldScore = scoreMaxAbs()
      if oldScore <= options.conditionalScoreTolerance then
        return Right(ConditionalModeResult(iter, oldScore, options.conditionalScoreTolerance, halvings))
      refresh() match
        case Left(error) => return Left(error)
        case Right(_)    => ()
      val previous = u.map(_.clone())
      val proposed = pls.ranefUGivenBeta(beta)
      u.indices.foreach(t => System.arraycopy(proposed(t), 0, u(t), 0, u(t).length))
      updateEta()
      var objective = conditionalObjective
      val rounding =
        64.0 * math.ulp(1.0) * math.sqrt(n.toDouble + u.map(_.length).sum) * math.max(1.0, math.abs(accepted))
      def acceptable: Boolean = objective.isFinite &&
        (objective <= accepted || objective - accepted <= rounding && scoreMaxAbs() < oldScore)
      var half = 0
      while !acceptable && half < 20 do
        u.indices.foreach: t =>
          u(t).indices.foreach(j => u(t)(j) = 0.5 * (u(t)(j) + previous(t)(j)))
        updateEta()
        objective = conditionalObjective
        half += 1
      halvings += half
      if !acceptable then
        u.indices.foreach(t => System.arraycopy(previous(t), 0, u(t), 0, u(t).length))
        updateEta()
        return Left(MixedModelError.PirlsStepFailed(iter + 1, half))
      accepted = objective
      iter += 1
    val score = scoreMaxAbs()
    if score <= options.conditionalScoreTolerance then
      Right(ConditionalModeResult(iter, score, options.conditionalScoreTolerance, halvings))
    else Left(MixedModelError.PirlsNotConverged(iter))

  private def conditionalObjective: Double =
    response.minusTwoLogDensity(eta) + u.iterator.flatMap(_.iterator).map(x => x * x).sum

  private def refresh(forSolve: Boolean = true): FitResult[Unit] =
    var i = 0
    while i < n do
      val variance = link.muEta(eta(i))
      val weight = response.workingWeights(i) * variance
      if !eta(i).isFinite || !mu(i).isFinite || !weight.isFinite || weight < 0.0 then
        return Left(MixedModelError.Optimization("invalid joint conditional mean or curvature"))
      // Canonical Newton working response; zero curvature at a rounded endpoint
      // cannot yield a certified solve unless its actual score has already passed.
      val residual = response.scoreResidual(i, eta(i), mu(i))
      if forSolve && weight == 0.0 && residual != 0.0 then
        return Left(MixedModelError.Optimization("unresolvable zero joint conditional curvature"))
      sqrtWeights(i) = math.sqrt(weight)
      workingY(i) = if forSolve && weight > 0.0 then eta(i) + residual / weight - response.offset(i) else 0.0
      i += 1
    pls.updateIrlsWeights(sqrtWeights, workingY).flatMap(_ => pls.updateRandomL())

  private def scoreMaxAbs(): Double =
    u.indices.foreach(t => System.arraycopy(u(t), 0, scores(t), 0, u(t).length))
    val terms = pls.fittedReterms
    var i = 0
    while i < n do
      if !eta(i).isFinite || !mu(i).isFinite then return Double.PositiveInfinity
      val residual = response.scoreResidual(i, eta(i), mu(i))
      var t = 0
      while t < terms.length do
        val term = terms(t)
        var col = 0
        while col < term.vsize do
          var zLambda = 0.0
          var row = 0
          while row < term.vsize do
            zLambda += term.z(row, i) * term.lambda(row, col)
            row += 1
          scores(t)(term.refs(i) * term.vsize + col) -= zLambda * residual
          col += 1
        t += 1
      i += 1
    var maximum = 0.0
    scores.foreach: block =>
      block.foreach: score =>
        if !score.isFinite then maximum = Double.PositiveInfinity
        else maximum = math.max(maximum, math.abs(score))
    maximum

  private def updateEta(): Unit =
    var i = 0
    while i < n do
      var value = response.offset(i)
      var j = 0
      while j < p do
        value += pls.design.xy.xy(i, j) * beta(j)
        j += 1
      val terms = pls.fittedReterms
      var t = 0
      while t < terms.length do
        val term = terms(t)
        val block = u(t)
        var col = 0
        while col < term.vsize do
          var coefficient = 0.0
          var row = 0
          while row < term.vsize do
            coefficient += term.z(row, i) * term.lambda(row, col)
            row += 1
          value += coefficient * block(term.refs(i) * term.vsize + col)
          col += 1
        t += 1
      eta(i) = value
      mu(i) = link.linkinv(value)
      i += 1

private[glmm] object JointWorkspace:
  def create(
      design: CompiledDesign,
      family: Family,
      link: Link,
      response: GlmmResponse,
      options: GlmmOptions
  ): FitResult[JointWorkspace] =
    for
      quadrature <- AdaptiveQuadratureWorkspace.create(design, link, options)
      pls <- PlsWorkspace(design, reml = false)
    yield new JointWorkspace(pls, family, link, response, options, quadrature)
