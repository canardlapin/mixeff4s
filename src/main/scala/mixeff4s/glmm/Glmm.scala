package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.design.{CompiledDesign, Design}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.Formula
import mixeff4s.lmm.{Pls, PlsWorkspace}
import mixeff4s.model.{Family, Link}
import mixeff4s.optimizer.TrustBq

/** Profiled fast-PIRLS options. Joint Laplace/AGQ is refused until a later phase. */
final case class GlmmOptions(
    fast: Boolean = true,
    nAgq: Int = 1,
    weights: Option[Vector[Double]] = None,
    offset: Option[Vector[Double]] = None
)

object GlmmOptions:
  def fastLaplace: GlmmOptions = GlmmOptions()

/** How the GLMM objective was approximated. Never labelled as `lme4::glmer`. */
enum Approximation:
  case FastPirls

  def label: String =
    this match
      case Approximation.FastPirls => "fast-PIRLS"

/** Certified profiled fast-PIRLS fit. */
final case class GlmmFit(
    formula: Formula,
    family: Family,
    link: Link,
    approximation: Approximation,
    theta: Vector[Double],
    beta: Vector[Double],
    deviance: Double,
    feNames: Vector[String],
    nAgq: Int
):
  def algorithmLabel: String = approximation.label

/** Generalized mixed-model front door: compile, then labelled fast-PIRLS. */
object Glmm:
  def compile(formula: Formula, frame: ModelFrame): FitResult[CompiledDesign] =
    Design.compile(formula, frame)

  def compile(source: String, frame: ModelFrame): FitResult[CompiledDesign] =
    Formula.parse(source) match
      case Left(err)      => Left(MixedModelError.Formula(err))
      case Right(formula) => compile(formula, frame)

  def fit(
      formula: Formula,
      frame: ModelFrame,
      family: Family,
      link: Option[Link] = None,
      options: GlmmOptions = GlmmOptions.fastLaplace
  ): FitResult[GlmmFit] =
    val resolvedLink = link.getOrElse(family.canonicalLink)
    for
      _ <- Pirls.refuseFamilyLink(family, resolvedLink)
      _ <- refuseUnimplementedPath(options)
      design <- compile(formula, frame)
      y <- response(design)
      _ <- Pirls.validateResponse(family, y)
      _ <- refuseConstantResponse(y)
      weights <- options.weights match
        case None    => Right(Vector.fill(design.n)(1.0))
        case Some(w) =>
          Pirls.validateCaseWeights(w, design.n).map(_ => w)
      offset <- options.offset match
        case None    => Right(Vector.fill(design.n)(0.0))
        case Some(o) =>
          Pirls.validateOffset(o, design.n).map(_ => o)
      pls <- PlsWorkspace(design, reml = false)
      workspace = GlmmWorkspace(pls, family, resolvedLink, y, weights, offset)
      _ = workspace.initializeBeta()
      result <- TrustBq.minimize(
        workspace.theta,
        workspace.lowerBounds,
        workspace.upperBounds,
        Pls.smallFamilyOptions(workspace.nTheta)
      )(th => Right(workspace.devianceAt(th)))
      deviance <-
        val value = workspace.devianceAt(result.x)
        if value.isFinite then Right(value)
        else Left(MixedModelError.Optimization("fast-PIRLS deviance is not finite at the TrustBQ minimizer"))
    yield GlmmFit(
      formula,
      family,
      resolvedLink,
      Approximation.FastPirls,
      workspace.theta,
      workspace.currentBeta,
      deviance,
      workspace.feNames,
      options.nAgq
    )

  /** Laplace deviance at a fixed θ after inner PIRLS. Used by the scorecard gate. */
  def profiledDeviance(
      formula: Formula,
      frame: ModelFrame,
      family: Family,
      theta: Vector[Double],
      link: Option[Link] = None,
      options: GlmmOptions = GlmmOptions.fastLaplace
  ): FitResult[Double] =
    val resolvedLink = link.getOrElse(family.canonicalLink)
    for
      _ <- Pirls.refuseFamilyLink(family, resolvedLink)
      _ <- refuseUnimplementedPath(options)
      design <- compile(formula, frame)
      y <- response(design)
      _ <- Pirls.validateResponse(family, y)
      weights <- options.weights match
        case None    => Right(Vector.fill(design.n)(1.0))
        case Some(w) =>
          Pirls.validateCaseWeights(w, design.n).map(_ => w)
      offset <- options.offset match
        case None    => Right(Vector.fill(design.n)(0.0))
        case Some(o) =>
          Pirls.validateOffset(o, design.n).map(_ => o)
      pls <- PlsWorkspace(design, reml = false)
    yield
      val workspace = GlmmWorkspace(pls, family, resolvedLink, y, weights, offset)
      workspace.initializeBeta()
      workspace.devianceAt(theta)

  def fit(
      source: String,
      frame: ModelFrame,
      family: Family,
      link: Option[Link],
      options: GlmmOptions
  ): FitResult[GlmmFit] =
    Formula.parse(source) match
      case Left(err)      => Left(MixedModelError.Formula(err))
      case Right(formula) => fit(formula, frame, family, link, options)

  def fit(source: String, frame: ModelFrame, family: Family): FitResult[GlmmFit] =
    fit(source, frame, family, None, GlmmOptions.fastLaplace)

  private def refuseUnimplementedPath(options: GlmmOptions): FitResult[Unit] =
    if !options.fast || options.nAgq != 1 then
      Left(
        MixedModelError.Unsupported(
          "joint Laplace/AGQ is not implemented; Glmm.fit is labelled fast-PIRLS"
        )
      )
    else Right(())

  private def response(design: CompiledDesign): FitResult[Vector[Double]] =
    Right(Vector.tabulate(design.n)(i => design.xy.xy(i, design.p)))

  private def refuseConstantResponse(y: Vector[Double]): FitResult[Unit] =
    if y.nonEmpty && y.forall(_ == y.head) then Left(MixedModelError.ConstantResponse)
    else Right(())

private final class GlmmWorkspace(
    val pls: PlsWorkspace,
    val family: Family,
    val link: Link,
    val y: Vector[Double],
    val weights: Vector[Double],
    val offset: Vector[Double]
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

  def devianceAt(values: Vector[Double]): Double =
    (for
      _ <- pls.setTheta(values)
      _ <- pls.updateL()
      _ <- pirls()
    yield laplaceObjective) match
      case Right(obj) if obj.isFinite => obj
      case _                          => 1e300

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

  def pirls(): FitResult[Boolean] =
    var t = 0
    while t < u.length do
      java.util.Arrays.fill(u(t), 0.0)
      t += 1
    updateBAndEta()
    val uPrev = u.map(_.clone())
    val betaPrev = beta.clone()
    var obj0 = laplaceObjective
    var halvingBound = obj0 * 1.0001
    var converged = false
    var iter = 0
    while iter < Pirls.MaxIter do
      var obs = 0
      while obs < n do
        val (sw, wy) =
          Pirls.workingObservation(family, link, y(obs), eta(obs), mu(obs), weights(obs), offset(obs))
        sqrtwts(obs) = sw
        workingY(obs) = wy
        obs += 1
      pls.updateIrlsWeights(sqrtwts, workingY) match
        case Left(err) => return Left(err)
        case Right(_)  => ()
      pls.updateL() match
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
      if Pirls.converged(obj, obj0) then
        converged = true
        iter = Pirls.MaxIter
      else
        t = 0
        while t < u.length do
          System.arraycopy(u(t), 0, uPrev(t), 0, u(t).length)
          t += 1
        System.arraycopy(beta, 0, betaPrev, 0, p)
        obj0 = obj
        halvingBound = obj
        iter += 1
    Right(converged)

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
