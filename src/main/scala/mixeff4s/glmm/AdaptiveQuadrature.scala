package mixeff4s.glmm

import gale.linalg.Matrix
import gale.spectral.{Eigen, EigenSelection, EigenVectors}
import mixeff4s.design.CompiledDesign
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.model.Link

/** Normalized Gauss-Hermite rule for a standard-normal measure. */
private[glmm] final case class NormalQuadrature(nodes: Vector[Double], logWeights: Vector[Double])

private[glmm] object AdaptiveQuadrature:
  val MaxNodes: Int = 31
  private val rules = scala.collection.mutable.Map.empty[Int, FitResult[NormalQuadrature]]

  def validate(design: CompiledDesign, options: GlmmOptions): FitResult[Unit] =
    if options.nAgq > 1 && (options.fast || design.reterms.length != 1 || design.reterms.head.vsize != 1) then
      Left(MixedModelError.Unsupported("AGQ requires fast=false and exactly one scalar random-effect term"))
    else Right(())

  def rule(count: Int): FitResult[NormalQuadrature] = synchronized:
    if count < 1 || count > MaxNodes then
      Left(MixedModelError.Unsupported(s"quadrature order must be between 1 and $MaxNodes"))
    else rules.getOrElseUpdate(count, build(count))

  private def build(count: Int): FitResult[NormalQuadrature] =
    if count == 1 then Right(NormalQuadrature(Vector(0.0), Vector(0.0)))
    else
      // Golub-Welsch: probabilists' orthonormal Hermite Jacobi matrix.
      val jacobi = Matrix.tabulate(count, count): (i, j) =>
        if math.abs(i - j) == 1 then math.sqrt(math.max(i, j).toDouble) else 0.0
      Eigen
        .eigSymmetric(jacobi, EigenSelection.All, EigenVectors.Right)
        .left
        .map(error => MixedModelError.Optimization(s"quadrature rule eigendecomposition: $error"))
        .flatMap: decomposition =>
          val nodes = Vector.tabulate(count)(decomposition.eigenvalues(_))
          val weights = Vector.tabulate(count)(j => math.pow(decomposition.eigenvectors(0, j), 2))
          val sum = weights.sum
          if nodes.exists(v => !v.isFinite) || weights
              .exists(w => !w.isFinite || w <= 0.0) || !sum.isFinite || sum <= 0.0
          then Left(MixedModelError.Optimization("quadrature rule has invalid nodes or weights"))
          else Right(NormalQuadrature(nodes, weights.map(w => math.log(w / sum))))

/** Cached grouping topology; the sweep never mutates accepted modes or numerical factors. */
private[glmm] final class AdaptiveQuadrature private[glmm] (
    design: CompiledDesign,
    link: Link,
    rule: NormalQuadrature
):
  private val term = design.reterms.head
  private val rows =
    val buffers = Array.fill(term.nLevels)(scala.collection.mutable.ArrayBuffer.empty[Int])
    term.refs.indices.foreach(i => buffers(term.refs(i)) += i)
    buffers.map(_.toArray)

  def correct(base: JointEvaluation, theta: Double, response: GlmmResponse): FitResult[JointEvaluation] =
    if theta == 0.0 then Right(base)
    else
      var totalLogCorrection = 0.0
      var group = 0
      while group < rows.length do
        val mode = base.sphericalModes.head(group)
        var curvature = 1.0
        var i = 0
        while i < rows(group).length do
          val observation = rows(group)(i)
          val zTheta = term.z(0, observation) * theta
          curvature += response.workingWeights(observation) * link.muEta(base.eta(observation)) * zTheta * zTheta
          i += 1
        if !curvature.isFinite || curvature < 1.0 then
          return Left(MixedModelError.Optimization("invalid adaptive quadrature curvature"))
        val scale = 1.0 / math.sqrt(curvature)
        val logTerms = new Array[Double](rule.nodes.length)
        var node = 0
        while node < rule.nodes.length do
          val z = rule.nodes(node)
          val delta = z * scale
          var change = delta * (2.0 * mode + delta)
          i = 0
          while i < rows(group).length do
            val observation = rows(group)(i)
            val eta = base.eta(observation)
            val shifted = eta + term.z(0, observation) * theta * delta
            // Response-only constants cancel exactly in the node difference.
            change += response.parameterDeviance(observation, shifted) - response.parameterDeviance(observation, eta)
            i += 1
          logTerms(node) = rule.logWeights(node) - 0.5 * change + 0.5 * z * z
          if logTerms(node).isNaN || logTerms(node) == Double.PositiveInfinity then
            return Left(MixedModelError.Optimization("nonfinite adaptive quadrature node"))
          node += 1
        val maximum = logTerms.max
        if !maximum.isFinite then return Left(MixedModelError.Optimization("all adaptive quadrature nodes vanished"))
        totalLogCorrection += maximum + math.log(logTerms.iterator.map(value => math.exp(value - maximum)).sum)
        group += 1
      val correction = -2.0 * totalLogCorrection
      val objective = base.objective + correction
      if !correction.isFinite || !objective.isFinite then
        Left(MixedModelError.Optimization("nonfinite adaptive quadrature objective"))
      else Right(base.copy(objective = objective, quadratureCorrection = correction))

private[glmm] object AdaptiveQuadratureWorkspace:
  def create(design: CompiledDesign, link: Link, options: GlmmOptions): FitResult[Option[AdaptiveQuadrature]] =
    AdaptiveQuadrature
      .validate(design, options)
      .flatMap: _ =>
        if options.nAgq == 1 then Right(None)
        else AdaptiveQuadrature.rule(options.nAgq).map(rule => Some(new AdaptiveQuadrature(design, link, rule)))
