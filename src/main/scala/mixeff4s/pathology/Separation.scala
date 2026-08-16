package mixeff4s.pathology

import mixeff4s.data.ModelFrame
import mixeff4s.model.Family

/** Two-tier binomial separation report. The MLE does not exist when either tier fires. */
final case class SeparationReport(
    feKind: Option[FeSeparationKind],
    conditionalGroups: Vector[Int],
    hyperplaneDirection: Option[Vector[Double]]
):
  def isSeparated: Boolean = feKind.isDefined || conditionalGroups.nonEmpty
  def nConditionalGroups: Int = conditionalGroups.length
  def kind: Option[SeparationKind] =
    (feKind, conditionalGroups.isEmpty) match
      case (Some(fe), true)  => Some(SeparationKind.FixedEffect(fe))
      case (None, false)     => Some(SeparationKind.Conditional(conditionalGroups.length))
      case (Some(fe), false) => Some(SeparationKind.Both(fe, conditionalGroups.length))
      case (None, true)      => None

object SeparationReport:
  def empty: SeparationReport = SeparationReport(None, Vector.empty, None)

/** Konis (2007) FE trichotomy for p ≤ 2, plus a per-group all-zero / all-one scan. */
object Separation:
  val MarginTol = 1e-6

  def detect(spec: GeneratorSpec): SeparationReport =
    if spec.family != Family.Bernoulli then SeparationReport.empty
    else
      Generate(spec) match
        case Left(_)          => SeparationReport.empty
        case Right(generated) => detect(spec, generated.frame)

  def detect(spec: GeneratorSpec, frame: ModelFrame): SeparationReport =
    frame.numeric(spec.responseName) match
      case None => SeparationReport.empty
      case Some(y) if y.isEmpty => SeparationReport.empty
      case Some(y) =>
        val x = designMatrix(spec, frame, y.length)
        val (feKind, beta) = detectFe(x, y)
        val groups = frame.factor(spec.groupName).map(_.refs).getOrElse(Vector.empty)
        SeparationReport(feKind, detectConditional(y, groups), beta)

  /** Rows of `x` are observations; column 0 is the intercept when the spec has one. */
  def detectFe(x: Vector[Vector[Double]], y: Vector[Double]): (Option[FeSeparationKind], Option[Vector[Double]]) =
    val n = x.length
    val p = x.headOption.map(_.length).getOrElse(0)
    if n == 0 || p == 0 || y.length != n || p > 2 then (None, None)
    else
      val z = Vector.tabulate(n, p): (i, j) =>
        val sign = if y(i) > 0.5 then 1.0 else -1.0
        sign * x(i)(j)
      val (epsStar, betaA) = maximinBox(z)
      if epsStar > MarginTol then (Some(FeSeparationKind.Complete), Some(betaA))
      else
        val (obj, betaB) = residualCone(z)
        if obj > MarginTol then (Some(FeSeparationKind.QuasiComplete), Some(betaB))
        else (None, None)

  def detectConditional(y: Vector[Double], groups: Vector[Int]): Vector[Int] =
    if y.isEmpty || groups.length != y.length then Vector.empty
    else
      val nGroups = groups.maxOption.map(_ + 1).getOrElse(0)
      if nGroups == 0 then Vector.empty
      else
        val sums = Array.fill(nGroups)(0)
        val counts = Array.fill(nGroups)(0)
        groups.zipWithIndex.foreach: (g, i) =>
          if g >= 0 && g < nGroups then
            counts(g) += 1
            if y(i) > 0.5 then sums(g) += 1
        (0 until nGroups).filter(g => counts(g) > 0 && (sums(g) == 0 || sums(g) == counts(g))).toVector

  private def designMatrix(spec: GeneratorSpec, frame: ModelFrame, n: Int): Vector[Vector[Double]] =
    val intercept = if spec.hasIntercept then 1 else 0
    Vector.tabulate(n, intercept + spec.nFePredictors): (i, j) =>
      if spec.hasIntercept && j == 0 then 1.0
      else
        val name = s"x${j + 1 - intercept}"
        frame.numeric(name).map(_(i)).getOrElse(0.0)

  private def maximinBox(z: Vector[Vector[Double]]): (Double, Vector[Double]) =
    z.headOption.map(_.length).getOrElse(0) match
      case 1 =>
        val col = z.map(_(0))
        val atPos = col.min
        val atNeg = -col.max
        if atPos >= atNeg then (atPos, Vector(1.0)) else (atNeg, Vector(-1.0))
      case 2 =>
        val corners = Vector(Vector(-1.0, -1.0), Vector(-1.0, 1.0), Vector(1.0, -1.0), Vector(1.0, 1.0))
        val edges = Vector(
          (corners(0), corners(1)),
          (corners(1), corners(3)),
          (corners(3), corners(2)),
          (corners(2), corners(0))
        )
        edges.map((a, b) => maximinSegment(z, a, b)).maxBy(_._1)
      case _ =>
        (0.0, Vector.fill(z.headOption.map(_.length).getOrElse(0))(0.0))

  private def maximinSegment(
      z: Vector[Vector[Double]],
      a: Vector[Double],
      b: Vector[Double]
  ): (Double, Vector[Double]) =
    val d = b.zip(a).map(_ - _)
    val alpha = z.map(dot(_, a))
    val gamma = z.map(dot(_, d))
    var bestT = 0.0
    var bestF = alpha.min
    val atOne = alpha.zip(gamma).map(_ + _).min
    if atOne > bestF then
      bestF = atOne
      bestT = 1.0
    val n = z.length
    var i = 0
    while i < n do
      var j = i + 1
      while j < n do
        val denom = gamma(i) - gamma(j)
        if math.abs(denom) > 1e-15 then
          val t = (alpha(j) - alpha(i)) / denom
          if t > 0.0 && t < 1.0 then
            val f = alpha(i) + t * gamma(i)
            var k = 0
            var isMin = true
            while isMin && k < n do
              if alpha(k) + t * gamma(k) < f - 1e-12 then isMin = false
              k += 1
            if isMin && f > bestF then
              bestF = f
              bestT = t
        j += 1
      i += 1
    (bestF, a.zip(d).map((ai, di) => ai + bestT * di))

  private def residualCone(z: Vector[Vector[Double]]): (Double, Vector[Double]) =
    z.headOption.map(_.length).getOrElse(0) match
      case 1 =>
        val col = z.map(_(0))
        if col.forall(_ >= -MarginTol) && col.exists(_ > MarginTol) then (col.sum, Vector(1.0))
        else if col.forall(_ <= MarginTol) && col.exists(_ < -MarginTol) then (-col.sum, Vector(-1.0))
        else (0.0, Vector(0.0))
      case 2 =>
        val corners = Vector(Vector(-1.0, -1.0), Vector(-1.0, 1.0), Vector(1.0, -1.0), Vector(1.0, 1.0))
        val edgeHits = Vector(
          (corners(0), corners(1)),
          (corners(1), corners(3)),
          (corners(3), corners(2)),
          (corners(2), corners(0))
        ).flatMap: (a, b) =>
          edgeZeros(z, a, b)
        val candidates = corners ++ edgeHits
        val feasible = candidates.flatMap: beta =>
          val zb = z.map(dot(_, beta))
          if zb.forall(_ >= -MarginTol) then Some((zb.sum, beta)) else None
        if feasible.isEmpty then (0.0, Vector(0.0, 0.0))
        else feasible.maxBy(_._1)
      case _ =>
        (0.0, Vector.fill(z.headOption.map(_.length).getOrElse(0))(0.0))

  private def edgeZeros(
      z: Vector[Vector[Double]],
      a: Vector[Double],
      b: Vector[Double]
  ): Vector[Vector[Double]] =
    val d = b.zip(a).map(_ - _)
    z.flatMap: zi =>
      val denom = dot(zi, d)
      if math.abs(denom) <= 1e-15 then None
      else
        val t = -dot(zi, a) / denom
        if t > 1e-12 && t < 1.0 - 1e-12 then Some(a.zip(d).map((ai, di) => ai + t * di))
        else None

  private def dot(u: Vector[Double], v: Vector[Double]): Double =
    u.zip(v).map(_ * _).sum
