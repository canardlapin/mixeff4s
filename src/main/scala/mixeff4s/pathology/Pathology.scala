package mixeff4s.pathology

import mixeff4s.design.{CompiledDesign, ReMat}
import mixeff4s.error.{FitResult, LinAlgError, MixedModelError}
import mixeff4s.model.Family

/** Design-time pathology front door. */
object Pathology:
  val ContractVersion = "v0.1"
  val CorpusContractVersion = "v0.6"
  val WeakIdThreshold = Fisher.WeakIdThreshold
  val BoundaryTol = 1e-8
  private val ZeroVarianceTol = 1e-10
  private val UnitCorrelationTol = 1e-6
  private val RankRelTol = 1e-12

  def generate(spec: GeneratorSpec): FitResult[Generated] = Generate(spec)

  def detectSeparation(spec: GeneratorSpec): SeparationReport = Separation.detect(spec)

  def certify(spec: GeneratorSpec): Certificate =
    val q = spec.reDim
    val cov = spec.reCovTruth
    val shapeProblem =
      if cov.length != q || cov.exists(_.length != q) then
        Some(s"re_cov_truth is ${cov.length}×${cov.headOption.map(_.length).getOrElse(0)} but re_dim = $q")
      else None
    val eigvals = shapeProblem match
      case Some(_) => Right(Vector.empty)
      case None    => SymmetricPsd.eigvals(cov)
    val feCorrEigs =
      if spec.nFePredictors == 0 then Vector.empty
      else SymmetricPsd.eigvals(spec.predictorCorr).getOrElse(Vector.empty)
    val feRankTruth = effectiveRank(feCorrEigs)
    val issue = shapeProblem
      .orElse(eigvals.left.toOption.map(_.message))
      .map(StructuralIssue.MalformedSpec(_))
      .orElse:
        if spec.nFePredictors > 0 && feRankTruth < spec.nFePredictors then
          Some(StructuralIssue.CollinearFixedEffects(feRankTruth, spec.nFePredictors))
        else None
      .orElse:
        spec.crossedSummary
          .filter(_.nComponents > 1)
          .map(summary => StructuralIssue.DisconnectedCrossings(summary.nComponents))
      .orElse:
        if spec.family == Family.Bernoulli then
          Separation.detect(spec).kind.map(StructuralIssue.Separation(_))
        else None
      .orElse:
        if spec.nReSlopes > 0 && spec.minGroupSize == 1 then
          Some(StructuralIssue.SingletonsWithSlope(spec.label, spec.minGroupSize))
        else None
    val vals = eigvals.getOrElse(Vector.empty)
    val rankTruth = effectiveRank(vals)
    val boundaries = if issue.exists(_.code == "malformed_spec") then Vector.empty else boundaryDirections(cov)
    val fisher = Fisher.correlationEigvals(spec)
    val weakScore = Fisher.weakIdScore(spec.n, fisher)
    val weakId = weakScore.isFinite && weakScore < WeakIdThreshold
    val (stratum, expected) = expectedFromTruth(issue, rankTruth, q, boundaries, weakId)
    Certificate(
      CorpusContractVersion,
      FitStatus.NotAssessed,
      expected,
      stratum,
      spec.n,
      spec.nParams,
      spec.minGroupSize,
      spec.maxGroupSize,
      spec.feRank,
      spec.nTheta,
      issue,
      Vector(spec.label),
      rankTruth,
      q,
      boundaries,
      fisher,
      weakScore,
      WeakIdThreshold,
      weakId,
      spec.crossedSummary
    )

  /** Classify a fitted θ. Design-time expected statuses are left unchanged. */
  def assessFit(
      cert: Certificate,
      theta: Vector[Double],
      parmap: Vector[(Int, Int, Int)]
  ): Certificate =
    if cert.stratum == Stratum.Refusal then
      cert.copy(notes = cert.notes :+ "fit status is not claimed for a refusal design")
    else classifyTheta(cert, theta, parmap)

  /** Map an engine outcome. Unlike `assessFit`, a refusal still receives a status. */
  def assessOutcome(
      cert: Certificate,
      outcome: FitResult[Vector[Double]],
      parmap: Vector[(Int, Int, Int)]
  ): Certificate =
    outcome match
      case Left(err)    => cert.copy(fitStatus = mapError(err), notes = cert.notes :+ err.code)
      case Right(theta) => classifyTheta(cert, theta, parmap)

  def certify(design: CompiledDesign): Certificate =
    val sizes = design.reterms.flatMap(groupSizes)
    val minG = if sizes.isEmpty then 0 else sizes.min
    val maxG = if sizes.isEmpty then 0 else sizes.max
    val nParams = design.p + design.nTheta + 1
    val issue = structuralIssue(design, nParams)
    val (stratum, expected) = issue match
      case Some(_) =>
        (Stratum.Refusal, Vector(FitStatus.NotIdentifiable, FitStatus.NotOptimized))
      case None =>
        (Stratum.Easy, Vector(FitStatus.ConvergedInterior))
    Certificate(
      ContractVersion,
      FitStatus.NotAssessed,
      expected,
      stratum,
      design.n,
      nParams,
      minG,
      maxG,
      design.p,
      design.nTheta,
      issue,
      Vector.empty
    )

  def fromError(error: MixedModelError): Certificate =
    val status = mapError(error)
    Certificate(
      ContractVersion,
      status,
      Vector(status),
      Stratum.Refusal,
      0,
      0,
      0,
      0,
      0,
      0,
      None,
      Vector(error.code, error.message)
    )

  def mapError(error: MixedModelError): FitStatus =
    error match
      case MixedModelError.Singular(_) | MixedModelError.RankSaturatedFixedEffects(_, _) |
          MixedModelError.PosDefException | MixedModelError.ConstantResponse | MixedModelError.NoRandomEffects =>
        FitStatus.NotIdentifiable
      case MixedModelError.LinAlg(
            LinAlgError.RankDeficient(_, _) | LinAlgError.Singular | LinAlgError.NotPositiveDefinite
          ) =>
        FitStatus.NotIdentifiable
      case _ =>
        FitStatus.NotOptimized

  private def classifyTheta(
      cert: Certificate,
      theta: Vector[Double],
      parmap: Vector[(Int, Int, Int)]
  ): Certificate =
    if theta.length != cert.nTheta || theta.length != parmap.length then
      cert.copy(
        fitStatus = FitStatus.NotAssessed,
        notes = cert.notes :+ s"theta length ${theta.length} does not match n_theta ${cert.nTheta}"
      )
    else if theta.exists(v => !v.isFinite) then
      cert.copy(
        fitStatus = FitStatus.NotAssessed,
        notes = cert.notes :+ "theta contains a non-finite value"
      )
    else
      val onBound = parmap
        .zip(theta)
        .exists:
          case ((_, row, col), value) =>
            row == col && math.abs(value) <= BoundaryTol
      val status =
        if onBound then FitStatus.ConvergedBoundary else FitStatus.ConvergedInterior
      cert.copy(fitStatus = status)

  private def structuralIssue(design: CompiledDesign, nParams: Int): Option[StructuralIssue] =
    if design.p >= design.n then Some(StructuralIssue.RankSaturated(design.p, design.n))
    else
      design.reterms
        .collectFirst:
          case rt if rt.nLevels < 2 =>
            StructuralIssue.FewLevels(rt.groupingName, rt.nLevels)
        .orElse:
          design.reterms.collectFirst:
            case rt if rt.vsize > 1 && groupSizes(rt).exists(_ == 1) =>
              StructuralIssue.SingletonsWithSlope(rt.groupingName, groupSizes(rt).min)
        .orElse:
          if nParams >= design.n then Some(StructuralIssue.InformationSaturated(nParams, design.n))
          else None

  private def groupSizes(rt: ReMat): Vector[Int] =
    rt.levels.indices.map(level => rt.refs.count(_ == level)).toVector

  private def expectedFromTruth(
      issue: Option[StructuralIssue],
      rankTruth: Int,
      rankRequested: Int,
      boundaries: Vector[BoundaryKind],
      weakIdentification: Boolean
  ): (Stratum, Vector[FitStatus]) =
    if issue.exists(_.code == "separation") then
      (
        Stratum.Refusal,
        Vector(FitStatus.NotIdentifiable, FitStatus.NotOptimized, FitStatus.ConvergedPenalised)
      )
    else if issue.exists(_.code == "disconnected_crossings") then
      (
        Stratum.Refusal,
        Vector(
          FitStatus.NotIdentifiable,
          FitStatus.NotOptimized,
          FitStatus.ConvergedBoundary,
          FitStatus.ConvergedReducedRank,
          FitStatus.ConvergedInterior
        )
      )
    else if issue.isDefined then
      (
        Stratum.Refusal,
        Vector(
          FitStatus.NotIdentifiable,
          FitStatus.NotOptimized,
          FitStatus.ConvergedBoundary
        )
      )
    else if rankTruth < rankRequested then
      (
        Stratum.ReducedRank,
        Vector(FitStatus.ConvergedReducedRank, FitStatus.ConvergedBoundary)
      )
    else if boundaries.nonEmpty then
      (Stratum.Boundary, Vector(FitStatus.ConvergedBoundary, FitStatus.ConvergedInterior))
    else if weakIdentification then
      (
        Stratum.Easy,
        Vector(FitStatus.ConvergedInterior, FitStatus.ConvergedBoundary, FitStatus.ConvergedReducedRank)
      )
    else (Stratum.Easy, Vector(FitStatus.ConvergedInterior))

  private def boundaryDirections(cov: Vector[Vector[Double]]): Vector[BoundaryKind] =
    val q = cov.length
    val zeros = (0 until q).collect:
      case i if math.abs(cov(i)(i)) <= ZeroVarianceTol =>
        BoundaryKind.ZeroVariance(i)
    val corrs = (0 until q).flatMap: i =>
      ((i + 1) until q).collect:
        case j if unitCorrelation(cov(i)(i), cov(j)(j), cov(i)(j)) =>
          BoundaryKind.UnitCorrelation(i, j)
    zeros.toVector ++ corrs.toVector

  private def unitCorrelation(a: Double, b: Double, c: Double): Boolean =
    if a <= ZeroVarianceTol || b <= ZeroVarianceTol then false
    else
      val rho = c / math.sqrt(a * b)
      math.abs(math.abs(rho) - 1.0) <= UnitCorrelationTol

  private def effectiveRank(eigvals: Vector[Double]): Int =
    if eigvals.isEmpty then 0
    else
      val trace = eigvals.map(math.abs).sum
      val cutoff = RankRelTol * math.max(trace, 1e-15)
      eigvals.count(_ > cutoff)

