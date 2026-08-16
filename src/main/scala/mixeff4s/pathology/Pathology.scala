package mixeff4s.pathology

import mixeff4s.design.{CompiledDesign, ReMat}
import mixeff4s.error.{LinAlgError, MixedModelError}

/** Design-time pathology front door. */
object Pathology:
  val ContractVersion = "v0.1"
  val BoundaryTol = 1e-8

  /** Classify a fitted θ. Design-time expected statuses are left unchanged. */
  def assessFit(
      cert: Certificate,
      theta: Vector[Double],
      parmap: Vector[(Int, Int, Int)]
  ): Certificate =
    if cert.stratum == Stratum.Refusal then
      cert.copy(notes = cert.notes :+ "fit status is not claimed for a refusal design")
    else if theta.length != cert.nTheta || theta.length != parmap.length then
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
      val onBound = parmap.zip(theta).exists:
        case ((_, row, col), value) =>
          row == col && math.abs(value) <= BoundaryTol
      val status =
        if onBound then FitStatus.ConvergedBoundary else FitStatus.ConvergedInterior
      cert.copy(fitStatus = status)

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
          MixedModelError.PosDefException | MixedModelError.ConstantResponse |
          MixedModelError.NoRandomEffects =>
        FitStatus.NotIdentifiable
      case MixedModelError.LinAlg(LinAlgError.RankDeficient(_, _) | LinAlgError.Singular |
            LinAlgError.NotPositiveDefinite) =>
        FitStatus.NotIdentifiable
      case _ =>
        FitStatus.NotOptimized

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
