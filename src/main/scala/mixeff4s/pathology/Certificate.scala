package mixeff4s.pathology

/** Engine-free identifiability certificate. Computed from the design or a refusal, never from optimizer folklore. */
final case class Certificate(
    contractVersion: String,
    fitStatus: FitStatus,
    expectedStatuses: Vector[FitStatus],
    stratum: Stratum,
    n: Int,
    nParams: Int,
    minGroupSize: Int,
    maxGroupSize: Int,
    feRank: Int,
    nTheta: Int,
    structuralIssue: Option[StructuralIssue],
    notes: Vector[String]
)
