package mixeff4s.comparison

/** Compare a mixeff4s numeric fit to a frozen reference. */
object FitCompare:
  val DefaultObjectiveTol = 1e-2
  val DefaultBetaTol = 1e-3
  val DefaultThetaTol = 1e-3
  val DefaultSigmaTol = 1e-3

  def row(dataset: String, estimator: String): ScorecardRow =
    Scorecard.loadEmbedded.rows
      .find(r => r.key.dataset == dataset && r.key.estimator == estimator)
      .getOrElse(throw IllegalArgumentException(s"missing scorecard row $dataset $estimator"))

  def claimed(row: ScorecardRow): FrozenResult =
    FrozenCatalog
      .loadEmbedded
      .claimed(row)
      .getOrElse(throw IllegalArgumentException(s"missing frozen claim for ${row.key.id}"))

  def objectiveTol(row: ScorecardRow, ref: FrozenResult): Double =
    row.objectiveAbsTol.orElse(ref.tolerances.objective).getOrElse(DefaultObjectiveTol)

  def betaTol(row: ScorecardRow, ref: FrozenResult): Double =
    row.betaAbsTol.orElse(ref.tolerances.beta).getOrElse(DefaultBetaTol)

  def thetaTol(row: ScorecardRow, ref: FrozenResult): Double =
    row.thetaAbsTol.orElse(ref.tolerances.theta).getOrElse(DefaultThetaTol)

  def sigmaTol(row: ScorecardRow, ref: FrozenResult): Double =
    row.sigmaAbsTol.orElse(ref.tolerances.sigma).getOrElse(DefaultSigmaTol)
