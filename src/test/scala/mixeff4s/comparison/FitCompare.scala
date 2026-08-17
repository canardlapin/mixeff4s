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

  /** Map `recipe: B:temperature: 185` onto lme4's `recipeB:temperature185`. */
  def contrastName(name: String): String =
    if name == "(Intercept)" then name
    else
      val pieces = raw"([A-Za-z_][\w.]*): ([^:]+)".r.findAllMatchIn(name).toVector
      if pieces.isEmpty then name
      else pieces.map(m => m.group(1) + m.group(2).trim).mkString(":")

  def alignedBeta(
      obtainedNames: Vector[String],
      obtained: Vector[Double],
      referenceNames: Vector[String],
      reference: Vector[Double]
  ): Vector[(String, Double, Double)] =
    if obtained.length != obtainedNames.length then
      throw IllegalArgumentException("obtained beta and names differ in length")
    if reference.length != referenceNames.length then
      throw IllegalArgumentException("reference beta and names differ in length")
    val byName = obtainedNames.map(contrastName).zip(obtained).toMap
    val missing = referenceNames.filterNot(byName.contains)
    if missing.nonEmpty then
      throw IllegalArgumentException(
        s"fit is missing reference contrasts $missing; have ${obtainedNames.map(contrastName)}"
      )
    referenceNames.zip(reference).map: (name, value) =>
      (name, byName(name), value)
