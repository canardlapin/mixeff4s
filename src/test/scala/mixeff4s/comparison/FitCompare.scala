package mixeff4s.comparison

/** Compare a mixeff4s numeric fit to a frozen reference. */
object FitCompare:
  val DefaultObjectiveTol = 1e-2
  val DefaultBetaTol = 1e-3
  val DefaultThetaTol = 1e-3
  val DefaultSigmaTol = 1e-3

  def row(dataset: String, estimator: String, formula: String): ScorecardRow =
    Scorecard.loadEmbedded.rows
      .find: r =>
        r.key.dataset == dataset && r.key.estimator == estimator && r.key.formula == formula
      .getOrElse(throw IllegalArgumentException(s"missing scorecard row $dataset $estimator $formula"))

  def row(dataset: String, estimator: String): ScorecardRow =
    val matches = Scorecard.loadEmbedded.rows.filter: r =>
      r.key.dataset == dataset && r.key.estimator == estimator
    matches match
      case Vector(single) => single
      case Vector() =>
        throw IllegalArgumentException(s"missing scorecard row $dataset $estimator")
      case many =>
        throw IllegalArgumentException(s"ambiguous scorecard row $dataset $estimator: ${many.map(_.key.formula)}")

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
    else raw"([A-Za-z_][\w.]*): ([^:]+)".r.replaceAllIn(name, m => m.group(1) + m.group(2).trim)

  /** Scalar multi-term θ is compared sorted so nRanef reordering is not a mismatch. */
  def alignedTheta(
      obtained: Vector[Double],
      reference: Vector[Double],
      scalarTerms: Boolean
  ): Vector[(Double, Double)] =
    if obtained.length != reference.length then
      throw IllegalArgumentException(s"theta length ${obtained.length} != ${reference.length}")
    val (got, expected) =
      if scalarTerms then (obtained.sorted, reference.sorted) else (obtained, reference)
    got.zip(expected)

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
