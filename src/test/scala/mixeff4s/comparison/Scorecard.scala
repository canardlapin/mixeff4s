package mixeff4s.comparison

enum ScorecardClass:
  case ReleaseBlockingParity, DocumentedDivergence, UnsupportedWithContract, StressOptIn, PerformanceKnownSlow

  def id: String =
    this match
      case ScorecardClass.ReleaseBlockingParity   => "release_blocking_parity"
      case ScorecardClass.DocumentedDivergence    => "documented_divergence"
      case ScorecardClass.UnsupportedWithContract => "unsupported_with_contract"
      case ScorecardClass.StressOptIn             => "stress_opt_in"
      case ScorecardClass.PerformanceKnownSlow    => "performance_known_slow"

object ScorecardClass:
  val all: Vector[ScorecardClass] = ScorecardClass.values.toVector

  def parse(id: String): ScorecardClass =
    all
      .find(_.id == id)
      .getOrElse(throw IllegalArgumentException(s"unknown scorecard class `$id`"))

final case class FitKey(
    dataset: String,
    formula: String,
    family: String,
    link: String,
    estimator: String
):
  def id: String = s"$dataset | $estimator | $formula"

final case class ScorecardRow(
    key: FitKey,
    classification: ScorecardClass,
    reference: String,
    reason: Option[String],
    objectiveAbsTol: Option[Double],
    betaAbsTol: Option[Double],
    thetaAbsTol: Option[Double],
    sigmaAbsTol: Option[Double]
)

final case class Scorecard(
    schemaVersion: String,
    classes: Vector[ScorecardClass],
    rows: Vector[ScorecardRow]
):
  def keys: Vector[FitKey] = rows.map(_.key)

object Scorecard:
  val embeddedToml: String =
    """# Machine-readable parity classification for every checked-in
      |# dataset/formula/estimator triple. Keep this in sync with
      |# comparison/fixtures.toml and comparison/frozen/references.json.
      |#
      |# This is not AGENTS.md Phase 8. Classes are rust vocabulary.
      |# Do not mark a row release_blocking_parity against lme4 until a
      |# frozen lme4 number exists and a mixeff4s fit is compared to it.
      |
      |schema_version = "1.0.0"
      |
      |classes = [
      |  "release_blocking_parity",
      |  "documented_divergence",
      |  "unsupported_with_contract",
      |  "stress_opt_in",
      |  "performance_known_slow",
      |]
      |
      |[[row]]
      |dataset = "sleepstudy"
      |formula = "reaction ~ 1 + days + (1 + days | subj)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "ML"
      |class = "release_blocking_parity"
      |reference = "mixedmodels.jl"
      |reason = "pinned in SleepstudyFitSuite against MixedModels.jl; lme4 is not claimed until a later ticket compares the frozen lme4 row"
      |
      |[[row]]
      |dataset = "sleepstudy"
      |formula = "reaction ~ 1 + days + (1 + days | subj)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "REML"
      |class = "release_blocking_parity"
      |reference = "mixedmodels.jl"
      |reason = "pinned in SleepstudyFitSuite against MixedModels.jl; lme4 is not claimed until a later ticket compares the frozen lme4 row"
      |
      |[[row]]
      |dataset = "penicillin"
      |formula = "diameter ~ 1 + (1 | plate) + (1 | sample)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "ML"
      |class = "release_blocking_parity"
      |reference = "mixedmodels.jl"
      |reason = "pinned in CrossedFitSuite against MixedModels.jl"
      |
      |[[row]]
      |dataset = "penicillin"
      |formula = "diameter ~ 1 + (1 | plate) + (1 | sample)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "REML"
      |class = "release_blocking_parity"
      |reference = "mixeff-rs"
      |reason = "only rust-labelled numeric pin in CrossedFitSuite; frozen lme4 numbers exist but are not yet a mixeff4s claim"
      |
      |[[row]]
      |dataset = "pastes"
      |formula = "strength ~ 1 + (1 | batch / cask)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "ML"
      |class = "release_blocking_parity"
      |reference = "mixedmodels.jl"
      |theta_abs_tol = 0.09
      |reason = "pinned in CrossedFitSuite against MixedModels.jl with a slack theta tolerance"
      |
      |[[row]]
      |dataset = "dyestuff"
      |formula = "Yield ~ 1 + (1 | Batch)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "ML"
      |class = "release_blocking_parity"
      |reference = "lme4"
      |reason = "scalar random-intercept floor; mixeff4s fit compared to frozen lme4"
      |
      |[[row]]
      |dataset = "dyestuff"
      |formula = "Yield ~ 1 + (1 | Batch)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "REML"
      |class = "release_blocking_parity"
      |reference = "lme4"
      |reason = "scalar random-intercept floor; mixeff4s fit compared to frozen lme4"
      |
      |[[row]]
      |dataset = "dyestuff2"
      |formula = "Yield ~ 1 + (1 | Batch)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "ML"
      |class = "release_blocking_parity"
      |reference = "lme4_boundary"
      |reason = "singular boundary is expected; assessFit must report converged_boundary"
      |
      |[[row]]
      |dataset = "dyestuff2"
      |formula = "Yield ~ 1 + (1 | Batch)"
      |family = "Gaussian"
      |link = "Identity"
      |estimator = "REML"
      |class = "release_blocking_parity"
      |reference = "lme4_boundary"
      |reason = "singular boundary is expected; assessFit must report converged_boundary"
      |
      |[[row]]
      |dataset = "contraception"
      |formula = "use_num ~ 1 + age + age2 + urban + livch + (1 | urban_dist)"
      |family = "Bernoulli"
      |link = "Logit"
      |estimator = "fast_pirls"
      |class = "documented_divergence"
      |reference = "mixedmodels.jl_fast_pirls"
      |objective_abs_tol = 1.0
      |reason = "labelled fast-PIRLS; not lme4::glmer. Objective constants are comparable only to MixedModels.jl fast=true"
      |""".stripMargin

  def loadEmbedded: Scorecard = parse(embeddedToml)

  def parse(text: String): Scorecard =
    val doc = TomlTables.parse(text)
    val schemaVersion = doc.root("schema_version").asString
    val classes = doc.root("classes").asStrings.map(ScorecardClass.parse)
    if classes != ScorecardClass.all then
      throw IllegalArgumentException(s"scorecard classes must be ${ScorecardClass.all.map(_.id)}")
    val rows = doc.rows.map(parseRow)
    val dupes = rows.groupBy(_.key).collect { case (key, group) if group.length > 1 => key.id }
    if dupes.nonEmpty then throw IllegalArgumentException(s"duplicate scorecard keys: $dupes")
    Scorecard(schemaVersion, classes, rows)

  private def parseRow(fields: Map[String, TomlTables.Value]): ScorecardRow =
    ScorecardRow(
      key = FitKey(
        dataset = fields("dataset").asString,
        formula = fields("formula").asString,
        family = fields("family").asString,
        link = fields("link").asString,
        estimator = fields("estimator").asString
      ),
      classification = ScorecardClass.parse(fields("class").asString),
      reference = fields("reference").asString,
      reason = fields.get("reason").map(_.asString),
      objectiveAbsTol = fields.get("objective_abs_tol").map(_.asDouble),
      betaAbsTol = fields.get("beta_abs_tol").map(_.asDouble),
      thetaAbsTol = fields.get("theta_abs_tol").map(_.asDouble),
      sigmaAbsTol = fields.get("sigma_abs_tol").map(_.asDouble)
    )
