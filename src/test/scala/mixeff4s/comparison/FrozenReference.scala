package mixeff4s.comparison

enum ReferenceEngine:
  case Lme4, MixedModelsJl, MixeffRs, Mixeff4s

  def id: String =
    this match
      case ReferenceEngine.Lme4          => "lme4"
      case ReferenceEngine.MixedModelsJl => "mixedmodels.jl"
      case ReferenceEngine.MixeffRs      => "mixeff-rs"
      case ReferenceEngine.Mixeff4s      => "mixeff4s"

object ReferenceEngine:
  def parse(id: String): ReferenceEngine =
    id match
      case "lme4"           => ReferenceEngine.Lme4
      case "mixedmodels.jl" => ReferenceEngine.MixedModelsJl
      case "mixeff-rs"      => ReferenceEngine.MixeffRs
      case "mixeff4s"       => ReferenceEngine.Mixeff4s
      case other            => throw IllegalArgumentException(s"unknown engine `$other`")

final case class Tolerances(
    objective: Option[Double],
    beta: Option[Double],
    theta: Option[Double],
    sigma: Option[Double]
)

final case class FrozenResult(
    key: FitKey,
    engine: ReferenceEngine,
    status: String,
    nObs: Option[Int],
    objective: Option[Double],
    objectiveDefinition: Option[String],
    objectiveComparable: Boolean,
    beta: Option[Vector[Double]],
    coefNames: Option[Vector[String]],
    theta: Option[Vector[Double]],
    sigma: Option[Double],
    isSingular: Option[Boolean],
    tolerances: Tolerances
)

final case class FrozenCatalog(
    mixeffRsRevision: Option[String],
    lme4Tool: Option[String],
    note: String,
    results: Vector[FrozenResult]
):
  def forKey(key: FitKey): Vector[FrozenResult] =
    results.filter(_.key == key)

  def claimed(row: ScorecardRow): Option[FrozenResult] =
    val engineId = row.reference match
      case "mixedmodels.jl_fast_pirls" => "mixedmodels.jl"
      case "lme4_boundary" | "lme4_joint_laplace" | "lme4_joint_agq" |
          "lme4_numeric_without_objective_constants" =>
        "lme4"
      case other => other
    forKey(row.key).find(_.engine.id == engineId)

object FrozenCatalog:
  def loadEmbedded: FrozenCatalog = parse(EmbeddedFrozen.json)

  def parse(text: String): FrozenCatalog =
    val root = JsonValue.parse(text)
    val schema = root.req("schema")
    if schema.req("name").asString != "mixeff4s.frozen_reference" then
      throw IllegalArgumentException("frozen catalog schema name mismatch")
    if schema.req("version").asInt != 1 then throw IllegalArgumentException("frozen catalog schema version mismatch")
    val source = root.req("source")
    val results = root.req("results").asArray.map(parseResult)
    val forbidden = results.exists: result =>
      result.key.estimator.toLowerCase.contains("pvalue")
    if forbidden then throw IllegalArgumentException("p-values are not admitted")
    FrozenCatalog(
      mixeffRsRevision = source.field("mixeff_rs_revision").map(_.asString),
      lme4Tool = source.field("lme4_tool").map(_.asString),
      note = source.req("note").asString,
      results = results
    )

  private def parseResult(value: JsonValue): FrozenResult =
    if value.field("p_value").isDefined || value.field("p_values").isDefined || value.field("pvalue").isDefined then
      throw IllegalArgumentException("frozen references must not contain p-values")
    val tols = value.field("tolerances")
    FrozenResult(
      key = FitKey(
        dataset = value.req("dataset").asString,
        formula = value.req("formula").asString,
        family = value.req("family").asString,
        link = value.req("link").asString,
        estimator = value.req("estimator").asString
      ),
      engine = ReferenceEngine.parse(value.req("engine").asString),
      status = value.req("status").asString,
      nObs = value.field("n_obs").map(_.asInt),
      objective = optionalNumber(value, "objective"),
      objectiveDefinition = value.field("objective_definition").map(_.asString),
      objectiveComparable = value.field("objective_comparable").map(_.asBoolean).getOrElse(true),
      beta = optionalNumbers(value, "beta"),
      coefNames = value.field("coef_names").map(_.asArray.map(_.asString)),
      theta = optionalNumbers(value, "theta"),
      sigma = optionalNumber(value, "sigma"),
      isSingular = value.field("is_singular").map(_.asBoolean),
      tolerances = Tolerances(
        objective = tols.flatMap(_.field("objective").map(_.asDouble)),
        beta = tols.flatMap(_.field("beta").map(_.asDouble)),
        theta = tols.flatMap(_.field("theta").map(_.asDouble)),
        sigma = tols.flatMap(_.field("sigma").map(_.asDouble))
      )
    )

  private def optionalNumber(value: JsonValue, name: String): Option[Double] =
    value
      .field(name)
      .flatMap:
        case JsonValue.Null => None
        case other          => Some(other.asDouble)

  private def optionalNumbers(value: JsonValue, name: String): Option[Vector[Double]] =
    value
      .field(name)
      .flatMap:
        case JsonValue.Null => None
        case other          => Some(other.asArray.map(_.asDouble))
