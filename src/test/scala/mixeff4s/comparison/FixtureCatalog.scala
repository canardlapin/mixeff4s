package mixeff4s.comparison

final case class FixtureCatalog(
    schemaVersion: String,
    builtin: Vector[String],
    vendored: Vector[String],
    mixeffRsRevision: String
):
  def available: Vector[String] = builtin ++ vendored

object FixtureCatalog:
  val embeddedToml: String =
    """schema_version = "1.0.0"
      |
      |builtin = [
      |  "sleepstudy",
      |  "penicillin",
      |  "pastes",
      |  "contraception",
      |]
      |
      |vendored = [
      |  "dyestuff",
      |  "dyestuff2",
      |  "rail",
      |  "cake",
      |  "ergostool",
      |  "oats",
      |  "orthodont",
      |  "oxide",
      |  "machines",
      |  "station_season_duration",
      |  "singular",
      |]
      |
      |mixeff_rs_revision = "d6b81ce8b5974b2f0e42fcf650ad6811e5202a4f"
      |""".stripMargin

  def loadEmbedded: FixtureCatalog = parse(embeddedToml)

  def parse(text: String): FixtureCatalog =
    val doc = TomlTables.parse(text)
    FixtureCatalog(
      schemaVersion = doc.root("schema_version").asString,
      builtin = doc.root("builtin").asStrings,
      vendored = doc.root("vendored").asStrings,
      mixeffRsRevision = doc.root("mixeff_rs_revision").asString
    )
