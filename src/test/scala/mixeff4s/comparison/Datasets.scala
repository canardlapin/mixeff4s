package mixeff4s.comparison

import mixeff4s.data.ModelFrame
import mixeff4s.fixtures.{Contra, Pastes, Penicillin, Sleepstudy}

enum ColumnKind:
  case Numeric, Categorical

final case class ColumnMeta(name: String, kind: ColumnKind, levels: Option[Vector[String]])

final case class DatasetMeta(
    name: String,
    nRows: Int,
    source: String,
    license: String,
    columns: Vector[ColumnMeta]
)

final case class Dataset(name: String, meta: DatasetMeta, frame: ModelFrame)

object Datasets:
  def names: Vector[String] = FixtureCatalog.loadEmbedded.available

  def load(name: String): Dataset =
    name match
      case "sleepstudy" =>
        builtin(
          name,
          nRows = 180,
          source = "lme4::sleepstudy / MixedModels.jl",
          license = "GPL (>= 2) — see lme4 DESCRIPTION",
          columns = Vector(
            ColumnMeta("reaction", ColumnKind.Numeric, None),
            ColumnMeta("days", ColumnKind.Numeric, None),
            ColumnMeta("subj", ColumnKind.Categorical, Some(Sleepstudy.subjects))
          ),
          Sleepstudy.frame
        )
      case "penicillin" =>
        builtin(
          name,
          nRows = 144,
          source = "lme4::Penicillin",
          license = "GPL (>= 2) — see lme4 DESCRIPTION",
          columns = Vector(
            ColumnMeta("diameter", ColumnKind.Numeric, None),
            ColumnMeta("plate", ColumnKind.Categorical, None),
            ColumnMeta("sample", ColumnKind.Categorical, None)
          ),
          Penicillin.frame
        )
      case "pastes" =>
        builtin(
          name,
          nRows = 60,
          source = "lme4::Pastes",
          license = "GPL (>= 2) — see lme4 DESCRIPTION",
          columns = Vector(
            ColumnMeta("strength", ColumnKind.Numeric, None),
            ColumnMeta("batch", ColumnKind.Categorical, None),
            ColumnMeta("cask", ColumnKind.Categorical, None)
          ),
          Pastes.frame
        )
      case "contraception" =>
        builtin(
          name,
          nRows = 1934,
          source = "MixedModels.jl contra",
          license = "Per MixedModels.jl test data",
          columns = Vector(
            ColumnMeta("use_num", ColumnKind.Numeric, None),
            ColumnMeta("age", ColumnKind.Numeric, None),
            ColumnMeta("age2", ColumnKind.Numeric, None),
            ColumnMeta("urban", ColumnKind.Categorical, None),
            ColumnMeta("livch", ColumnKind.Categorical, None),
            ColumnMeta("urban_dist", ColumnKind.Categorical, None)
          ),
          Contra.frame
        )
      case other =>
        val csv = EmbeddedFrames.csv.getOrElse(other, throw IllegalArgumentException(s"unknown dataset `$other`"))
        val meta = parseMeta(EmbeddedFrames.meta(other))
        Dataset(other, meta, frameFromCsv(meta, csv))

  def loadAll: Vector[Dataset] = names.map(load)

  def parseMeta(json: String): DatasetMeta =
    val root = JsonValue.parse(json)
    DatasetMeta(
      name = root.req("name").asString,
      nRows = root.req("n_rows").asInt,
      source = root.req("source").asString,
      license = root.req("license").asString,
      columns = root
        .req("columns")
        .asArray
        .map: col =>
          val kind = col.req("type").asString match
            case "numeric"     => ColumnKind.Numeric
            case "categorical" => ColumnKind.Categorical
            case other         => throw IllegalArgumentException(s"unknown column type `$other`")
          ColumnMeta(
            name = col.req("name").asString,
            kind = kind,
            levels = col.field("levels").map(_.asArray.map(_.asString))
          )
    )

  def frameFromCsv(meta: DatasetMeta, csv: String): ModelFrame =
    val (header, rows) = Csv.parse(csv)
    if header != meta.columns.map(_.name) then
      throw IllegalArgumentException(s"${meta.name}: CSV header $header != ${meta.columns.map(_.name)}")
    if rows.length != meta.nRows then
      throw IllegalArgumentException(s"${meta.name}: expected ${meta.nRows} rows, got ${rows.length}")
    val cols = meta.columns.zipWithIndex.map: (spec, idx) =>
      val values = rows.map(_(idx))
      spec.kind match
        case ColumnKind.Numeric =>
          spec.name -> ModelFrame.numeric(values.map(_.toDouble))
        case ColumnKind.Categorical =>
          val factor =
            spec.levels match
              case Some(levels) =>
                ModelFrame.factor(values, levels).fold(err => throw IllegalArgumentException(err.message), identity)
              case None => ModelFrame.factor(values)
          spec.name -> factor
    ModelFrame.of(cols*).fold(err => throw IllegalStateException(err.message), identity)

  private def builtin(
      name: String,
      nRows: Int,
      source: String,
      license: String,
      columns: Vector[ColumnMeta],
      frame: ModelFrame
  ): Dataset =
    Dataset(name, DatasetMeta(name, nRows, source, license, columns), frame)
