package mixeff4s.compiler

import mixeff4s.design.{CompiledDesign, ReMat}
import mixeff4s.pathology.{Certificate, Pathology}

/** Versioned compiled-design certificate. Unstable; schema may change. */
final case class Schema(name: String, version: Int, libraryVersion: String)

final case class ModelBoundary(
    modelKind: String,
    responseDistribution: String,
    link: String,
    objectiveApproximation: String,
    inferenceAvailability: String
)

final case class RandomTermCard(
    index: Int,
    grouping: String,
    nLevels: Int,
    basis: Vector[String],
    nRanef: Int,
    nTheta: Int,
    covariance: String
)

final case class ThetaSlot(
    globalIndex: Int,
    termIndex: Int,
    row: Int,
    col: Int,
    name: String,
    constraint: String
)

final case class CompiledArtifact(
    schema: Schema,
    requestedFormula: String,
    n: Int,
    p: Int,
    nReTerms: Int,
    nRanef: Int,
    nTheta: Int,
    feNames: Vector[String],
    randomTerms: Vector[RandomTermCard],
    thetaSlots: Vector[ThetaSlot],
    parmap: Vector[(Int, Int, Int)],
    modelBoundary: ModelBoundary,
    pathology: Certificate
):
  def toJson: String = CompiledArtifact.encode(this)

object CompiledArtifact:
  val SchemaName = "mixeff4s.compiled_design_artifact"
  val SchemaVersion = 1
  val LibraryVersion = "0.1.0-SNAPSHOT"

  def fromDesign(design: CompiledDesign): CompiledArtifact =
    val slots = design.parmap.zipWithIndex.map:
      case ((term, row, col), global) =>
        val basis = design.reterms(term).cnames
        val rowName = basis.lift(row).getOrElse(row.toString)
        val colName = basis.lift(col).getOrElse(col.toString)
        ThetaSlot(
          global,
          term,
          row,
          col,
          s"theta[$term:$rowName,$colName]",
          if row == col then "lower_bound_0" else "unconstrained"
        )
    CompiledArtifact(
      Schema(SchemaName, SchemaVersion, LibraryVersion),
      design.formula.toString,
      design.n,
      design.p,
      design.nReTerms,
      design.nRanef,
      design.nTheta,
      design.fe.fullRankNames,
      design.reterms.zipWithIndex.map(randomCard),
      slots,
      design.parmap,
      ModelBoundary(
        "linear_mixed_model",
        "gaussian",
        "identity",
        "exact_gaussian",
        "not_assessed"
      ),
      Pathology.certify(design)
    )

  private def encode(artifact: CompiledArtifact): String =
    Json.pretty(
      Vector(
        "schema" -> Json.pretty(
          Vector(
            "name" -> Json.str(artifact.schema.name),
            "version" -> Json.num(artifact.schema.version),
            "library_version" -> Json.str(artifact.schema.libraryVersion)
          ),
          indent = 1
        ),
        "requested_formula" -> Json.str(artifact.requestedFormula),
        "n" -> Json.num(artifact.n),
        "p" -> Json.num(artifact.p),
        "n_re_terms" -> Json.num(artifact.nReTerms),
        "n_ranef" -> Json.num(artifact.nRanef),
        "n_theta" -> Json.num(artifact.nTheta),
        "fe_names" -> Json.arr(artifact.feNames.map(Json.str)),
        "random_terms" -> Json.arr(artifact.randomTerms.map(cardJson)),
        "theta_slots" -> Json.arr(artifact.thetaSlots.map(slotJson)),
        "parmap" -> Json.arr(artifact.parmap.map(tripleJson)),
        "model_boundary" -> Json.pretty(
          Vector(
            "model_kind" -> Json.str(artifact.modelBoundary.modelKind),
            "response_distribution" -> Json.str(artifact.modelBoundary.responseDistribution),
            "link" -> Json.str(artifact.modelBoundary.link),
            "objective_approximation" -> Json.str(artifact.modelBoundary.objectiveApproximation),
            "inference_availability" -> Json.str(artifact.modelBoundary.inferenceAvailability)
          ),
          indent = 1
        ),
        "pathology" -> pathologyJson(artifact.pathology)
      )
    )

  private def pathologyJson(cert: Certificate): String =
    Json.pretty(
      Vector(
        "contract_version" -> Json.str(cert.contractVersion),
        "fit_status" -> Json.str(cert.fitStatus.code),
        "expected_statuses" -> Json.arr(cert.expectedStatuses.map(s => Json.str(s.code))),
        "stratum" -> Json.str(cert.stratum.code),
        "n" -> Json.num(cert.n),
        "n_params" -> Json.num(cert.nParams),
        "min_group_size" -> Json.num(cert.minGroupSize),
        "max_group_size" -> Json.num(cert.maxGroupSize),
        "fe_rank" -> Json.num(cert.feRank),
        "n_theta" -> Json.num(cert.nTheta),
        "structural_issue" -> cert.structuralIssue.fold(Json.nul)(issue =>
          Json.obj("code" -> Json.str(issue.code), "details" -> Json.str(issue.details))
        ),
        "notes" -> Json.arr(cert.notes.map(Json.str))
      ),
      indent = 1
    )

  private def randomCard(rt: ReMat, index: Int): RandomTermCard =
    RandomTermCard(index, rt.groupingName, rt.nLevels, rt.cnames, rt.nRanef, rt.nTheta, family(rt))

  private def family(rt: ReMat): String =
    val diagonal = rt.inds.forall: idx =>
      val row = idx % rt.vsize
      val col = idx / rt.vsize
      row == col
    if rt.vsize == 1 then "scalar"
    else if diagonal then "diagonal"
    else "full_cholesky"

  private def cardJson(card: RandomTermCard): String =
    Json.obj(
      "index" -> Json.num(card.index),
      "grouping" -> Json.str(card.grouping),
      "n_levels" -> Json.num(card.nLevels),
      "basis" -> Json.arr(card.basis.map(Json.str)),
      "n_ranef" -> Json.num(card.nRanef),
      "n_theta" -> Json.num(card.nTheta),
      "covariance" -> Json.str(card.covariance)
    )

  private def slotJson(slot: ThetaSlot): String =
    Json.obj(
      "global_index" -> Json.num(slot.globalIndex),
      "term_index" -> Json.num(slot.termIndex),
      "row" -> Json.num(slot.row),
      "col" -> Json.num(slot.col),
      "name" -> Json.str(slot.name),
      "constraint" -> Json.str(slot.constraint)
    )

  private def tripleJson(triple: (Int, Int, Int)): String =
    Json.arr(Vector(Json.num(triple._1), Json.num(triple._2), Json.num(triple._3)))
