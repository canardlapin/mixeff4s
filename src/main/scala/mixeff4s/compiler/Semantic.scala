package mixeff4s.compiler

import mixeff4s.formula.{FixedTerm, Formula, GroupingFactor, RandomCovariance, RandomTerm}

enum CovarianceForm:
  case Scalar, Diagonal, Full
  case Structured(kind: String)

  def code: String =
    this match
      case CovarianceForm.Scalar         => "scalar"
      case CovarianceForm.Diagonal       => "diagonal"
      case CovarianceForm.Full           => "full"
      case CovarianceForm.Structured(k)  => k

  def support: CovarianceSupport =
    this match
      case CovarianceForm.Scalar | CovarianceForm.Diagonal | CovarianceForm.Full =>
        CovarianceSupport.Supported
      case CovarianceForm.Structured(_) =>
        CovarianceSupport.ParsedRefused

enum CovarianceSupport:
  case Supported, ParsedRefused

  def code: String =
    this match
      case CovarianceSupport.Supported     => "supported"
      case CovarianceSupport.ParsedRefused => "parsed_refused"

enum InterceptPolicy:
  case Included, Omitted

  def code: String =
    this match
      case InterceptPolicy.Included => "included"
      case InterceptPolicy.Omitted  => "omitted"

enum CoefficientKind:
  case Intercept, Slope, Interaction

  def code: String =
    this match
      case CoefficientKind.Intercept   => "intercept"
      case CoefficientKind.Slope       => "slope"
      case CoefficientKind.Interaction => "interaction"

final case class RandomCoefficient(name: String, kind: CoefficientKind, source: String)

final case class SemanticRandomTerm(
    id: String,
    grouping: GroupingFactor,
    basis: Vector[RandomCoefficient],
    covariance: CovarianceForm,
    support: CovarianceSupport,
    intercept: InterceptPolicy,
    sourceText: String,
    written: Option[String]
)

/** Formula-level semantic IR. No design matrices, no optimizer folklore. */
final case class SemanticModel(
    response: String,
    fixedTerms: Vector[String],
    randomTerms: Vector[SemanticRandomTerm]
)

object SemanticModel:
  val SchemaName = "mixeff4s.semantic_model"
  val SchemaVersion = 1

  def from(formula: Formula): SemanticModel =
    SemanticModel(
      formula.response,
      formula.fixedTerms.map(_.toString),
      formula.randomTerms.zipWithIndex.map(compileTerm)
    )

  private def compileTerm(term: RandomTerm, index: Int): SemanticRandomTerm =
    val intercept =
      if term.terms.contains(FixedTerm.NoIntercept) then InterceptPolicy.Omitted
      else InterceptPolicy.Included
    val slopes = term.terms.flatMap:
      case FixedTerm.Column(name) =>
        Vector(RandomCoefficient(name, CoefficientKind.Slope, name))
      case FixedTerm.Interaction(names) =>
        val src = names.mkString(":")
        Vector(RandomCoefficient(src, CoefficientKind.Interaction, src))
      case FixedTerm.Intercept | FixedTerm.NoIntercept =>
        Vector.empty
    val basis =
      if intercept == InterceptPolicy.Included then
        Vector(RandomCoefficient("intercept", CoefficientKind.Intercept, "1")) ++ slopes
      else slopes
    val form =
      if term.zerocorr && basis.length > 1 then CovarianceForm.Diagonal
      else covarianceForm(term.covariance, basis.length)
    val written = term.source.map(_.written).filter(_ != term.toString)
    SemanticRandomTerm(
      s"r$index",
      term.grouping,
      basis,
      form,
      form.support,
      intercept,
      term.toString,
      written
    )

  private def covarianceForm(covariance: RandomCovariance, basisLen: Int): CovarianceForm =
    (covariance, basisLen) match
      case (_, 0) =>
        CovarianceForm.Structured("empty_basis")
      case (RandomCovariance.Full, 1) =>
        CovarianceForm.Scalar
      case (RandomCovariance.Full, _) =>
        CovarianceForm.Full
      case (RandomCovariance.Diagonal, _) =>
        CovarianceForm.Diagonal
      case (RandomCovariance.CompoundSymmetry, _) =>
        CovarianceForm.Structured("compound_symmetry")
      case (RandomCovariance.Ar1, _) =>
        CovarianceForm.Structured("ar1")

  private[compiler] def encode(model: SemanticModel, indent: Int): String =
    Json.pretty(
      Vector(
        "schema" -> Json.str(SchemaName),
        "version" -> Json.num(SchemaVersion),
        "response" -> Json.str(model.response),
        "fixed_terms" -> Json.arr(model.fixedTerms.map(Json.str)),
        "random_terms" -> Json.arr(model.randomTerms.map(termJson))
      ),
      indent = indent
    )

  private def termJson(term: SemanticRandomTerm): String =
    val fields = Vector.newBuilder[(String, String)]
    fields += "id" -> Json.str(term.id)
    fields += "grouping" -> groupingJson(term.grouping)
    fields += "basis" -> Json.arr(term.basis.map(coefJson))
    fields += "covariance" -> Json.str(term.covariance.code)
    fields += "support" -> Json.str(term.support.code)
    fields += "intercept" -> Json.str(term.intercept.code)
    fields += "source_text" -> Json.str(term.sourceText)
    term.written.foreach(text => fields += "written" -> Json.str(text))
    Json.obj(fields.result()*)

  private def coefJson(coef: RandomCoefficient): String =
    Json.obj(
      "name" -> Json.str(coef.name),
      "kind" -> Json.str(coef.kind.code),
      "source" -> Json.str(coef.source)
    )

  private def groupingJson(grouping: GroupingFactor): String =
    grouping match
      case GroupingFactor.Single(name) =>
        Json.obj("kind" -> Json.str("single"), "name" -> Json.str(name))
      case GroupingFactor.Interaction(names) =>
        Json.obj("kind" -> Json.str("interaction"), "names" -> Json.arr(names.map(Json.str)))
      case GroupingFactor.Cell(names) =>
        Json.obj("kind" -> Json.str("cell"), "names" -> Json.arr(names.map(Json.str)))
