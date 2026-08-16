package mixeff4s.formula.dsl

import mixeff4s.formula.{FixedTerm, Formula, GroupingFactor, RandomCovariance, RandomTerm}

final case class Response(name: String)

final case class FixedExpr(terms: Vector[FixedTerm]):
  def +(term: FixedTerm): FixedExpr = FixedExpr(terms :+ term)
  def +(random: RandomTerm): FormulaRhs = FormulaRhs(terms, Vector(random))
  def +(rhs: FormulaRhs): FormulaRhs = rhs.copy(fixed = terms ++ rhs.fixed)
  def |(grouping: GroupingFactor): RandomTerm =
    RandomTerm(terms, grouping, zerocorr = false, RandomCovariance.Full)
  def ||(grouping: GroupingFactor): RandomTerm =
    RandomTerm(terms, grouping, zerocorr = true, RandomCovariance.Diagonal)

final case class FormulaRhs(fixed: Vector[FixedTerm], random: Vector[RandomTerm]):
  def +(term: FixedTerm): FormulaRhs = copy(fixed = fixed :+ term)
  def +(randomTerm: RandomTerm): FormulaRhs = copy(random = random :+ randomTerm)
  def +(other: FormulaRhs): FormulaRhs =
    FormulaRhs(fixed ++ other.fixed, random ++ other.random)

  def toFormula(response: String): Formula =
    val hasExplicit = fixed.exists(t => t == FixedTerm.Intercept || t == FixedTerm.NoIntercept)
    val resolved =
      if hasExplicit then
        if fixed.exists(_ == FixedTerm.NoIntercept) then
          fixed.filterNot(t => t == FixedTerm.Intercept || t == FixedTerm.NoIntercept)
        else fixed
      else FixedTerm.Intercept +: fixed
    Formula(response, resolved, random)

object dsl:
  def response(name: String): Response = Response(name)
  def intercept: FixedTerm = FixedTerm.Intercept
  def noIntercept: FixedTerm = FixedTerm.NoIntercept
  def col(name: String): FixedTerm = FixedTerm.Column(name)
  def factor(name: String): GroupingFactor = GroupingFactor.Single(name)
  def interaction(names: String*): FixedTerm = FixedTerm.Interaction(names.toVector)

  extension (lhs: Response)
    def ~(rhs: FormulaRhs): Formula = rhs.toFormula(lhs.name)
    def ~(term: FixedTerm): Formula = FormulaRhs(Vector(term), Vector.empty).toFormula(lhs.name)
    def ~(random: RandomTerm): Formula = FormulaRhs(Vector.empty, Vector(random)).toFormula(lhs.name)
    def ~(expr: FixedExpr): Formula = FormulaRhs(expr.terms, Vector.empty).toFormula(lhs.name)

  extension (term: FixedTerm)
    def +(other: FixedTerm): FixedExpr = FixedExpr(Vector(term, other))
    def +(random: RandomTerm): FormulaRhs = FormulaRhs(Vector(term), Vector(random))
    def +(rhs: FormulaRhs): FormulaRhs = rhs.copy(fixed = term +: rhs.fixed)
    def |(grouping: GroupingFactor): RandomTerm =
      RandomTerm(Vector(term), grouping, zerocorr = false, RandomCovariance.Full)
    def ||(grouping: GroupingFactor): RandomTerm =
      RandomTerm(Vector(term), grouping, zerocorr = true, RandomCovariance.Diagonal)

  extension (sc: StringContext)
    def formula(args: Any*): Formula =
      require(args.isEmpty, "formula interpolator does not take arguments")
      Formula.parse(sc.parts.mkString) match
        case Right(parsed) => parsed
        case Left(error)   => throw IllegalArgumentException(error.message)
