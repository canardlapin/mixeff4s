package mixeff4s.design

import gale.linalg.*
import mixeff4s.data.{CategoricalCoding, Column, EncodedColumn, ModelFrame, encodedColumns}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.{FixedTerm, Formula, GroupingFactor, RandomCovariance, RandomTerm}

/** Compiled mixed-model design: rank-truncated X, ReMats, and the live θ map. */
final case class CompiledDesign(
    formula: Formula,
    fe: FeTerm,
    xy: FeMat,
    reterms: Vector[ReMat],
    parmap: Vector[(Int, Int, Int)],
    prediction: Option[PredictionDesign] = None
):
  def n: Int = xy.n
  def p: Int = fe.rank
  def nReTerms: Int = reterms.length
  def nRanef: Int = reterms.map(_.nRanef).sum
  def nTheta: Int = parmap.length
  def theta: Vector[Double] = reterms.flatMap(_.theta)

object Design:
  def compile(formula: Formula, frame: ModelFrame): FitResult[CompiledDesign] =
    for
      materialized <- formula.materialize(frame)
      _ <- refuseUnsupportedCovariance(formula)
      y <- response(formula, materialized)
      rawX <- buildFixedDesign(formula, materialized)
      fe <- FeTerm(rawX._1, rawX._2)
      xy <- FeMat(fe, y)
      reterms <- buildReTerms(formula, materialized)
    yield
      val paired = reterms.zip(formula.randomTerms).sortBy(pair => -pair._1.nRanef)
      val ordered = paired.map(_._1)
      CompiledDesign(
        formula,
        fe,
        xy,
        ordered,
        ReMat.buildParmap(ordered),
        Some(PredictionDesign.capture(formula, materialized, fe, xy, paired.map(_._2), ordered))
      )

  private def refuseUnsupportedCovariance(formula: Formula): FitResult[Unit] =
    formula.randomTerms.collectFirst:
      case term if !term.covariance.isSupportedForFit =>
        val label = term.covariance match
          case RandomCovariance.CompoundSymmetry => "cs(...)"
          case RandomCovariance.Ar1              => "ar1(...)"
          case other                             => other.label
        MixedModelError.Unsupported(s"$label random-effect covariance is parsed and refused for fitting")
    match
      case Some(err) => Left(err)
      case None      => Right(())

  private def response(formula: Formula, frame: ModelFrame): FitResult[Vector[Double]] =
    frame.column(formula.response) match
      case Some(Column.Numeric(values)) => Right(values)
      case Some(_)                      =>
        Left(MixedModelError.InvalidArgument(s"response `${formula.response}` must be numeric"))
      case None =>
        Left(MixedModelError.InvalidArgument(s"column `${formula.response}` is not present in the model frame"))

  private[design] def buildFixedDesign(
      formula: Formula,
      frame: ModelFrame
  ): FitResult[(DMat, Vector[String])] =
    val n = frame.nRows
    val intercept =
      if formula.hasIntercept then Vector(EncodedColumn("(Intercept)", Vector.fill(n)(1.0)))
      else Vector.empty
    buildBasis(formula.fixedTerms, formula.hasIntercept, frame, intercept)
      .map: built =>
        if built.isEmpty then (Matrix.zeros(n, 0), Vector.empty)
        else (Matrix.tabulate(n, built.length)((i, j) => built(j).values(i)), built.map(_.name))

  private def buildReTerms(formula: Formula, frame: ModelFrame): FitResult[Vector[ReMat]] =
    if formula.randomTerms.isEmpty then Left(MixedModelError.NoRandomEffects)
    else
      formula.randomTerms.foldLeft[FitResult[Vector[ReMat]]](Right(Vector.empty)): (acc, term) =>
        acc.flatMap(soFar => buildReMat(term, frame).map(soFar :+ _))

  private def buildReMat(rt: RandomTerm, frame: ModelFrame): FitResult[ReMat] =
    for
      grouping <- groupingLevels(rt.grouping, frame)
      (groupName, refs, levels) = grouping
      basis <- reBasis(rt, frame)
      z = Matrix.tabulate(basis.length, frame.nRows)((row, obs) => basis(row).values(obs))
      remat <- ReMat(groupName, refs, levels, basis.map(_.name), z)
    yield
      if rt.zerocorr || rt.covariance == RandomCovariance.Diagonal then remat.zerocorr
      else remat

  private[design] def reBasis(rt: RandomTerm, frame: ModelFrame): FitResult[Vector[EncodedColumn]] =
    val hasIntercept = rt.terms.exists(_ == FixedTerm.Intercept) || rt.terms.isEmpty
    val start =
      if hasIntercept then Vector(EncodedColumn("(Intercept)", Vector.fill(frame.nRows)(1.0)))
      else Vector.empty
    buildBasis(rt.terms, hasIntercept, frame, start)

  /** Apply R's contrast assignment to either fixed terms or a random-slope basis. */
  private def buildBasis(
      terms: Vector[FixedTerm],
      hasIntercept: Boolean,
      frame: ModelFrame,
      start: Vector[EncodedColumn]
  ): FitResult[Vector[EncodedColumn]] =
    val fullCoded = noInterceptFullCodingTarget(terms, hasIntercept, frame)
    val order = variableOrder(terms)
    terms.foldLeft[FitResult[Vector[EncodedColumn]]](Right(start)): (acc, term) =>
      acc.flatMap: soFar =>
        term match
          case FixedTerm.Intercept | FixedTerm.NoIntercept =>
            Right(soFar)
          case FixedTerm.Column(name) =>
            val coding =
              if fullCoded.contains(name) then CategoricalCoding.CellMeans
              else CategoricalCoding.Treatment
            expandFactor(name, frame, coding).map(soFar ++ _)
          case FixedTerm.Interaction(vars) =>
            val treatment = interactionTreatmentVariables(terms, vars) -- fullCoded
            val ordered = order.filter(vars.contains)
            expandInteraction(ordered, frame, frame.nRows, treatment).map(soFar ++ _)

  private def groupingLevels(
      grouping: GroupingFactor,
      frame: ModelFrame
  ): FitResult[(String, Vector[Int], Vector[String])] =
    grouping match
      case GroupingFactor.Single(name) =>
        frame.factor(name) match
          case Some(factor) => Right(compactGrouping(name, factor))
          case None         =>
            Left(
              MixedModelError.InvalidArgument(s"Grouping factor '$name' not found or not categorical")
            )
      case GroupingFactor.Interaction(names) =>
        combinedGrouping(names, frame)
      case GroupingFactor.Cell(names) =>
        combinedGrouping(names, frame)

  private def combinedGrouping(
      names: Vector[String],
      frame: ModelFrame
  ): FitResult[(String, Vector[Int], Vector[String])] =
    val cats = names.map(frame.factor)
    if cats.exists(_.isEmpty) then
      val missing = names.zip(cats).collect { case (name, None) => name }.mkString(", ")
      Left(MixedModelError.InvalidArgument(s"Grouping factor '$missing' not found"))
    else
      val factors = cats.flatten
      // Keep the tuple structured: labels such as ("a_b", "c") and ("a", "b_c")
      // must not collapse to one grouping level.
      val keys = Vector.tabulate(frame.nRows): obs =>
        factors.map(f => f.levels(f.refs(obs)))
      val order = scala.collection.mutable.LinkedHashMap.empty[Vector[String], Int]
      val refs = keys.map: key =>
        order.getOrElseUpdate(key, order.size)
      val labels = order.keysIterator.map(tuple => tuple.map(value => s"${value.length}:$value").mkString("|"))
      Right((names.mkString(" & "), refs, labels.toVector))

  private def compactGrouping(
      name: String,
      factor: Column.Factor
  ): (String, Vector[Int], Vector[String]) =
    val used = factor.refs.toSet
    val retained = factor.levels.zipWithIndex.collect { case (level, old) if used(old) => (level, old) }
    val oldToNew = retained.zipWithIndex.map { case ((_, old), next) => old -> next }.toMap
    (name, factor.refs.map(oldToNew), retained.map(_._1))

  private def expandFactor(
      name: String,
      frame: ModelFrame,
      coding: CategoricalCoding
  ): FitResult[Vector[EncodedColumn]] =
    frame.column(name) match
      case Some(Column.Numeric(values)) =>
        Right(Vector(EncodedColumn(name, values)))
      case Some(factor: Column.Factor) =>
        Right(factor.encodedColumns(name, coding))
      case None =>
        Left(MixedModelError.InvalidArgument(s"Column '$name' not found in data"))

  private def expandInteraction(
      vars: Vector[String],
      frame: ModelFrame,
      n: Int,
      treatmentVariables: Set[String]
  ): FitResult[Vector[EncodedColumn]] =
    vars
      .foldLeft[FitResult[Vector[Vector[EncodedColumn]]]](Right(Vector.empty)): (acc, name) =>
        val coding =
          if treatmentVariables(name) then CategoricalCoding.Treatment
          else CategoricalCoding.CellMeans
        acc.flatMap(soFar => expandFactor(name, frame, coding).map(soFar :+ _))
      .map(cartesianInteraction(_, n))

  /** Variables in their formula first-appearance order. */
  private def variableOrder(terms: Vector[FixedTerm]): Vector[String] =
    terms.foldLeft(Vector.empty[String]): (order, term) =>
      term.columnNames.foldLeft(order): (soFar, name) =>
        if soFar.contains(name) then soFar else soFar :+ name

  /** The first categorical variable in R's no-intercept term order gets full indicators. */
  private def noInterceptFullCodingTarget(
      terms: Vector[FixedTerm],
      hasIntercept: Boolean,
      frame: ModelFrame
  ): Option[String] =
    if hasIntercept then None
    else
      val order = variableOrder(terms)
      terms.zipWithIndex
        .collect { case (term, index) if term.columnNames.nonEmpty => (index, term.columnNames) }
        .sortBy { case (index, names) => (names.length, index) }
        .iterator
        .flatMap: (_, names) =>
          order.iterator.filter(name => names.contains(name) && frame.factor(name).isDefined)
        .toSeq
        .headOption

  /** R's marginality rule: only variables present in every not-yet-covered ANOVA component of an interaction receive
    * ordinary treatment contrasts.
    */
  private def interactionTreatmentVariables(
      terms: Vector[FixedTerm],
      vars: Vector[String]
  ): Set[String] =
    val lower = terms.flatMap: term =>
      val names = term.columnNames
      if names.nonEmpty && names.length < vars.length && names.forall(vars.contains) then Some(names.toSet)
      else None
    val components = nonEmptySubsets(vars.distinct).filterNot: component =>
      lower.exists(lowerTerm => component.subsetOf(lowerTerm))
    components.reduceOption(_ intersect _).getOrElse(Set.empty)

  private def nonEmptySubsets(names: Vector[String]): Vector[Set[String]] =
    (1 until (1 << names.length)).toVector.map: mask =>
      names.indices.collect { case index if (mask & (1 << index)) != 0 => names(index) }.toSet

  private def cartesianInteraction(
      perVar: Vector[Vector[EncodedColumn]],
      n: Int
  ): Vector[EncodedColumn] =
    if perVar.isEmpty then Vector.empty
    else
      perVar.foldRight(Vector(EncodedColumn("", Vector.fill(n)(1.0)))): (cols, acc) =>
        acc.flatMap: right =>
          cols.map: left =>
            val name = if right.name.isEmpty then left.name else s"${left.name}:${right.name}"
            EncodedColumn(name, left.values.zip(right.values).map(_ * _))
