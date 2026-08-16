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
    parmap: Vector[(Int, Int, Int)]
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
      val ordered = reterms.sortBy(rt => -rt.nRanef)
      CompiledDesign(formula, fe, xy, ordered, ReMat.buildParmap(ordered))

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
      case Some(_) =>
        Left(MixedModelError.InvalidArgument(s"response `${formula.response}` must be numeric"))
      case None =>
        Left(MixedModelError.InvalidArgument(s"column `${formula.response}` is not present in the model frame"))

  private def buildFixedDesign(
      formula: Formula,
      frame: ModelFrame
  ): FitResult[(DMat, Vector[String])] =
    val n = frame.nRows
    val intercept =
      if formula.hasIntercept then Vector(EncodedColumn("(Intercept)", Vector.fill(n)(1.0)))
      else Vector.empty
    formula.fixedTerms
      .foldLeft[FitResult[Vector[EncodedColumn]]](Right(intercept)): (acc, term) =>
        acc.flatMap: soFar =>
          term match
            case FixedTerm.Intercept | FixedTerm.NoIntercept =>
              Right(soFar)
            case FixedTerm.Column(name) =>
              expandFactor(name, frame, CategoricalCoding.Treatment).map(soFar ++ _)
            case FixedTerm.Interaction(vars) =>
              expandInteraction(vars, frame, n, CategoricalCoding.Treatment).map(soFar ++ _)
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

  private def reBasis(rt: RandomTerm, frame: ModelFrame): FitResult[Vector[EncodedColumn]] =
    val coding =
      if rt.terms.exists(_ == FixedTerm.NoIntercept) then CategoricalCoding.CellMeans
      else CategoricalCoding.Treatment
    val hasIntercept = rt.terms.exists(_ == FixedTerm.Intercept) || rt.terms.isEmpty
    val start =
      if hasIntercept then Vector(EncodedColumn("(Intercept)", Vector.fill(frame.nRows)(1.0)))
      else Vector.empty
    rt.terms.foldLeft[FitResult[Vector[EncodedColumn]]](Right(start)): (acc, term) =>
      acc.flatMap: soFar =>
        term match
          case FixedTerm.Intercept | FixedTerm.NoIntercept =>
            Right(soFar)
          case FixedTerm.Column(name) =>
            expandFactor(name, frame, coding).map(soFar ++ _)
          case FixedTerm.Interaction(vars) =>
            expandInteraction(vars, frame, frame.nRows, coding).map(soFar ++ _)

  private def groupingLevels(
      grouping: GroupingFactor,
      frame: ModelFrame
  ): FitResult[(String, Vector[Int], Vector[String])] =
    grouping match
      case GroupingFactor.Single(name) =>
        frame.factor(name) match
          case Some(factor) => Right((name, factor.refs, factor.levels))
          case None =>
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
      val keys = Vector.tabulate(frame.nRows): obs =>
        factors.map(f => f.levels(f.refs(obs))).mkString("_")
      val order = scala.collection.mutable.LinkedHashMap.empty[String, Int]
      val refs = keys.map: key =>
        order.getOrElseUpdate(key, order.size)
      Right((names.mkString(" & "), refs, order.keys.toVector))

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
      coding: CategoricalCoding
  ): FitResult[Vector[EncodedColumn]] =
    vars
      .foldLeft[FitResult[Vector[Vector[EncodedColumn]]]](Right(Vector.empty)): (acc, name) =>
        acc.flatMap(soFar => expandFactor(name, frame, coding).map(soFar :+ _))
      .map(cartesianInteraction(_, n))

  private def cartesianInteraction(
      perVar: Vector[Vector[EncodedColumn]],
      n: Int
  ): Vector[EncodedColumn] =
    if perVar.isEmpty then Vector.empty
    else
      perVar.foldLeft(Vector(EncodedColumn("", Vector.fill(n)(1.0)))): (acc, cols) =>
        acc.flatMap: left =>
          cols.map: right =>
            val name = if left.name.isEmpty then right.name else s"${left.name}:${right.name}"
            EncodedColumn(name, left.values.zip(right.values).map(_ * _))
