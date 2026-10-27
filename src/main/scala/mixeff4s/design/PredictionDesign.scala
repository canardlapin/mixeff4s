package mixeff4s.design

import gale.linalg.*
import mixeff4s.data.{Column, ModelFrame}
import mixeff4s.error.{FitResult, MixedModelError}
import mixeff4s.formula.{Expr, Formula, GroupingFactor, RandomTerm}

final case class PredictionRandomRows(refs: Vector[Int], basis: DMat)
final case class PredictionRows(fixed: DMat, random: Vector[PredictionRandomRows])

/** Frozen training encoding and alias relations. New rows never select another rank or reorder random terms. */
final class PredictionDesign private (
    formula: Formula,
    schema: Map[String, Option[Vector[String]]],
    layout: FixedEffectLayout,
    randomTerms: Vector[RandomTerm],
    randomLevels: Vector[Vector[String]],
    randomNames: Vector[Vector[String]],
    aliases: FitResult[Vector[(Int, Vector[Double], Double)]],
    expectedFixed: FeTerm,
    expectedXY: FeMat,
    expectedRandom: Vector[ReMat]
):
  val aliasTolerance: Double = 1e-7
  private val expectedParmap = ReMat.buildParmap(expectedRandom)

  private[mixeff4s] def matches(design: CompiledDesign): Boolean =
    if design.formula != formula || design.fe != expectedFixed || design.reterms != expectedRandom ||
      design.parmap != expectedParmap || design.xy.n != expectedFixed.n || design.xy.rank != layout.rank
    then false
    else if design.xy eq expectedXY then true
    else
      var i = 0
      while i < design.n do
        var j = 0
        while j < layout.rank do
          if design.xy.xy(i, j) != expectedFixed.x(i, j) then return false
          j += 1
        i += 1
      true

  def encode(frame: ModelFrame, conditional: Boolean, allowNewGroups: Boolean): FitResult[PredictionRows] =
    val usedRandom = if conditional then randomTerms else Vector.empty
    val direct = (formula.fixedTerms.flatMap(_.columnNames) ++ usedRandom.flatMap(_.terms.flatMap(_.columnNames))).toSet
    val needed = requiredColumns(direct)
    val selected = formula.copy(randomTerms = usedRandom, derived = formula.derived.filter(d => needed(d.label)))
    for
      aligned <- align(frame, needed)
      materialized <- selected.materialize(aligned)
      raw <- Design.buildFixedDesign(selected, materialized)
      _ <-
        if raw._2 == layout.requestedNames then Right(())
        else refusal("prediction_encoding_mismatch", "Prediction columns differ from the fitted encoding.")
      relations <- aliases
      _ <- checkAliases(raw._1, relations)
      random <- if conditional then encodeRandom(materialized, allowNewGroups) else Right(Vector.empty)
    yield PredictionRows(
      Matrix.tabulate(frame.nRows, layout.rank)((i, j) => raw._1(i, layout.retainedIndices(j))),
      random
    )

  private def requiredColumns(initial: Set[String]): Set[String] =
    def variables(expr: Expr): Set[String] = expr match
      case Expr.Lit(_)              => Set.empty
      case Expr.Col(name)           => Set(name)
      case Expr.Neg(inner)          => variables(inner)
      case Expr.Bin(_, left, right) => variables(left) ++ variables(right)
      case Expr.Call(_, arg)        => variables(arg)
    var required = initial
    var changed = true
    while changed do
      val expanded =
        formula.derived.filter(d => required(d.label)).foldLeft(required)((names, d) => names ++ variables(d.expr))
      changed = expanded != required
      required = expanded
    required

  private def align(frame: ModelFrame, needed: Set[String]): FitResult[ModelFrame] =
    frame.columns
      .foldLeft[FitResult[Vector[(String, Column)]]](Right(Vector.empty)): (acc, entry) =>
        acc.flatMap: columns =>
          val (name, column) = entry
          if !needed(name) then Right(columns :+ entry)
          else
            schema.get(name) match
              case Some(None) =>
                column match
                  case _: Column.Numeric => Right(columns :+ entry)
                  case _ => refusal("prediction_column_type", s"Fitted numeric predictor `$name` must remain numeric.")
              case Some(Some(levels)) =>
                column match
                  case factor: Column.Factor =>
                    ModelFrame
                      .factor(factor.values, levels)
                      .left
                      .map(_ =>
                        MixedModelError.InferenceUnavailable(
                          "prediction_unknown_fixed_level",
                          s"Predictor `$name` contains a level absent from the training encoding."
                        )
                      )
                      .map(value => columns :+ (name -> value))
                  case _ =>
                    refusal("prediction_column_type", s"Fitted categorical predictor `$name` must remain categorical.")
              case None => Right(columns :+ entry)
      .flatMap(columns => ModelFrame.of(columns*))

  private def checkAliases(raw: DMat, relations: Vector[(Int, Vector[Double], Double)]): FitResult[Unit] =
    var i = 0
    while i < raw.rows do
      var j = 0
      while j < raw.cols do
        if !raw(i, j).isFinite then
          return refusal("prediction_nonfinite_design", s"Expanded predictor at row $i, column $j is non-finite.")
        j += 1
      i += 1
    var relation = 0
    while relation < relations.length do
      val (dropped, weights, trainingScale) = relations(relation)
      if weights.exists(w => !w.isFinite) then
        return refusal("prediction_alias_unavailable", "Fitted alias coefficients are non-finite.")
      var row = 0
      while row < raw.rows do
        var predicted = 0.0
        var scale = 0.0
        var j = 0
        while j < weights.length do
          val value = raw(row, layout.retainedIndices(j)) * weights(j)
          predicted += value
          scale += math.abs(value)
          j += 1
        val actual = raw(row, dropped)
        if !predicted.isFinite || !scale.isFinite then
          return refusal("prediction_nonfinite_design", s"Alias evaluation overflows at row $row.")
        if math.abs(actual - predicted) > aliasTolerance * math.max(trainingScale, math.max(math.abs(actual), scale))
        then
          return refusal(
            "prediction_nonestimable",
            s"Row $row violates the fitted alias relation for `${layout.requestedNames(dropped)}`."
          )
        row += 1
      relation += 1
    Right(())

  private def encodeRandom(frame: ModelFrame, allowNew: Boolean): FitResult[Vector[PredictionRandomRows]] =
    randomTerms.zipWithIndex.foldLeft[FitResult[Vector[PredictionRandomRows]]](Right(Vector.empty)): (acc, entry) =>
      acc.flatMap: built =>
        val (term, index) = entry
        for
          basis <- Design.reBasis(term, frame)
          _ <-
            if basis.map(_.name) == randomNames(index) then Right(())
            else refusal("prediction_encoding_mismatch", "Random-effect basis differs from training.")
          labels <- groupLabels(term.grouping, frame)
          refs <-
            val lookup = randomLevels(index).zipWithIndex.toMap
            val unknown = labels.find(label => !lookup.contains(label))
            if unknown.nonEmpty && !allowNew then
              refusal("prediction_new_group", s"Unknown level `${unknown.get}` in grouping `${term.grouping}`.")
            else Right(labels.map(label => lookup.getOrElse(label, -1)))
        yield built :+ PredictionRandomRows(
          refs,
          Matrix.tabulate(basis.length, frame.nRows)((r, c) => basis(r).values(c))
        )

  private def groupLabels(grouping: GroupingFactor, frame: ModelFrame): FitResult[Vector[String]] =
    val names = grouping.factorNames
    if names.exists(name => frame.factor(name).isEmpty) then
      refusal("prediction_group_type", s"Grouping `$grouping` requires categorical columns ${names.mkString(", ")}.")
    else
      val factors = names.map(name => frame.factor(name).get)
      grouping match
        case GroupingFactor.Single(_) => Right(factors.head.values)
        case _                        =>
          Right(
            Vector.tabulate(frame.nRows)(i =>
              factors.map(_.values(i)).map(value => s"${value.length}:$value").mkString("|")
            )
          )

  private def refusal[A](code: String, message: String): FitResult[A] = Left(
    MixedModelError.InferenceUnavailable(code, message)
  )

object PredictionDesign:
  private[design] def capture(
      formula: Formula,
      frame: ModelFrame,
      fe: FeTerm,
      xy: FeMat,
      terms: Vector[RandomTerm],
      random: Vector[ReMat]
  ): PredictionDesign =
    val schema = frame.columns.map:
      case (name, Column.Numeric(_))           => name -> None
      case (name, Column.Factor(levels, _, _)) => name -> Some(levels)
    new PredictionDesign(
      formula,
      schema.toMap,
      fe.fixedEffects,
      terms,
      random.map(_.levels),
      random.map(_.cnames),
      aliasMap(fe),
      fe,
      xy,
      random
    )

  /** Two-pass thin QR of the retained training basis; no normal-equation squaring of its condition number. */
  private def aliasMap(fe: FeTerm): FitResult[Vector[(Int, Vector[Double], Double)]] =
    val layout = fe.fixedEffects
    if layout.droppedIndices.isEmpty then Right(Vector.empty)
    else
      val p = fe.rank
      val q = Array.ofDim[Double](p, fe.n)
      val r = Array.ofDim[Double](p, p)
      var j = 0
      while j < p do
        val residual = Array.tabulate(fe.n)(i => fe.x(i, j))
        var pass = 0
        while pass < 2 do
          var k = 0
          while k < j do
            var projection = 0.0
            var i = 0
            while i < fe.n do
              projection += q(k)(i) * residual(i)
              i += 1
            r(k)(j) += projection
            i = 0
            while i < fe.n do
              residual(i) -= projection * q(k)(i)
              i += 1
            k += 1
          pass += 1
        val norm = residual.foldLeft(0.0)(math.hypot)
        if !norm.isFinite || norm <= 0.0 then
          return Left(
            MixedModelError.InferenceUnavailable(
              "prediction_alias_unavailable",
              "Cannot establish the fitted alias relations."
            )
          )
        r(j)(j) = norm
        var i = 0
        while i < fe.n do
          q(j)(i) = residual(i) / norm
          i += 1
        j += 1
      Right(layout.droppedIndices.map: dropped =>
        val col = fe.piv.indexOf(dropped)
        val coefficients = Array.tabulate(p)(k => (0 until fe.n).iterator.map(i => q(k)(i) * fe.x(i, col)).sum)
        var k = p - 1
        while k >= 0 do
          var value = coefficients(k)
          var m = k + 1
          while m < p do
            value -= r(k)(m) * coefficients(m)
            m += 1
          coefficients(k) = value / r(k)(k)
          k -= 1
        val trainingScale = (0 until fe.n).iterator.map(i => math.abs(fe.x(i, col))).foldLeft(0.0)(math.max)
        (dropped, coefficients.toVector, trainingScale))
