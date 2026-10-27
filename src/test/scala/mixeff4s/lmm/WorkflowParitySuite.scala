package mixeff4s.lmm

import mixeff4s.comparison.{JsonValue, WorkflowFrozen}
import mixeff4s.data.ModelFrame
import mixeff4s.fixtures.{Sleepstudy, Pastes, Penicillin}

class WorkflowParitySuite extends munit.FunSuite:
  private def vector(value: JsonValue): Vector[Double] = value match
    case JsonValue.Num(n) => Vector(n)
    case _                => value.asArray.map(_.asDouble)
  private def strings(value: JsonValue): Vector[String] = value match
    case JsonValue.Str(s) => Vector(s)
    case _                => value.asArray.map(_.asString)
  private val references = JsonValue.parse(WorkflowFrozen.json).req("results").asArray
  private def frame(name: String): ModelFrame = name match
    case "sleepstudy" => Sleepstudy.frame
    case "pastes"     => Pastes.frame
    case "penicillin" => Penicillin.frame
    case other        => fail(other)
  private def close(actual: Vector[Double], expected: Vector[Double], tolerance: Double): Unit =
    assertEquals(actual.length, expected.length)
    actual.zip(expected).zipWithIndex.foreach { case ((a, b), i) => assertEqualsDouble(a, b, tolerance, clues(i)) }

  for reference <- references; strategy <- LmmOptimizer.values do
    val dataset = reference.req("dataset").asString
    val criterion = if reference.req("criterion").asString == "ML" then Criterion.ML else Criterion.REML
    test(s"frozen lme4 fit/prediction/residual/modes: $dataset $criterion $strategy"):
      val data = frame(dataset)
      val fit = Lmm
        .fit(reference.req("formula").asString, data, FitOptions(criterion, optimizer = strategy))
        .fold(e => fail(e.message), identity)
      close(fit.beta, vector(reference.req("beta")), 1e-3)
      assertEqualsDouble(fit.objective, reference.req("objective").asDouble, 1e-3)
      assertEqualsDouble(fit.sigma, reference.req("sigma").asDouble, 1e-3)
      close(fit.theta, vector(reference.req("theta")), if dataset == "pastes" then 5e-3 else 1e-3)
      close(Prediction.fitted(fit).toOption.get, vector(reference.req("conditional")), 1e-2)
      close(Prediction.fitted(fit, PredictionMode.Population).toOption.get, vector(reference.req("population")), 1e-3)
      close(Prediction.predict(fit, data).toOption.get, vector(reference.req("conditional")), 1e-2)
      close(Prediction.residuals(fit).toOption.get, vector(reference.req("residuals")), 1e-2)
      val covariance = reference.req("vcov").asArray.map(vector)
      fit.vcov.zip(covariance).foreach((a, b) => close(a, b, 1e-2))
      val context = fit.validatedContext.toOption.get
      context.design.reterms.zipWithIndex.foreach: (term, block) =>
        val names = term.groupingName.split(" & ").toVector.sorted
        val matched = reference.req("modes").asObject.find((name, _) => name.split(":").toVector.sorted == names).get
        val groupVariables = matched._1.split(":").toVector
        val modes = matched._2
        val levels = strings(modes.req("levels")).zipWithIndex.toMap
        val coefficientNames = strings(modes.req("names"))
        val values = modes.req("values").asArray.map(vector)
        var i = 0
        while i < fit.n do
          val label = groupVariables.map(name => data.factor(name).get.values(i)).mkString(":")
          val expected = values(levels(label))
          term.cnames.zipWithIndex.foreach: (name, j) =>
            assertEqualsDouble(
              context.modes(block)(term.refs(i) * term.vsize + j),
              expected(coefficientNames.indexOf(name)),
              1e-2
            )
          i += 1
