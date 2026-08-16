package mixeff4s.data

import mixeff4s.error.MixedModelError
import mixeff4s.formula.Formula

class ModelFrameSuite extends munit.FunSuite:
  test("builds aligned numeric and factor columns"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(1.0, 2.0, 3.0)),
        "g" -> ModelFrame.factor(Vector("b", "a", "b"))
      )
      .getOrElse(fail("frame"))
    assertEquals(frame.nRows, 3)
    assertEquals(frame.factor("g").map(_.levels), Some(Vector("b", "a")))
    assertEquals(frame.factor("g").map(_.refs), Some(Vector(0, 1, 0)))

  test("rejects length mismatch"):
    ModelFrame.of(
      "y" -> ModelFrame.numeric(Vector(1.0, 2.0)),
      "x" -> ModelFrame.numeric(Vector(1.0))
    ) match
      case Left(MixedModelError.DimensionMismatch(_)) => ()
      case other                                      => fail(s"expected dimension mismatch, got $other")

  test("materializes I() onto the frame"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(1.0, 4.0, 9.0)),
        "x" -> ModelFrame.numeric(Vector(1.0, 2.0, 3.0)),
        "g" -> ModelFrame.factor(Vector("a", "b", "c"))
      )
      .getOrElse(fail("frame"))
    val formula = Formula.parse("y ~ I(x^2) + (1 | g)").getOrElse(fail("parse"))
    val materialized = formula.materialize(frame).getOrElse(fail("materialize"))
    assertEquals(materialized.numeric("I(x^2)"), Some(Vector(1.0, 4.0, 9.0)))
