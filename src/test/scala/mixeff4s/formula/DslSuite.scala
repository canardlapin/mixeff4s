package mixeff4s.formula

import mixeff4s.formula.dsl.dsl.*

class DslSuite extends munit.FunSuite:
  test("typed DSL lowers to the same IR as the string parser"):
    val fromDsl =
      response("y") ~ (FixedTerm.Intercept + col("x") + (FixedTerm.Intercept | factor("g")))
    val fromString = Formula.parse("y ~ 1 + x + (1 | g)").getOrElse(fail("parse failed"))
    assertEquals(fromDsl.toString, fromString.toString)

  test("string interpolator"):
    val parsed = formula"y ~ 1 + x + (1 | g)"
    assertEquals(parsed.toString, "y ~ 1 + x + (1 | g)")
