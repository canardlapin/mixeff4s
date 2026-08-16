package mixeff4s

import mixeff4s.prelude.*

class ReadmeExampleSuite extends munit.FunSuite:
  test("README quick start compiles against the public API"):
    val formula = Formula.parse("y ~ 1 + x + (1 | g)").getOrElse(fail("parse"))
    assertEquals(formula.toString, "y ~ 1 + x + (1 | g)")

    val frame = ModelFrame
      .of(
        "y" -> numeric(Vector(1.0, 2.1, 3.0, 4.2, 5.1, 6.0)),
        "x" -> numeric(Vector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)),
        "g" -> factorCol(Vector("a", "a", "b", "b", "c", "c"))
      )
      .getOrElse(fail("frame"))

    Lmm.fit(formula, frame, FitOptions.reml) match
      case Right(fit) =>
        assertEquals(fit.beta.length, 2)
        assert(fit.sigma > 0.0, clues(fit.sigma))
      case other =>
        fail(s"expected a certified fit, got $other")
