package mixeff4s.design

import mixeff4s.data.ModelFrame
import mixeff4s.formula.Formula

class UpstreamDesignSuite extends munit.FunSuite:
  private val frame = ModelFrame
    .of(
      "y" -> ModelFrame.numeric(Vector(1.0, 2.0, 3.0, 4.0)),
      "f" -> ModelFrame.factor(Vector("a", "a", "b", "b")),
      "h" -> ModelFrame.factor(Vector("u", "v", "u", "v")),
      "x" -> ModelFrame.numeric(Vector(1.0, 2.0, 3.0, 4.0)),
      "g" -> ModelFrame.factor(Vector("g1", "g1", "g2", "g2"))
    )
    .getOrElse(fail("frame"))

  test("no-intercept main effects full-code only the first categorical variable"):
    val design = compile("y ~ 0 + f + h + (1 | g)")
    assertEquals(design.fe.fixedEffects.requestedNames, Vector("f: a", "f: b", "h: v"))

  test("no-intercept factor interaction full-codes both factors"):
    val design = compile("y ~ 0 + f:h + (1 | g)")
    assertEquals(
      design.fe.fixedEffects.requestedNames,
      Vector("f: a:h: u", "f: b:h: u", "f: a:h: v", "f: b:h: v")
    )

  test("interaction contrasts follow marginality and formula variable order"):
    val design = compile("y ~ h + f:h + (1 | g)")
    assertEquals(
      design.fe.fixedEffects.requestedNames,
      Vector("(Intercept)", "h: v", "h: u:f: b", "h: v:f: b")
    )

  test("numeric by factor interaction uses full factor indicators when uncovered"):
    val design = compile("y ~ 0 + x:f + (1 | g)")
    assertEquals(design.fe.fixedEffects.requestedNames, Vector("x:f: a", "x:f: b"))

  test("random-slope bases use the same no-intercept contrast assignment"):
    val design = compile("y ~ 1 + (0 + f + h | g)")
    assertEquals(design.reterms.head.cnames, Vector("f: a", "f: b", "h: v"))

  test("random no-intercept interaction full-codes both factor bases"):
    val design = compile("y ~ 1 + (0 + f:h | g)")
    assertEquals(
      design.reterms.head.cnames,
      Vector("f: a:h: u", "f: b:h: u", "f: a:h: v", "f: b:h: v")
    )

  test("grouping drops unused declared levels in their declared order"):
    val grouped = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(1.0, 2.0)),
        "x" -> ModelFrame.numeric(Vector(1.0, 2.0)),
        "g" -> ModelFrame
          .factor(Vector("b", "d"), Vector("a", "b", "c", "d"))
          .getOrElse(fail("factor"))
      )
      .getOrElse(fail("frame"))
    val design = compile("y ~ x + (1 | g)", grouped)
    assertEquals(design.reterms.head.levels, Vector("b", "d"))
    assertEquals(design.reterms.head.refs, Vector(0, 1))

  test("combined grouping retains distinct tuples with underscore-containing levels"):
    val grouped = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(1.0, 2.0)),
        "x" -> ModelFrame.numeric(Vector(1.0, 2.0)),
        "g1" -> ModelFrame.factor(Vector("a_b", "a")),
        "g2" -> ModelFrame.factor(Vector("c", "b_c"))
      )
      .getOrElse(fail("frame"))
    val design = compile("y ~ x + (1 | g1:g2)", grouped)
    assertEquals(design.reterms.head.nLevels, 2)

  test("expanded fixed columns agree with independent R model.matrix numeric witnesses"):
    val witnesses = Vector(
      "y ~ 0 + f + h + (1|g)" -> Vector(
        Vector(1.0, 1.0, 0.0, 0.0),
        Vector(0.0, 0.0, 1.0, 1.0),
        Vector(0.0, 1.0, 0.0, 1.0)
      ),
      "y ~ 0 + f:h + (1|g)" -> Vector(
        Vector(1.0, 0.0, 0.0, 0.0),
        Vector(0.0, 0.0, 1.0, 0.0),
        Vector(0.0, 1.0, 0.0, 0.0),
        Vector(0.0, 0.0, 0.0, 1.0)
      ),
      "y ~ h + f:h + (1|g)" -> Vector(
        Vector(1.0, 1.0, 1.0, 1.0),
        Vector(0.0, 1.0, 0.0, 1.0),
        Vector(0.0, 0.0, 1.0, 0.0),
        Vector(0.0, 0.0, 0.0, 1.0)
      ),
      "y ~ 0 + x:f + (1|g)" -> Vector(Vector(1.0, 2.0, 0.0, 0.0), Vector(0.0, 0.0, 3.0, 4.0))
    )
    witnesses.foreach: (formula, expected) =>
      val fe = compile(formula).fe
      val actual = fe.piv.indices.map: original =>
        val working = fe.piv.indexOf(original)
        Vector.tabulate(fe.n)(i => fe.x(i, working))
      assertEquals(actual.toVector, expected)

  private def compile(source: String, data: ModelFrame = frame): CompiledDesign =
    val formula = Formula.parse(source).getOrElse(fail(s"parse: $source"))
    Design.compile(formula, data).getOrElse(fail(s"compile: $source"))
