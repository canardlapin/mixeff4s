package mixeff4s.design

import mixeff4s.fixtures.Sleepstudy
import mixeff4s.formula.Formula
import mixeff4s.lmm.Lmm

class DesignSuite extends munit.FunSuite:
  test("sleepstudy design shapes match mixeff-rs"):
    val formula = Formula.parse("reaction ~ 1 + days + (1 + days | subj)").getOrElse(fail("parse"))
    val design = Lmm.compile(formula, Sleepstudy.frame).getOrElse(fail("compile"))
    assertEquals(design.n, 180)
    assertEquals(design.p, 2)
    assertEquals(design.nReTerms, 1)
    assertEquals(design.nRanef, 36)
    assertEquals(design.nTheta, 3)
    assertEquals(design.fe.rank, 2)
    assertEquals(design.fe.piv, Vector(0, 1))
    assertEquals(design.fe.fullRankNames, Vector("(Intercept)", "days"))
    assertEquals(design.xy.xy.rows, 180)
    assertEquals(design.xy.xy.cols, 3)
    val re = design.reterms.head
    assertEquals(re.groupingName, "subj")
    assertEquals(re.nLevels, 18)
    assertEquals(re.vsize, 2)
    assertEquals(re.cnames, Vector("(Intercept)", "days"))
    assertEquals(re.z.rows, 2)
    assertEquals(re.z.cols, 180)
    assertEquals(re.adjA.rows, 36)
    assertEquals(re.adjA.cols, 180)
    assertEquals(design.parmap, Vector((0, 0, 0), (0, 1, 0), (0, 1, 1)))
    assertEquals(design.theta, Vector(1.0, 0.0, 1.0))
    assert(Sleepstudy.frame.factor("subj").exists(_.levels == Sleepstudy.subjects))
    val counts = re.levels.map(level => re.refs.count(_ == re.levels.indexOf(level)))
    assert(counts.forall(_ == 10), clues(counts))

  test("zerocorr sleepstudy has two theta slots"):
    val design = Lmm
      .compile("reaction ~ 1 + days + (1 + days || subj)", Sleepstudy.frame)
      .getOrElse(fail("compile"))
    assertEquals(design.nTheta, 2)
    assertEquals(design.parmap, Vector((0, 0, 0), (0, 1, 1)))
    assertEquals(design.theta, Vector(1.0, 1.0))
