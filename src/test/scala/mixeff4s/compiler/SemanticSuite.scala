package mixeff4s.compiler

import mixeff4s.formula.GroupingFactor

class SemanticSuite extends munit.FunSuite:
  test("sleepstudy IR is a full Cholesky random slope"):
    val model = Compiler.semantic("reaction ~ 1 + days + (1 + days | subj)").getOrElse(fail("parse"))
    assertEquals(model.response, "reaction")
    assertEquals(model.fixedTerms, Vector("1", "days"))
    assertEquals(model.randomTerms.length, 1)
    val term = model.randomTerms.head
    assertEquals(term.id, "r0")
    assertEquals(term.grouping, GroupingFactor.Single("subj"))
    assertEquals(term.basis.map(_.kind), Vector(CoefficientKind.Intercept, CoefficientKind.Slope))
    assertEquals(term.covariance, CovarianceForm.Full)
    assertEquals(term.support, CovarianceSupport.Supported)
    assertEquals(term.intercept, InterceptPolicy.Included)
    assertEquals(term.sourceText, "(1 + days | subj)")
    assertEquals(term.written, None)

  test("zerocorr is one diagonal term, not a silent split"):
    val model = Compiler.semantic("reaction ~ 1 + days + (1 + days || subj)").getOrElse(fail("parse"))
    assertEquals(model.randomTerms.length, 1)
    assertEquals(model.randomTerms.head.covariance, CovarianceForm.Diagonal)
    assertEquals(model.randomTerms.head.support, CovarianceSupport.Supported)
    assertEquals(model.randomTerms.head.sourceText, "(1 + days || subj)")

  test("cs is typed IR and marked parsed_refused"):
    val model = Compiler.semantic("y ~ 1 + cs(1 + x | g)").getOrElse(fail("parse"))
    val term = model.randomTerms.head
    assertEquals(term.covariance, CovarianceForm.Structured("compound_symmetry"))
    assertEquals(term.support, CovarianceSupport.ParsedRefused)

  test("nested grouping keeps the written slash"):
    val model = Compiler.semantic("strength ~ 1 + (1 | batch / cask)").getOrElse(fail("parse"))
    assertEquals(model.randomTerms.map(_.grouping), Vector(GroupingFactor.Single("batch"), GroupingFactor.Cell(Vector("batch", "cask"))))
    assert(model.randomTerms.forall(_.written.contains("(1 | batch / cask)")), clues(model.randomTerms.map(_.written)))
    assertEquals(model.randomTerms.map(_.covariance), Vector(CovarianceForm.Scalar, CovarianceForm.Scalar))
