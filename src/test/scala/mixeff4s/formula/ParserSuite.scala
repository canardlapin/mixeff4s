package mixeff4s.formula

class ParserSuite extends munit.FunSuite:
  test("fixed effects with explicit intercept"):
    val f = parse("y ~ 1 + x1 + x2")
    assertEquals(f.response, "y")
    assertEquals(
      f.fixedTerms,
      Vector(FixedTerm.Intercept, FixedTerm.Column("x1"), FixedTerm.Column("x2"))
    )
    assert(f.randomTerms.isEmpty)

  test("implicit intercept"):
    val f = parse("y ~ x1 + x2")
    assertEquals(f.fixedTerms.head, FixedTerm.Intercept)
    assertEquals(f.fixedTerms.length, 3)

  test("empty RHS is rejected"):
    assertEquals(Formula.parse("y ~"), Left(FormulaError.EmptyRhs))
    assertEquals(Formula.parse("y ~   "), Left(FormulaError.EmptyRhs))

  test("no intercept with zero"):
    val f = parse("y ~ 0 + x1 + x2")
    assert(!f.fixedTerms.contains(FixedTerm.Intercept))
    assert(!f.fixedTerms.contains(FixedTerm.NoIntercept))
    assertEquals(f.fixedTerms, Vector(FixedTerm.Column("x1"), FixedTerm.Column("x2")))

  test("no intercept with minus one"):
    val f = parse("y ~ -1 + x1")
    assert(!f.fixedTerms.contains(FixedTerm.Intercept))
    assertEquals(f.fixedTerms, Vector(FixedTerm.Column("x1")))

  test("colon interaction"):
    val f = parse("y ~ x1 + x2:x3")
    assert(f.fixedTerms.contains(FixedTerm.Interaction(Vector("x2", "x3"))))

  test("star expansion"):
    val f = parse("y ~ x1 * x2")
    assert(f.fixedTerms.contains(FixedTerm.Column("x1")))
    assert(f.fixedTerms.contains(FixedTerm.Column("x2")))
    assert(f.fixedTerms.contains(FixedTerm.Interaction(Vector("x1", "x2"))))

  test("three-way star expansion"):
    val f = parse("y ~ A * B * C")
    assertEquals(
      f.fixedTerms,
      Vector(
        FixedTerm.Intercept,
        FixedTerm.Column("A"),
        FixedTerm.Column("B"),
        FixedTerm.Column("C"),
        FixedTerm.Interaction(Vector("A", "B")),
        FixedTerm.Interaction(Vector("A", "C")),
        FixedTerm.Interaction(Vector("B", "C")),
        FixedTerm.Interaction(Vector("A", "B", "C"))
      )
    )

  test("nesting slash"):
    val f = parse("y ~ x1 / x2")
    assert(f.fixedTerms.contains(FixedTerm.Column("x1")))
    assert(f.fixedTerms.contains(FixedTerm.Interaction(Vector("x1", "x2"))))
    assert(!f.fixedTerms.contains(FixedTerm.Column("x2")))

  test("random intercept"):
    val f = parse("y ~ x1 + (1 | group)")
    assertEquals(f.randomTerms.length, 1)
    val rt = f.randomTerms.head
    assertEquals(rt.terms, Vector(FixedTerm.Intercept))
    assertEquals(rt.grouping, GroupingFactor.Single("group"))
    assert(!rt.zerocorr)

  test("random intercept and slope"):
    val f = parse("y ~ x1 + x2 + (1 + x1 | group)")
    assertEquals(
      f.randomTerms.head.terms,
      Vector(FixedTerm.Intercept, FixedTerm.Column("x1"))
    )

  test("random slope syntax has implicit intercept"):
    val f = parse("y ~ x1 + (x1 | group)")
    assertEquals(
      f.randomTerms.head.terms,
      Vector(FixedTerm.Intercept, FixedTerm.Column("x1"))
    )
    assertEquals(f.randomTerms.head.source.map(_.written), Some("(x1 | group)"))

  test("zerocorr"):
    val f = parse("y ~ x1 + (1 + x1 || group)")
    val rt = f.randomTerms.head
    assert(rt.zerocorr)
    assertEquals(rt.covariance, RandomCovariance.Diagonal)

  test("covariance wrappers parse as random terms"):
    val cases = Vector(
      ("us(1 + x1 | group)", RandomCovariance.Full, false),
      ("diag(1 + x1 | group)", RandomCovariance.Diagonal, false),
      ("cs(1 + x1 | group)", RandomCovariance.CompoundSymmetry, false),
      ("ar1(0 + x1 | group)", RandomCovariance.Ar1, false)
    )
    cases.foreach: (term, covariance, zerocorr) =>
      val f = parse(s"y ~ x1 + $term")
      assertEquals(f.randomTerms.length, 1)
      val rt = f.randomTerms.head
      assertEquals(rt.covariance, covariance)
      assertEquals(rt.zerocorr, zerocorr)
      assertEquals(rt.source.map(_.written), Some(term))

  test("crossed random effects"):
    val f = parse("y ~ x1 + (1 | g1) + (1 | g2)")
    assertEquals(f.randomTerms.map(_.grouping), Vector(GroupingFactor.Single("g1"), GroupingFactor.Single("g2")))

  test("interaction grouping"):
    val f = parse("y ~ x1 + (1 | g1 & g2)")
    assertEquals(f.randomTerms.head.grouping, GroupingFactor.Interaction(Vector("g1", "g2")))

  test("cell grouping with colon"):
    val f = parse("y ~ x1 + (1 | subject:item)")
    assertEquals(f.randomTerms.head.grouping, GroupingFactor.Cell(Vector("subject", "item")))
    assertEquals(f.toString, "y ~ 1 + x1 + (1 | subject:item)")

  test("nested grouping expands to main and cell"):
    val f = parse("y ~ x1 + (1 | school/class)")
    assertEquals(f.randomTerms.length, 2)
    assertEquals(f.randomTerms(0).grouping, GroupingFactor.Single("school"))
    assertEquals(f.randomTerms(1).grouping, GroupingFactor.Cell(Vector("school", "class")))
    assertEquals(f.toString, "y ~ 1 + x1 + (1 | school) + (1 | school:class)")
    assert(f.randomTerms.forall: term =>
      term.source.exists: source =>
        source.written == "(1 | school/class)" && source.expansion.contains(RandomTermExpansion.NestedGrouping)
    )

  test("crossed star grouping expands to main effects and cell"):
    val f = parse("y ~ x1 + (1 | subject*item)")
    assertEquals(f.randomTerms.length, 3)
    assertEquals(f.randomTerms(0).grouping, GroupingFactor.Single("subject"))
    assertEquals(f.randomTerms(1).grouping, GroupingFactor.Single("item"))
    assertEquals(f.randomTerms(2).grouping, GroupingFactor.Cell(Vector("subject", "item")))
    assertEquals(f.toString, "y ~ 1 + x1 + (1 | subject) + (1 | item) + (1 | subject:item)")

  test("only random term still has implicit intercept"):
    val f = parse("y ~ (1 | g)")
    assertEquals(f.fixedTerms, Vector(FixedTerm.Intercept))
    assertEquals(f.randomTerms.length, 1)

  test("whitespace handling"):
    val f1 = parse("y~1+x1+(1|g)")
    val f2 = parse("  y  ~  1 +  x1  + ( 1 | g )  ")
    assertEquals(f1.response, f2.response)
    assertEquals(f1.fixedTerms, f2.fixedTerms)
    assertEquals(f1.randomTerms.length, f2.randomTerms.length)

  test("dotted and underscored identifiers"):
    val dotted = parse("y.resp ~ x.pred + (1 | g.group)")
    assertEquals(dotted.response, "y.resp")
    assert(dotted.fixedTerms.contains(FixedTerm.Column("x.pred")))
    val underscored = parse("y_resp ~ x_pred + (1 | g_group)")
    assertEquals(underscored.response, "y_resp")

  test("backtick identifiers"):
    val f = parse("`reaction time` ~ `day of study` + (1 | `subject id`)")
    assertEquals(f.response, "reaction time")
    assert(f.fixedTerms.contains(FixedTerm.Column("day of study")))
    assertEquals(f.randomTerms.head.grouping, GroupingFactor.Single("subject id"))

  test("unterminated backtick is actionable"):
    Formula.parse("y ~ `oops + (1 | g)") match
      case Left(FormulaError.Other(msg)) =>
        assert(msg.contains("unterminated backtick"), clues(msg))
      case other =>
        fail(s"expected Other, got $other")

  test("empty backtick is rejected"):
    Formula.parse("y ~ `` + x") match
      case Left(FormulaError.Other(msg)) =>
        assert(msg.contains("empty backtick"), clues(msg))
      case other =>
        fail(s"expected Other, got $other")

  test("bare arithmetic is refused"):
    Vector("y ~ x^2", "y ~ x %in% g", "y ~ x > 0").foreach: src =>
      Formula.parse(src) match
        case Left(FormulaError.Other(msg)) =>
          assert(msg.contains("not supported") && msg.contains("precompute"), clues(src, msg))
        case other =>
          fail(s"expected Other for `$src`, got $other")

  test("stateful transforms are refused"):
    Vector(
      "y ~ poly(x, 2) + (1 | g)",
      "y ~ scale(x) + (1 | g)",
      "y ~ ns(x, 3) + (1 | g)",
      "y ~ factor(g2) + (1 | g)",
      "log(reaction, 2) ~ x + (1 | g)"
    ).foreach: src =>
      Formula.parse(src) match
        case Left(FormulaError.Other(msg)) =>
          assert(
            (msg.contains("stateless") || msg.contains("stateful")) &&
              (msg.contains("precompute") || msg.contains("host wrapper") || msg.contains("out of scope")),
            clues(src, msg)
          )
        case other =>
          fail(s"expected Other for `$src`, got $other")

  test("stateless subset parses with canonical labels"):
    val f = parse("log(reaction) ~ days + I(days^2) + (1 | subj)")
    assertEquals(f.response, "log(reaction)")
    assert(f.fixedTerms.contains(FixedTerm.Column("I(days^2)")))
    assert(f.fixedTerms.contains(FixedTerm.Column("days")))
    val labels = f.derived.map(_.label)
    assert(labels.contains("log(reaction)"))
    assert(labels.contains("I(days^2)"))
    assertEquals(f.derived.length, 2)
    Vector(
      "y ~ I(a*b) + (1|g)" -> "I(a*b)",
      "y ~ I(1/x) + (1|g)" -> "I(1/x)",
      "y ~ I(-x) + (1|g)" -> "I(-x)",
      "y ~ I( x  +  1 ) + (1|g)" -> "I(x+1)",
      "y ~ sqrt(I(x+1)) + (1|g)" -> "sqrt(I(x+1))",
      "y ~ I((a+b)*x) + (1|g)" -> "I((a+b)*x)"
    ).foreach: (src, label) =>
      val parsed = parse(src)
      assert(
        parsed.fixedTerms.contains(FixedTerm.Column(label)),
        clues(src, label, parsed.fixedTerms)
      )

  test("formula display canonicalizes terms"):
    val formula = Formula(
      "reaction",
      Vector(FixedTerm.Intercept, FixedTerm.Column("days")),
      Vector(
        RandomTerm(
          Vector(FixedTerm.Intercept, FixedTerm.Column("days")),
          GroupingFactor.Single("subj"),
          zerocorr = false,
          RandomCovariance.Full
        )
      )
    )
    assertEquals(formula.toString, "reaction ~ 1 + days + (1 + days | subj)")

  private def parse(src: String): Formula =
    Formula.parse(src) match
      case Right(formula) => formula
      case Left(error)    => fail(s"parse failed for `$src`: ${error.message}")
