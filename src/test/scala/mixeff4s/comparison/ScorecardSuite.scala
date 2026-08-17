package mixeff4s.comparison

class ScorecardSuite extends munit.FunSuite:
  test("embedded scorecard uses the rust class vocabulary"):
    val scorecard = Scorecard.loadEmbedded
    assertEquals(scorecard.schemaVersion, "1.0.0")
    assertEquals(scorecard.classes, ScorecardClass.all)
    assertEquals(scorecard.rows.length, 20)
    assertEquals(
      scorecard.rows.map(_.classification).toSet,
      Set(ScorecardClass.ReleaseBlockingParity, ScorecardClass.DocumentedDivergence)
    )
    assert(scorecard.rows.exists(_.reference == "lme4"), clues(scorecard.rows.map(_.reference)))
    assert(
      !scorecard.rows.exists: row =>
        row.key.estimator == "fast_pirls" &&
          row.classification == ScorecardClass.ReleaseBlockingParity &&
          row.reference.startsWith("lme4"),
      clues(scorecard.rows.map(r => (r.key.estimator, r.classification, r.reference)))
    )

  test("dyestuff2 is a frozen lme4 boundary claim"):
    val rows = Scorecard.loadEmbedded.rows.filter(_.key.dataset == "dyestuff2")
    assertEquals(rows.map(_.reference).distinct, Vector("lme4_boundary"))
    assert(rows.forall(_.classification == ScorecardClass.ReleaseBlockingParity))

  test("contrast names align mixeff4s dummy labels to lme4"):
    assertEquals(FitCompare.contrastName("(Intercept)"), "(Intercept)")
    assertEquals(FitCompare.contrastName("Type: T2"), "TypeT2")
    assertEquals(FitCompare.contrastName("recipe: B:temperature: 185"), "recipeB:temperature185")
    assertEquals(FitCompare.contrastName("age:Sex: Female"), "age:SexFemale")

  test("sleepstudy zerocorr is a frozen MixedModels.jl claim"):
    val row = FitCompare.row("sleepstudy", "ML", "reaction ~ 1 + days + (1 + days || subj)")
    assertEquals(row.classification, ScorecardClass.ReleaseBlockingParity)
    assertEquals(row.reference, "mixedmodels.jl")
    val ref = FitCompare.claimed(row)
    assertEquals(ref.theta.map(_.length), Some(2))

  test("fast-PIRLS contraception is documented divergence"):
    val row = Scorecard.loadEmbedded.rows.find(_.key.dataset == "contraception").get
    assertEquals(row.classification, ScorecardClass.DocumentedDivergence)
    assertEquals(row.key.estimator, "fast_pirls")
    assertEquals(row.reference, "mixedmodels.jl_fast_pirls")

  test("every scorecard row has a claimed frozen reference"):
    val scorecard = Scorecard.loadEmbedded
    val frozen = FrozenCatalog.loadEmbedded
    val missing = scorecard.rows.filter(frozen.claimed(_).isEmpty).map(_.key.id)
    assertEquals(missing, Vector.empty)

  test("frozen catalog refuses p-values"):
    intercept[IllegalArgumentException]:
      FrozenCatalog.parse(
        """{"schema":{"name":"mixeff4s.frozen_reference","version":1},"source":{"note":"x"},"results":[{"dataset":"sleepstudy","formula":"y ~ 1","family":"Gaussian","link":"Identity","estimator":"ML","engine":"lme4","status":"ok","p_values":[0.01]}]}"""
      )
