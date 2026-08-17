package mixeff4s.comparison

class ScorecardSuite extends munit.FunSuite:
  test("embedded scorecard uses the rust class vocabulary"):
    val scorecard = Scorecard.loadEmbedded
    assertEquals(scorecard.schemaVersion, "1.0.0")
    assertEquals(scorecard.classes, ScorecardClass.all)
    assertEquals(scorecard.rows.length, 6)
    assertEquals(
      scorecard.rows.map(_.classification).toSet,
      Set(ScorecardClass.ReleaseBlockingParity, ScorecardClass.DocumentedDivergence)
    )
    assert(!scorecard.rows.exists(_.reference == "lme4"), clues(scorecard.rows.map(_.reference)))

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
