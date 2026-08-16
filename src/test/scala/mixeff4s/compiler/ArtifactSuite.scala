package mixeff4s.compiler

import mixeff4s.fixtures.Sleepstudy

class ArtifactSuite extends munit.FunSuite:
  test("sleepstudy compiled-design snapshot"):
    val artifact =
      Compiler.compile("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame).getOrElse(fail("compile"))
    assertEquals(artifact.n, 180)
    assertEquals(artifact.p, 2)
    assertEquals(artifact.nReTerms, 1)
    assertEquals(artifact.nRanef, 36)
    assertEquals(artifact.nTheta, 3)
    assertEquals(artifact.feNames, Vector("(Intercept)", "days"))
    assertEquals(artifact.randomTerms.head.covariance, "full_cholesky")
    assertEquals(artifact.parmap, Vector((0, 0, 0), (0, 1, 0), (0, 1, 1)))
    assertEquals(artifact.pathology, "not_assessed")
    assertEquals(artifact.toJson, SleepstudySnapshot.json)

  test("zerocorr sleepstudy is a diagonal theta map"):
    val artifact =
      Compiler.compile("reaction ~ 1 + days + (1 + days || subj)", Sleepstudy.frame).getOrElse(fail("compile"))
    assertEquals(artifact.nTheta, 2)
    assertEquals(artifact.randomTerms.head.covariance, "diagonal")
    assertEquals(artifact.parmap, Vector((0, 0, 0), (0, 1, 1)))
    assertEquals(artifact.thetaSlots.map(_.constraint), Vector("lower_bound_0", "lower_bound_0"))

private object SleepstudySnapshot:
  val json: String =
    """{
      |  "schema": {
      |    "name": "mixeff4s.compiled_design_artifact",
      |    "version": 1,
      |    "library_version": "0.1.0-SNAPSHOT"
      |  },
      |  "requested_formula": "reaction ~ 1 + days + (1 + days | subj)",
      |  "n": 180,
      |  "p": 2,
      |  "n_re_terms": 1,
      |  "n_ranef": 36,
      |  "n_theta": 3,
      |  "fe_names": ["(Intercept)", "days"],
      |  "random_terms": [{"index": 0, "grouping": "subj", "n_levels": 18, "basis": ["(Intercept)", "days"], "n_ranef": 36, "n_theta": 3, "covariance": "full_cholesky"}],
      |  "theta_slots": [{"global_index": 0, "term_index": 0, "row": 0, "col": 0, "name": "theta[0:(Intercept),(Intercept)]", "constraint": "lower_bound_0"}, {"global_index": 1, "term_index": 0, "row": 1, "col": 0, "name": "theta[0:days,(Intercept)]", "constraint": "unconstrained"}, {"global_index": 2, "term_index": 0, "row": 1, "col": 1, "name": "theta[0:days,days]", "constraint": "lower_bound_0"}],
      |  "parmap": [[0, 0, 0], [0, 1, 0], [0, 1, 1]],
      |  "model_boundary": {
      |    "model_kind": "linear_mixed_model",
      |    "response_distribution": "gaussian",
      |    "link": "identity",
      |    "objective_approximation": "exact_gaussian",
      |    "inference_availability": "not_assessed"
      |  },
      |  "pathology": "not_assessed"
      |}""".stripMargin
