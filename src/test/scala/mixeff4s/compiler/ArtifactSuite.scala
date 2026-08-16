package mixeff4s.compiler

import mixeff4s.fixtures.{Pastes, Penicillin, Sleepstudy}

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
    assertEquals(artifact.pathology.stratum.code, "easy")
    assertEquals(artifact.pathology.fitStatus.code, "not_assessed")
    assertEquals(artifact.toJson, SleepstudySnapshot.json)

  test("zerocorr sleepstudy is a diagonal theta map"):
    val artifact =
      Compiler.compile("reaction ~ 1 + days + (1 + days || subj)", Sleepstudy.frame).getOrElse(fail("compile"))
    assertEquals(artifact.nTheta, 2)
    assertEquals(artifact.randomTerms.head.covariance, "diagonal")
    assertEquals(artifact.parmap, Vector((0, 0, 0), (0, 1, 1)))
    assertEquals(artifact.thetaSlots.map(_.constraint), Vector("lower_bound_0", "lower_bound_0"))

  test("assessing sleepstudy ML theta is converged_interior"):
    val compiled =
      Compiler.compile("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame).getOrElse(fail("compile"))
    val assessed = Compiler.assess(compiled, Vector(0.9292297167514472, 0.01816466496782548, 0.22264601131030412))
    assertEquals(assessed.pathology.fitStatus.code, "converged_interior")
    assertEquals(compiled.pathology.fitStatus.code, "not_assessed")

  test("penicillin compiled-design snapshot"):
    val artifact =
      Compiler.compile("diameter ~ 1 + (1 | plate) + (1 | sample)", Penicillin.frame).getOrElse(fail("compile"))
    assertEquals(artifact.n, 144)
    assertEquals(artifact.p, 1)
    assertEquals(artifact.nReTerms, 2)
    assertEquals(artifact.nRanef, 30)
    assertEquals(artifact.nTheta, 2)
    assertEquals(artifact.randomTerms.map(_.grouping), Vector("plate", "sample"))
    assertEquals(artifact.parmap, Vector((0, 0, 0), (1, 0, 0)))
    assertEquals(artifact.pathology.stratum.code, "easy")
    assertEquals(artifact.pathology.minGroupSize, 6)
    assertEquals(artifact.pathology.maxGroupSize, 24)
    assertEquals(artifact.toJson, PenicillinSnapshot.json)
    val assessed = Compiler.assess(artifact, Vector(1.5375939045981573, 3.219792193110907))
    assertEquals(assessed.pathology.fitStatus.code, "converged_interior")

  test("pastes compiled-design snapshot"):
    val artifact =
      Compiler.compile("strength ~ 1 + (1 | batch / cask)", Pastes.frame).getOrElse(fail("compile"))
    assertEquals(artifact.n, 60)
    assertEquals(artifact.p, 1)
    assertEquals(artifact.nReTerms, 2)
    assertEquals(artifact.nRanef, 40)
    assertEquals(artifact.nTheta, 2)
    assertEquals(artifact.requestedFormula, "strength ~ 1 + (1 | batch) + (1 | batch:cask)")
    assertEquals(artifact.randomTerms.map(_.grouping), Vector("batch & cask", "batch"))
    assertEquals(artifact.parmap, Vector((0, 0, 0), (1, 0, 0)))
    assertEquals(artifact.pathology.stratum.code, "easy")
    assertEquals(artifact.pathology.minGroupSize, 2)
    assertEquals(artifact.pathology.maxGroupSize, 6)
    assertEquals(artifact.toJson, PastesSnapshot.json)
    val assessed = Compiler.assess(artifact, Vector(3.5269029347766856, 1.3299137410046242))
    assertEquals(assessed.pathology.fitStatus.code, "converged_interior")

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
      |  "pathology": {
      |    "contract_version": "v0.1",
      |    "fit_status": "not_assessed",
      |    "expected_statuses": ["converged_interior"],
      |    "stratum": "easy",
      |    "n": 180,
      |    "n_params": 6,
      |    "min_group_size": 10,
      |    "max_group_size": 10,
      |    "fe_rank": 2,
      |    "n_theta": 3,
      |    "structural_issue": null,
      |    "notes": []
      |  }
      |}""".stripMargin

private object PenicillinSnapshot:
  val json: String =
    """{
      |  "schema": {
      |    "name": "mixeff4s.compiled_design_artifact",
      |    "version": 1,
      |    "library_version": "0.1.0-SNAPSHOT"
      |  },
      |  "requested_formula": "diameter ~ 1 + (1 | plate) + (1 | sample)",
      |  "n": 144,
      |  "p": 1,
      |  "n_re_terms": 2,
      |  "n_ranef": 30,
      |  "n_theta": 2,
      |  "fe_names": ["(Intercept)"],
      |  "random_terms": [{"index": 0, "grouping": "plate", "n_levels": 24, "basis": ["(Intercept)"], "n_ranef": 24, "n_theta": 1, "covariance": "scalar"}, {"index": 1, "grouping": "sample", "n_levels": 6, "basis": ["(Intercept)"], "n_ranef": 6, "n_theta": 1, "covariance": "scalar"}],
      |  "theta_slots": [{"global_index": 0, "term_index": 0, "row": 0, "col": 0, "name": "theta[0:(Intercept),(Intercept)]", "constraint": "lower_bound_0"}, {"global_index": 1, "term_index": 1, "row": 0, "col": 0, "name": "theta[1:(Intercept),(Intercept)]", "constraint": "lower_bound_0"}],
      |  "parmap": [[0, 0, 0], [1, 0, 0]],
      |  "model_boundary": {
      |    "model_kind": "linear_mixed_model",
      |    "response_distribution": "gaussian",
      |    "link": "identity",
      |    "objective_approximation": "exact_gaussian",
      |    "inference_availability": "not_assessed"
      |  },
      |  "pathology": {
      |    "contract_version": "v0.1",
      |    "fit_status": "not_assessed",
      |    "expected_statuses": ["converged_interior"],
      |    "stratum": "easy",
      |    "n": 144,
      |    "n_params": 4,
      |    "min_group_size": 6,
      |    "max_group_size": 24,
      |    "fe_rank": 1,
      |    "n_theta": 2,
      |    "structural_issue": null,
      |    "notes": []
      |  }
      |}""".stripMargin

private object PastesSnapshot:
  val json: String =
    """{
      |  "schema": {
      |    "name": "mixeff4s.compiled_design_artifact",
      |    "version": 1,
      |    "library_version": "0.1.0-SNAPSHOT"
      |  },
      |  "requested_formula": "strength ~ 1 + (1 | batch) + (1 | batch:cask)",
      |  "n": 60,
      |  "p": 1,
      |  "n_re_terms": 2,
      |  "n_ranef": 40,
      |  "n_theta": 2,
      |  "fe_names": ["(Intercept)"],
      |  "random_terms": [{"index": 0, "grouping": "batch & cask", "n_levels": 30, "basis": ["(Intercept)"], "n_ranef": 30, "n_theta": 1, "covariance": "scalar"}, {"index": 1, "grouping": "batch", "n_levels": 10, "basis": ["(Intercept)"], "n_ranef": 10, "n_theta": 1, "covariance": "scalar"}],
      |  "theta_slots": [{"global_index": 0, "term_index": 0, "row": 0, "col": 0, "name": "theta[0:(Intercept),(Intercept)]", "constraint": "lower_bound_0"}, {"global_index": 1, "term_index": 1, "row": 0, "col": 0, "name": "theta[1:(Intercept),(Intercept)]", "constraint": "lower_bound_0"}],
      |  "parmap": [[0, 0, 0], [1, 0, 0]],
      |  "model_boundary": {
      |    "model_kind": "linear_mixed_model",
      |    "response_distribution": "gaussian",
      |    "link": "identity",
      |    "objective_approximation": "exact_gaussian",
      |    "inference_availability": "not_assessed"
      |  },
      |  "pathology": {
      |    "contract_version": "v0.1",
      |    "fit_status": "not_assessed",
      |    "expected_statuses": ["converged_interior"],
      |    "stratum": "easy",
      |    "n": 60,
      |    "n_params": 4,
      |    "min_group_size": 2,
      |    "max_group_size": 6,
      |    "fe_rank": 1,
      |    "n_theta": 2,
      |    "structural_issue": null,
      |    "notes": []
      |  }
      |}""".stripMargin
