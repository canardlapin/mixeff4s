package mixeff4s.stats

import mixeff4s.examples.LikelihoodInference

class InferenceWorkflowSuite extends munit.FunSuite:
  private lazy val result = LikelihoodInference.run().fold(e => fail(s"${e.code}: ${e.message}"), identity)
  test("documented inference workflow keeps coefficient information, null pairs and qualified profile endpoints"):
    assertEquals(result.comparison.attempts.length, 39)
    println(
      s"NULL_BOOTSTRAP qualified=${result.comparison.attempts.count(_.qualified)} calibration=${result.comparison.calibration}"
    )
    result.comparison.attempts.foreach: a =>
      println(
        s"NULL_PAIR index=${a.index} seed=${a.seed} null=${a.nullFit.left.toOption.map(_.code)} alternative=${a.alternativeFit.left.toOption.map(_.code)} statistic=${a.statistic.map(_.likelihoodRatio)}"
      )
    assert(result.comparison.calibration.isRight, clues(result.comparison.failedOrUnqualifiedIndices))
    assert(result.coefficients.intervals.forall(_.standardError > 0))
    assertEquals(result.profiles.length, 4)
    result.profiles.foreach: profile =>
      println(s"PROFILE_INTERVAL target=${profile.target} lower=${profile.lower} upper=${profile.upper}")
      assert(profile.lower.value.nonEmpty && profile.upper.value.nonEmpty, clues(profile))
