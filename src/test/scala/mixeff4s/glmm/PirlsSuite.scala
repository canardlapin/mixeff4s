package mixeff4s.glmm

import mixeff4s.model.{Family, Link}

class PirlsSuite extends munit.FunSuite:
  test("Bernoulli logit working observation at eta=0"):
    val mu = Link.Logit.linkinv(0.0)
    val (sw, wy) = Pirls.workingObservation(Family.Bernoulli, Link.Logit, 1.0, 0.0, mu, 1.0)
    assertEqualsDouble(sw, 0.5, 1e-12)
    assertEqualsDouble(wy, 2.0, 1e-12)

  test("working observation subtracts a supplied offset"):
    val mu = Link.Logit.linkinv(0.0)
    val (_, wy) = Pirls.workingObservation(Family.Bernoulli, Link.Logit, 1.0, 0.0, mu, 1.0, 0.25)
    assertEqualsDouble(wy, 1.75, 1e-12)

  test("bounded Bernoulli mean stays off the 0/1 boundary"):
    val (mu, eta) = Pirls.boundedMeanAndEta(Family.Bernoulli, Link.Logit, 0.0, -40.0)
    assert(mu > 0.0 && mu < 1.0, clues(mu))
    assert(eta.isFinite, clues(eta))
