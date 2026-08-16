package mixeff4s.model

class FamilySuite extends munit.FunSuite:
  test("canonical links"):
    assertEquals(Family.Bernoulli.canonicalLink, Link.Logit)
    assertEquals(Family.Poisson.canonicalLink, Link.Log)
    assertEquals(Family.Normal.canonicalLink, Link.Identity)

  test("Normal+Identity is not a supported GLMM pair"):
    assert(!Family.Normal.allows(Link.Identity))
    assert(Family.Bernoulli.allows(Link.Logit))

  test("logit is inverse of logistic"):
    val eta = 0.5
    val mu = Link.Logit.linkinv(eta)
    assertEqualsDouble(Link.Logit.link(mu), eta, 1e-12)

  test("probit is centered at zero"):
    assertEqualsDouble(Link.Probit.linkinv(0.0), 0.5, 1e-8)
    assertEqualsDouble(Link.Probit.link(0.5), 0.0, 1e-8)
