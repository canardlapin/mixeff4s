package mixeff4s.lmm

import mixeff4s.fixtures.{Pastes, Penicillin}

class CrossedFitSuite extends munit.FunSuite:
  test("penicillin crossed intercepts match MixedModels.jl"):
    val fit = Lmm
      .fit("diameter ~ 1 + (1 | plate) + (1 | sample)", Penicillin.frame, FitOptions.ml)
      .fold(err => fail(err.message), identity)
    assertEquals(fit.beta.length, 1)
    assertEquals(fit.theta.length, 2)
    assertEqualsDouble(fit.objective, 332.1883486700085, 1e-2)
    assertEqualsDouble(fit.beta(0), 22.97222222222222, 1e-3)
    assertEqualsDouble(fit.theta(0), 1.5375939045981573, 1e-2)
    assertEqualsDouble(fit.theta(1), 3.219792193110907, 1e-2)
    assertEqualsDouble(fit.stderror(0), 0.7446037806555799, 0.01)
    assertEquals(fit.varcorr.components.length, 2)

  test("penicillin REML matches mixeff-rs"):
    val fit = Lmm
      .fit("diameter ~ 1 + (1 | plate) + (1 | sample)", Penicillin.frame, FitOptions.reml)
      .fold(err => fail(err.message), identity)
    assertEqualsDouble(fit.objective, 330.86058899126897, 1e-2)
    assertEqualsDouble(fit.beta(0), 22.97222222222248, 1e-3)
    assertEqualsDouble(fit.theta(0), 1.5396773350745998, 1e-2)
    assertEqualsDouble(fit.theta(1), 3.51241122181154, 1e-2)
    assertEqualsDouble(fit.sigma, 0.549923173720829, 1e-3)

  test("pastes nested intercepts match MixedModels.jl"):
    val fit = Lmm
      .fit("strength ~ 1 + (1 | batch / cask)", Pastes.frame, FitOptions.ml)
      .fold(err => fail(err.message), identity)
    assertEquals(fit.theta.length, 2)
    assertEqualsDouble(fit.objective, 247.9944658624955, 1e-2)
    assertEqualsDouble(fit.beta(0), 60.0533333333333, 1e-3)
    assertEqualsDouble(fit.theta(0), 3.5269029347766856, 0.09)
    assertEqualsDouble(fit.theta(1), 1.3299137410046242, 0.09)
