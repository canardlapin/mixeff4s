# mixeff4s

[Repository](https://github.com/canardlapin/mixeff4s)

mixeff4s is a mixed-effects modeling library for Scala 3 that turns an lme4-style
formula and a column store into either a certified fit or a typed refusal. Use
it when you want Bates/MixedModels PLS semantics without inventing a p-value.

> This is pre-release development software. APIs, semantics, and package
> boundaries may change. No stable support or binary-compatibility promise
> should be inferred unless the package explicitly states one.

**Status:** `prototype`. The formula language, model frame, family/link types,
Gale-backed design compilation, and a blocked-Cholesky LMM kernel are
implemented for single-term, nested, and crossed random intercepts.

## Quick start

The library is source-only. From this checkout:

```scala
import mixeff4s.prelude.*

val formula = Formula.parse("y ~ 1 + x + (1 | g)").toOption.get
formula.toString
// y ~ 1 + x + (1 | g)

val frame = ModelFrame.of(
  "y" -> numeric(Vector(1.0, 2.1, 3.0, 4.2, 5.1, 6.0)),
  "x" -> numeric(Vector(0.0, 1.0, 0.0, 1.0, 0.0, 1.0)),
  "g" -> factorCol(Vector("a", "a", "b", "b", "c", "c"))
).toOption.get

val fit = Lmm.fit(formula, frame, FitOptions.reml).toOption.get
(fit.beta, fit.sigma, fit.objective)
```

The same formula can be built without a string:

```scala
import mixeff4s.formula.dsl.dsl.*

val spec = response("y") ~ (intercept + col("x") + (intercept | factor("g")))
```

`cs(...)` and `ar1(...)` parse, then `Lmm.fit` refuses them. That is intentional:
a number means what it says, otherwise you get a matchable reason code.

## What it covers

- Parse the supported lme4 subset: intercepts, `:` / `*` / `/`, `(re | g)`,
  `(re || g)`, nested and crossed grouping, `us` / `diag` / `cs` / `ar1`
- Lower a typed DSL into the same formula IR
- Evaluate the stateless `I(...)` / `log` / `sqrt` transform subset onto a frame
- Name GLMM families and links, including the rule that Normal+Identity is an LMM
- Compile a formula and frame into `FeTerm`, `ReMat`s, and the live θ `parmap`
- Fit profiled (RE)ML with blocked Cholesky + TrustBQ, including crossed intercepts
- Report VarCorr, Wald SEs, and z-statistics; p-values stay a typed refusal

## Fit and boundaries

Good fit for describing mixed models and for fitting the current LMM kernel
(sleepstudy, penicillin, pastes). Structured RE covariance and GLMMs are not
implemented yet.

Out of scope for this slice, and for the Rust 1.0 line we are porting:

- Multivariate responses
- Residual AR/spatial correlation
- Structured RE covariance as a fitted family
- GLMM profile likelihood
- A published Maven coordinate

Design compilation depends on [Gale](https://github.com/canardlapin/gale)
through an exact source pin in `build.sbt`. A local checkout is admitted only
via `-Dmixeff4s.gale.build=/path/to/gale`.

## Maturity and verification

| Question | Answer |
| --- | --- |
| What does it do? | Mixed-model language, frame, design compilation, profiled LMM fitting, and honest summaries |
| Smallest useful example? | Parse `y ~ 1 + x + (1 \| g)` and inspect the IR |
| Maturity / publication? | `prototype`, source-only, unpublished at [canardlapin/mixeff4s](https://github.com/canardlapin/mixeff4s) |
| How to verify? | `sbt -Dmixeff4s.gale.build=/path/to/gale test` |

Platform: JVM. Scala 3.7.4, sbt 1.12.14. Scala.js is deferred until the kernel
is Gale-only.

This is a Scala 3 re-expression of
[mixeff-rs](https://github.com/bbuchsbaum/mixeff-rs), itself an independent
implementation of MixedModels.jl's PLS/PIRLS formulation. It is not a
line-for-line Rust translation.

## Documentation

- [AGENTS.md](AGENTS.md) — package layout, layering, and port phases
- mixeff-rs guide — [what is supported](https://github.com/bbuchsbaum/mixeff-rs/blob/main/docs/guide/05_what_is_supported.md)

## License

Apache-2.0. The Rust reference crate is MIT.
