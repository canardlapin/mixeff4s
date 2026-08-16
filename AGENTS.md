# mixeff4s

Scala 3 port of [mixeff-rs](https://github.com/bbuchsbaum/mixeff-rs): linear and
generalized linear mixed-effects models with an explicit inference contract.
The numerical reference is MixedModels.jl; lme4 is the later release scorecard.

This is a re-expression, not a line-for-line translation. Fit is a value.
Unavailable inference is a refusal with a stable code. Do not invent p-values.

## Layout

Single JVM project, many packages. Split artifacts when GLMM or the compiler
has a reason to exist independently.

```
mixeff4s.error      MixedModelError, LinAlgError, FitResult
mixeff4s.formula    AST, string parser, I() transforms, typed DSL
mixeff4s.data       ModelFrame — fitting substrate, not a dataframe library
mixeff4s.model      Family, Link
mixeff4s.design     FeTerm, FeMat, ReMat, parmap (Gale)
mixeff4s.linalg     WorkMat, MatrixBlock, blocked Cholesky
mixeff4s.optimizer  TrustBQ
mixeff4s.lmm        FitOptions, Lmm.compile, Lmm.fit (single-term PLS)
```

Layer rules, enforced by `ArchitectureSuite`:

- `error` does not import `data`, `lmm`, or `model`
- `formula` does not import `lmm` or `model`
- `data` does not import `formula`, `lmm`, or `model`
- `model` does not import `lmm`, `formula`, or `data`
- `linalg` does not import `lmm`, `formula`, `data`, `model`, `design`, or `optimizer`
- `optimizer` does not import `lmm`, `design`, `formula`, `data`, or `model`
- `lmm` does not import a future `compiler` or `stats` package

## Build

```sh
sbt -Dmixeff4s.gale.build=/path/to/gale test
sbt scalafmtAll
```

Scala 3.7.4, sbt 1.12.14. Gale is pinned at `galeRevision` in `build.sbt`. A local checkout is admitted
only through `-Dmixeff4s.gale.build=/path/to/gale`. Do not add Cats Effect in
core.

Issue tracking is **mote**. Check `mote ready`, reserve paths, and publish
through the CLI. Do not hand-edit `.mote/ops`.

## Port phases

| Phase | Scope | Gate |
| --- | --- | --- |
| 0–1 | Language, frame, stub `Lmm.fit` | Parser + architecture tests |
| 2 | `FeTerm` / `ReMat` / `parmap` on Gale | Design shapes on sleepstudy |
| 3 | Blocked-Cholesky PLS + TrustBQ | sleepstudy β, θ, σ, objective (current) |
| 4 | Public LMM API, VarCorr, Wald | Honest summaries |
| 5 | GLMM fast-PIRLS, labelled | Scorecard fast-PIRLS rows |
| 6 | LRT / profile / bootstrap | Refusal contracts |
| 7 | Compiler / pathology | JSON snapshots |

Keep index-parallel PLS loops close to mixeff-rs until parity is green.

## Invariants

1. No fake p-values.
2. Certified fit or precise diagnostic.
3. lme4 formula semantics; Normal+Identity is an LMM, refused as a GLMM.
4. θ is Cholesky of relative covariance with an explicit parmap.
5. Profiled (RE)ML: the optimizer sees θ only.
6. GLMM `fast=true` is labelled fast-PIRLS, never sold as `lme4::glmer`.
7. Unsupported optimizer or covariance → error, not silent fallback.
8. Layer acyclicity.
