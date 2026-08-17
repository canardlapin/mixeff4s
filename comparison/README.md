# Parity evidence harness

Frozen-first comparison files for mixeff4s. This is not AGENTS.md Phase 8.
Default `sbt test` reads these files. It does not call R or Julia.

| File | Role |
| --- | --- |
| `parity_scorecard.toml` | Class for every checked-in dataset/formula/estimator triple |
| `fixtures.toml` | Builtin and vendored dataset names the loader must materialize |
| `datasets/*/data.csv` | Small LMM frames copied from mixeff-rs |
| `datasets/*/meta.json` | Column types and factor level order |
| `schema/frozen_reference.schema.json` | Frozen numeric catalog |
| `frozen/references.json` | Pinned MixedModels.jl / mixeff-rs / lme4 numbers |

Classes are the rust vocabulary: `release_blocking_parity`,
`documented_divergence`, `unsupported_with_contract`, `stress_opt_in`,
`performance_known_slow`.

A scorecard `reference` of `mixedmodels.jl` or `mixeff-rs` is a claim
against that engine. Frozen `lme4` rows may exist without that claim.
Fast-PIRLS rows stay `documented_divergence` and are never sold as
`lme4::glmer`. GLMM objective constants are comparable only when the row
says so. There are no p-values.

The JVM harness gate fails if the scorecard, fixture catalog, and
on-disk datasets drift.
