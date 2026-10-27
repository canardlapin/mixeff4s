**mixeff4s recovery assessment — 26 October 2026**

Continue the Scala implementation, with correctness and evidence contracts as the next milestone. The core is useful and worth preserving: both platforms pass the current suite, including the unfinished station work. It is not yet a certified counterpart of current Rust. Some important gaps predate the recent Rust improvements; adding faster optimization alone will not close them.

This assessment changes no production code or existing tests. The sequence below is a proposal for implementation. It does not authorize publication, replace the existing phases, or claim completion of the parity program.

**Assessed state and executed evidence**

| Component | Observed revision/state |
| --- | --- |
| Scala | `ec1937e2d2611c021519c4049c7eda49249f6311`, branch `cursor/semantic-model-ir`; seven modified parity files plus untracked `StationFitSuite.scala`, all preserved |
| Rust | Clean local `696e5e74623301280142e73cbb6bf7f2650647d7`, RC.5; 40 commits after the Scala frozen catalog's `d6b81ce8b5974b2f0e42fcf650ad6811e5202a4f` reference |
| R working directory | `bb352ff430eea4674e2c748eb40fa998f7ceb41a`, old August corpus branch; vendored engine `1f3f6892226e65dc96fb540e6b3806db0105389e`; unrelated untracked files preserved |
| Newer R branch inspected through Git | `codex/vignette-review-remediation` at `5557a756a6ccf25b4bc1471eda147ae19350242c`; vendor manifest pins Rust RC.5 `696e5e7`. Its worktree record points to a prunable `/private/tmp` checkout. No checkout switched |
| Tested dependency | Exact Gale pin `f869613cec0a89e57b6c995b0a02cf471ac7127c`; the newer local Gale checkout was deliberately not substituted |
| Runtime | sbt 1.12.14, Scala 3.7.4, sbt-selected Homebrew Java 25.0.1, Node 26.7.0 |

`sbt -batch test` passed **163 JVM tests and 159 Scala.js tests**, including architecture and frozen comparisons. The first launch was blocked by sandbox access to `~/.sbt/boot/sbt.boot.lock`; the authorized retry passed. The shell's `java -version` reported Java 22, so the sbt log, not the shell default, identifies the tested JVM.

Five additional JVM observational probes confirmed the defects below. These probes intentionally assert current bad behavior: a passing probe is evidence of a defect, not an acceptance test. The initial four-probe run had one failed hypothesis because a separate RE-nesting bug intercepted the intended test; the revised five-probe run isolates both cases. Both logs are retained.

Evidence is in [assessment-2026-10-26](assessment-2026-10-26/): [baseline log](assessment-2026-10-26/baseline-authorized.log), [probe log](assessment-2026-10-26/probes-v2.log), [probe source](assessment-2026-10-26/AssessmentProbeSuite.scala), and [revision/hash/exit receipt](assessment-2026-10-26/receipt.json). Replay from the repository root:

```sh
sbt -batch test
sbt -batch \
  'set mixeff4sJVM / Test / unmanagedSources += file("comparison/assessment-2026-10-26/AssessmentProbeSuite.scala")' \
  'mixeff4sJVM/testOnly mixeff4s.assessment.AssessmentProbeSuite'
```

**What works, and what the evidence does not establish**

| Claim | Implementation and evidence | Implication |
| --- | --- | --- |
| Formula/frame/design → profiled LMM | Connected; frozen ML/REML, nested/crossed intercepts, slopes, categorical slopes, diagonal covariance and boundary examples pass on JVM/JS | Preserve parser, DSL, Gale design, explicit parmap, blocked PLS and fixtures |
| Labelled fast-PIRLS | Connected; contraception fit and fixed-θ reference tests pass; unsupported joint/AGQ paths refuse | Preserve the approximation label; this is not joint-lme4 parity or broad GLMM family qualification |
| Certified fit or precise diagnostic | Partial and contradicted by probes and fit-path inspection | Highest-priority repair; successful numerical return currently carries too little evidence |
| Coefficient inference alignment | Ordinary full-rank examples pass; rank-deficient alignment is contradicted | Repair before expanding inference APIs |
| Nested LRT | Some reference examples pass; semantic guards fail in both directions | Refuse unverified comparisons until actual nesting can be established |
| Profile/bootstrap and broader R workflows | Profile/bootstrap are explicit refusal stubs; no corresponding public prediction/refit/simulation APIs were found in the inspected Scala surfaces | Original phase completion does not imply functional parity with the R package |
| Efficiency | Array-backed mutable kernels exist; sparse nested off-diagonal blocks, analytic gradient and reusable refit path are missing | Plausible optimization opportunities, not measured speed claims |

**Correctness findings that change the next milestone**

1. **RC-01 — convergence claims lack convergence evidence.** `Pls.fit` (`lmm/Pls.scala:636`) and `Glmm.fit` (`glmm/Glmm.scala:80`) discard TrustBQ's stop reason, including exhausted evaluation budgets. `GlmmWorkspace.devianceAt` (`Glmm.scala:197`) ignores the Boolean returned by `pirls()` and converts errors to finite `1e300`; final GLMM acceptance only checks finiteness. The public `profiledDeviance` probe with a wrong-length θ returns `Right(1e300)`. Separately, `Pathology.classifyTheta` (`pathology/Pathology.scala:159`) labels an arbitrary unfitted θ of `12345` as `ConvergedInterior`. Shape/boundary classification cannot establish convergence. The ignored stop/PIRLS results are source-confirmed; a public fit exhausting its budget was not reproduced in this assessment.

2. **RC-02 — active coefficients and covariance use different layouts.** `PlsWorkspace.beta` and `feNames` expose retained columns, while `vcov` calls `unpivotVcov` and expands to the original columns (`Pls.scala:397–458`). For sleepstudy with an all-zero middle predictor, names are `(Intercept), days` but SEs are `[9.5062, NaN, 0.801735]`. `CoefTable.wald` zips them, assigning `NaN` to the retained days coefficient. Choose one public layout and carry explicit alias metadata; do not let zip silently truncate mismatched arrays.

3. **RC-03 — LRT checks names instead of design spaces and mishandles RE intercepts.** `stats/Lrt.scala:42–83` checks response equality and FE name subsets, not predictor identity/column spaces. A fit with `days` containing squared values and a fit with ordinary `days` plus `cube` are accepted as nested and produce `p = 0.0131547`, though their fixed-effect spaces are not nested. Conversely, the genuinely nested `(1 | subj)` → `(1 + days | subj)` comparison is refused: `columnsOf` drops the intercept when other columns exist. Retain enough design provenance to verify nesting or return a precise refusal. Verifying structural nesting still does not by itself justify the chi-square reference law for every variance-boundary test.

**Rust changes worth porting, in dependency order**

| Work | Rust evidence | Scala decision and acceptance gate |
| --- | --- | --- |
| Fit/derivative evidence lifecycle | `5967b20`, `977ef11`; `compiler/audit.rs`, `generalized/certify.rs`, `generalized/optimizer.rs::finalize_theta_after_optimizer` | Port distinctions between not assessed, unavailable and failed. Budget exhaustion, unfinished PIRLS and invalid probes must not produce certified convergence; Hessian failure must not erase valid gradient evidence. Final-PIRLS checking already existed before August: this is also a missed port contract |
| Sparse nested blocks | `44e534a`; `linear/blocks.rs::subtract_sparse_sparse_t_in_pattern`, block matvecs | Preserve sparsity where valid, with complete pattern-coverage preflight before mutation and dense fallback for fill-in. Check blocked factors/objectives against independently assembled dense Cholesky for nested, partially nested, crossed and vector terms |
| Weighted cross-product pattern reuse | `6f06369`; `ScalarCrossPattern` | Cache topology, refresh values. Compare to full rebuild under changing weights, including zero→nonzero transitions |
| Analytic profiled LMM gradient | `a303a32`, `80e767e`; `linear/gradient.rs`, `docs/profiled_deviance_gradient.md` | Qualify as a separate kernel before making it an optimizer dependency. Compare with independent dense trace/envelope calculations and finite differences for ML/REML, boundaries and multiple RE terms. Avoid division through singular Λ; bound selected-inverse memory |
| Gradient-driven TrustBQ | `f636282` plus `977ef11`; `optimizer/trust_bq.rs` | Add an explicit strategy with projected-gradient stopping and the shifted-Newton original-quadratic ray correction. Gate nested `pw2_nested`, crossed, boundary and maximal-singular cases. Retain a specified derivative-free path |
| Refit and response workspace | `84c87c4`, `3096ee4`, `cfd1710`, `897dd25` | Exact response-only block refresh first; explicit warm starts second. Compare refresh/full rebuild and warm/cold objectives; preserve optimizer controls; bootstrap starts must not depend accidentally on replicate order |
| Deferred certification | `217e841`, `3f93ffd`, `15d374f`, `5967b20` | Only after evidence semantics work. Check first/cached inspection equality, invalidation after refit, and no fitted-mode mutation during inspection |
| Joint Laplace GLMM | `25db276`, `0fbfd83`, `7fd1360`, `977ef11`; `generalized/joint.rs`, `pirls.rs` | Separate estimator milestone: optimize β and θ, verify conditional modes with final factors, preserve objective constants and fallback provenance. Require an independent scalar conditional-mode oracle and frozen joint-lme4 rows before changing labels. AGQ remains separate |

Do not copy every Rust diff. Scala already returns the original dense matrix from `MatrixBlock.asDense` and already uses observation-outer FE×RE rebuild loops, so those Rust optimizations are not missing here. The recent shifted-Newton ray correction is enabled in Rust's gradient-driven path, not its derivative-free path. RC.5's conditional-score PIRLS repair is guarded by fixed β in the joint estimator; applying it wholesale to fast-PIRLS changes the current approximation. The normal Wald-tail fix is inapplicable while Scala deliberately refuses Wald p-values.

**Idiomatic, typed, efficient Scala direction**

Keep immutable public results, algebraic error/option types, one formula IR, and tightly scoped mutable array workspaces. The current enums, `Either` results, DSL and package separation are good foundations. A functional collection rewrite of numerical loops would add allocation without improving their contract.

Strengthen types at the boundaries that currently lose information: a fitted result should carry optimization/evidence status; θ should be validated against its parameter map; coefficients and covariance should share a validated layout; diagnostic availability should be an enum with reasons. The present DSL types syntax, not column existence or matrix dimensions. Public case-class constructors and raw vectors still admit inconsistent values. Hide unchecked construction where it protects an actual invariant; a wholesale type-level-dimension rewrite is not required.

Profile before changing storage. Concrete allocation candidates include `ReMat.withTheta` rebuilding immutable λ matrices, `lambdaWork` recreating work matrices inside factor updates, and off-diagonal RE blocks being dense. Keep Gale pinned during recovery. Compare end-to-end construction, cold/first/warm fit, repeated response fit, inspection and peak allocation separately on JVM and JS; include matched numerical tolerances and estimator settings. No speed ratio was measured here.

**Recovery sequence and existing Mote alignment**

| Order | Deliverable | Exit gate / existing work |
| --- | --- | --- |
| 0 | Preserve and review the station checkpoint; record source baselines; reconcile stale board descriptions | Current JVM/JS baseline is green. `parity-lmm-station` is closed but its files are uncommitted. `parity-zerocorr` still says compile-only although its fit/row children are closed and tests pass. Do not redo those implementations |
| 1 | Repair RC-01, RC-02 and RC-03 with targeted regression tests | Invalid/unfinished work cannot become a certified fit; retained coefficient/covariance layouts agree; nonnested comparisons refuse; valid supported comparisons work. Preserve all frozen rows, revising only labels that lacked evidence |
| 2 | Establish matched performance baselines, then qualify sparse/rebuild kernels and analytic gradients | Independent dense numerical witnesses, boundary tests, memory bounds; JVM/JS regression gates. Enable gradient-driven optimization only after gradient qualification |
| 3 | Expand useful LMM operations through an immutable compiled plan and private reusable workspace | Cold/warm refit equivalence and control preservation before prediction/simulation/bootstrap/profile breadth. Keep unavailable inference explicit |
| 4 | Broaden fast-PIRLS evidence, then implement joint Laplace separately | Reuse `parity-glmm-fast`, family children and `parity-glmm-joint` / `parity-joint-design`; align frames, links, weights, constants and estimator. AGQ follows its own gate |
| 5 | Broaden stress, live comparisons and release integration | Reuse pathology/stress/live epics; Maven/hosted CI is a later publication decision, not inferred from local tests |

The first implementation slice should be RC-01: preserve optimizer and final PIRLS outcomes, separate covariance geometry from convergence, and make final acceptance reject invalid penalties. Then repair the coefficient layout and LRT contract before expanding numerical claims. The existing parity backlog alone does not cover these prerequisites.

The strongest reason to start performance work sooner would be a demonstrated representative workload that is unusably slow. There is no such current Scala measurement in this assessment. Independent kernel work could proceed alongside contract repairs, but should not be promoted as a certified public fit path until the repairs pass.

This is a local implementation assessment. Rust/R suites, fresh Julia/lme4 reference generation, hosted CI, remote branch parity, publishing, JVM/JS performance and statistical calibration were not executed or claimed. Frozen station boundary agreement includes the pre-existing `sigma_abs_tol = 0.03` allowance; it does not establish optimizer certification. The scalar/categorical examples do not establish universal formula compatibility. Production changes remain to be implemented.
