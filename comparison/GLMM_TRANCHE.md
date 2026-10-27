# Qualified canonical GLMM tranche

The local GLMM tranche is complete: explicit response/trial semantics, deterministic fast-PIRLS evaluation with accepted-mean factor refresh, a separately labelled canonical joint-Laplace estimator, independent numerical references, a runnable workflow, and full JVM/Scala.js regression qualification. Mote authority is `workflow-glmm`; this closes the binomial-trials, joint-design, cbpp-joint and culcita-Laplace slices. The broader GLMM corpus and AGQ epic remain open.

## Delivered contract

Joint v1 supports Bernoulli/logit, binomial/logit and Poisson/log with `GlmmOptions.jointLaplace` (`fast=false`, `nAgq=1`). For each outer `[beta; theta]`, beta remains fixed while the inner solver minimizes conditional minus twice log likelihood plus `u'u`. Only afterward is the Laplace determinant added. A dedicated private workspace reuses the blocked PLS random-effect factor and supplied-beta solve, without factoring or reading the unnecessary fixed/response Schur block.

Each objective probe starts with `u=0`. It uses actual canonical-link score and curvature, conditional-only line search, and an explicit maximum absolute spherical score tolerance of `1e-8`. A failed or exhausted solve returns its typed diagnostic. There are no arbitrary finite penalties, automatic retries, or silent fast-PIRLS fallback. Final curvature and determinant are recomputed at the accepted modes, including after step halving. Fast-PIRLS retains its separate historical objective-change rule, and now resets fixed-effect starts for deterministic objective evaluation and refreshes final accepted-mean factors.

Binomial y is a proportion. `trials: Option[Vector[Int]]` supplies positive integer trials; without it only binary responses are admissible. Successes must be integral within absolute count tolerance `1e-8`; accepted near-integer values are canonicalized once. Bernoulli and Poisson's existing count tolerances are canonicalized too, so likelihood, score, initialization and working response use the same response. Poisson counts must be nonnegative integers. Separate positive case weights multiply the entire log density, score and curvature; noninteger weights define a power likelihood. Offsets are explicit finite link-scale vectors in `GlmmOptions.offset`, including `log(exposure)` for Poisson.

The joint objective retains binomial combinatorial and Poisson factorial constants. Stable softplus, logistic inverse/derivative and endpoint residual formulas avoid overflow or loss of representable curvature. Degenerate combinatorial constants and log-Gamma at one/two are exactly zero, preventing constant approximation error from being amplified by large case weights. No response constant is used to scale assessment tolerances.

`GlmmFit` has read-only fields, package-private construction and no public copy. `GlmmInnerResult` distinguishes fast objective-change results from conditional-score results. Optimizer numerical stopping is separate from `JointAssessment`, whose finite-difference projected gradient has an explicit absolute default tolerance `1e-4` in retained parameter units. Feasible positive-theta probes establish boundary descent when found; otherwise a boundary remains unresolved. A fit can be returned with failed/unavailable outer evidence. No Hessian, global optimum, finite-MLE existence, or statistical coverage certificate is implied.

Training `conditionalEta`, `conditionalMean` and link-scale `randomEffects` are immutable snapshots, with term, level and coefficient identities. New-data GLMM prediction, GLMM intervals/bootstrap, profile likelihood, KR/Satterthwaite inference, AGQ, negative binomial, noncanonical joint links and estimated-dispersion families are excluded. Dispersion families now refuse explicitly rather than imply an unqualified unit-dispersion fit. Existing allowed noncanonical fast links have not acquired a new qualification claim.

## Independent numerical evidence

The scalar oracle uses raw rows and score bisection, independently of engine design/factors/PIRLS. It checks modes, final-mean determinant and retained-constant objective, including theta zero, offsets, grouped binomial and weighted Poisson. A separate raw dense Newton/Cholesky oracle checks crossed/vector effects and reduced-rank covariance. Tests also cover `A -> B -> A` probe order, malformed counts/trials, budgets and invalid paths, high-weight near-integer responses, extreme logit endpoints, grouped-versus-expanded equivalence, and a zero-theta saddle that passes projected gradient but has feasible variance descent.

Live R 4.5.1 / lme4 2.0.1 generated eight references: four datasets at `nAGQ=1` and `nAGQ=0`. The controls are bobyqa, maxfun 100000, rhoend `1e-9`, and explicit `tolPwrss=1e-10`. Raw frames, categorical level order, generator hashes, optimizer receipts, full objectives, fixed effects, theta, eta, means and conditional modes are frozen. R reported that lme4 was built under R 4.5.2; that loading warning and locale warnings are retained with the executing R 4.5.1 session. The 300-row synthetic Poisson case includes correlated random slopes, crossed effects, zero counts and exposure offsets; its generating R seed is 90210 and the exact generated frame is retained.

| Dataset | lme4 full objective | JVM absolute difference | JS absolute difference | Largest observed projected gradient |
| --- | ---: | ---: | ---: | ---: |
| cbpp | 184.052563742994 | 1.35e-08 | 1.35e-08 | 4.7e-08 |
| culcitalogreg | 60.705890386427 | 4.92e-05 | 4.92e-05 | 3.34e-08 |
| synthetic_poisson_offset | 68.333634742643 | 1.78e-12 | 1.78e-12 | 3.51e-08 |
| synthetic_poisson_crossed | 1293.855308058347 | 4.79e-06 | 4.79e-06 | 1.08e-06 |

All four joint fits on both platforms pass the `1e-4` first-order gate. Frozen comparison gates are objective `1e-4`, beta/theta `1e-3`, eta and modes `2e-3`; means use the corresponding log-link relative or logit Lipschitz bound. These are finite tolerance comparisons, not exact arithmetic identities. [Final fit results](glmm-tranche/fit-results.json) retain all platform-specific values and differences. The scorecard adds four joint lme4 parity rows and four documented fast-divergence rows; the policy against fast-PIRLS/lme4 release-parity claims is preserved.

The small Poisson fast case differs from R nAGQ=0 by approximately `0.0141` in theta, `0.00311` in full-constant objective, `0.0155` in eta and `0.0905` in fitted mean. Its explicit divergence guards are theta `0.02`, objective `0.01`, and eta/modes `0.02`, with means bounded by the corresponding multiplicative log-link change. Those guards document bounded divergence; they do not establish estimator parity. The other fast rows retain tighter sanity comparisons, and remain labelled divergence too. Existing contra/MixedModels.jl tests retain their previous thresholds.

## Adverse evidence and review

Default R `tolPwrss=1e-7` initially reported cbpp objective `184.05313277509163`, about `0.000569` above the independent conditional-mode result. Tightening the explicit reference solve, rather than weakening the comparison, yields `184.05256374299438`. The original references, provenance and log remain under `glmm-tranche/prior-default-pwrss/`. Culcita and the crossed Poisson retain smaller differences attributable to numerical fitting and are compared under their recorded gates.

Initial joint numerical stops left projected gradient norms up to about `0.00471`. Tighter existing TrustBQ controls (final radius `1e-8`, zero objective/stall tolerances, stall count 20) now pass the same absolute assessment threshold. This is a stopping-policy choice, not a new optimizer algorithm. It is not a guarantee on other data. Main evaluation budgets exclude independent post-fit gradient/boundary probes; final evaluations and errors remain visible. Intermediate compile failures, reference discrepancies and failed checks are retained in `glmm-tranche/checks/` with actual exits.

An independent read-only scientific review examined the conditional formulation and implementation, found four concrete issues that were corrected, then returned no remaining blocker at the source hashes in `review.json`. The reviewer inspected test receipts but did not execute an additional test run. Review does not extend the qualification beyond these fixtures and contracts.

## Final checks and preservation

**357 JVM and 353 Scala.js tests pass**, including 27 additional tests per platform; two intentionally opt-in JS benchmark suites remain ignored. The documented CLI workflow also succeeds. Default tests call neither R nor Julia. No Scala compiler warnings were emitted in final closure; an earlier build emitted sbt's informational multiple-main-classes notice after adding the second runnable example. No performance or allocation claim is made for this tranche.

All 54 prior frozen records are preserved, including tolerances and duplicate multiplicity; eight were added. All starting files outside the 19 explicit existing-file changes retain their exact hashes. The formatter touched 13 unrelated files; their exact original bytes were restored from hash-verified Git/source-archive content before the final full suite. Earlier dirty work is preserved. The prior workflow receipt was fully verified at entry: 134 source/input hashes and 86 evidence hashes matched.

```sh
sbt -batch test
sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.GlmmWorkflow'
Rscript comparison/scripts/freeze_glmm_references.R
python3 comparison/scripts/embed_glmm_references.py
python3 comparison/scripts/update_glmm_catalog.py
```

The reference generators use checked-in frames. `vendor_harness.py` preserves the new pins on regeneration. [receipt.json](glmm-tranche/receipt.json) binds the final source archive, input/evidence hashes, actual exits and preservation checks. The tested dependency remains Gale `f869613cec0a89e57b6c995b0a02cf471ac7127c`; Rust reviewed through `33c8639d4e4eaf055dfa09ac36794c5b11912332`. Source Git HEAD remains `ec1937e2d2611c021519c4049c7eda49249f6311`, so the archive and hashes identify the uncommitted implementation.

The larger fast corpus (arabidopsis, gopherdat2, contraception slopes, verbagg and rare-event rows), AGQ and estimated dispersion still need their own evidence. This closes the bounded canonical GLMM workflow, not the entire parity program. No commit, push, hosted CI, merge, deployment or publication was performed.
