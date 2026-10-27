# Efficient reusable LMM engine tranche

Outcome: implement the six-part tranche authorized after RC-03: upstream correctness alignment, reproducible baseline, sparse kernels, analytic ML/REML gradients, gradient-driven TrustBQ, and compile-once/response-only refitting.

Mote authority: `recovery-performance` followed by `recovery-refit`. Dashboard session `codex-a5afe69b` is a supplementary progress view.

Completed locally with the qualification limits below. The [receipt](engine-tranche/receipt.json) binds source and evidence hashes; the [measurement comparison](engine-tranche/comparison.json) retains every phase and objective check.

Preserve prior uncommitted RC-01/02/03 and parity work; keep Gale pinned. Immutable typed public results, private mutable numerical workspaces, precise refusals, and no fake p-values remain mandatory. No publication requested. GLMM joint-Laplace, inference breadth and prediction/simulation breadth are deferred.

## Acceptance and order

1. Audit Rust `696e5e7..33c8639`, port applicable current-LMM fixes and reproduce with independent witnesses. Resolve discrepancies against actual R model matrices, not only source-to-source matching.
2. Measure unchanged baseline before production edits. Reproducible workload/criterion settings; separate construction, first fit, warm fits, repeated-response rebuild baseline, inspection, and memory. JVM and JS reported separately. Bind inputs and runtime metadata.
3. Sparse off-diagonal blocks, coverage-before-write, dense fill-in fallback, and weighted topology reuse. Dense factor/objective witnesses for nested, partial nesting, crossed and vector terms; zero-to-nonzero weights.
4. Analytic profiled ML/REML gradients, valid at theta diagonal zero without inverting Lambda. Independent dense derivative and finite-difference checks. Bound complete working scratch memory.
5. Explicit gradient TrustBQ strategy with projected-gradient and curvature checks, original-quadratic shifted-Newton ray correction; keep derivative-free path. Qualify nested, crossed, boundary and maximal-singular cases before speed claims.
6. Immutable compiled plan and safe reusable response-only refit session. Fit from compiled design without a second compile; response refresh versus rebuild, cold/warm agreement, options and independent batch starts, immutable prior snapshots.
7. Full JVM and JS tests, before/after measurements, no weakened frozen tolerances; requirement-by-requirement final audit.

The baseline does not have a public refit API; its repeated-response operation rebuilds and fits each frame. Later results must label this difference rather than pretending the baseline reused a workspace. Inspection means existing coefficient/variance summaries and must not be mislabeled Hessian certification. Heap deltas are observations, not peak allocation bounds; per-thread allocations are available only on JVM.

## Delivered implementation

- **Upstream alignment:** R-verified categorical and interaction coding, order-independent interaction subtraction, finite numeric/raw factor validation, observed grouping levels, early-column alias retention when ranks agree, and invalid REML pivot refusal. Independent R scripts and output are in `engine-tranche/upstream/`.
- **Sparse factorization:** CSC off-diagonal blocks retain structural zero slots across weight changes. Product, rank-K and triangular-solve destinations promote before mutation when fill-in requires it. Dense Gram/Cholesky witnesses cover nested, crossed, crossed-within-site, partial nesting and vector random terms across weight, theta and response refreshes.
- **Analytic derivatives:** exact profiled ML/REML gradients use per-level solves for one random term and streamed random-space solves for multiple terms. They do not invert Lambda or construct observation-space covariance. Independent observation-space derivative calculations and finite differences include singular Lambda and zero variance. Checked scratch sizing caps additional numeric gradient buffers; it excludes existing factors and runtime object overhead.
- **Gradient optimization:** an explicit analytic-gradient TrustBQ strategy preserves the existing derivative-free default. It combines gradient-difference curvature, secant updates, feasible steps, shifted-Newton ray correction, and local model refresh. Oracle probes count against the budget and returned gradient evidence belongs to the returned best point. Numerical stopping is not certification of stationarity or global optimality.
- **Reusable fitting:** `Lmm.prepare` creates an immutable compiled plan; `plan.fit`, `plan.session` and `plan.fitResponses` avoid repeated compilation. Sessions refresh response cross-products, own their mutable workspace, serialize calls and return immutable snapshots. Cold/previous starts and independent/warm-sequential batches are explicit. A failed batch retains the last successful prefix start; it is not transactional. The README example is exercised on both platforms.

## Upstream disposition

Applicable current-LMM correctness changes were reviewed through Rust `33c8639d4e4eaf055dfa09ac36794c5b11912332`, from the earlier assessment at `696e5e7`.

| Rust change | Scala disposition |
| --- | --- |
| `8795083`: design/parser, grouping levels, rank selection, invalid REML pivots | Applied with R and numerical regression witnesses. |
| `8795083`: comparison geometry and boundary calibration | Prior RC-03 checks actual fitted-design nesting; covariance/boundary p-values remain conservatively refused. |
| `5969b36`: fit compiled design without recompilation | Implemented through compiled fit/plan/session APIs; lazy host-certificate hooks remain deferred. |
| `8795083`: prior weights, weighted simulation/bootstrap, prediction and KR derivatives; `5db36cd`: KR F scaling | Deferred with those unsupported public APIs. Positive finite prior weights and residual simulation SD `sigma/sqrt(weight)` are mandatory when added. |
| `10f1768`: blocked-factor beta profiling and leverage/prediction-variance row solves | Deferred with profile and prediction/influence breadth. |
| `8795083`, `5db36cd`, `43e4c51`: GLMM dispersion/joint-start work | Deferred to GLMM qualification. The accepted-mean factor refresh also affects fast-PIRLS and is explicitly recorded on `parity-glmm-fast`. |
| `28798b3`: native-worker FP environment; `004a5d9`/`33c8639`: restored host callbacks | Inapplicable to the current Scala execution/API surface. |

The audit found no further current-LMM correctness blocker within these boundaries. It does not claim that every upstream feature or Rust performance optimization has been ported.

## Qualification status

The final full test run passes **286 JVM tests and 282 Scala.js tests**, with one deliberately disabled opt-in benchmark suite on JS. `engine-tranche/checks/full-density-policy.log` and its exit metadata retain the complete output. Both runtimes agree on the difficult-fit results. The README compiled/refit example runs in the public-API test.

The optimizer policy is Rust-derived: explicit dimension-dependent LMM plateau settings and a `1e-9` gradient tolerance, with a tighter Scala gradient-path final radius of `1e-6`. Generic TrustBQ options remain caller-controlled. A positive observed secant curvature now contributes its update even when the previous model is indefinite, matching Rust. Objective-based stops require the projected-gradient guard; step/radius stops do not certify stationarity.

**Maximal-singular limitation:** objective **760.6570831961285**, **593 main-loop evaluations**, `StepBelowTolerance`, projected-gradient norm **0.5367257607808256**. This passes the unchanged broad `<770` regression gate, but fails both the strict first-order band (~`7.61e-7`) and loosened objective-stop guard (~`7.61e-4`). It is about 0.334 worse than the earlier stopped endpoint (760.323269), and an intermediate budget-exhausted experiment observed a lower objective. This is retained adverse optimization evidence, not endpoint parity, a stationary solution, or an accuracy-preserving speedup. The strategy remains explicit and opt-in. The difficult singular fixture is excluded from timing ratios.

The final receipt binds the delivered source. Earlier intermediate logs retain failures and stop reasons but are not substitutes for the final source-bound qualification. These include the sparse vector-solve gap, default-budget failures and earlier passing endpoint. An independent R/lme4 diagnostic on the singular fixture also exhausted its budget; that output is a diagnostic, not a frozen optimum reference. The pw2 nested fit meets its existing frozen lme4 beta/sigma/log-likelihood gates.

The initial repeated timing exposed severe sparse-kernel regressions. The initial measured source is archived under `engine-tranche/after/source.tar.gz`, with all 105 input hashes verified, and its adverse ratios retained in `comparison-initial.json`. The repair replaces dense coverage scans with CSC shared-column checks, uses sparse outer-product kernels for dense destinations, and loads streamed gradient columns directly. A second measured pass showed dense crossings still paying CSC overhead. Storage now chooses dense blocks at at least 50% structural incidence (when representable), with sparse topology reserved for sparser designs. This choice is independent of weights and numeric zeros. Scalar sparse scaling, solves and diagonal downdates reuse arrays directly. The final full test run above includes all repairs; all final timing and objective checks have completed. The first repaired source and measurements remain under `engine-tranche/final/`; the delivered measurements are under `engine-tranche/qualified/`.

## Measured results

Three independent processes per side produced 2,340 baseline samples and 3,888 delivered-source samples. All 72 cross-version/strategy objective checks passed the existing absolute 0.002 comparison band. The full 156 phase comparisons, process-median ranges, allocations and raw data are in [comparison.json](engine-tranche/comparison.json). Host: Apple M3 Max, macOS 14.3 arm64; Homebrew Java 25.0.1, Node 26.7.0, sbt 1.12.14, Scala 3.7.4, pinned Gale `f869613cec0a89e57b6c995b0a02cf471ac7127c`.

Times below are milliseconds, median of three process medians. Ordinary fit uses the same derivative-free strategy on both sides. Compiled gradient fitting changes both API and optimizer. Batch comparisons use five intercept-shifted responses; independent refits keep cold starts.

### JVM

| Workload | Ordinary fit before → after | Compiled gradient fit | Five rebuilds → independent refits | Five warm gradient refits |
| --- | ---: | ---: | ---: | ---: |
| sleepstudy | 3.343 → 3.639 | 1.974 | 12.449 → 6.634 | 2.875 |
| pastes | 2.072 → 3.087 | 8.425 | 9.290 → 4.457 | 3.712 |
| penicillin | 3.240 → 2.699 | 2.544 | 5.086 → 2.540 | 1.674 |
| station_boundary | 1.944 → 2.392 | 1.137 | 8.920 → 4.671 | 1.632 |
| nested1200 | 3.060 → 3.269 | 5.373 | 15.601 → 5.196 | 9.756 |
| crossed1200 | 3.637 → 3.978 | 8.888 | 17.557 → 11.151 | 13.258 |

### Scala.js

| Workload | Ordinary fit before → after | Compiled gradient fit | Five rebuilds → independent refits | Five warm gradient refits |
| --- | ---: | ---: | ---: | ---: |
| sleepstudy | 3.162 → 3.255 | 1.306 | 10.445 → 6.758 | 2.261 |
| pastes | 3.414 → 2.482 | 2.031 | 14.592 → 6.938 | 3.259 |
| penicillin | 2.087 → 2.387 | 1.956 | 8.046 → 5.588 | 3.205 |
| station_boundary | 3.746 → 3.497 | 1.725 | 16.971 → 15.599 | 2.835 |
| nested1200 | 7.114 → 4.253 | 12.915 | 34.340 → 11.961 | 21.722 |
| crossed1200 | 9.717 → 10.951 | 39.138 | 47.687 → 43.092 | 74.272 |

Independent refit batches improve every measured row (1.09–3.00×), with lower JVM current-thread allocation in all six rows (about 6–33% less). Ordinary fits are mixed: Scala.js nested/pastes improve, while several other rows regress, including JVM pastes. The JVM process ranges overlap for several short-fit comparisons; retain them rather than generalizing a small ratio. Gradient fits help sleepstudy/station but the streamed backend remains expensive for larger nested/crossed models, and JVM pastes is also slower. A selected-trace acceleration would need new evidence before broadening the gradient performance claim. The default remains derivative-free.

| JVM five-response workload | Rebuild allocation (MiB) | Independent refit allocation (MiB) |
| --- | ---: | ---: |
| sleepstudy | 5.342 | 4.664 |
| pastes | 8.199 | 7.395 |
| penicillin | 4.105 | 3.411 |
| station_boundary | 10.023 | 9.392 |
| nested1200 | 17.636 | 11.891 |
| crossed1200 | 16.921 | 11.352 |

The measurements support a reusable fitting substrate and workload-specific optimizer choices. They do not establish a universal speedup, peak-memory guarantee, or singular-endpoint accuracy improvement.

## Measurement protocol

The initial baseline contains 780 samples and 12 fit witnesses, with an exact source archive and verified hashes. Its three replay processes each completed successfully; the source-only replay emitted sbt Git-version lookup messages because the archive intentionally has no `.git` directory. The measured code and pinned dependency were hash-verified, and these messages did not prevent compilation or either benchmark suite. Repeated runs use separate sbt/Node processes and the same pinned Gale dependency. Six ML workloads cover sleepstudy, pastes, penicillin, the station boundary example, and 1,200-row nested/crossed synthetic designs.

Phases distinguish compile, workspace construction, first use, warm fit, fixed-theta objective, summary inspection and batches of five responses. New APIs add compiled fitting, objective-plus-gradient, independent refits and explicit warm-sequential refits. The five benchmark responses differ by intercept shifts, so warm-start timings represent this favorable case; they are not a claim about arbitrary bootstrap responses. Functional refit tests also exercise nonconstant response changes.

Timings are process-level medians with all raw samples retained. First fit means the workload's first use in that process, not process startup. JVM allocation covers the current thread only; JS allocation is unavailable. Heap deltas include garbage collection and are not peak memory bounds. Host timings do not imply universal speed guarantees, and comparisons that change API or optimizer are labelled separately.

Reproduce from the repository root:

```sh
python3 comparison/scripts/run_engine_benchmarks.py --baseline comparison/engine-tranche/baseline --output /tmp/mixeff4s-before --repeats 3
python3 comparison/scripts/run_engine_benchmarks.py --output /tmp/mixeff4s-after --repeats 3
python3 comparison/scripts/compare_engine_benchmarks.py /tmp/mixeff4s-before /tmp/mixeff4s-after /tmp/mixeff4s-comparison.json
```

## Boundaries

This tranche qualifies the existing LMM engine and reusable fitting substrate. It does not add joint-Laplace GLMMs, AGQ, prediction/simulation, profile/bootstrap inference, KR/Satterthwaite inference, or new statistical calibration. Fast-PIRLS retains its label. Selected-inverse trace acceleration is not required for the exact streamed gradient and remains a possible measured optimization. No hosted CI, commit, push, merge, publication or release is claimed.
