# GLMM analysis workflow and scalar AGQ

This tranche completes a chosen-estimator GLMM workflow: compile, fit, assess,
predict, compute residuals, simulate, and bootstrap specified coefficients or
contrasts. Scalar adaptive Gauss–Hermite quadrature (AGQ) is independently
qualified. The implementation remains a source-only prototype; these are local
checks, not hosted CI, a release or publication.

The reviewed starting checkpoint is
`d2ac2a1115a6c3ad1b3026f9e6eb9e210dd26ba2`. Gale remains pinned at
`f869613cec0a89e57b6c995b0a02cf471ac7127c`. Rust `33c8639d4e4eaf055dfa09ac36794c5b11912332`
provides the vendored rare-event diagnostic frame. Numerical reference versions
are R 4.5.1, lme4 2.0.1 and jsonlite 2.0.0. Raw R session information and its
package-build/locale warnings are retained. The final source and evidence are
bound by [receipt.json](glmm-analysis-tranche/receipt.json), not by an assumption
that the working checkout will remain unchanged.

## Delivered contracts

| Work package | Behavior and boundary |
| --- | --- |
| Scalar AGQ | `GlmmOptions.jointAgq(nodes)` accepts orders 1–31. One node dispatches to joint Laplace exactly. Higher orders require one scalar random-effect term, including a scalar random slope. Crossed/vector terms refuse higher-order AGQ. |
| Reusable fitting | `Glmm.prepare` freezes the coding, retained coefficient basis and family/link. `session(options)` owns private serialized scratch; every response starts independently. Trials, offsets, weights, estimator and numerical controls are fixed. Every batch attempt is retained. Direct/prepared defaults both remain labelled fast-PIRLS; the qualified example explicitly chooses AGQ5. |
| Prediction/residuals | Frozen training coding supplies population/conditional predictions on link/response scales. Population means set random effects to zero. Unknown groups refuse unless zero effect is explicit. New offsets are required when fitting used offsets. Binomial predictions are proportions. Raw and Pearson residuals are available. |
| Seeded simulation | Canonical Bernoulli/logit, binomial/logit and Poisson/log responses; draw `b = Lambda z` or explicitly condition on fitted modes. No Gaussian residual scale is introduced. Random-effect and observation streams are separate. Non-unit case weights refuse. Poisson means above `1e6`, invalid parameters and exhausted sampler budgets have precise errors. |
| Assessment/bootstrap | Fresh conditional and outer checks evaluate the selected quadrature objective. Qualification additionally requires accepted numerical termination, interior/full-rank covariance and a sufficient fixed-effect recession/overlap witness. Every seed, simulated response, fit/error, assessment and target statistic is retained. Any failed or unqualified replicate prevents every interval; attempts are never replaced or dropped. |
| Qualification | Full JVM/JS suites, independent integration, frozen lme4 fits, R refits of the actual bootstrap response vectors, calibrated samplers, analytic/integrated response moments, and repeated matched rebuild/refit measurements. |

Fit snapshots are immutable and remain usable after later session refits. Separate
sessions have independent scratch. Tests interleave successful, invalid and
successful responses and compare against independent fresh fits, including
weights, trials, offsets and objectives. Large-trial binomial `K/n` roundtrips use
an ULP-aware success-count tolerance before canonicalization. Unsupported
approximations never silently fall back.

`GlmmAssessment` separates numerical stopping, conditional score, projected outer
gradient, covariance rank and fixed-effect recession evidence. Its defaults are
conditional score `1e-8`, absolute projected gradient `1e-4` in retained parameter
units, and fresh objective relative agreement `1e-9`. The sufficient witnesses are
deliberately incomplete. Mixed binomial rows, positive Poisson rows or exact
mixed-outcome cells spanning the retained design can exclude fixed-effect
recession; limited one/two-coordinate geometry has separate checks. Unavailable
or numerically uncertain evidence refuses bootstrap qualification. Pure-outcome
groups are risk flags, not proof that Gaussian-penalized conditional modes fail
to exist. `finite_joint_mle_not_assessed` remains visible even for an eligible fit.
There is no finite-joint-MLE, global-optimality, Hessian or inferential-coverage
certificate.

## Independent numerical evidence

[references.json](glmm-analysis-tranche/references.json) retains raw lme4 objectives
at nAGQ 1, 5, 9 and 15 for cbpp, culcitalogreg and a synthetic Poisson-offset
frame. [integrals.json](glmm-analysis-tranche/integrals.json) uses stabilized R
adaptive integration of the full conditional response density times the standard
normal density per group, independently of the fitting engine's PIRLS equations.
Four additional fixed-parameter oracles cover scalar random slopes and zero
variance for binomial and Poisson responses.

High-order lme4 raw objectives use response-only saturated-density constants.
The measured correction is `83.956707537871154` for cbpp and
`55.837423853260` for the Poisson fixture; binary culcita needs zero correction.
The raw value and correction are both preserved. Comparisons use the common full
response-density convention. Joint Laplace/AGQ and fast-PIRLS retain distinct
labels and objective contracts.

All nine fitted AGQ5/9/15 rows pass the frozen objective (`1e-4`), coefficient and
theta (`1e-3`), conditional predictor (`2e-3`) and independent first-order gates
on JVM and JS. At frozen nAGQ15 parameters, 31-node evaluation differs from R
integration by:

| Frame | Absolute full-density objective difference |
| --- | ---: |
| cbpp | `1.14e-13` |
| culcitalogreg | `3.18e-6` |
| Poisson with offset | `1.42e-14` |

Finite quadrature order is still an approximation. Culcita's AGQ15 objective is
about `1.71e-4` from the independent integral at those parameters. Neither a
monotone order-error claim nor exact integration is implied. Gaussian-rule
moments cover every admitted order, one-node identity is exact, and zero theta
has zero correction. Scalar-slope oracle objective tolerance is `1e-7`.

The catalog grows from 62 to 71 records, and the scorecard from 32 to 41 rows.
Every prior frozen record is preserved exactly as a multiset. The nine actual
AGQ comparisons support their `release_blocking_parity / lme4_joint_agq` rows.
The imported 120-row rare Bernoulli frame is used for diagnostic timing only;
this tranche adds no parity claim for that corpus. Fast-PIRLS remains documented
divergence.

## Complete executable example and bootstrap witness

Run:

```sh
sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.GlmmAnalysis'
```

The example uses a frozen 90-row binomial frame (30 groups, three predictor
values, 20 trials per row), R generator seed `2026102701`, with true slope `0.5`.
It chooses AGQ5, produces predictions at zero random effects and residuals,
simulates, and requests percentile intervals for slope `x` and contrast `2*x`.
Bootstrap seed is `20261027`, with 39 independent cold refits. The successful
example is not the failure policy: tests inject failed middle attempts, budget
failures, constant rare responses, foreign design/estimator results, boundary
fits, unsupported weights and unavailable witnesses. Later attempts still run;
intervals refuse if any attempt is bad.

The executable JVM fit reports full-density objective `453.1936138884049`, fresh
conditional score about `1.24e-11` and maximum projected gradient about `5.08e-7`.
Its zero-random-effect predictions at x = -1, 0, 1 are approximately
`0.284921558, 0.417059489, 0.562291147`.

All 39 example attempts qualify on each runtime. The actual 90-value response
vectors and indexed seeds are identical on JVM and JS and are frozen in
[bootstrap-responses.csv](glmm-analysis-tranche/bootstrap-responses.csv).
Independent R/lme4 AGQ5 refits read those vectors unchanged; no R/Scala seed
identity is assumed. Base/replicate coefficients, theta, objectives and base
predictions are compared. The x percentile interval is:

| Producer | Lower | Upper |
| --- | ---: | ---: |
| JVM | 0.5051113020066174 | 0.7101344329487466 |
| Scala.js | 0.5051113022195416 | 0.7101344332508962 |
| R/lme4 | 0.5051113026090784 | 0.7101344272287792 |

The largest endpoint difference is below `6.1e-9`. The generating slope `0.5`
lies outside this particular interval. This bounded 39-replicate demonstration
qualifies implementation consistency, not coverage or Monte Carlo precision.
Intervals use empirical type-7 percentiles. Conditioning on estimated modes is
explicitly labelled and has a different sampling interpretation; it does not
sample posterior random-effect uncertainty.

## Sampling evidence and retained defects

The exact NumPy distribution source and license are pinned at commit
`b833c1b9bf6712cc8a03bb301b77c48984a2373d`, with upstream byte hashes verified.
The license is included in both locally built artifacts. The Scala implementation
uses finite inversion/BTPE for binomial and inversion/PTRS for Poisson; it does
not use a Gaussian approximation or silently clamp/restart tails.

Review found and fixed three substantive sampler defects: missing `log(k)` in
a small-count Stirling error, an incorrect near-mean deviance-series seed, and a
BTPE rejection branch that could enter an enormous loop. The retained regression
seed 116 with `Int.MaxValue` binomial trials exercises the last case. Raw failed
and passing runs remain in `samplers/`; failures are not rewritten as successes.

Each runtime checks 20,000 draws in each of nine distribution regimes (180,000
primitive draws), with mean and known-mean second moments calibrated using true
fourth moments, plus independent R CDF bins around 0.01, 0.5 and 0.99. Nine
8-seed vectors match across runtimes in this environment. Direct Poisson PMF
references cover small counts and near the admitted mean limit. Independently
computed response moments cover Poisson-lognormal crossed/sloped covariance,
logistic-normal binomial and Bernoulli mixtures, and conditional Poisson
independence, using 4,000 simulated response vectors per experiment.

These bounded distribution checks are not a proof of exact RNG behavior at every
parameter or across all future runtimes. The sampler budget and mean limit remain
public refusals. Power-likelihood fitting does not invent a sampling law.

## Runtime checks and evidence custody

`checks/full-final.log` and its sidecar record `sbt -batch test`, the runnable
example and local JVM/JS `packageBin`, actual exit 0:

- JVM: 411 passed, no failures/errors.
- Scala.js: 407 passed, no failures/errors; three ignored entries remain visible in the raw log. Opt-in benchmark
  execution is recorded separately below.
- The log retains one sbt warning that multiple executable main classes exist;
  there are no compiler warnings.
- Architecture and harness/scorecard gates execute locally. No R or Julia is
  needed by the default frozen-first test suite.
- Exact upstream license bytes are verified inside both local jars in
  `package-license-check.json`. No artifact was published.

Formatting was applied to tranche sources. Twelve unrelated formatter-only
changes were restored to their verified starting bytes. All 387 tracked files
outside the 17-file existing-file allowlist and all 22 pre-existing untracked
names are preserved. See `format-preservation.json` and `preservation.json`.

The [independent source review](glmm-analysis-tranche/review.md) is bound to the
85-input manifest SHA-256
`0119c281037610932ee5915acd374d0276bf6057985d1776e14027ff0070c97b`.
The reviewer recomputed all hashes and found no unresolved consequential defect
within the stated scope. Runtime and benchmark claims come from the parent's
executed checks, not that source-only review.

## Performance qualification

The three final processes completed with actual exit 0 in about 487 seconds
including startup, warm-up and assessments. There are **960 measured fit attempts
and 480 matched pairs**, using 44 distinct newly simulated responses: 12 each
binomial, small Poisson-offset and rare Bernoulli, plus eight crossed Poisson.
The same bank is reused across processes and estimator orders; repetitions are
not additional independent response draws. Response values agree across runtime
and process runs within `1e-12`. Each path has three unmeasured response warm-ups.

Host: Apple M3 Max, macOS 14.3 arm64; Java 25.0.1 and Node 26.7.0. Each entry
below is the median of three process medians in milliseconds, over every attempt
including failures. Allocation is the corresponding current-thread JVM median
in MiB; JS allocation is unavailable. The fit timer includes its built-in outer
diagnostic. An additional fresh inference assessment is timed separately in the
machine-readable results. These are neither peak-memory nor cold-start numbers.

| Runtime | Workload | Nodes | Rebuild ms | Refit ms | Rebuild → refit MiB | Fit failures | Unqualified successful fits | Attempts per path |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| js | binomial | 1 | 167.25 | 166.43 | unavailable | 0 | 0 | 36 |
| js | binomial | 5 | 209.65 | 208.86 | unavailable | 0 | 0 | 36 |
| js | crossed_poisson | 1 | 3986.49 | 3959.23 | unavailable | 3 | 0 | 24 |
| js | poisson_offset | 1 | 101.75 | 101.54 | unavailable | 0 | 18 | 36 |
| js | poisson_offset | 5 | 137.37 | 137.21 | unavailable | 0 | 18 | 36 |
| js | rare_bernoulli | 1 | 141.62 | 140.79 | unavailable | 0 | 15 | 36 |
| js | rare_bernoulli | 5 | 188.73 | 187.91 | unavailable | 0 | 15 | 36 |
| jvm | binomial | 1 | 48.90 | 47.65 | 29.20 → 29.07 | 0 | 0 | 36 |
| jvm | binomial | 5 | 68.01 | 67.58 | 38.10 → 37.95 | 0 | 0 | 36 |
| jvm | crossed_poisson | 1 | 1032.33 | 1038.42 | 148.57 → 147.98 | 3 | 0 | 24 |
| jvm | poisson_offset | 1 | 27.83 | 27.39 | 24.30 → 24.25 | 0 | 18 | 36 |
| jvm | poisson_offset | 5 | 40.60 | 39.83 | 32.27 → 32.21 | 0 | 18 | 36 |
| jvm | rare_bernoulli | 1 | 44.55 | 43.90 | 14.79 → 14.63 | 0 | 15 | 36 |
| jvm | rare_bernoulli | 5 | 65.81 | 65.28 | 19.67 → 19.52 | 0 | 15 | 36 |

All matched successes agree on objective within `1e-6` and beta/theta within
`1e-5`; estimator, options and quadrature identity match exactly. Failed pairs
match error codes. The source hash guard passes after every process completes;
its exact source archive and hashes are retained beside the raw logs.

Crossed measurements use an explicit 2,000-evaluation cap on both paths. The
response-generating fit and small workloads retain 10,000 evaluations. Production
defaults were not changed, and failed attempts were not retried. The one distinct
crossed response that exhausts its budget is counted three times per runtime/path
because it appears in all three process runs. Six of twelve small-Poisson and
five of twelve rare-Bernoulli responses have successful numerical fits but fail
strict bootstrap qualification on both orders/runtimes. Boundary covariance and
outer-check refusals remain visible in each row's ledger.

Session reuse modestly reduces most small-workload timings and allocations; fit
work dominates compilation. The crossed JVM refit median is slightly slower in
this run, while the JS difference is small. There is no universal speedup claim.
The earlier pilot used three responses and a 10,000-evaluation crossed cap; its
adverse failures and timings are retained separately and are not pooled with the
final measurements.

## Remaining boundaries

AGQ supports one scalar random-effect term only. Conservative bootstrap witnesses
can decline usable models; overriding tolerances is explicit and cannot certify
MLE existence. Dispersion families, negative binomial, noncanonical joint links,
profile intervals and KR/Satterthwaite inference remain outside this tranche.
The crossed optimizer cost and failures remain measured limitations; this work
does not introduce another optimizer project or claim universal session speedup.
