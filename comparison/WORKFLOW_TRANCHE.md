# Trustworthy LMM analysis workflow

The authorized six-package tranche is implemented locally: fit assessment, frozen-design prediction and residuals, seeded Gaussian simulation, conservative parametric bootstrap, classic lme4 parity, and end-to-end JVM/Scala.js qualification. The repository remains an unpublished prototype. Earlier engine and parity work is preserved.

Mote authority is `workflow-lmm` and its assessment, prediction, simulation, bootstrap and qualification children; existing `parity-lmm-repin-classic` and `parity-pathology-reduced-rank` are included. The source/evidence receipt is [workflow-tranche/receipt.json](workflow-tranche/receipt.json).

## Public workflow and evidence

| Package | Delivered behavior | Evidence |
| --- | --- | --- |
| Fit assessment | Independent numerical-stop, recomputed exact projected-gradient, boundary-adequacy and relative covariance-rank fields. Failed and unavailable evidence remain distinct. | Fresh-response checks; rank-one covariance with positive marginals; weak full-rank covariance; unavailable gradient retains rank; maximal-singular first-order failure. |
| Prediction and residuals | Frozen training types, factor levels, contrasts, pointwise transforms, retained columns, alias relations and random-term order. Population and conditional point predictions; explicit unknown-group refusal or zero effect. | Full fitted-state and copied-design guards, reordered/subset factor levels, missing/wrong inputs, alias violations including tiny units and overflow, old-fit stability after refits, independent lme4 predictions and modes. |
| Seeded simulation | Fresh effects `b = sigma Lambda z`, or fixed fitted modes `bHat = Lambda u`; Gaussian residual noise. Indexed streams preserve response identity after failures. | Analytic Gaussian means/covariances for scalar, vector, crossed, conditional and zero-variance cases; 14,000 generated responses per platform, with six-standard-error moment bounds. |
| Parametric bootstrap | Explicit seeds and retained-basis coefficient/contrast targets, one private reusable session, independent cold starts, all attempts retained, R type-7 percentile intervals. | Same 39-response experiment qualified under explicit tighter controls; injected middle failure, unchanged later seeds/statistics, cap precedence, invalid-control refusal, arbitrary confidence level and finite quantile endpoints. |
| LMM parity | Six classic ML/REML references including predictions, residuals, conditional modes and fixed-effect covariance. Classic scorecard rows now name lme4; pastes REML added. | Live R 4.5.1/lme4 2.0.1, frozen frames/provenance, twelve fit/strategy comparisons on each platform. Existing MixedModels.jl numerical gates retained. |
| End-to-end | Runnable public example, full frozen-first suites, matched rebuild/refit measurements on fresh simulated responses, separate assessment cost and failure accounting. | Full test logs, real command-line example, three independent benchmark processes per platform, independent R refits of all 39 actual response draws. |

The complete [example](../src/test/scala/mixeff4s/examples/LmmWorkflow.scala) uses the repository's sleepstudy frame by default and accepts a replacement frame with the same variables:

```sh
sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.LmmWorkflow'
```

It compiles, fits, assesses, predicts, computes residuals, simulates and bootstraps. It prints qualified counts and failed/unqualified indices; its output retains the bootstrap record even when interval extraction refuses. On the qualified reference experiment, the 95% retained-basis days interval is **[9.066105347530183, 11.74643118383828]**. R independently refitted the actual 39 responses and reproduced the endpoints within `1e-7`. Matching R and Scala integer seeds was not used as a substitute for matching response data.

## Assessment and inference boundaries

The default first-order band is `1e-9 * max(1, abs(objective))`. `StepBelowTolerance` and other numerical stops do not imply a passed first-order check. The inherited maximal-singular diagnostic still records objective 760.6570831961285, 593 main-loop evaluations and projected norm about 0.536726; the new assessment reports `GradientCheck.Failed` and refuses bootstrap eligibility. Its broad `<770` regression gate is retained, with no stationary-solution or endpoint-parity claim.

At boundary Cholesky directions, two feasible objective probes can establish descent. Failure to find descent remains unresolved; initial bootstrap intervals conservatively refuse such cases. Covariance rank comes from each relative `Lambda Lambda'` spectrum with explicit relative/absolute thresholds. Rank-unavailable is not rank zero. Legacy pathology `converged_*` labels describe numerical termination and geometry; their notes now include assessment failures and reduced covariance rank has its own label.

An early controlled experiment using ordinary analytic-gradient stopping qualified **23/39** ML replicates and correctly refused intervals; the REML original failed its first-order check. These adverse receipts are retained. Explicit `FitOptions.optimizerControls` expose the existing algorithm's controls, with finite-value validation and authoritative `FitOptions.maxEvaluations`. The diagonal advisory start is unchanged. Smaller final radius, zero objective stop tolerances and a longer stable-position stall allowance qualified all 39 replicates for the same ML experiment and for the REML workflow. The assessment threshold was unchanged, and no optimizer strategy or automatic retry was added.

Tighter controls are requests, not guarantees. An initial synthetic demonstration also produced an interval refusal; its source and log are retained. The broader fresh-response benchmark includes unqualified fits even under tighter controls, and qualification can differ near the threshold across runtimes. Bootstrap returns all attempts and refuses intervals containing any failed or unqualified replicate; it never silently drops or replaces one. The recorded qualification is numerical and first-order evidence, not global optimality, a Hessian certificate, nominal coverage validation, or proof of small Monte Carlo error. The 39-replicate demonstration is deliberately bounded.

Retained-basis targets describe the selected identifiable basis. They do not identify a dropped coefficient in the original redundant parameterization. New-data prediction checks the fitted dropped-column relations with a recorded relative tolerance scaled to the column's units. Fixed-level or alias incompatibility refuses; an explicit zero-effect policy applies only to new random groups. These APIs provide point predictions. Conditional simulation holds estimated modes fixed and does not sample their posterior uncertainty.

## Frozen references and reproducibility

`workflow-tranche/frames/` contains exact Scala fixture values and `frame-provenance.json` binds their source and CSV hashes. `freeze_workflow_lme4.R` creates `lme4-workflow.json`; the default tests use its portable Scala embedding. The catalog update preserves existing numerical pins. `freeze_workflow_bootstrap.R` fits the actual generated responses in `bootstrap-responses.csv`; response seed strings, CSV hash and originating log are retained in `bootstrap-response-provenance.json`. The package loaded successfully with a warning that lme4 was built under R 4.5.2; runtime R 4.5.1, package version and session output are retained. R's locale-startup warnings are also preserved.

```sh
Rscript comparison/scripts/freeze_workflow_lme4.R
python3 comparison/scripts/embed_workflow_references.py
Rscript comparison/scripts/freeze_workflow_bootstrap.R
python3 comparison/scripts/run_workflow_benchmarks.py --output /tmp/mixeff4s-workflow-bench --repeats 3
```

Default `sbt test` invokes neither R nor Julia. Independent checks use model identities and externally generated reference results; bootstrap refit comparisons use identical responses. No GLMM estimator changes, profile intervals, KR/Satterthwaite inference, bootstrap LRT calibration, publication or hosted-CI claim is included.

## Final checks and measurements

**330 JVM tests and 326 Scala.js tests pass**, with two intentionally opt-in benchmark suites ignored on JS. The runnable CLI example succeeds. Raw logs and actual exit metadata are under `workflow-tranche/checks/`; no new compiler warnings were recorded. All 49 prior frozen reference records remain unchanged, including tolerances and duplicate multiplicity; five records were added.

The final benchmark records **1,440 samples and 720 matched rebuild/refit pairs** across four workloads, two optimizer configurations and both platforms. Each process uses 15 distinct freshly simulated responses per workload; the same responses are repeated in three independent processes. Repeated-process counts are not 45 independent bootstrap draws. Rebuild/refit objectives, coefficients, options and stopping outcomes agree within the stated gates. All numerical fits completed. Strict assessment failures and unresolved boundary/rank checks are retained separately with their reason codes.

Times below are milliseconds, median of three process medians. Assessment cost is measured separately. AnalyticGradient selects explicit tighter controls; DerivativeFree retains its existing default policy. Each ratio compares identical strategy and controls.

### JVM

| Workload / strategy | Rebuild → cold refit (ms) | Ratio | Unqualified / repeated attempts |
| --- | ---: | ---: | ---: |
| pastes / AnalyticGradient | 1.656 → 1.246 | 1.33× | 21/45 |
| pastes / DerivativeFree | 1.658 → 1.037 | 1.60× | 45/45 |
| penicillin / AnalyticGradient | 1.786 → 1.225 | 1.46× | 0/45 |
| penicillin / DerivativeFree | 1.739 → 1.327 | 1.31× | 45/45 |
| sleepstudy_intercept / AnalyticGradient | 1.311 → 0.664 | 1.98× | 3/45 |
| sleepstudy_intercept / DerivativeFree | 1.720 → 0.756 | 2.28× | 45/45 |
| sleepstudy_slopes / AnalyticGradient | 1.439 → 0.948 | 1.52× | 6/45 |
| sleepstudy_slopes / DerivativeFree | 2.355 → 1.638 | 1.44× | 45/45 |

### Scala.js

| Workload / strategy | Rebuild → cold refit (ms) | Ratio | Unqualified / repeated attempts |
| --- | ---: | ---: | ---: |
| pastes / AnalyticGradient | 2.307 → 1.837 | 1.26× | 21/45 |
| pastes / DerivativeFree | 1.440 → 1.000 | 1.44× | 45/45 |
| penicillin / AnalyticGradient | 1.998 → 1.585 | 1.26× | 0/45 |
| penicillin / DerivativeFree | 1.934 → 1.477 | 1.31× | 45/45 |
| sleepstudy_intercept / AnalyticGradient | 0.974 → 0.461 | 2.11× | 0/45 |
| sleepstudy_intercept / DerivativeFree | 1.223 → 0.652 | 1.88× | 45/45 |
| sleepstudy_slopes / AnalyticGradient | 1.535 → 1.023 | 1.50× | 3/45 |
| sleepstudy_slopes / DerivativeFree | 2.226 → 1.684 | 1.32× | 45/45 |

Cold-refit fit-time ratios span **1.26–2.28×** on this host. [Full results](workflow-tranche/benchmarks-final/results.json) retain process medians, allocation, assessment costs and all outcomes. These timings include successful numerical fits whose strict assessment fails; they are not all qualified inference runs. Current-thread JVM allocation is measured; JS allocation is unavailable. No peak-memory bound or universal speedup is inferred.

The 39 exported indexed responses agree across runtimes and repeated processes to a maximum absolute difference of **1.14e-13**. Near-threshold fit qualification differs on some other benchmark responses, so runtime-specific outcomes remain visible. Independent R refits use the original exported values. The final-run JVM payload matches that original payload exactly.

Final and earlier benchmark source archives are verified against all recorded input hashes. The receipt binds source, scripts, frozen inputs and checks. Earlier refusals and intermediate failures are retained. Host: Apple M3 Max, macOS 14.3 arm64; Homebrew Java 25.0.1, Node 26.7.0, sbt 1.12.14, Scala 3.7.4 and pinned Gale f869613cec0a89e57b6c995b0a02cf471ac7127c.

Changes remain uncommitted. No hosted CI, push, merge, deployment or publication was performed.
