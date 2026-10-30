# Likelihood inference and model comparison

This tranche provides selected-objective GLMM coefficient uncertainty, explicit
ML LMM profile targets, verified GLMM likelihood comparisons, and fitted-null
parametric-bootstrap calibration for LMM/GLMM model pairs. Results carry numerical
evidence or precise refusals. Local checks and simulation observations are
separate from scientific calibration, hosted CI and publication.

Starting checkpoint: `74239f432790d9f2b90d1864280b66a3bfdc173d`. Gale stays pinned
at `f869613cec0a89e57b6c995b0a02cf471ac7127c`. Independent R witnesses use R 4.5.1,
lme4 2.0.1 and jsonlite 2.0.0. The complete source/evidence receipt is
[inference-tranche/receipt.json](inference-tranche/receipt.json). The formulation
contract is [CONTRACT.md](inference-tranche/CONTRACT.md).

## Delivered behavior

| Package | Contract |
| --- | --- |
| GLMM coefficient inference | `GlmmInference.coefficients` evaluates the chosen joint-Laplace/AGQ full-density objective and returns retained-basis covariance, standard errors, z statistics and approximate normal intervals for coefficients or specified contrasts. No coefficient p-value is supplied. |
| Nuisance adjustment | With D = -2 log likelihood, covariance is the beta block of **2 H(D)^-1**, including theta nuisance coordinates. Holding theta fixed would use the wrong covariance. Central Hessians at two step sizes must pass shared-diagonal scaled Hessian/covariance and coordinatewise SE stability checks. Positive scaled Cholesky only; no ridge or pseudo-inverse. |
| LMM profiles | `Profile.confint(fit,target,options)` supports retained coefficients, residual SD, and actual SD/variance for an independent scalar random term. ML is required; REML needs explicit `refitWithML`, which retains criterion identity and preserves the actual response even after response-only refitting. Vector/correlated component and GLMM profiles remain unavailable. |
| Profile diagnostics | Center likelihood, constrained nuisance termination/stationarity, brackets, final LR values, failed evaluations and controls are retained. Better constrained likelihood invalidates the baseline. Exhausted searches and unresolved roots have no fabricated endpoint. Nuisance boundary evidence is conservatively refused. |
| GLMM comparisons | `Lrt.statistic` requires identical canonical response, family/link, trials, offsets, unit case weights, approximation/order, plus exact fixed-column and random-partition/covariance inclusion. Equivalent bases without an exact witness can refuse. `Lrt.compare` has ordinary asymptotic chi-square calibration only for interior fixed-effect changes with unchanged random structure. |
| Null bootstrap | `LikelihoodBootstrap.compare` draws new random effects and responses under the fitted null, independently cold-refits both models to every same response, and retains every pair, response, seed and statistic/error. Any failed pair prevents calibration; no filtering, retries or replacement. ML LMMs and explicitly joint GLMMs are supported. |

For each observed and simulated model, the comparison gate checks fresh objective
identity, numerical termination, stationarity and positive step-stable curvature.
Boundary qualification is method-specific. Independent scalar boundary terms use
strict positive, stable one-sided scores in **variance = theta squared**, stationary
free coordinates, positive tangent curvature and feasible single/combined probes.
The theta gradient at zero is automatically zero and is insufficient. Weak scores,
near-boundary central stencils and vector/full-covariance boundaries can refuse.
A covariance-rank deficiency expected under a null does not automatically inherit
the coefficient-bootstrap interior prohibition.

These are local numerical checks. Finite joint-MLE existence, global optimality,
exact finite-sample reference laws and nominal coverage remain unassessed. The
existing `finite_joint_mle_not_assessed` risk is preserved.

## Profile likelihood criteria

All profiles use ML. At theta, the existing PLS kernel supplies q = pwrss,
log|V0|, beta-hat, and A^-1 = vcov / sigma-hat squared. For beta_j=t, the exact
constrained Gaussian quadratic is q_t=q+(t-beta-hat_j)^2/(A^-1)_jj. Theta is
re-optimized using log|V0|+n[1+log(2 pi q_t/n)]. For fixed residual SD s, use
log|V0|+n log(2 pi s squared)+q/s squared. For actual random SD tau, impose
theta=tau/sigma and optimize residual scale and other theta. A relative-theta
curve is not reported as an absolute variance profile.

A support endpoint crosses the nominal chi-square-one cutoff. This is a computed
likelihood support set with an explicitly approximate reference, not a coverage
certificate. The shared cutoff search expands its bracket, including extreme
admitted confidence levels. The existing chi-square tail uses an error-function
approximation: at 95% the cutoff differs from R by about `4.7e-6`. Endpoint tests
separately check parameter brackets, LR residuals and independent R agreement.
Residual-scale lower searches approach zero from positive values; actual RE zero
is evaluated on its permitted constrained face.

## Independent witnesses

[references.json](inference-tranche/references.json) contains both raw and normalized
likelihoods, full R/lme4 Hessians/covariances, order-specific null/alternative LR
values, scalar-boundary likelihoods, lme4 profile limits, and independently
optimized dense Gaussian constrained likelihoods. The complete JSON is embedded
for frozen-first default tests; R is not called by `sbt test`.

The binomial witness has 90 rows, 30 groups and 20 trials per row. Distinct
normalizations are preserved:

| Estimator | Raw R objective | Added response-only constant | Full objective |
| --- | ---: | ---: | ---: |
| nAGQ1 | 453.4246865863183 | 0 | 453.4246865863183 |
| nAGQ5 | 174.1079847148045 | 279.0856291736002 | 453.1936138884047 |

The nAGQ1 raw scale is already full density; adding the high-order constant to it
would be wrong. The corrected metadata explicitly guards this distinction. The
LR statistics are approximately 82.3713536049 (Laplace) and 82.3584320361 (AGQ5).
Scala agrees within `1e-5` (Laplace discrepancy about `1.1e-6`), and nuisance-adjusted
beta covariance agrees within `1e-7` for each selected estimator. AGQ5 covariance
is approximately `[[0.032700988,-0.000232224],[-0.000232224,0.004365722]]`.

The Gaussian witness uses ML sleepstudy with a scalar subject intercept. Its
independent dense-normal objective equals lme4 to `2.3e-13`. Seven constrained
points (coefficient days at 8/10/12, random SD at 20/40, residual SD at 20/30)
agree with Scala within `1.1e-7`, against a test tolerance of `1e-5`. The nominal
95% profile limits are:

| Target | Scala lower | Scala upper |
| --- | ---: | ---: |
| days coefficient | 8.8865541562 | 12.0480177629 |
| Subject SD | 26.0070846126 | 52.9359464510 |
| Subject variance | 676.3683459568 | 2802.2142029461 |
| Residual SD | 27.8138546521 | 34.5910525306 |

Endpoints agree with frozen lme4 within `1e-3`; final LR residuals are checked to
`2e-4`. The variance support set is also checked against the squared SD axis.

The scalar binomial boundary witness repeats successes 5/8/11 over x=-1/0/1 in
12 groups, with 20 trials per row. Both fitted variances are zero; independent
R LR is `45.79397009115837`. Fresh variance scores qualify the observed scalar
fits while coefficient uncertainty refuses their boundary covariance. A separate
zero-theta saddle has a negative variance derivative despite zero theta gradient.

## Executable workflow and complete attempt ledgers

```sh
sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.LikelihoodInference'
```

The example computes selected-AGQ5 coefficient uncertainty, compares an intercept
null with the x alternative, runs 39 paired null simulations with seed 20261029,
and obtains the four supported ML profile target kinds above. All 39 GLMM pairs
qualify in this environment. Observed LR is `82.35843203693366`, exceedances are
zero, and the plus-one estimate is `1/40 = 0.025`. The 95% Wilson interval for the
conditional exceedance probability is approximately `[0,0.08967]`. The estimated
binomial Monte Carlo SE is zero at this observed extreme; the Wilson interval
preserves the remaining uncertainty. This is an approximate fitted-null plug-in
calibration, not an exact-size test or evidence that the true p-value is 0.025.

A separate Gaussian fixed-effect comparison has seven qualified pairs with
explicit tight analytic-gradient fitting controls. Its plus-one estimate is
`1/8 = 0.125`; this is a bounded execution witness, not precision qualification.
Default derivative-free controls can fail the stricter gradient gate; the adverse
run is retained, rather than weakening that gate.

The boundary-null seven-pair run retains three near-boundary curvature refusals
(indices 2,4,5); its calibrated result is unavailable. The crossed scalar-boundary
fit terminates numerically but fails a fresh free-tangent gradient (about 0.00754),
so its comparison refuses. Neither becomes a successful calibration. Tests also
retain deliberately exhausted budgets, middle-pair failures, mismatched responses/
offsets/weights, wrong targets, zero-information contrasts and unsupported shapes.
Both models are attempted even when the other fails. A material negative LR
(below `-1e-7` on GLMM objective scale) is a failed pair; only documented numerical
roundoff can be clamped to zero. Numerical ties count as exceedances under the
recorded `LikelihoodBootstrapOptions.statisticTolerance`, default `1e-7`; setting
zero requests exact comparison. The predicate is T_rep >= T_observed - tolerance.
It can only increase the estimated tail relative to exact comparison. Wilson
intervals/MCSE refer to that tolerance-adjusted exceedance probability.

The initial study (retained in `study/`) exposed an effectively zero observed LR
with roundoff-scale runtime differences: exact comparisons gave p=1 on JVM versus
p=0.95 on JS. A regression fixes the declared numerical-tie convention. The final
study reruns the identical predetermined seeds under `study-final/`; initial
results are neither erased nor pooled with the corrected run.

## Bounded calibration, coverage and allocation observations

The executed opt-in study and its final numerical table are recorded in
`inference-tranche/study-final/results.json`. It uses predetermined outer seeds, every
requested response and all attempted pairs. The same 16 outer responses are
replayed on JVM and JS; these are not 32 independent scientific replications.
Normal and profile intervals have nominal 90% levels, and null rejections are
counted at 0.1. Thirty-nine bootstrap pairs imply p-value resolution 0.025.

Refusals remain in requested denominators. Conditional rates are labelled as
such, with possible-outcome bounds treating unqualified cases as either outcome.
This bounded study cannot establish nominal coverage or calibration accuracy.
The corrected run completed with exit 0 in about 215 seconds, including sbt
startup and both runtimes. Both runtimes have 624 attempted inner pairs, all
qualified, from 16 distinct outer binomial response vectors. Replay differences
are zero in this environment. Every final bootstrap tail estimate now agrees
across runtimes, including the numerical-zero case (p=1 on both).

| Measure | JVM | Scala.js |
| --- | ---: | ---: |
| Qualified null comparisons | 16/16 | 16/16 |
| Rejections at nominal 0.1 | 1/16 | 1/16 |
| Qualified approximate-normal slope intervals | 16/16 | 16/16 |
| Normal intervals containing generating slope zero | 15/16 | 15/16 |
| Qualified Gaussian coefficient profiles | 11/16 | 11/16 |
| Gaussian profiles containing generating slope, among qualified | 9/11 | 10/11 |
| Requested-cohort profile coverage bounds, including refusals | [9/16,14/16] | [10/16,15/16] |
| Median full GLMM inference workflow | 2253.21 ms | 8623.31 ms |
| Median cumulative current-thread allocation per workflow | 2.09 GiB | unavailable |
| Median Gaussian profile attempt, including refusals | 8.93 ms | 30.36 ms |

The Gaussian qualification sets differ: JVM refuses indices 3,4,7,11,14; JS
refuses 3,4,7,13,14. Thus equal refusal counts do not mean identical availability.
This evidence does not establish platform-identical profile qualification. The
reference profiles and all qualified endpoint checks pass on both runtimes;
strict baseline assessment can still differ on newly simulated fits. The
possible-outcome bounds above describe this finite requested cohort, not a
population confidence interval. Sixteen outer draws are too few to qualify
nominal coverage or test size.

Measurements use Java 25.0.1 and Node 26.7.0 on macOS 14.3 arm64. The GLMM timer
includes both outer fits, coefficient information, all 39 inner response draws,
paired fits and fresh comparison assessments. Outer response generation and
printed inner-pair ledgers are outside that timer. JVM allocation is cumulative
on the current thread, about 2.09 GiB per full workflow, not resident or peak
memory. It exposes a material allocation cost rather than an efficiency claim.
JS allocation is unavailable. There is one process per runtime and these
attempt medians include refusals; no universal speedup claim is made. Exact
source hashes and archive accompany each study run.

## Qualification and evidence custody

The full local suite before the opt-in study passes **440 JVM and 436 Scala.js
checks**, with four ignored JS entries visible in the raw log. Architecture,
frozen-reference and existing workflow gates execute. The executable inference
workflow also completes. Raw logs, actual exits, failed runs and R warnings remain
in `inference-tranche/`; only regenerated build outputs were cleaned after an
overlapping-build class-cache failure. Checks are serialized thereafter.

The existing numerical catalog is unchanged. Unrelated tracked/untracked work
is preserved; twelve formatter-only paths were restored to verified starting
bytes. The receipt binds final source, fixtures, helpers, references, logs and
review evidence. No hosted CI, push, release or publication is claimed.

GLMM profiles, correlated/vector variance profiles, KR/Satterthwaite, estimated
dispersion and optimizer redesign remain outside this tranche. Strict local
qualification can refuse usable models, particularly near covariance boundaries;
those limitations and any study failures remain explicit.
