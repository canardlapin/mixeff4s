# mixeff4s

[Repository](https://github.com/canardlapin/mixeff4s)

mixeff4s is a mixed-effects modeling library for Scala 3 that turns an lme4-style
formula and a column store into either a fit with recorded numerical termination or a typed refusal. Use
it when you want Bates/MixedModels PLS semantics without inventing a p-value.

> This is pre-release development software. APIs, semantics, and package
> boundaries may change. No stable support or binary-compatibility promise
> should be inferred unless the package explicitly states one.

**Status:** `prototype`. The formula language, model frame, family/link types,
Gale-backed design compilation, a blocked-Cholesky LMM kernel, and labelled
fast-PIRLS, canonical-link joint Laplace, and scalar adaptive-quadrature GLMMs
are implemented, with prediction, seeded simulation and qualified bootstrap.

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

## Compile once and refit responses

`Lmm.prepare` returns an immutable LMM plan. It fixes the row order, predictors,
grouping partitions, coding and retained columns. Fitting that plan does not
parse or compile the formula again:

```scala
val options = FitOptions.reml.copy(optimizer = LmmOptimizer.AnalyticGradient)
val responses = Vector(frame.numeric("y").get, frame.numeric("y").get.map(_ + 0.25))
val fittedResponses = for
  plan <- Lmm.prepare(formula, frame)
  fits <- plan.fitResponses(responses, options)
yield fits
```

For repeated calls, `plan.session(options)` returns a session with a private
reusable workspace. `session.refit(y)` refreshes only response cross-products
and starts the optimizer cold. Returned fits, including their response and LRT
evidence, remain immutable. Separate sessions own separate workspaces; calls
on one session are serialized. Responses must be finite and have exactly the
compiled row count; changing predictors or groups requires another plan.

Warm starts are explicit: `session.refit(y, RefitStart.Previous)` requires a
previous successful fit. Batches default to `BatchStart.Independent`, which
resets the optimizer start for each response while reusing design work. Choose
`BatchStart.WarmSequential` explicitly for order-dependent starts. A failed
call cannot replace the last successful start. A batch that fails after some
successful fits retains that successful prefix's last start; batches are not
transactional. These are fitting operations;
they do not supply bootstrap calibration.

The default strategy remains `LmmOptimizer.DerivativeFree`.
`LmmOptimizer.AnalyticGradient` uses exact profiled ML/REML derivatives and a
gradient-driven TrustBQ search. It records gradient evaluation counts and the
projected-gradient norm at the returned point. Its numerical stopping rules
and feasible curvature checks do not certify stationarity or a global optimum. No optimizer
strategy silently switches to the other.

Random-effect cross-products preserve sparse structural slots through changes
in weights and covariance parameters, promoting to dense storage when fill-in
requires it. Analytic gradients use per-level solves for a single random term
and streamed solves for multiple terms, without an observation-by-observation
covariance matrix. `ProfiledGradientOptions.maxWorkingBytes` bounds additional
numeric gradient buffers (64 MiB by default), excluding existing model/factor
storage and runtime object overhead. A smaller cap is a refusal, not a
finite-difference fallback. Measurements and qualification are documented in
[the engine tranche](comparison/ENGINE_TRANCHE.md).

## LMM analysis workflow

Run the complete [public-API example](src/test/scala/mixeff4s/examples/LmmWorkflow.scala)
against the repository's sleepstudy frame:

```sh
sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.LmmWorkflow'
```

It compiles once, fits, assesses, predicts population means at new rows, obtains
conditional residuals, simulates a response, and performs 39 cold-start bootstrap
refits. It prints the qualified count and failed/unqualified indices before
reporting a percentile interval or its refusal. The example's `run(frame)` accepts
a replacement `ModelFrame` with the same named variables. Its source is exercised
on JVM and Scala.js; the command runs on JVM.

`LmmAssessment.assess(fit)` records numerical termination, a recomputed exact
projected-gradient check, and each random-effect covariance's eigenvalues and
rank independently. Its default gradient band is `1e-9 * max(1, abs(objective))`.
Failed and unavailable checks remain distinct. A boundary probe can detect
feasible descent; an unresolved boundary check does not qualify a bootstrap
interval. The existing pathology `converged_*` labels still describe numerical
termination and geometry. Reduced-rank covariance now has its own label.

`Prediction.fitted(fit)` and `Prediction.residuals(fit)` use the captured conditional
modes. `Prediction.predict(fit, newData, PredictionMode.Population)` needs only
fixed-effect inputs. Conditional prediction also needs grouping variables and
random-slope inputs; unknown groups refuse by default. Explicit
`NewGroupPolicy.ZeroEffect` assigns zero to an unknown group's effect. Fixed
factor levels, contrasts, transforms and retained columns come from training;
new rows violating fitted alias relations refuse. These are point predictions,
not prediction intervals.

`Simulation.response(fit, seed, SimulationMode.DrawRandomEffects)` draws fresh
group effects with covariance `sigma² Lambda Lambda'` and independent residual
noise. `ConditionOnModes` holds the estimated modes fixed and varies residual
noise. It does not draw posterior uncertainty in those modes. Indexed seeds
make replicates independent of execution order and earlier failures. JVM/JS
transcendental rounding can differ slightly; measured draw agreement and outcome
qualification are recorded in the [workflow report](comparison/WORKFLOW_TRANCHE.md).

`Bootstrap.parametric(fit, BootstrapOptions(...))` returns every requested attempt,
including seed, response, fit/error, assessment and statistic. Targets are named
retained-basis coefficients or contrasts. A dropped original coefficient refuses.
Percentile endpoints use R's type-7 interpolation and an explicit confidence
level; extraction refuses any failed or unqualified replicate. The two-argument
legacy call requests an explicit seed and targets. The demonstration's 39 draws
keep it short; they do not establish nominal coverage or small Monte Carlo error.

The example selects tighter `FitOptions.optimizerControls` because ordinary
numerical stops can be too coarse for the strict post-fit check. Explicit controls
replace the main-loop family policy; `FitOptions.maxEvaluations` remains the
budget, and the advisory diagonal start retains its existing policy. Controls
propagate through sessions and bootstrap. They are requests, not guarantees of
qualification. No refit is retried automatically or switched to another strategy.

Prediction, simulation, assessment and model comparison validate the complete
fitted snapshot. Copying coefficients, variance, options or optimizer evidence
invalidates that seal and yields a precise refusal. Workspaces remain private;
old fit values stay valid after later refits.

## GLMM analysis workflow

Select an estimator explicitly. The direct and reusable front doors both default
to labelled fast-PIRLS. `GlmmOptions.jointLaplace` uses one-node Laplace;
`GlmmOptions.jointAgq(5)` selects labelled adaptive Gauss-Hermite quadrature.
Orders 1–31 are admitted. More than one node requires exactly one scalar random
term; scalar random slopes are supported. Use explicit Laplace options for crossed/vector terms.
Unsupported shapes refuse rather than switch approximations.

This function accepts a frame containing binomial proportions `y`, predictor
`x`, group `g`, and positive integer trials aligned to its rows:

```scala
import mixeff4s.prelude.*

def analyse(frame: ModelFrame, trials: Vector[Int]) =
  val options = GlmmOptions.jointAgq(5).copy(trials = Some(trials))
  for
    plan <- Glmm.prepare("y ~ 1 + x + (1 | g)", frame, Family.Binomial)
    fit <- plan.fit(options)
    newData <- ModelFrame.of("x" -> numeric(Vector(-1.0, 0.0, 1.0)))
    predicted <- GlmmPrediction.predict(fit, newData, PredictionMode.Population)
    simulated <- GlmmSimulation.response(fit, seed = 42L)
    bootstrap <- Bootstrap.parametric(fit, GlmmBootstrapOptions(
      nsim = 39, seed = 20261027L,
      targets = Vector(BootstrapTarget.Coefficient("x"))
    ))
  yield (fit, GlmmAssessment.assess(fit), predicted, simulated, bootstrap)
```

A plan's `session(options)` reuses private numerical scratch for response-only
`refit(y)` calls. Every refit resets beta/theta to independent cold starts;
trials, offsets, weights, coding, estimator and controls remain fixed. Previous
fit snapshots remain valid. `refitAll` retains a result for every response,
including failures. Use separate sessions for independent concurrent work.

Binomial y is a proportion in `[0, 1]`. `trials` supplies positive integer counts;
without trials only binary responses are admissible. Valid floating `K/n`
roundtrips are canonicalized to integer successes. `weights` are separate positive
case weights multiplying the entire log density, score and curvature; noninteger
weights define a power likelihood. Poisson responses are integer counts including
zero. Offsets are explicit finite link-scale vectors, such as `log(exposure)`.

`GlmmPrediction` uses frozen training coding on link or response scales, with
explicit new-group refusal or `NewGroupPolicy.ZeroEffect`. New data must supply an
offset if the fit used one. Binomial response predictions are proportions.
`PredictionMode.Population` sets random effects to zero; it does **not** integrate
over their distribution. Raw and Pearson residuals are available.

`GlmmSimulation` draws fresh `b = Lambda z`, or conditions on the estimated modes
when `SimulationMode.ConditionOnModes` is explicit. Bernoulli, binomial and Poisson
count samplers use independent indexed streams. Simulated binomial responses are
`K/trials`. Simulation and bootstrap currently require unit case weights;
power likelihood alone does not define a sampling distribution. Poisson simulation
admits means through `1e6`; larger or nonfinite draws return explicit diagnostics.

`fit.inner` records conditional-score evidence; `fit.optimizer` records numerical
termination. `GlmmAssessment.assess` recomputes evidence for the actual quadrature
objective (default score tolerance `1e-8`, absolute projected-gradient tolerance
`1e-4` in retained parameter units). Bootstrap requires passed conditional and outer checks, interior/full-rank
covariance and a sufficient fixed-effect overlap/recession witness. Unavailable
witnesses refuse. Pure-outcome groups remain risk flags, not proof that penalized
conditional modes do not exist. Finite joint-MLE existence and global optimality
remain unassessed.

`GlmmBootstrapResult` retains every response, seed, fit/error, assessment and target
statistic. Any failed or unqualified attempt prevents interval extraction; none are
replaced or silently dropped. Numerical controls may be overridden explicitly,
but estimator, quadrature, trials, offsets and weights must match the fitted model.
Intervals are empirical percentile intervals for specified retained-basis coefficients
or contrasts. The 39-replicate example is a bounded demonstration, not a coverage
or Monte Carlo precision qualification. Conditional-mode experiments are labelled
in the result and do not sample posterior random-effect uncertainty.

Joint Laplace and AGQ report full-density minus twice log likelihood, including
response constants. High-order lme4 AGQ references are normalized explicitly before
comparison; raw values remain in the evidence. Fast-PIRLS keeps its historical
objective and divergence label. GLMM profile intervals, estimated dispersion,
negative binomial, noncanonical joint links and KR/Satterthwaite remain unavailable.

Run the complete compile/fit/assess/predict/simulate/bootstrap example:

```sh
sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.GlmmAnalysis'
```

Independent integration, R refits of the actual simulated responses, cross-runtime
sampler checks and matched fresh-response measurements are documented in
[the analysis qualification report](comparison/GLMM_ANALYSIS_TRANCHE.md).
The [earlier joint-Laplace report](comparison/GLMM_TRANCHE.md) retains its original evidence.

## What it covers

- Parse the supported lme4 subset: intercepts, `:` / `*` / `/`, `(re | g)`,
  `(re || g)`, nested and crossed grouping, `us` / `diag` / `cs` / `ar1`
- Lower a typed DSL into the same formula IR
- Evaluate the stateless `I(...)` / `log` / `sqrt` transform subset onto a frame
- Name GLMM families and links, including the rule that Normal+Identity is an LMM
- Compile a formula and frame into `FeTerm`, `ReMat`s, and the live θ `parmap`
- Fit profiled (RE)ML with blocked Cholesky + TrustBQ, including crossed intercepts
- Report VarCorr, Wald SEs, and z-statistics; p-values stay a typed refusal
- Fit separately labelled fast-PIRLS, canonical joint-Laplace and scalar-AGQ GLMMs (`Glmm.fit`); Normal+Identity is refused as a GLMM
- Predict with frozen coding, compute residuals, simulate Gaussian or supported GLMM responses, and obtain qualified coefficient/contrast bootstrap intervals
- Compare verified nested LMMs or joint GLMMs; obtain explicit LMM profiles and fitted-null bootstrap LRT calibration with full failure ledgers
- Emit a versioned compiled-design JSON artifact with a design-time pathology certificate; raw θ assessment reports geometry, while assessment of an `LmmFit` also checks its optimizer evidence

## Fit and boundaries

Good fit for describing mixed models, fitting the current LMM kernel
(sleepstudy, penicillin, pastes), labelled fast-PIRLS GLMMs, and canonical-link
joint Laplace for Bernoulli/binomial and Poisson. Structured RE covariance,
multidimensional AGQ and estimated-dispersion GLMMs are not implemented.

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
| What does it do? | Mixed-model language, frame, design compilation, profiled LMM fitting, prediction, simulation, qualified coefficient bootstrap, labelled fast-PIRLS, joint-Laplace and scalar-AGQ GLMMs with analysis workflows, nested LRT and compiled-design artifacts |
| Smallest useful example? | Parse `y ~ 1 + x + (1 \| g)` and inspect the IR |
| Maturity / publication? | `prototype`, source-only, unpublished at [canardlapin/mixeff4s](https://github.com/canardlapin/mixeff4s) |
| How to verify? | `sbt -Dmixeff4s.gale.build=/path/to/gale test` (JVM and Scala.js). Default tests read frozen numbers; they do not call R or Julia. |

Platform: JVM and Scala.js (Node). Scala 3.7.4, sbt 1.12.14. The numeric kernel
is Gale, which already cross-builds; mixeff4s follows it. The source-tree
architecture walk is JVM-only.

Fits retain the TrustBQ result as `fit.optimizer`; GLMM fits also retain the
final inner-loop criterion as `fit.inner`: `GlmmInnerResult.FastPirls` carries
an objective-change result, while `ConditionalMode` carries a score result. Exhausted optimizer budgets and
unfinished PIRLS return typed errors (`optimizer_budget_exhausted`,
`pirls_not_converged`, or `pirls_step_failed`). `FitOptions.maxEvaluations`
limits the main LMM search, excluding the advisory diagonal warm start;
`GlmmOptions` exposes `maxEvaluations`, `pirlsMaxIterations` (fast) and
`jointMaxIterations` (conditional). The main GLMM budget excludes independent
post-fit gradient/boundary probes. Explicit `optimizerControls` are supported;
`maxEvaluations` remains authoritative.

These are numerical stopping criteria, not gradient/Hessian or global-optimum
certificates. `Compiler.assess(artifact, fit)` uses the LMM's recorded evidence.
`Compiler.assess(artifact, theta)` leaves `fitStatus` as `NotAssessed`, even
for an interior θ. Compiled-design schema version 2 reports `theta_geometry`
separately from `fit_status`. The existing `converged_*` codes on fit-backed
assessments mean the optimizer met its numerical stopping criterion.
Fast-PIRLS retains its objective-change convergence rule after finite step
halving; the fast result's `halvingLimitHits` records exhausted halving bounds. It does
not use the conditional-score criterion of the joint-Laplace estimator.

Fixed-effect outputs use retained design-column order. For an LMM, `beta`,
`feNames`, `stderror`, coefficient-table rows, and both axes of `vcov` describe
the same retained columns; covariance is not padded with aliased rows.
Both LMM and GLMM fits expose `fixedEffects: FixedEffectLayout`, containing
`requestedNames`, `retainedIndices`, `retainedNames`, `droppedIndices`, and
`droppedNames`. Indices are zero-based positions in the expanded design, not
formula-term indices. Dropped coefficients are not separately identified by
the selected reduced model. Retained estimates and covariance describe that
selected basis; they do not establish separate estimability of coefficients
in the original redundant parameterization.

Reading `fit.feNames` remains supported through `fixedEffects`. LMM case-class
construction/copy checks layout and covariance dimensions; subsequent workflow
operations additionally validate its immutable snapshot seal. GLMM fits have
package-private construction and no public `copy`, so fitted evidence cannot be
retained while replacing parameters. Fitting methods return `FitResult` for
recoverable model errors.

Likelihood comparisons retain immutable, row-aligned design evidence in
`LmmFit.designEvidence`. `Lrt.statistic(smaller, larger)` returns the likelihood
ratio and nominal parameter-count difference after checking response, criterion,
actual retained fixed columns, grouping partitions, random basis columns, and
covariance constraints. Intercepts are ordinary basis columns in this check.
The check proves nesting by exact column inclusion; it conservatively refuses
other equivalent parameterizations. REML requires the same retained fixed basis
up to permutation, since rescaling can change its likelihood normalization.

`Lrt.compare` supplies an **asymptotic chi-square approximation** only for nested
fixed-effect ML comparisons with unchanged random covariance design. It returns
`lrt_reference_unavailable` for covariance changes and for either fit having a
Cholesky diagonal at or below `1e-8`. This boundary screen is conservative and
is not a condition-number or asymptotic-adequacy certificate. Covariance-changing
comparisons can still use `Lrt.statistic`; no p-value is implied. For example,
random intercept → random intercept plus slope has a verified statistic, but
its reference law is not implemented. Other unsupported or invalid comparisons
return `lrt_unavailable`. Manually constructed fits without matching design and
numerical termination evidence are refused. Method-specific fitted-null bootstrap calibration
and explicit ML profile targets are described below.

Numerical evidence lives in `comparison/`. Classic sleepstudy, penicillin and
pastes ML/REML rows now compare against frozen lme4 fits, predictions, residuals,
conditional modes and summaries; existing MixedModels.jl checks remain. Other frozen
lme4 numbers are claims only when the scorecard names them. Fast-PIRLS
rows stay labelled and are never sold as `lme4::glmer`. This evidence
program is not a new port phase.

This is a Scala 3 re-expression of
[mixeff-rs](https://github.com/bbuchsbaum/mixeff-rs), itself an independent
implementation of MixedModels.jl's PLS/PIRLS formulation. It is not a
line-for-line Rust translation.

## Documentation

- [AGENTS.md](AGENTS.md) — package layout, layering, and port phases
- mixeff-rs guide — [what is supported](https://github.com/bbuchsbaum/mixeff-rs/blob/main/docs/guide/05_what_is_supported.md)

## License

Apache-2.0. The Rust reference crate is MIT.

## Likelihood inference and model comparison

`GlmmInference.coefficients(fit)` returns observed-information covariance, standard
errors, z statistics and approximate normal intervals for the selected joint
Laplace/AGQ objective. It includes covariance-parameter nuisance uncertainty;
it does not hold theta fixed. `GlmmInferenceOptions.targets` accepts retained
coefficients or contrasts. Step-halving, scale-aware curvature checks must pass;
boundary, indefinite or unstable curvature refuses without regularization. No
coefficient p-values are supplied.

`Profile.confint(fit, LmmProfileTarget.Coefficient("days"), options)` computes
an LMM likelihood support interval with nominal chi-square calibration. Targets
also include `ResidualStdDev`, `RandomEffectStdDev(group)` and
`RandomEffectVariance(group)` for scalar independent random terms. These are
actual SD/variance profiles, with nuisance parameters re-optimized. Use an ML
fit, or explicitly set `LmmProfileOptions(refitWithML = true)`; the result retains
the criterion and actual ML baseline, including response-only refits. Inspect
both endpoints and the complete evaluation ledger. Failed fits, unqualified
nuisance boundaries and unresolved roots do not become interval endpoints.
Tighter fit/profile controls may be needed; they are explicit in the runnable
example. Vector/correlated RE component profiles and GLMM profiles remain refused.

`Lrt.statistic(nullFit, alternativeFit)` verifies nesting and compatible data and
likelihoods. GLMM comparisons require identical family/link, trials, offsets,
weights, approximation and quadrature order. `Lrt.compare` supplies an ordinary
asymptotic chi-square reference only for qualified interior fixed-effect changes
with unchanged random structure. Covariance changes or boundaries require an
explicit calibration.

`LikelihoodBootstrap.compare(nullFit, alternativeFit,
LikelihoodBootstrapOptions(nsim = 39, seed = 20261029L))` simulates under the fitted
null, draws new random effects, and cold-refits both models to every identical
response. It supports ML LMMs and qualified joint GLMMs. All attempts remain in
the result; any failed or unqualified pair refuses calibration. The estimate is
`(exceedances + 1) / (nsim + 1)`, including numerical ties under the recorded
`statisticTolerance` (default `1e-7`; zero requests exact comparison). Wilson bounds
quantify Monte Carlo uncertainty in that fitted-null, tolerance-adjusted tail probability, not nuisance-estimation
uncertainty or exact test size. A zero exceedance count does not produce a zero
p-value.

Boundary LRT qualification is distinct from coefficient bootstrap. Independent
scalar boundary terms require stable one-sided variance scores, positive tangent
curvature and feasible-direction checks. Weak/near-boundary evidence and vector
boundaries refuse. Finite-MLE existence and global optimality remain unassessed.

Run the complete workflow, including all four supported profile target kinds:

```sh
sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.LikelihoodInference'
```

The [inference qualification report](comparison/INFERENCE_TRANCHE.md) records
independent R curvature, exact Gaussian constrained likelihoods, profile endpoints,
bounded calibration/coverage observations, refusal rates and allocation limits.
