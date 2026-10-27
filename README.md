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
fast-PIRLS and canonical-link joint-Laplace GLMMs are implemented. AGQ is refused.

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

## Joint-Laplace GLMM workflow

For Bernoulli/logit, binomial/logit or Poisson/log, select the estimator explicitly:

```scala
import mixeff4s.prelude.*

val joint: FitResult[GlmmFit] = Glmm.fit(
  "y ~ 1 + x + (1 | g)", frame, Family.Bernoulli, None,
  GlmmOptions.jointLaplace
)
```

Binomial responses are proportions in `[0, 1]`; supply positive integer counts
in `GlmmOptions.trials`. Counts greater than one in the response, or fractional
successes without matching trials, refuse. `weights` are separate positive case
weights multiplying the entire conditional log likelihood, score and curvature.
Noninteger case weights define a power likelihood. Poisson responses are integer
counts, including zero. Offsets are finite link-scale vectors in
`GlmmOptions.offset`; supply `log(exposure)` for a Poisson exposure offset.

`fit.algorithmLabel` is `joint-Laplace`. Its `deviance` is minus twice the Laplace
log likelihood **including response constants**. Fast-PIRLS keeps its historical
deviance-residual convention; consult `fit.objectiveDefinition` before comparing
objectives. Fast-PIRLS remains a documented divergence from joint fitting.

Joint fitting holds beta fixed in each conditional solve, checks the actual
spherical random-effect score (default maximum absolute value `1e-8`), and
refreshes the determinant at the accepted modes. `fit.inner`, `fit.optimizer`
and `fit.jointAssessment` report conditional evidence, numerical termination,
and outer finite-difference evidence separately. The outer threshold is an
explicit absolute `stationarityTolerance` (default `1e-4`) in retained parameter
units; it does not depend on response-only likelihood constants. A returned fit
may fail that assessment. At a covariance boundary, feasible descent probes can
establish failure; no detected descent remains unresolved. There is no Hessian,
global-optimum or finite-MLE existence certificate.

Fits are immutable snapshots with package-private construction. `conditionalEta`,
`conditionalMean` and `randomEffects` expose training predictions and modes; the
latter includes fitted term, level and coefficient order. GLMM new-data prediction,
intervals and advanced inference remain outside this workflow. Joint noncanonical
links, estimated dispersion, negative binomial and AGQ return explicit refusals.

Run the complete grouped-binomial example:

```sh
sbt 'mixeff4sJVM/Test/runMain mixeff4s.examples.GlmmWorkflow'
```

Four joint cases (cbpp, culcita, offset/zero Poisson, and correlated/crossed
Poisson) are compared to frozen lme4 fits on JVM and Scala.js, including means
and modes. The small Poisson fast path has a measured theta difference of about
`0.014` from lme4 nAGQ=0; that row stays documented divergence. Details and
reproducibility receipts are in [the GLMM qualification report](comparison/GLMM_TRANCHE.md).

## What it covers

- Parse the supported lme4 subset: intercepts, `:` / `*` / `/`, `(re | g)`,
  `(re || g)`, nested and crossed grouping, `us` / `diag` / `cs` / `ar1`
- Lower a typed DSL into the same formula IR
- Evaluate the stateless `I(...)` / `log` / `sqrt` transform subset onto a frame
- Name GLMM families and links, including the rule that Normal+Identity is an LMM
- Compile a formula and frame into `FeTerm`, `ReMat`s, and the live θ `parmap`
- Fit profiled (RE)ML with blocked Cholesky + TrustBQ, including crossed intercepts
- Report VarCorr, Wald SEs, and z-statistics; p-values stay a typed refusal
- Fit separately labelled fast-PIRLS and canonical-link joint-Laplace GLMMs (`Glmm.fit`); Normal+Identity is refused as a GLMM
- Predict with frozen training coding, obtain residuals, simulate Gaussian responses, and bootstrap retained-basis coefficients or contrasts with explicit qualification
- Compare nested LMMs with `Lrt.compare`; profile intervals and bootstrap LRT calibration stay typed refusals
- Emit a versioned compiled-design JSON artifact with a design-time pathology certificate; raw θ assessment reports geometry, while assessment of an `LmmFit` also checks its optimizer evidence

## Fit and boundaries

Good fit for describing mixed models, fitting the current LMM kernel
(sleepstudy, penicillin, pastes), labelled fast-PIRLS GLMMs, and canonical-link
joint Laplace for Bernoulli/binomial and Poisson. Structured RE covariance,
AGQ and estimated-dispersion GLMMs are not implemented.

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
| What does it do? | Mixed-model language, frame, design compilation, profiled LMM fitting, prediction, simulation, qualified coefficient bootstrap, labelled fast-PIRLS and joint-Laplace GLMMs, nested LRT and compiled-design artifacts |
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
numerical termination evidence are refused. Profile intervals and bootstrap LRT calibration
remain unavailable; coefficient/contrast bootstrap is described above.

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
