# Authorized likelihood inference tranche

Outcome: selected-objective GLMM coefficient uncertainty; LMM likelihood support
intervals with nominal chi-square calibration; nested GLMM statistics and
explicitly qualified asymptotic comparisons; fitted-null parametric-bootstrap
calibration for LMM/GLMM pairs, including bounded scalar boundary nulls; independent
references, bounded studies and a documented JVM/JS workflow.

Non-goals: GLMM profiles, KR/Satterthwaite, dispersion families, general optimizer
redesign, hosted CI, push or publication.

Strong formulation review (read-only `/root/glmm_workflow_contract`) established:

- For D=-2logL, coefficient covariance is 2(H_D^-1)beta,beta. All covariance
  nuisance coordinates are included. Central Hessian probes use the selected
  approximation/order and tighter cold conditional solves. Step-halving checks
  must compare shared-diagonal-equilibrated Hessian and covariance plus every SE;
  raw Frobenius norms alone can conceal unstable coordinates. Positive-definite
  scaled Cholesky only; no ridge or pseudoinverse.
- Profiles use an ML baseline. REML requires explicit ML refit; retain both
  criterion identity and baseline. Coefficient constraints update q via the
  exact Gaussian conditional quadratic. Residual scale uses the unprofiled-sigma
  likelihood. Scalar actual RE SD tau requires theta=tau/sigma, with sigma and
  other theta re-optimized. Relative-theta curves are not variance profiles.
- Constrained endpoints retain nuisance solutions, errors and brackets. Better
  constrained likelihood invalidates baseline. Failure is not an infinite or
  fabricated finite endpoint. A computed support set and its nominal chi-square
  reference are separate claims.
- Boundary-null qualification differs from coefficient bootstrap. Theta gradients
  vanish automatically at zero and are insufficient. Independent scalar terms
  require strict positive, stable one-sided derivatives in variance theta²,
  stationary free tangent coordinates, stable positive tangent curvature and
  feasible single/combined direction probes. Weak scores/vector boundaries
  refuse. This is local evidence, not finite-MLE/global-optimum certification.
- Nesting requires identical rowresponse/family/link/trials/offset/weights/
  approximation/order plus exact fixed basis and random partition/covariance
  inclusion. Covariance changes do not inherit ordinary chi-square calibration.
- Null bootstrap draws new random effects and responses under fitted null, then
  cold-refits BOTH models to each same response. Observed and all replicate fits
  pass the method-specific gate. Materially negative LR and every failed pair
  prevent calibrated extraction. No filtering, replacement or retries.
- Tail includes ties; (b+1)/(B+1) is a fitted-null approximate calibration, not an
  exact-size guarantee. Wilson uncertainty refers to conditional tail b/B only;
  not nuisance-estimation uncertainty or scientific interval coverage.

References consulted: https://lme4.github.io/lme4/reference/vcov.merMod.html and
https://hojsgaard.github.io/pbkrtest/reference/pb-modcomp.html. The latter documents
negative numerical LR in bootstrap samples; this implementation retains and
refuses material negatives instead of reducing its denominator.

Acceptance: independent full-objective/Hessian reference, Gaussian constrained
profile references and endpoint LR checks, scale-dominance regression, zero-theta
saddle and strict boundary example, injected failed pair with later attempts,
full default JVM/JS tests and executed opt-in simulation/performance studies.
All old numerical evidence and unrelated work remain preserved.
