# Bounded independent source review

Reviewer: `/root/glmm_workflow_contract` (lean_expert, read-only). Verdict: accepted within the reviewed numerical/API scope; no unresolved consequential defect identified.

Bound to `review-inputs.json`, SHA-256 `0119c281037610932ee5915acd374d0276bf6057985d1776e14027ff0070c97b`. The reviewer independently recomputed all 85 listed input hashes with zero mismatches. This verifies source identity, not exhaustive review of every listed file.

Reviewed scope:

- ULP-aware binomial success-count canonicalization accepts valid large-trial roundtrips while rejecting substantive fractional successes.
- Prepared fit/session/batch defaults match the direct fast-PIRLS default.
- AGQ standard-normal normalization, response constants, final conditional curvature, scalar slopes, zero-theta identity, shape refusals and immutable accepted state follow the agreed contract.
- Previously identified Poisson density and BTPE rejection defects are corrected; random-effect and observation streams remain separate.
- Sufficient fixed-effect witnesses remain distinct from joint-MLE existence. Fresh assessment evaluates the selected estimator. Cold refits preserve the frozen contract. Failed/unqualified replicates cannot enter intervals by filtering or replacement.
- Inspected regressions cover large-trial/default consistency, touching hulls and near-dependent witnesses, foreign-design/estimator rejection, failed-middle-replicate accounting, direct sampler densities and independent distribution moments/tails.

The reviewer did not execute sbt, benchmarks or R integrations. Parent runtime checks are separately recorded. This review does not establish finite joint-MLE existence, global optimality, Hessian adequacy, interval coverage, universal cross-platform reproducibility or performance. Fixed-effect witnesses are deliberately incomplete and unavailable evidence refuses bootstrap qualification. AGQ is restricted to one scalar random-effect term. Simulation retains its family, unit-weight and Poisson-mean limits; conditioning on fitted modes has a distinct recorded interpretation.
