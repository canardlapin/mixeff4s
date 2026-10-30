#!/usr/bin/env Rscript
# Independent R/lme4 likelihood-inference witnesses.  This script does not call
# mixeff4s.  It retains lme4's raw likelihood scale and an explicitly measured
# response-only normalization shift for the binomial witness.
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(jsonlite))

root <- "comparison/inference-tranche"
binomial_path <- "comparison/glmm-analysis-tranche/frames/binomial_workflow.csv"
sleepstudy_path <- "comparison/workflow-tranche/frames/sleepstudy.csv"
dir.create(root, recursive = TRUE, showWarnings = FALSE)

array <- function(x) as.list(unname(as.vector(x)))
sha256 <- function(path) {
  line <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub(" .*", "", line[[1]])
}
controls <- list(optimizer = "bobyqa", maxfun = 100000L, rhoend = 1e-9, tolPwrss = 1e-10)
glmer_controls <- glmerControl(optimizer = controls$optimizer,
  optCtrl = list(maxfun = controls$maxfun, rhoend = controls$rhoend), tolPwrss = controls$tolPwrss)
lmer_controls <- lmerControl(optimizer = controls$optimizer,
  optCtrl = list(maxfun = controls$maxfun, rhoend = controls$rhoend))

# This is deliberately a plain central-difference calculation of the lme4
# deviance closure.  Its parameter order is exactly lme4's c(theta, beta).
central_hessian <- function(f, x) {
  n <- length(x); h <- pmax(1e-4, abs(x) * 1e-4); out <- matrix(NA_real_, n, n)
  f0 <- f(x)
  for (i in seq_len(n)) {
    ei <- rep(0, n); ei[[i]] <- h[[i]]
    out[i, i] <- (f(x + ei) - 2 * f0 + f(x - ei)) / h[[i]]^2
    if (i > 1) for (j in seq_len(i - 1L)) {
      ej <- rep(0, n); ej[[j]] <- h[[j]]
      value <- (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) /
        (4 * h[[i]] * h[[j]])
      out[i, j] <- value; out[j, i] <- value
    }
  }
  list(hessian = out, steps = h)
}

binomial_witness <- function() {
  d <- read.csv(binomial_path, stringsAsFactors = FALSE)
  d$g <- factor(d$g, levels = sprintf("g%d", sort(unique(as.integer(sub("g", "", d$g))))))
  if (nrow(d) != 90L) stop("binomial workflow frame must have 90 rows")
  formula <- successes / trials ~ 1 + x + (1 | g)
  fit <- function(n_agq, rhs = formula) do.call(glmer, list(formula = rhs, data = d,
    family = binomial("logit"), weights = d$trials, nAGQ = n_agq, control = glmer_controls))
  agq1 <- fit(1L); agq5 <- fit(5L)
  null1 <- fit(1L, successes / trials ~ 1 + (1 | g))
  null5 <- fit(5L, successes / trials ~ 1 + (1 | g))
  if (agq1@optinfo$conv$opt != 0L || agq5@optinfo$conv$opt != 0L || null1@optinfo$conv$opt != 0L || null5@optinfo$conv$opt != 0L)
    stop("one or more binomial fits have nonzero optimizer code")
  devfun <- update(agq1, devFunOnly = TRUE)
  devfun5 <- update(agq5, devFunOnly = TRUE)
  parameter <- c(getME(agq1, "theta"), fixef(agq1))
  parameter5 <- c(getME(agq5, "theta"), fixef(agq5))
  h <- central_hessian(devfun, parameter)
  h5 <- central_hessian(devfun5, parameter5)
  covariance <- 2 * solve(h$hessian)
  covariance5 <- 2 * solve(h5$hessian)
  saturated_shift <- -2 * sum(dbinom(d$successes, d$trials, d$successes / d$trials, log = TRUE))
  one_fit <- function(fit, n_agq, closure) list(
    n_agq = n_agq, beta = array(fixef(fit)), beta_names = as.list(names(fixef(fit))),
    theta = array(getME(fit, "theta")), theta_names = as.list(names(getME(fit, "theta"))),
    minus2log_lme4_raw = -2 * as.numeric(logLik(fit)),
    deviance = as.numeric(deviance(fit)),
    objective_at_estimate = as.numeric(closure(c(getME(fit, "theta"), fixef(fit)))),
    optimizer = list(code = fit@optinfo$conv$opt, evaluations = fit@optinfo$feval,
      messages = fit@optinfo$conv$lme4$messages)
  )
  raw1 <- -2 * as.numeric(logLik(agq1)); raw5 <- -2 * as.numeric(logLik(agq5))
  full5 <- raw5 + saturated_shift
  expected_full5 <- 453.1936138884049
  if (abs(full5 - expected_full5) > 1e-8) stop("nAGQ5 full-density normalization disagrees with retained independent baseline")
  list(
    definition = "Binomial logit GLMM: successes/trials ~ 1 + x + (1|g), proportion response with weights=trials.",
    factor_coding = "g is an ordered explicit factor g1,...,g30; x is numeric.",
    parameter_order = list(lme4_devfun = c(as.list(names(getME(agq1, "theta"))), as.list(names(fixef(agq1)))),
      stored_beta_then_theta = c(as.list(names(fixef(agq1))), as.list(names(getME(agq1, "theta"))))),
    nAGQ1 = one_fit(agq1, 1L, devfun), nAGQ5 = one_fit(agq5, 5L, devfun5),
    full_objective_by_estimator = list(
      nAGQ1 = list(definition = "lme4 nAGQ=1 devFunOnly and logLik use the retained full objective scale",
        raw_minus2loglik = raw1, normalization_shift = 0, full_normalized_minus2loglik = raw1,
        devfun_equals_raw = abs(devfun(parameter) - raw1) < 1e-8),
      nAGQ5 = list(definition = "lme4 high-order AGQ raw -2 logLik plus measured response-only saturated binomial constant",
        raw_minus2loglik = raw5, response_constant_shift = saturated_shift, full_normalized_minus2loglik = full5,
        retained_independent_full_baseline = expected_full5, baseline_residual = full5 - expected_full5)
    ),
    central_difference = list(method = "second-order central finite differences of lme4 nAGQ=1 devFunOnly closure",
      step = array(h$steps), hessian_theta_then_beta = unname(split(h$hessian, row(h$hessian))),
      covariance_2_hessian_inverse_theta_then_beta = unname(split(covariance, row(covariance))),
      covariance_2_hessian_inverse_beta_block = unname(split(covariance[2:3, 2:3, drop = FALSE], row(covariance[2:3, 2:3, drop = FALSE]))),
      lme4_vcov_use_hessian_true_beta_block = unname(split(vcov(agq1, use.hessian = TRUE), row(vcov(agq1, use.hessian = TRUE))))
    ),
    selected_nAGQ5_central_difference = list(method = "second-order central finite differences of lme4 nAGQ=5 devFunOnly closure",
      parameter_order = c(as.list(names(getME(agq5, "theta"))), as.list(names(fixef(agq5)))), step = array(h5$steps),
      hessian_theta_then_beta = unname(split(h5$hessian, row(h5$hessian))),
      covariance_2_hessian_inverse_theta_then_beta = unname(split(covariance5, row(covariance5))),
      covariance_2_hessian_inverse_beta_block = unname(split(covariance5[2:3, 2:3, drop = FALSE], row(covariance5[2:3, 2:3, drop = FALSE]))),
      lme4_vcov_use_hessian_true_beta_block = unname(split(vcov(agq5, use.hessian = TRUE), row(vcov(agq5, use.hessian = TRUE)))),
      objective_at_estimate = devfun5(parameter5)),
    likelihood_ratio_small_null = list(
      null_formula = "successes/trials ~ 1 + (1|g)", alternative_formula = "successes/trials ~ 1 + x + (1|g)",
      nAGQ1 = list(null_minus2loglik = -2 * as.numeric(logLik(null1)), alternative_minus2loglik = raw1,
        statistic = 2 * (as.numeric(logLik(agq1)) - as.numeric(logLik(null1))), df = 1L),
      nAGQ5 = list(null_minus2loglik = -2 * as.numeric(logLik(null5)), alternative_minus2loglik = raw5,
        statistic = 2 * (as.numeric(logLik(agq5)) - as.numeric(logLik(null5))), df = 1L)
    )
  )
}

# Exact marginal normal ML objective.  It is independent of lme4's profiler and
# makes each retained constrained trace point explicit.
normal_ml <- function(y, x, group, beta0, beta1, tau, sigma) {
  if (!is.finite(tau) || !is.finite(sigma) || tau <= 0 || sigma <= 0 || tau > 1e6 || sigma > 1e6) return(1e300)
  n <- length(y); Z <- model.matrix(~ 0 + group)
  V <- tcrossprod(Z) * tau^2 + diag(sigma^2, n)
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(1e300)
  residual <- y - beta0 - beta1 * x
  sum(2 * log(diag(R))) + sum(forwardsolve(t(R), residual)^2) + n * log(2 * pi)
}

boundary_witness <- function() {
  d <- expand.grid(g = sprintf("g%d", 1:12), x = c(-1, 0, 1), KEEP.OUT.ATTRS = FALSE)
  d <- d[order(d$g, d$x), ]; d$g <- factor(d$g, levels = sprintf("g%d", 1:12))
  d$trials <- 20L; d$successes <- rep(c(5L, 8L, 11L), 12L)
  fit <- function(rhs) do.call(glmer, list(formula = rhs, data = d, family = binomial("logit"),
    weights = d$trials, nAGQ = 1L, control = glmer_controls))
  large <- fit(successes / trials ~ 1 + x + (1 | g)); small <- fit(successes / trials ~ 1 + (1 | g))
  list(definition = "Balanced synthetic binomial boundary witness: 12 groups x (-1,0,1), 20 trials, successes (5,8,11) in every group.",
    large = list(formula = "successes/trials ~ 1 + x + (1|g)", beta = array(fixef(large)), theta = array(getME(large, "theta")),
      minus2log_lme4_raw = -2 * as.numeric(logLik(large)), optimizer_code = large@optinfo$conv$opt),
    small = list(formula = "successes/trials ~ 1 + (1|g)", beta = array(fixef(small)), theta = array(getME(small, "theta")),
      minus2log_lme4_raw = -2 * as.numeric(logLik(small)), optimizer_code = small@optinfo$conv$opt),
    likelihood_ratio = 2 * (as.numeric(logLik(large)) - as.numeric(logLik(small))),
    limitation = "Boundary-only likelihood witness; no covariance, standard-error, p-value, or mixture-reference claim is made.")
}

sleepstudy_witness <- function() {
  d <- read.csv(sleepstudy_path, stringsAsFactors = FALSE)
  d$subj <- factor(d$subj, levels = unique(d$subj))
  fit <- lmer(reaction ~ 1 + days + (1 | subj), d, REML = FALSE, control = lmer_controls)
  if (fit@optinfo$conv$opt != 0L) stop("sleepstudy ML lmer optimizer code is nonzero")
  beta <- fixef(fit); sigma <- sigma(fit); tau <- as.numeric(VarCorr(fit)$subj[1, 1])^0.5
  direct_at_fit <- normal_ml(d$reaction, d$days, d$subj, beta[[1]], beta[[2]], tau, sigma)
  profile_result <- tryCatch(profile(fit, which = c("days", ".sig01", ".sigma"), maxpts = 100L), error = identity)
  profile_info <- if (inherits(profile_result, "error")) list(status = "failed", message = conditionMessage(profile_result), intervals = NULL, trace = NULL) else list(
    status = "ok", message = NULL, intervals = unname(split(confint(profile_result, level = 0.95), row(confint(profile_result, level = 0.95)))),
    interval_names = rownames(confint(profile_result, level = 0.95)), trace = profile_result)
  constrain <- function(kind, value) {
    start <- c(beta[[1]], beta[[2]], log(tau), log(sigma))
    if (kind == "beta_days") {
      objective <- function(p) normal_ml(d$reaction, d$days, d$subj, p[[1]], value, exp(p[[2]]), exp(p[[3]]))
      ans <- optim(c(start[[1]], start[[3]], start[[4]]), objective, method = "L-BFGS-B",
        lower = c(-1e4, -10, -10), upper = c(1e4, 10, 10))
      list(parameter = "days", constrained_value = value, minus2loglik = ans$value,
        beta = list(intercept = ans$par[[1]], days = value), random_intercept_sd = exp(ans$par[[2]]), residual_sd = exp(ans$par[[3]]), convergence = ans$convergence)
    } else if (kind == "tau") {
      objective <- function(p) normal_ml(d$reaction, d$days, d$subj, p[[1]], p[[2]], value, exp(p[[3]]))
      ans <- optim(c(start[[1]], start[[2]], start[[4]]), objective, method = "L-BFGS-B",
        lower = c(-1e4, -1e4, -10), upper = c(1e4, 1e4, 10))
      list(parameter = "random_intercept_sd", constrained_value = value, minus2loglik = ans$value,
        beta = list(intercept = ans$par[[1]], days = ans$par[[2]]), random_intercept_sd = value, residual_sd = exp(ans$par[[3]]), convergence = ans$convergence)
    } else {
      objective <- function(p) normal_ml(d$reaction, d$days, d$subj, p[[1]], p[[2]], exp(p[[3]]), value)
      ans <- optim(c(start[[1]], start[[2]], start[[3]]), objective, method = "L-BFGS-B",
        lower = c(-1e4, -1e4, -10), upper = c(1e4, 1e4, 10))
      list(parameter = "residual_sd", constrained_value = value, minus2loglik = ans$value,
        beta = list(intercept = ans$par[[1]], days = ans$par[[2]]), random_intercept_sd = exp(ans$par[[3]]), residual_sd = value, convergence = ans$convergence)
    }
  }
  list(
    definition = "Gaussian ML LMM: reaction ~ 1 + days + (1|subj), scalar random intercept.",
    factor_coding = "subj is an explicit factor in first-appearance CSV order.",
    fit = list(beta = array(beta), beta_names = as.list(names(beta)), random_intercept_sd = tau,
      random_intercept_variance = tau^2, residual_sd = sigma, residual_variance = sigma^2,
      minus2log_lme4 = -2 * as.numeric(logLik(fit)), exact_marginal_normal_minus2loglik = direct_at_fit,
      exact_minus_lme4 = direct_at_fit + 2 * as.numeric(logLik(fit))),
    profile = profile_info,
    variance_transformed_limits = if (identical(profile_info$status, "ok")) list(
      random_intercept_variance = array(profile_info$intervals[[which(profile_info$interval_names == ".sig01")]]^2),
      residual_variance = array(profile_info$intervals[[which(profile_info$interval_names == ".sigma")]]^2)) else NULL,
    constrained_ml_trace = c(lapply(c(8, 10, 12), function(v) constrain("beta_days", v)),
      lapply(c(20, 40), function(v) constrain("tau", v)), lapply(c(20, 30), function(v) constrain("sigma", v)))
  )
}

bin <- binomial_witness(); lmm <- sleepstudy_witness(); boundary <- boundary_witness()
references <- list(schema_version = 1L, generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source = list(R = as.character(getRversion()), lme4 = as.character(packageVersion("lme4")), jsonlite = as.character(packageVersion("jsonlite"))),
  controls = controls, binomial_workflow = bin, sleepstudy_ml_profile = lmm, scalar_boundary_null = boundary)
write_json(references, file.path(root, "references.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17)
provenance <- list(schema_version = 1L, checkpoint = "74239f432790d9f2b90d1864280b66a3bfdc173d",
  command = "Rscript comparison/scripts/freeze_inference.R", frames = list(
    binomial_workflow = list(path = binomial_path, sha256 = sha256(binomial_path), rows = 90L),
    sleepstudy = list(path = sleepstudy_path, sha256 = sha256(sleepstudy_path), rows = 180L)),
  artifacts = list(references_json_sha256 = sha256(file.path(root, "references.json"))),
  generator = list(path = "comparison/scripts/freeze_inference.R", sha256 = sha256("comparison/scripts/freeze_inference.R")),
  estimator_contract = "R/lme4 fits are independent witnesses. GLMM curvature retains separate lme4 deviance Hessians for nAGQ=1 and selected nAGQ=5 in theta,beta order; their likelihood normalizations are order-specific. LMM trace uses an exact marginal Gaussian ML likelihood, not Scala code.")
write_json(provenance, file.path(root, "reference-provenance.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17)
cat("Frozen R likelihood-inference witnesses\n")
