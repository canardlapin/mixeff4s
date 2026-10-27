#!/usr/bin/env Rscript
# Independent one-dimensional marginal-likelihood witness for the GLMM analysis
# tranche.  This intentionally does not call any mixeff4s code or reuse PIRLS.
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(jsonlite))

root <- "comparison/glmm-analysis-tranche"
frames <- "comparison/glmm-tranche/frames"
dir.create(root, recursive = TRUE, showWarnings = FALSE)

as_array <- function(x) as.list(unname(as.vector(x)))

read_frame <- function(name, levels) {
  d <- read.csv(file.path(frames, paste0(name, ".csv")), check.names = FALSE,
                stringsAsFactors = FALSE)
  for (nm in names(levels)) d[[nm]] <- factor(d[[nm]], levels = levels[[nm]])
  d
}

specs <- list(
  list(dataset = "cbpp", formula = "incidence / size ~ period + (1 | herd)",
       family = binomial("logit"), group = "herd", response = "incidence", trials = "size",
       weights = "size", levels = list(herd = as.character(1:15), period = as.character(1:4))),
  list(dataset = "culcitalogreg", formula = "predation ~ ttt + (1 | block)",
       family = binomial("logit"), group = "block", response = "predation", trials = NULL,
       levels = list(block = as.character(1:10), ttt = c("none", "crabs", "shrimp", "both"))),
  list(dataset = "synthetic_poisson_offset", formula = "count ~ x + offset(log(exposure)) + (1 | group)",
       family = poisson("log"), group = "group", response = "count", trials = NULL,
       levels = list(group = sprintf("g%d", 1:6)))
)

controls <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L, rhoend = 1e-9),
                         tolPwrss = 1e-10)

# Integrate exp(log conditional density + standard-normal RE density) per
# independent group.  The scale is found from the conditional mode so underflow
# cannot turn a valid small integral into zero.
group_integral <- function(log_density, mode) {
  scale <- log_density(mode)
  # integrate supplies vector abscissae; evaluate one scalar conditional
  # likelihood per abscissa rather than recycling observations over quadrature.
  f <- Vectorize(function(u) exp(log_density(u) - scale), "u")
  z <- integrate(f, lower = -Inf, upper = Inf, rel.tol = 1e-11, abs.tol = 1e-13,
                 subdivisions = 2000L, stop.on.error = FALSE)
  if (!is.finite(z$value) || z$value <= 0 || z$message != "OK")
    stop(sprintf("independent integration failed: value=%s abs.error=%s message=%s",
                 z$value, z$abs.error, z$message))
  list(log_value = scale + log(z$value), abs_error = z$abs.error * exp(scale),
       relative_error_bound = z$abs.error / z$value, subdivisions = z$subdivisions,
       message = z$message, scale = scale)
}

full_integral <- function(spec, d, beta, theta, modes) {
  group <- as.character(d[[spec$group]])
  groups <- levels(d[[spec$group]])
  # Model matrices deliberately use R's treatment coding with predeclared levels.
  # This is the same coding stated in the frozen reference, but not lme4's design.
  fixed_rhs <- if (spec$dataset == "cbpp") ~ period else if (spec$dataset == "culcitalogreg") ~ ttt else ~ x
  X <- model.matrix(fixed_rhs, d)
  eta_fixed <- drop(X %*% beta)
  if (spec$dataset == "synthetic_poisson_offset") eta_fixed <- eta_fixed + log(d$exposure)
  sd_re <- theta[[1]]
  rows <- lapply(seq_along(groups), function(i) {
    ix <- which(group == groups[[i]])
    log_density <- function(u) {
      eta <- eta_fixed[ix] + sd_re * u
      if (spec$family$family == "binomial") {
        n <- if (is.null(spec$trials)) rep(1, length(ix)) else d[[spec$trials]][ix]
        sum(dbinom(d[[spec$response]][ix], n, plogis(eta), log = TRUE)) + dnorm(u, log = TRUE)
      } else {
        sum(dpois(d[[spec$response]][ix], exp(eta), log = TRUE)) + dnorm(u, log = TRUE)
      }
    }
    # Do not trust a fitted conditional mode as an integration scale: at a
    # boundary theta can make lme4's stored spherical mode poorly conditioned.
    optimum <- optimize(function(u) -log_density(u), interval = c(-12, 12), tol = 1e-13)
    result <- group_integral(log_density, optimum$minimum)
    c(list(level = groups[[i]], mode_standard_normal = optimum$minimum,
           fitted_mode_standard_normal = modes[[i]]), result)
  })
  log_value <- sum(vapply(rows, `[[`, 0.0, "log_value"))
  list(minus2log_full_density = -2 * log_value, groups = rows,
       max_relative_error_bound = max(vapply(rows, `[[`, 0.0, "relative_error_bound")),
       sum_absolute_error_bound = sum(vapply(rows, `[[`, 0.0, "abs_error")),
       tolerance = list(relative = 1e-11, absolute = 1e-13, subdivisions = 2000L))
}

saturated_constant <- function(spec, d) {
  if (spec$family$family == "binomial") {
    n <- if (is.null(spec$trials)) rep(1, nrow(d)) else d[[spec$trials]]
    -2 * sum(dbinom(d[[spec$response]], n, d[[spec$response]] / n, log = TRUE))
  } else -2 * sum(dpois(d[[spec$response]], d[[spec$response]], log = TRUE))
}

fit_one <- function(spec, n_agq) {
  d <- read_frame(spec$dataset, spec$levels)
  args <- list(formula = as.formula(spec$formula), data = d, family = spec$family,
               nAGQ = n_agq, control = controls)
  if (!is.null(spec$weights)) args$weights <- d[[spec$weights]]
  fit <- do.call(glmer, args)
  if (fit@optinfo$conv$opt != 0L) stop(sprintf("%s nAGQ=%d optimizer code %s", spec$dataset, n_agq, fit@optinfo$conv$opt))
  eta <- predict(fit, type = "link")
  population_eta <- predict(fit, re.form = NA, type = "link")
  list(n_agq = n_agq, beta = as_array(fixef(fit)), beta_names = as.list(names(fixef(fit))),
       theta = as_array(getME(fit, "theta")), theta_names = as.list(names(getME(fit, "theta"))),
       random_effect_standard_normal_modes = as_array(getME(fit, "u")),
       conditional_eta = as_array(eta), conditional_mean = as_array(fitted(fit)),
       population_eta = as_array(population_eta), population_mean = as_array(spec$family$linkinv(population_eta)),
       logLik = as.numeric(logLik(fit)), minus2log_lme4_raw = as.numeric(-2 * logLik(fit)),
       deviance = as.numeric(deviance(fit)), optimizer = list(code = fit@optinfo$conv$opt,
       evaluations = fit@optinfo$feval, messages = fit@optinfo$conv$lme4$messages))
}

# Fixed-parameter scalar/slope oracles exercise the integrator without borrowing
# fitted parameters.  They use both response families and include zero outcomes.
oracle <- function(name, family, y, x, beta, theta, exposure = NULL, trials = NULL, z = rep(1, length(x))) {
  log_density <- function(u) {
    eta <- beta[[1]] + beta[[2]] * x + z * theta * u + if (is.null(exposure)) 0 else log(exposure)
    if (family == "binomial") sum(dbinom(y, trials, plogis(eta), log = TRUE)) + dnorm(u, log = TRUE)
    else sum(dpois(y, exp(eta), log = TRUE)) + dnorm(u, log = TRUE)
  }
  opt <- optimize(function(u) -log_density(u), c(-12, 12), tol = 1e-13)
  integral_result <- group_integral(log_density, opt$minimum)
  list(name = name, family = family, beta = as_array(beta), theta = theta, x = as_array(x), z = as_array(z),
       y = as_array(y), exposure = if (is.null(exposure)) NULL else as_array(exposure),
       trials = if (is.null(trials)) NULL else as_array(trials),
       mode_standard_normal = opt$minimum, minus2log_full_density = -2 * integral_result$log_value,
       relative_error_bound = integral_result$relative_error_bound, absolute_error_bound = integral_result$abs_error)
}

all_fits <- list(); integrals <- list()
for (spec in specs) {
  fits <- lapply(c(1L, 5L, 9L, 15L), function(n) fit_one(spec, n))
  names(fits) <- paste0("nAGQ", c(1, 5, 9, 15))
  d <- read_frame(spec$dataset, spec$levels)
  high <- fits[["nAGQ15"]]
  full <- full_integral(spec, d, unlist(high$beta), unlist(high$theta), unlist(high$random_effect_standard_normal_modes))
  sat <- saturated_constant(spec, d)
  # lme4's high-order AGQ logLik uses the deviance-residual normalization.
  # Store both scales; their difference must be the response-only saturated term.
  high$minus2log_full_density_normalized <- high$minus2log_lme4_raw + sat
  high$full_density_normalization_shift_expected <- sat
  high$full_density_normalization_shift_measured <- full$minus2log_full_density - high$minus2log_lme4_raw
  high$full_density_normalization_residual <- full$minus2log_full_density - (high$minus2log_lme4_raw + sat)
  fits[["nAGQ15"]] <- high
  all_fits[[spec$dataset]] <- fits
  integrals[[spec$dataset]] <- c(list(dataset = spec$dataset, family = spec$family$family,
    response_only_saturated_minus2loglik = sat, n = nrow(d)), full)
}

oracles <- list(
  oracle("binomial_scalar_slope_unequal_trials", "binomial", c(0, 2, 1, 4), c(-1, 0, .5, 1), c(-.4, .7), .9, trials = c(2, 3, 4, 5), z = c(-1, 0, .5, 1)),
  oracle("poisson_scalar_slope_offset_zero", "poisson", c(0, 1, 3, 0), c(-1, 0, .5, 1), c(.2, .35), .65, exposure = c(.5, 1, 1.5, 2), z = c(-1, 0, .5, 1)),
  oracle("binomial_zero_theta", "binomial", c(0, 2, 1, 4), c(-1, 0, .5, 1), c(-.4, .7), 0, trials = c(2, 3, 4, 5)),
  oracle("poisson_zero_theta_offset", "poisson", c(0, 1, 3, 0), c(-1, 0, .5, 1), c(.2, .35), 0, exposure = c(.5, 1, 1.5, 2))
)

write_json(list(schema_version = 1L, generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source = list(R = as.character(getRversion()), lme4 = as.character(packageVersion("lme4")),
                jsonlite = as.character(packageVersion("jsonlite"))),
  controls = list(optimizer = "bobyqa", maxfun = 100000L, rhoend = 1e-9, tolPwrss = 1e-10),
  factor_coding = "R treatment contrasts with explicit levels listed per dataset", fits = all_fits),
  file.path(root, "references.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17)
write_json(list(schema_version = 1L, integration_definition = "-2 log integral product_i p(y_i|u_g,beta,theta) phi(u_g) du_g; u is standard normal and theta scales it",
  integrals = integrals, fixed_parameter_oracles = oracles), file.path(root, "integrals.json"),
  auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17)
cat("Frozen independent AGQ and full-density integration witnesses\n")
