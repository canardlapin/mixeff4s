#!/usr/bin/env Rscript
# Generate independent lme4 GLMM pins.  Tests consume the generated Scala only.
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(jsonlite))

root <- "comparison/glmm-tranche"
frames <- file.path(root, "frames")
dir.create(frames, recursive = TRUE, showWarnings = FALSE)

# Compact offset/zero-count fixture.  It stays separate from the larger
# correlated/crossed case below so scalar-offset behavior has a short witness.
poisson_small <- data.frame(
  group = rep(sprintf("g%d", 1:6), each = 4),
  x = rep(c(-1, 0, 0.5, 1), 6),
  exposure = rep(c(0.5, 1, 1.5, 2), 6),
  count = c(0, 1, 1, 3, 0, 1, 2, 4, 1, 2, 3, 5, 0, 2, 2, 6, 1, 3, 4, 8, 0, 1, 3, 7)
)
write.csv(poisson_small, file.path(frames, "synthetic_poisson_offset.csv"), row.names = FALSE, quote = TRUE)

# This deliberately small frame has zero counts, a log exposure offset, a
# correlated intercept/slope term, and an independently crossed batch effect.
poisson <- expand.grid(site = sprintf("s%d", 1:12), batch = sprintf("b%d", 1:5), visit = 0:4,
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
poisson <- poisson[order(poisson$site, poisson$batch, poisson$visit), ]
poisson$x <- rep(c(-1, -0.5, 0, 0.5, 1), length.out = nrow(poisson))
poisson$exposure <- rep(c(0.5, 1, 1.5, 2, 0.75), length.out = nrow(poisson))
site_intercept <- seq(-0.8, 0.85, length.out = 12)
site_slope <- c(-0.6, 0.2, -0.25, 0.65, -0.1, 0.4, -0.45, 0.1, 0.55, -0.35, 0.7, 0)
batch_intercept <- c(-0.4, -0.15, 0.05, 0.25, 0.5)
set.seed(90210)
poisson$count <- rpois(nrow(poisson), poisson$exposure * exp(
  1.0 + 0.35 * poisson$x + site_intercept[match(poisson$site, sprintf("s%d", 1:12))] +
    site_slope[match(poisson$site, sprintf("s%d", 1:12))] * poisson$x +
    batch_intercept[match(poisson$batch, sprintf("b%d", 1:5))]
))
poisson$count[c(1, 17, 58, 121, 233)] <- 0L
write.csv(poisson, file.path(frames, "synthetic_poisson_crossed.csv"), row.names = FALSE, quote = TRUE)

read_frame <- function(name, categorical) {
  d <- read.csv(file.path(frames, paste0(name, ".csv")), check.names = FALSE, stringsAsFactors = FALSE)
  for (col in names(categorical)) d[[col]] <- factor(d[[col]], levels = categorical[[col]])
  d
}

specs <- list(
  list(
    dataset = "cbpp", formula = "incidence / size ~ period + (1 | herd)", family = binomial(link = "logit"),
    levels = list(herd = as.character(1:15), period = as.character(1:4)), weights = "size"
  ),
  list(
    dataset = "culcitalogreg", formula = "predation ~ ttt + (1 | block)", family = binomial(link = "logit"),
    levels = list(block = as.character(1:10), ttt = c("none", "crabs", "shrimp", "both"))
  ),
  list(
    dataset = "synthetic_poisson_offset", formula = "count ~ x + offset(log(exposure)) + (1 | group)",
    family = poisson(link = "log"), levels = list(group = sprintf("g%d", 1:6))
  ),
  list(
    dataset = "synthetic_poisson_crossed", formula = "count ~ x + offset(log(exposure)) + (1 + x | site) + (1 | batch)",
    family = poisson(link = "log"), levels = list(site = sprintf("s%d", 1:12), batch = sprintf("b%d", 1:5))
  )
)

array <- function(x) as.list(unname(as.vector(x)))
mode_payload <- function(fit) lapply(ranef(fit), function(x) list(
  levels = rownames(x), names = colnames(x), values = unname(as.matrix(x))
))
fit_payload <- function(spec, n_agq) {
  d <- read_frame(spec$dataset, spec$levels)
  control <- glmerControl(
    optimizer = "bobyqa", optCtrl = list(maxfun = 100000L, rhoend = 1e-9),
    tolPwrss = 1e-10
  )
  args <- list(formula = as.formula(spec$formula), data = d, family = spec$family, nAGQ = n_agq, control = control)
  if (!is.null(spec$weights)) args$weights <- d[[spec$weights]]
  fit <- do.call(glmer, args)
  # A boundary is a valid fitted reference when the optimizer itself reports
  # success; preserve its lme4 diagnostic below rather than discarding it.
  if (fit@optinfo$conv$opt != 0L)
    stop(sprintf("%s nAGQ=%d optimizer failed with code %d", spec$dataset, n_agq, fit@optinfo$conv$opt))
  theta <- getME(fit, "theta")
  eta <- predict(fit, type = "link")
  eta_population <- predict(fit, re.form = NA, type = "link")
  list(
    dataset = spec$dataset, formula = spec$formula,
    family = if (spec$family$family == "binomial") "Binomial" else "Poisson",
    link = if (spec$family$link == "logit") "Logit" else "Log",
    n_agq = n_agq,
    approximation = if (n_agq == 1L) "joint_laplace" else "documented_fast_divergence",
    status = "ok", n = nobs(fit),
    objective_minus2loglik_with_constants = as.numeric(-2 * logLik(fit)),
    # lme4's residual/conditional deviance at the final conditional modes.
    conditional_deviance_at_final_modes = as.numeric(deviance(fit)),
    beta = array(fixef(fit)), coef_names = as.list(names(fixef(fit))),
    theta = array(theta), theta_names = as.list(names(theta)),
    random_covariance = lapply(VarCorr(fit), function(x) unname(as.matrix(x))),
    conditional_eta = array(eta), conditional_mu = array(fitted(fit)),
    population_eta = array(eta_population), population_mu = array(spec$family$linkinv(eta_population)),
    modes = mode_payload(fit), singular = isSingular(fit),
    optimizer = list(name = "bobyqa", maxfun = 100000L, rhoend = 1e-9,
                     code = fit@optinfo$conv$opt, evaluations = fit@optinfo$feval,
                     lme4_messages = fit@optinfo$conv$lme4$messages)
  )
}

results <- unlist(lapply(specs, function(s) list(fit_payload(s, 1L), fit_payload(s, 0L))), recursive = FALSE)
out <- list(
  schema_version = 1L,
  source = list(R = as.character(getRversion()), lme4 = as.character(packageVersion("lme4")),
                jsonlite = as.character(packageVersion("jsonlite")),
                generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                controls = "glmerControl(optimizer=bobyqa, maxfun=100000, rhoend=1e-9, tolPwrss=1e-10)"),
  results = results
)
write_json(out, file.path(root, "references.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
cat(sprintf("Frozen %d GLMM fits: %s\n", length(results), paste(vapply(results, `[[`, "", "dataset"), collapse = ", ")))
