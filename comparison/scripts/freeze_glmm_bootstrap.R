#!/usr/bin/env Rscript
# Freeze lme4 refits for the supplied Scala bootstrap responses.  The response
# CSV is an input witness, never generated or reseeded here.
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(jsonlite))

root <- "comparison/glmm-analysis-tranche"
frame_path <- file.path(root, "frames", "binomial_workflow.csv")
responses_path <- file.path(root, "bootstrap-responses.csv")
out_path <- file.path(root, "bootstrap-lme4.json")
control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L, rhoend = 1e-9), tolPwrss = 1e-10)

sha256 <- function(path) tail(strsplit(system2("openssl", c("dgst", "-sha256", path), stdout = TRUE), " ")[[1]], 1)
as_array <- function(x) as.list(unname(as.vector(x)))

frame <- read.csv(frame_path, check.names = FALSE, stringsAsFactors = FALSE)
frame$group <- factor(frame$group, levels = sprintf("g%d", 1:30))
stopifnot(nrow(frame) == 90L, identical(levels(frame$group), sprintf("g%d", 1:30)))
response_classes <- rep("numeric", 92L)
names(response_classes) <- c("index", "seed", paste0("y", 0:89))
response_classes[c("index", "seed")] <- c("integer", "character")
responses <- read.csv(responses_path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = response_classes)
stopifnot(nrow(responses) == 39L, identical(names(responses)[3:92], paste0("y", 0:89)))

saturated_minus2loglik <- function(proportion, trials) {
  successes <- as.integer(round(proportion * trials))
  if (any(abs(successes / trials - proportion) > 1e-12)) stop("bootstrap response is not an integer/trials proportion")
  -2 * sum(dbinom(successes, trials, successes / trials, log = TRUE))
}

fit_one <- function(proportion, identity, seed) {
  d <- frame
  d$proportion <- proportion
  warnings <- character()
  result <- tryCatch(
    withCallingHandlers(
      glmer(proportion ~ 1 + x + (1 | group), data = d, weights = trials, family = binomial("logit"),
            nAGQ = 5L, control = control),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) return(list(index = identity, seed = seed, status = "failure",
    failure = list(class = class(result), message = conditionMessage(result)), warnings = as.list(warnings)))
  beta <- fixef(result)
  theta <- getME(result, "theta")
  raw <- as.numeric(-2 * logLik(result))
  correction <- saturated_minus2loglik(proportion, d$trials)
  list(index = identity, seed = seed, status = "ok", n = nobs(result),
       beta = as_array(beta), beta_names = as.list(names(beta)), theta = as_array(theta), theta_names = as.list(names(theta)),
       conditional_mean = as_array(fitted(result)),
       population_mean_x = list(x = list(-1, 0, 1), mean = as_array(plogis(beta[[1]] + beta[[2]] * c(-1, 0, 1)))),
       minus2log_lme4_raw = raw, response_only_saturated_minus2loglik = correction,
       minus2log_full_density_normalized = raw + correction,
       optimizer = list(code = result@optinfo$conv$opt, evaluations = result@optinfo$feval,
         lme4_messages = result@optinfo$conv$lme4$messages, warnings = as.list(warnings), singular = isSingular(result)))
}

base <- fit_one(frame$successes / frame$trials, "base", NA_character_)
replicates <- lapply(seq_len(nrow(responses)), function(i) fit_one(as.numeric(responses[i, 3:92]), responses$index[[i]], responses$seed[[i]]))
successes <- Filter(function(x) identical(x$status, "ok"), replicates)
if (length(successes) != 39L) warning(sprintf("retained %d bootstrap fitting failures", 39L - length(successes)))
coefficient <- vapply(successes, function(x) x$beta[[2]], 0.0)
contrast <- 2 * coefficient
interval <- function(x) as_array(quantile(x, c(.025, .975), type = 7, names = FALSE))

payload <- list(schema_version = 1L, source = list(R = as.character(getRversion()), lme4 = as.character(packageVersion("lme4")),
  jsonlite = as.character(packageVersion("jsonlite")), generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  controls = "glmerControl(optimizer=bobyqa, maxfun=100000, rhoend=1e-9, tolPwrss=1e-10)", nAGQ = 5L),
  factor_coding = list(group = sprintf("g%d", 1:30), contrasts = "R treatment contrasts; x numeric"),
  base_fit = base, bootstrap = list(requested = 39L, numerical_successes = length(successes), failures = 39L - length(successes),
    coefficient_x = as_array(coefficient), contrast_2x = as_array(contrast), percentile_type = 7L,
    percentile_95 = list(x = interval(coefficient), `2x` = interval(contrast)), replicates = replicates))
write_json(payload, out_path, auto_unbox = TRUE, pretty = TRUE, digits = 17, null = "null")

provenance <- list(schema_version = 1L, inputs_sha256 = list(
  binomial_workflow_csv = sha256(frame_path), bootstrap_responses_csv = sha256(responses_path),
  bootstrap_response_provenance_json = sha256(file.path(root, "bootstrap-response-provenance.json")),
  generator = sha256("comparison/scripts/freeze_glmm_bootstrap.R")),
  source = list(R = as.character(getRversion()), lme4 = as.character(packageVersion("lme4")),
                jsonlite = as.character(packageVersion("jsonlite"))),
  response_policy = "Read supplied Scala DrawRandomEffects rows unchanged; seed is character to prevent integer precision loss.")
write_json(provenance, file.path(root, "bootstrap-r-provenance.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("Frozen base fit and %d/%d bootstrap refits\n", length(successes), nrow(responses)))
