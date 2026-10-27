#!/usr/bin/env Rscript
# Independent R 4.5.1 calibration fixture; base stats only.
binom <- list(c(100, .5), c(10000, .7), c(2147483647, 5e-9), c(2147483647, .5))
pois <- c(9.9, 10, 10.1, 1000, 1e6)
row <- function(kind, a, b = NA_real_) {
  qs <- if (kind == "binomial") qbinom(c(.01,.5,.99), a, b) else qpois(c(.01,.5,.99), a)
  ps <- if (kind == "binomial") pbinom(qs, a, b) else ppois(qs, a)
  sprintf('{"kind":"%s","a":%.17g,"b":%s,"q":[%s],"cdf":[%s]}', kind, a,
    ifelse(is.na(b), "null", format(b, digits=17)), paste(format(qs, scientific=FALSE, trim=TRUE), collapse=","),
    paste(format(ps, digits=17, scientific=FALSE, trim=TRUE), collapse=","))
}
writeLines(paste0('{"r_version":"', R.version.string, '","cases":[',
  paste(c(vapply(binom, function(x) row("binomial",x[1],x[2]), ""), vapply(pois, function(x) row("poisson",x), "")), collapse=","), ']}'))
