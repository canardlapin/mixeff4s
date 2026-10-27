# Synthetic, reproducible grouped-binomial workflow fixture. No restricted data.
set.seed(2026102701)
frame <- expand.grid(time = 0:2, group = seq_len(30), KEEP.OUT.ATTRS = FALSE)
frame$x <- frame$time - 1
frame$trials <- 20L
random <- rnorm(30, sd = 1.1)
frame$successes <- rbinom(nrow(frame), frame$trials, plogis(-0.3 + 0.5 * frame$x + random[frame$group]))
frame$group <- paste0('g', frame$group)
dir.create('comparison/glmm-analysis-tranche/frames', recursive = TRUE, showWarnings = FALSE)
write.csv(frame, 'comparison/glmm-analysis-tranche/frames/binomial_workflow.csv', row.names = FALSE)
cat('Generated', nrow(frame), 'grouped-binomial rows with seed 2026102701; R', as.character(getRversion()), '\n')
