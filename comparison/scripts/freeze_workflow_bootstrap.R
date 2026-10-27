suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(jsonlite))
root <- 'comparison/workflow-tranche'
data <- read.csv(file.path(root,'frames','sleepstudy.csv'),check.names=FALSE)
data$subj <- factor(data$subj,levels=unique(data$subj))
responses <- read.csv(file.path(root,'bootstrap-responses.csv'),check.names=FALSE)
stopifnot(nrow(responses)==nrow(data),ncol(responses)==39L)
provenance <- read_json(file.path(root,'bootstrap-response-provenance.json'),simplifyVector=TRUE)
results <- lapply(seq_len(ncol(responses)),function(i) {
  data$reaction <- responses[[i]]
  fit <- lmer(reaction ~ 1 + days + (1 | subj),data,REML=FALSE,
    control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000L,rhoend=1e-9)))
  stopifnot(is.null(fit@optinfo$conv$lme4$messages),fit@optinfo$conv$opt==0L)
  list(index=i-1L,seed=provenance$responses$seed[i],beta=unname(fixef(fit)),sigma=sigma(fit),
    theta=unname(getME(fit,'theta')),objective=as.numeric(-2*logLik(fit)))
})
statistics <- vapply(results,function(r) r$beta[2],numeric(1))
result <- list(schema_version=1L,source=list(R=as.character(getRversion()),lme4=as.character(packageVersion('lme4')),
  generated_utc=format(Sys.time(),tz='UTC',usetz=TRUE),responses_sha256=provenance$csv_sha256,
  note='Refitted actual Scala-generated responses; R RNG seeds were not substituted.'),results=results,
  percentile=list(level=0.95,quantile_type=7L,days=unname(quantile(statistics,c(0.025,0.975),type=7)),
    two_days=unname(quantile(2*statistics,c(0.025,0.975),type=7))))
write_json(result,file.path(root,'bootstrap-lme4.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('Refitted all39 actual response draws and froze R type-7 intervals:',result$percentile$days,'\n')
