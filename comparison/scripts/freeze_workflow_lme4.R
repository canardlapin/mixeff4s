suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(jsonlite))
root <- 'comparison/workflow-tranche'
specs <- list(
  sleepstudy=list(formula='reaction ~ 1 + days + (1 + days | subj)', factors='subj'),
  penicillin=list(formula='diameter ~ 1 + (1 | plate) + (1 | sample)', factors=c('plate','sample')),
  pastes=list(formula='strength ~ 1 + (1 | batch / cask)', factors=c('batch','cask'))
)
rows <- list()
for (dataset in names(specs)) {
  spec <- specs[[dataset]]
  data <- read.csv(file.path(root, 'frames', paste0(dataset,'.csv')), check.names=FALSE)
  for (name in spec$factors) data[[name]] <- factor(data[[name]], levels=unique(data[[name]]))
  for (criterion in c('ML','REML')) {
    fit <- lmer(as.formula(spec$formula), data, REML=criterion=='REML',
      control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000L, rhoend=1e-9)))
    stopifnot(is.null(fit@optinfo$conv$lme4$messages), fit@optinfo$conv$opt==0L)
    covariance <- lapply(VarCorr(fit), function(x) unname(as.matrix(x)))
    modes <- lapply(ranef(fit), function(x) list(levels=rownames(x), names=colnames(x), values=unname(as.matrix(x))))
    rows[[length(rows)+1L]] <- list(dataset=dataset, formula=spec$formula, criterion=criterion,
      n=nobs(fit), objective=as.numeric(-2*logLik(fit)), beta=unname(fixef(fit)), coef_names=names(fixef(fit)),
      theta=unname(getME(fit,'theta')), sigma=sigma(fit), vcov=unname(as.matrix(vcov(fit))),
      conditional=unname(predict(fit)), population=unname(predict(fit,re.form=NA)), residuals=unname(residuals(fit)),
      random_covariance=covariance, modes=modes, singular=isSingular(fit),
      optimizer=list(name='bobyqa', code=fit@optinfo$conv$opt, evaluations=fit@optinfo$feval))
  }
}
result <- list(schema_version=1L, source=list(R=as.character(getRversion()),lme4=as.character(packageVersion('lme4')),
  jsonlite=as.character(packageVersion('jsonlite')), generated_utc=format(Sys.time(),tz='UTC',usetz=TRUE)), results=rows)
write_json(result, file.path(root,'lme4-workflow.json'), auto_unbox=TRUE, pretty=TRUE, digits=16)
writeLines(capture.output(sessionInfo()), file.path(root,'r-session.txt'))
cat('Frozen',length(rows),'classic LMM fits with predictions, covariance, residuals and conditional modes\n')
