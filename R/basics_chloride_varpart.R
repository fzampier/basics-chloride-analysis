# basics_chloride_varpart.R — Analysis 3: Bayesian longitudinal variance partition.
# Partitions the variance in serum chloride (level, days 1-3, ANCOVA on baseline) into
# fluid (postbaseline study-fluid load + total volume) / study day / patient
# (baseline Cl, age, creatinine, SOFA) / fixed-predictor covariance /
# between-patient / residual.
# Requires Dados_BASICS_MainPaper.RData in the working dir. Saves chl_varpart_fit.rds (git-ignored).
suppressMessages(library(brms))
if (!exists("BASICS_PATHS", inherits = TRUE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grepl("^--file=", .args)][[1]])
  source(file.path(dirname(normalizePath(.file)), "project_paths.R"))
}
P <- get("BASICS_PATHS", inherits = TRUE)
source(file.path(P$root, "R", "exposure_construction.R"))
load(P$data_file)
d<-dados
d$arm<-ifelse(d$grupo=="Saline Solution","Saline","Plasma-Lyte")
d$age<-as.numeric(difftime(d$dtarrol,d$dtnasc,units="days"))/365.25
d<-add_postbaseline_exposures(d)
long<-build_varpart_long(d)
long<-long[complete.cases(long),]
long<-scale_varpart_predictors(long)

fit<-brm(clev ~ load+vol+day+clobas+age+creabas+sofa + (1|id), data=long,
         chains=4, iter=4000, warmup=2000, cores=4, seed=1, refresh=0,
         prior=prior(normal(0,5),class="b"))
saveRDS(fit, P$varpart_fit)

diag_vars<-c(
  "b_load","b_vol","b_day2","b_day3","b_clobas","b_age","b_creabas","b_sofa",
  "sd_id__Intercept","sigma"
)
diag<-posterior::summarise_draws(
  posterior::as_draws_array(fit,variable=diag_vars),
  "rhat","ess_bulk","ess_tail"
)
nuts<-brms::nuts_params(fit)
cat(sprintf(
  "Diagnostics: max Rhat %.4f; min bulk ESS %.0f; min tail ESS %.0f; divergences %d\n",
  max(diag$rhat,na.rm=TRUE), min(diag$ess_bulk,na.rm=TRUE),
  min(diag$ess_tail,na.rm=TRUE),
  sum(nuts$Parameter=="divergent__" & nuts$Value==1)
))

# variance shares (% of total) with 95% credible intervals, computed from the posterior
post<-as_draws_df(fit)
Xf<-as.matrix(long[,c("load","vol")]); Xp<-as.matrix(long[,c("clobas","age","creabas","sofa")]); Xt<-model.matrix(~day,long)[,-1,drop=FALSE]
bf<-as.matrix(post[,c("b_load","b_vol")]); bp<-as.matrix(post[,c("b_clobas","b_age","b_creabas","b_sofa")]); bt<-as.matrix(post[,c("b_day2","b_day3")])
# Compute the variance of each posterior linear predictor from the predictor
# covariance matrix. This is algebraically identical to materializing every
# observation-by-draw prediction, but uses far less memory.
Vf<-rowSums((bf%*%cov(Xf))*bf); Vp<-rowSums((bp%*%cov(Xp))*bp); Vt<-rowSums((bt%*%cov(Xt))*bt)
Xfull<-cbind(Xf,Xp,Xt); bfull<-cbind(bf,bp,bt)
Vfull<-rowSums((bfull%*%cov(Xfull))*bfull)
Vcov<-Vfull-Vf-Vp-Vt
Vbp<-post$sd_id__Intercept^2; Vr<-post$sigma^2; Vtot<-Vfull+Vbp+Vr
qd<-function(x) sprintf("%.1f%% [%.1f, %.1f]",100*median(x),100*quantile(x,.025),100*quantile(x,.975))
cat("Fluid      :",qd(Vf/Vtot),"\nStudy day  :",qd(Vt/Vtot),"\nCovariance :",qd(Vcov/Vtot),"\nPatient    :",qd(Vp/Vtot),"\nBetween    :",qd(Vbp/Vtot),"\nResidual   :",qd(Vr/Vtot),"\n")

shares<-list(fluid=Vf/Vtot,study_day=Vt/Vtot,covariance=Vcov/Vtot,patient=Vp/Vtot,between_patient=Vbp/Vtot,residual=Vr/Vtot)
result<-do.call(rbind,lapply(names(shares),function(component){
  x<-100*shares[[component]]
  data.frame(analysis="Primary level model",outcome="Serum chloride level on days 1-3",component=component,observations=nrow(long),patients=length(unique(long$id)),estimate_percent=median(x),conf_low_percent=unname(quantile(x,.025)),conf_high_percent=unname(quantile(x,.975)),stringsAsFactors=FALSE)
}))
write.csv(result,file.path(P$reproduction,"varpart_primary_results.csv"),row.names=FALSE)
