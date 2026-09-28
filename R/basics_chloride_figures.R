# =====================================================================
# BaSICS chloride paper — Figures 1-4 (reproducible, single source)
#   Fig 1 = severity and chloride rise (descriptive, by APACHE quartile)
#   Fig 2 = chloride trajectory and randomized arm difference by eGFR (brms)
#   Fig 3 = randomized normalization and within-saline volume association
#   Fig 4 = Bayesian variance decomposition
# Run basics_chloride_analysis.R (.chl_cache.rds) and basics_chloride_brms.R
# (chl_brms_fit.rds) first.
# =====================================================================
suppressMessages(suppressWarnings(library(brms)))
if (!exists("BASICS_PATHS", inherits = TRUE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grepl("^--file=", .args)][[1]])
  source(file.path(dirname(normalizePath(.file)), "project_paths.R"))
}
P <- get("BASICS_PATHS", inherits = TRUE)
FIG_OUT <- Sys.getenv("BASICS_FIGURE_DIR", unset = P$figures)
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)
C<-readRDS(P$cache_file)
source(file.path(P$root, "R", "exposure_construction.R"))
d<-C$d; meas<-C$meas; m2<-C$m2; col_s<-C$col_s; col_p<-C$col_p; col_acc<-C$col_acc; z<-zero_if_missing; d$id<-1:nrow(d)
mse<-function(x){x<-x[!is.na(x)];c(m=mean(x),lo=mean(x)-1.96*sd(x)/sqrt(length(x)),hi=mean(x)+1.96*sd(x)/sqrt(length(x)))}
egfr<-function(scr,age,fem){k<-ifelse(fem==1,0.7,0.9);a<-ifelse(fem==1,-0.241,-0.302);142*pmin(scr/k,1)^a*pmax(scr/k,1)^(-1.2)*0.9938^age*ifelse(fem==1,1.012,1)}
d$eGFR<-egfr(d$creabas,d$age,d$female); d$eGFR[!is.finite(d$eGFR)|d$eGFR<5|d$eGFR>200]<-NA
for(k in 1:3) d[[paste0("dCl",k)]]<-d[[paste0("clod",k)]]-d$clobas

## ---- FIGURE 1: severity and chloride rise (descriptive) ----
# Day-2 paired cohort, by APACHE II quartile: mean rise by arm, with the
# saline-group median cumulative chloride load annotated for each quartile.
fdf<-m2[!is.na(m2$apache)&!is.na(m2$dCl2)&!is.na(m2$load2),]
fdf$q<-cut(fdf$apache,quantile(fdf$apache,c(0,.25,.5,.75,1)),include.lowest=TRUE,labels=FALSE)
# cut() is right-closed, so display the realized non-overlapping integer ranges.
fqlab<-c("Q1\n(0-9)","Q2\n(10-12)","Q3\n(13-17)","Q4\n(18-45)")
for(.k in 1:2){ if(.k==1) png(file.path(FIG_OUT,"Figure1.png"),width=950,height=700,res=160) else pdf(file.path(FIG_OUT,"Figure1.pdf"),width=950/160,height=700/160); par(mar=c(6,5,4,1),mgp=c(3,0.8,0))
rstat<-function(a,q){s<-fdf[fdf$arm==a&fdf$q==q,];m<-mean(s$dCl2);se<-sd(s$dCl2)/sqrt(nrow(s));c(m,m-1.96*se,m+1.96*se,median(s$load2))}
plot(NA,xlim=c(0.5,4.5),ylim=c(-1.6,3.6),xaxt="n",xlab="",ylab="Day-2 rise in serum chloride (mmol/L)",bty="l",las=1)
abline(h=0,col="grey70",lty=2); axis(1,at=1:4,labels=fqlab,tick=FALSE,line=0.6,cex.axis=0.85)
for(a in c("Saline","Plasma-Lyte")){cl<-if(a=="Saline")col_s else col_p; M<-sapply(1:4,function(q)rstat(a,q))
  off<-if(a=="Saline")-0.06 else 0.06; xx<-(1:4)+off
  segments(xx,M[2,],xx,M[3,],col=cl,lwd=2.2); points(xx,M[1,],pch=19,col=cl,cex=1.5)}
Ms<-sapply(1:4,function(q)rstat("Saline",q)); text(1:4,rep(-1.35,4),sprintf("saline study-fluid load\n%.0f mmol",Ms[4,]),cex=0.62,col="grey30")
legend("topleft",pch=19,col=c(col_s,col_p),legend=c("Saline","Plasma-Lyte"),bty="n",cex=0.9)
title("Chloride rise by severity and randomized arm",adj=0,font.main=2,cex.main=1.0); mtext("APACHE II quartile",1,line=4.2,cex=0.9)
dev.off() }

## ---- FIGURE 2: longitudinal randomized-arm contrasts ----
fit<-readRDS(P$brms_fit)
cells<-expand.grid(dayf=factor(1:3), eGFRc=c(-4,-2,0))
mk<-function(a) transform(cells, arm=factor(a,levels=c("Plasma-Lyte","Saline")), clobas=0, day=as.numeric(as.character(dayf)))
dif<-posterior_epred(fit,newdata=mk("Saline"),re_formula=NA)-posterior_epred(fit,newdata=mk("Plasma-Lyte"),re_formula=NA)
cells$est<-apply(dif,2,median); cells$lo<-apply(dif,2,quantile,.025); cells$hi<-apply(dif,2,quantile,.975)
cells$eGFR<-cells$eGFRc*15+90; cells$day<-as.numeric(as.character(cells$dayf))
for(.k in 1:2){ if(.k==1) png(file.path(FIG_OUT,"Figure2.png"),width=1500,height=620,res=150) else pdf(file.path(FIG_OUT,"Figure2.pdf"),width=1500/150,height=620/150); par(mfrow=c(1,2),mar=c(4.5,4.6,3,1.2),mgp=c(2.6,0.8,0))
# R2 (reviewer 1, comment 21): start both panels at baseline (day 0 = 0 by definition).
daylab<-c("Baseline","1","2","3")
plot(NA,xlim=c(-0.2,3.2),ylim=c(-2.5,3.5),xaxt="n",xlab="Day after randomization",ylab="Change from baseline (mmol/L)")
title("A   Serum chloride change by arm",adj=0,cex.main=1.0,font.main=2); axis(1,0:3,daylab); abline(h=0,col="grey75",lty=2)
for(a in c("Saline","Plasma-Lyte")){col<-ifelse(a=="Saline",col_s,col_p); M<-cbind(c(m=0,lo=0,hi=0),sapply(1:3,function(k) mse(d[[paste0("dCl",k)]][d$arm==a])))
  polygon(c(0:3,3:0),c(M["lo",],rev(M["hi",])),col=adjustcolor(col,0.15),border=NA); lines(0:3,M["m",],col=col,lwd=3,type="b",pch=19,cex=1.2)}
legend("topleft",c("0.9% saline","Plasma-Lyte"),col=c(col_s,col_p),lwd=3,pch=19,bty="n",cex=0.95)
eg<-c("30"="#922B21","60"="#B7950B","90"="#1F618D")
plot(NA,xlim=c(-0.2,3.65),ylim=c(0,4.2),xaxt="n",xlab="Day after randomization",ylab="Saline minus Plasma-Lyte (mmol/L)")
title("B   Randomized arm difference by eGFR",adj=0,cex.main=1.0,font.main=2); axis(1,0:3,daylab); abline(h=0,col="grey80",lty=3)
for(e in c(30,60,90)){ s<-cells[cells$eGFR==e,]; s<-s[order(s$day),]; s<-rbind(transform(s[1,],day=0,est=0,lo=0,hi=0),s); col<-eg[as.character(e)]
  polygon(c(s$day,rev(s$day)),c(s$lo,rev(s$hi)),col=adjustcolor(col,0.13),border=NA); lines(s$day,s$est,col=col,lwd=3,type="b",pch=19)
  text(3.12,s$est[s$day==3],paste("eGFR",e),pos=4,col=col,cex=0.78,font=2,xpd=NA) }
dev.off() }

## ---- FIGURE 3: randomized normalization and within-saline association ----
paired<-m2[complete.cases(m2[,c("dCl2","load2","study_vol2","nonstudy_vol2","arm")]),]
dd<-mean(paired$dCl2[paired$arm=="Saline"])-mean(paired$dCl2[paired$arm=="Plasma-Lyte"])
ld<-mean(paired$load2[paired$arm=="Saline"])-mean(paired$load2[paired$arm=="Plasma-Lyte"])
conc<-as.numeric(100*dd/ld)
set.seed(20240607)
conc_boot<-numeric(10000); ii_s<-which(paired$arm=="Saline"); ii_p<-which(paired$arm=="Plasma-Lyte")
for(b in seq_len(10000)){
  # Sample Plasma-Lyte then saline to preserve the documented stratified-bootstrap RNG order.
  ip<-sample(ii_p,length(ii_p),replace=TRUE); is<-sample(ii_s,length(ii_s),replace=TRUE)
  ddb<-mean(paired$dCl2[is])-mean(paired$dCl2[ip])
  ldb<-mean(paired$load2[is])-mean(paired$load2[ip])
  conc_boot[b]<-100*ddb/ldb
}
conc_ci<-as.numeric(quantile(conc_boot,c(.025,.975),na.rm=TRUE))

cc<-m2[complete.cases(m2[,c("dCl2","load2","clobas","creabas","peso","sofa","pambas","age")]),]
sal<-cc[cc$arm=="Saline",]
gps<-lm(load2~creabas+peso+sofa+pambas+age+clobas,sal)
sw<-dnorm(sal$load2,mean(sal$load2),sd(sal$load2))/dnorm(sal$load2,fitted(gps),sigma(gps))
sw<-pmin(sw,quantile(sw,.99))
fit_voln<-lm(dCl2~load2,sal)
fit_voli<-lm(dCl2~load2,sal,weights=sw)
voln<-as.numeric(100*coef(fit_voln)["load2"]); voli<-as.numeric(100*coef(fit_voli)["load2"])
set.seed(20240608)
voln_boot<-rep(NA_real_,10000); voli_boot<-rep(NA_real_,10000)
for(b in seq_len(10000)){
  ib<-sample(seq_len(nrow(sal)),nrow(sal),replace=TRUE); sb<-sal[ib,]
  fit_n_b<-try(lm(dCl2~load2,sb),silent=TRUE)
  gps_b<-try(lm(load2~creabas+peso+sofa+pambas+age+clobas,sb),silent=TRUE)
  if(inherits(fit_n_b,"try-error")||inherits(gps_b,"try-error")||!is.finite(sigma(gps_b))||sigma(gps_b)<=0) next
  wb<-dnorm(sb$load2,mean(sb$load2),sd(sb$load2))/dnorm(sb$load2,fitted(gps_b),sigma(gps_b))
  wb<-pmin(wb,quantile(wb,.99,na.rm=TRUE))
  fit_i_b<-try(lm(dCl2~load2,sb,weights=wb),silent=TRUE)
  if(inherits(fit_i_b,"try-error")) next
  voln_boot[b]<-100*coef(fit_n_b)["load2"]
  voli_boot[b]<-100*coef(fit_i_b)["load2"]
}
if(sum(is.finite(voln_boot))<9900||sum(is.finite(voli_boot))<9900) stop("Too many failed Analysis 2 bootstrap replicates.")
voln_ci<-as.numeric(quantile(voln_boot,c(.025,.975),na.rm=TRUE))
voli_ci<-as.numeric(quantile(voli_boot,c(.025,.975),na.rm=TRUE))
volume_results<-data.frame(
  analysis=c("Unweighted saline-only volume association","IPW-adjusted saline-only volume association"),
  outcome="Day-2 serum chloride change from baseline (mmol/L)",
  exposure_scale="Per additional 100 mmol post-randomization study-fluid chloride delivered as saline through day 2",
  outcome_model_n=nrow(sal),
  propensity_model_n=nrow(sal),
  estimate=c(voln,voli),
  conf_low=c(voln_ci[1],voli_ci[1]),
  conf_high=c(voln_ci[2],voli_ci[2]),
  ci_method="Participant bootstrap with propensity-score refitting; 10000 replicates; seed 20240608",
  stringsAsFactors=FALSE
)
write.csv(volume_results,file.path(P$reproduction,"analysis2_volume_results.csv"),row.names=FALSE)
concentration_results<-data.frame(
  analysis="Descriptive normalization of randomized between-group differences",
  outcome="Day-2 serum chloride change from baseline (mmol/L)",
  exposure_scale="Per 100 mmol between-group excess post-randomization study-fluid chloride associated with randomized fluid assignment",
  analysis_n=nrow(paired), saline_n=sum(paired$arm=="Saline"), plasma_lyte_n=sum(paired$arm=="Plasma-Lyte"),
  chloride_change_difference=dd, chloride_load_difference=ld,
  estimate=conc, conf_low=conc_ci[1], conf_high=conc_ci[2],
  ci_method="Stratified nonparametric participant bootstrap; 10000 replicates; seed 20240607",
  stringsAsFactors=FALSE
)
write.csv(concentration_results,file.path(P$reproduction,"analysis2_concentration_results.csv"),row.names=FALSE)
weight_diagnostics<-data.frame(
  analysis="IPW-adjusted saline-only volume association",
  n=nrow(sal),
  weight_mean=mean(sw), weight_sd=sd(sw), weight_min=min(sw),
  weight_median=median(sw), weight_p99=quantile(sw,.99), weight_max=max(sw),
  effective_sample_size=sum(sw)^2/sum(sw^2),
  stringsAsFactors=FALSE
)
write.csv(weight_diagnostics,file.path(P$reproduction,"analysis2_weight_diagnostics.csv"),row.names=FALSE)

mean_study_s<-mean(paired$study_vol2[paired$arm=="Saline"])
mean_study_p<-mean(paired$study_vol2[paired$arm=="Plasma-Lyte"])
mean_nonstudy_s<-mean(paired$nonstudy_vol2[paired$arm=="Saline"])
mean_nonstudy_p<-mean(paired$nonstudy_vol2[paired$arm=="Plasma-Lyte"])
conc_component<-(154-98)*(mean_study_s+mean_study_p)/2/1000
volume_component<-(mean_study_s-mean_study_p)*(154+98)/2/1000
apparent_volume_l<-ld/dd
apparent_volume_ci<-sort(100/conc_ci)
descriptive_checks<-data.frame(
  measure=c("Mean study-fluid volume, saline, days 1-2","Mean study-fluid volume, Plasma-Lyte, days 1-2","Between-group study-fluid volume difference","Mechanical concentration component of study-fluid load gap","Mechanical volume component of study-fluid load gap","Mechanical concentration percentage of study-fluid load gap","Mechanical volume percentage of study-fluid load gap","Study-fluid load gap","Median non-study-fluid volume, saline, days 1-2","Median non-study-fluid volume, Plasma-Lyte, days 1-2","Between-group mean non-study-fluid volume difference","Load gap after pricing mean non-study-volume difference at 98 mmol/L","Load gap after pricing mean non-study-volume difference at 154 mmol/L","Apparent distribution volume","Apparent distribution volume, lower 95% limit","Apparent distribution volume, upper 95% limit","Observed/ECF-predicted serum chloride rise for 70 kg","Observed/TBW-predicted serum chloride rise for 70 kg"),
  value=c(mean_study_s,mean_study_p,mean_study_s-mean_study_p,conc_component,volume_component,100*conc_component/ld,100*volume_component/ld,ld,median(paired$nonstudy_vol2[paired$arm=="Saline"]),median(paired$nonstudy_vol2[paired$arm=="Plasma-Lyte"]),mean_nonstudy_s-mean_nonstudy_p,ld+(mean_nonstudy_s-mean_nonstudy_p)*98/1000,ld+(mean_nonstudy_s-mean_nonstudy_p)*154/1000,apparent_volume_l,apparent_volume_ci[1],apparent_volume_ci[2],100*14/apparent_volume_l,100*42/apparent_volume_l),
  unit=c("mL","mL","mL","mmol","mmol","percent","percent","mmol","mL","mL","mL","mmol","mmol","L","L","L","percent","percent"),
  stringsAsFactors=FALSE
)
write.csv(descriptive_checks,file.path(P$reproduction,"analysis2_descriptive_checks.csv"),row.names=FALSE)
cat(sprintf("Analysis 2 saline-only volume association (n=%d): unweighted %.5f (95%% CI %.5f to %.5f); IPW-adjusted %.5f (95%% CI %.5f to %.5f) mmol/L per 100 mmol.\n",
            nrow(sal),voln,voln_ci[1],voln_ci[2],voli,voli_ci[1],voli_ci[2]))
cat(sprintf("Analysis 2 descriptive randomized normalization (n=%d): ΔCl=%+.5f mmol/L; Δstudy-load=%.5f mmol; normalized ratio %.5f (bootstrap 95%% CI %.5f to %.5f) mmol/L per 100 mmol between-group excess study-fluid chloride.\n",
            nrow(paired),dd,ld,conc,conc_ci[1],conc_ci[2]))
for(.k in 1:2){ if(.k==1) png(file.path(FIG_OUT,"Figure3.png"),width=850,height=650,res=150) else pdf(file.path(FIG_OUT,"Figure3.pdf"),width=850/150,height=650/150); par(mar=c(6.4,4.9,3,1),mgp=c(2.9,0.8,0))
bb<-c(conc,voln,voli); blo<-c(conc_ci[1],voln_ci[1],voli_ci[1]); bhi<-c(conc_ci[2],voln_ci[2],voli_ci[2]); names(bb)<-c("Randomized\nbetween-group ratio","Within saline\n(unweighted)","Within saline\n(IPW-adjusted)")
bp<-barplot(bb,col=c(col_s,"#AEB6BF","#5D6D7E"),ylab="Normalized estimate\n(mmol/L per 100 mmol)",ylim=c(0,2.15),las=1,cex.names=0.75); title("Randomized normalization and within-saline association",adj=0,cex.main=0.98,font.main=2)
segments(bp,blo,bp,bhi,lwd=1.6,col="grey20"); segments(bp-0.06,blo,bp+0.06,blo,lwd=1.6,col="grey20"); segments(bp-0.06,bhi,bp+0.06,bhi,lwd=1.6,col="grey20")
text(bp,bhi+0.10,sprintf("%.2f",bb),cex=0.92)
mtext("Randomized bar: descriptive ratio of between-group differences",side=1,line=4.35,cex=0.66,col="grey30")
mtext("Within-saline bars: observational study-fluid volume associations",side=1,line=5.15,cex=0.66,col="grey30"); dev.off() }

## ---- FIGURE 4: Bayesian variance partition (level, all days); fit from basics_chloride_varpart.R ----
long<-build_varpart_long(d)
long<-long[complete.cases(long),]; long<-scale_varpart_predictors(long)
post<-brms::as_draws_df(readRDS(P$varpart_fit))
Xf<-as.matrix(long[,c("load","vol")]); Xp<-as.matrix(long[,c("clobas","age","creabas","sofa")]); Xt<-model.matrix(~day,long)[,-1,drop=FALSE]
bf<-as.matrix(post[,c("b_load","b_vol")]); bp_<-as.matrix(post[,c("b_clobas","b_age","b_creabas","b_sofa")]); bt<-as.matrix(post[,c("b_day2","b_day3")])
Vf<-rowSums((bf%*%cov(Xf))*bf); Vp<-rowSums((bp_%*%cov(Xp))*bp_); Vt<-rowSums((bt%*%cov(Xt))*bt); Xfull<-cbind(Xf,Xp,Xt); bfull<-cbind(bf,bp_,bt); Vfull<-rowSums((bfull%*%cov(Xfull))*bfull); Vcov<-Vfull-Vf-Vp-Vt; Vbp<-post$sd_id__Intercept^2; Vr<-post$sigma^2; Vtot<-Vfull+Vbp+Vr
S<-cbind(Fluid=Vf/Vtot,StudyDay=Vt/Vtot,Covariance=Vcov/Vtot,Patient=Vp/Vtot,Between=Vbp/Vtot,Residual=Vr/Vtot); med<-100*apply(S,2,median); lo<-100*apply(S,2,quantile,.025); hi<-100*apply(S,2,quantile,.975)
for(.k in 1:2){ if(.k==1) png(file.path(FIG_OUT,"Figure4.png"),width=950,height=620,res=150) else pdf(file.path(FIG_OUT,"Figure4.pdf"),width=950/150,height=620/150); par(mar=c(5,9,3,2),mgp=c(2.6,0.8,0))
o<-c("Fluid","StudyDay","Covariance","Residual","Between","Patient"); labs<-c(Fluid="Fluid\n(study load + total volume)",StudyDay="Study day",Covariance="Shared predictor\nvariation (covariance)",Residual="Day-to-day\nwithin-patient",Between="Persistent\nbetween-patient",Patient="Measured baseline\n(mainly baseline Cl)")
bp2<-barplot(med[o],horiz=TRUE,las=1,names.arg=labs[o],col=c(col_acc,"#5499C7","#D4AC0D","grey75","grey55","grey45"),xlab="% of modeled total variance in serum chloride",xlim=c(-2,60),cex.names=0.73); abline(v=0,col="grey45"); title("Variation in serum chloride, by model component",adj=0,cex.main=0.98,font.main=2)
segments(lo[o],bp2,hi[o],bp2,lwd=1.6,col="grey25"); segments(lo[o],bp2-0.15,lo[o],bp2+0.15,lwd=1.6,col="grey25"); segments(hi[o],bp2-0.15,hi[o],bp2+0.15,lwd=1.6,col="grey25")
text(pmax(hi[o]+1.5,med[o]+1.5),bp2,sprintf("%.1f%%",med[o]),cex=0.8,pos=4,offset=0.1,xpd=NA)
mtext("Bayesian longitudinal partition (days 1-3); bars = posterior median, lines = 95% CrI",side=1,line=3.4,cex=0.74,col="grey30"); dev.off() }

cat("Figures 1-4 written to", FIG_OUT, "\n")
