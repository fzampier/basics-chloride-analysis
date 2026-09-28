# ============================================================
# BaSICS chloride — Bayesian longitudinal effect-modification model: PRIMARY
# serum Cl ~ arm*day*eGFR + baseline Cl + (1+day|patient)
# Partial pooling for repeated serum chloride measurements on days 1-3.
# Run basics_chloride_analysis.R first (writes .chl_cache.rds).
# ============================================================
suppressMessages(suppressWarnings(library(brms)))
if (!exists("BASICS_PATHS", inherits = TRUE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grepl("^--file=", .args)][[1]])
  source(file.path(dirname(normalizePath(.file)), "project_paths.R"))
}
P <- get("BASICS_PATHS", inherits = TRUE)
C <- readRDS(P$cache_file); d <- C$d; d$id <- 1:nrow(d)
egfr <- function(scr,age,fem){k<-ifelse(fem==1,0.7,0.9);a<-ifelse(fem==1,-0.241,-0.302)
  142*pmin(scr/k,1)^a*pmax(scr/k,1)^(-1.2)*0.9938^age*ifelse(fem==1,1.012,1)}
d$eGFR <- egfr(d$creabas,d$age,d$female); d$eGFR[!is.finite(d$eGFR)|d$eGFR<5|d$eGFR>200] <- NA
long <- do.call(rbind, lapply(1:3, function(k){ v<-d[[paste0("clod",k)]]; ok<-!is.na(v)&!is.na(d$clobas)&!is.na(d$eGFR)
  data.frame(id=d$id[ok], cl=v[ok], day=k, arm=factor(d$arm[ok],levels=c("Plasma-Lyte","Saline")),
             clobas=d$clobas[ok]-106, eGFR=d$eGFR[ok]) }))   # clobas centred at 106
long$dayf <- factor(long$day); long$eGFRc <- (long$eGFR-90)/15

pr <- set_prior("normal(0,10)", class="b")
fit <- brm(cl ~ arm*dayf*eGFRc + clobas + (1+day|id), data=long, family=gaussian,
           prior=pr, chains=4, cores=4, iter=4000, seed=1, refresh=0,
           control=list(adapt_delta=0.9))
saveRDS(fit, P$brms_fit)   # gitignored (*.rds)

sink(file.path(P$reproduction,"chl_brms_results.txt"))
cat("=== brms longitudinal model (n=",nrow(long)," obs, ",length(unique(long$id))," patients) ===\n",sep="")
fe <- fixef(fit)
cat("\n--- treatment-coded 3-way terms (arm x day x eGFRc, per +15 mL/min) ---\n")
print(round(fe[grep("armSaline:dayf.*:eGFRc", rownames(fe)), c("Estimate","Q2.5","Q97.5")],3))
cat("\n--- Arm effect on serum Cl (Saline - Plasma-Lyte) by day x eGFR [95% CrI] ---\n")
cells <- expand.grid(dayf=factor(1:3), eGFRc=c(-4,-2,0))
mk <- function(a) transform(cells, arm=factor(a,levels=c("Plasma-Lyte","Saline")), clobas=0, day=as.numeric(as.character(dayf)))
dif <- posterior_epred(fit, newdata=mk("Saline"), re_formula=NA) - posterior_epred(fit, newdata=mk("Plasma-Lyte"), re_formula=NA)
cells$est<-apply(dif,2,median); cells$lo<-apply(dif,2,quantile,.025); cells$hi<-apply(dif,2,quantile,.975); cells$eGFR<-cells$eGFRc*15+90
o <- cells[order(cells$eGFR,cells$dayf),]
for (i in 1:nrow(o)) cat(sprintf("  eGFR %2.0f day %s: %+.2f [%+.2f, %+.2f]\n", o$eGFR[i], as.character(o$dayf[i]), o$est[i], o$lo[i], o$hi[i]))
cat("\n--- Directional posterior probability of arm-by-eGFR effect modification ---\n")
mkone <- function(a, d, e) data.frame(
  dayf=factor(d, levels=1:3), eGFRc=e,
  arm=factor(a, levels=c("Plasma-Lyte", "Saline")),
  clobas=0, day=d
)
for (d in 2:3) {
  contrast0 <- posterior_epred(fit, newdata=mkone("Saline", d, 0), re_formula=NA) -
    posterior_epred(fit, newdata=mkone("Plasma-Lyte", d, 0), re_formula=NA)
  contrast1 <- posterior_epred(fit, newdata=mkone("Saline", d, 1), re_formula=NA) -
    posterior_epred(fit, newdata=mkone("Plasma-Lyte", d, 1), re_formula=NA)
  interaction <- as.numeric(contrast1 - contrast0)
  cat(sprintf(
    "  day %d: Pr(arm difference decreases per +15 mL/min/1.73 m2 eGFR) = %.4f%%\n",
    d, 100 * mean(interaction < 0)
  ))
}
cat("\n--- diagnostics ---\n")
cat("max Rhat:", round(max(rhat(fit), na.rm=TRUE),3), "| min Bulk_ESS:", round(min(summary(fit)$fixed[,"Bulk_ESS"]),0),
    "| divergences:", sum(subset(nuts_params(fit), Parameter=="divergent__")$Value), "\n")
sink()
cat("DONE\n")
