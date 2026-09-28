# =====================================================================
# BaSICS — Chloride load, serum chloride, and outcomes
# Master data derivation and aggregate checks.
# Run from revision/ or set BASICS_REVISION_ROOT.
# =====================================================================
if (!exists("BASICS_PATHS", inherits = TRUE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grepl("^--file=", .args)][[1]])
  source(file.path(dirname(normalizePath(.file)), "project_paths.R"))
}
P <- get("BASICS_PATHS", inherits = TRUE)
source(file.path(P$root, "R", "exposure_construction.R"))
load(P$data_file)
z <- zero_if_missing
col_s <- "#C0392B"; col_p <- "#2471A3"; col_acc <- "#884EA0"

# ---- derived variables ----
d <- dados
d$arm  <- ifelse(d$grupo=="Saline Solution","Saline","Plasma-Lyte")
d$armN <- as.integer(d$arm=="Saline")
d$age  <- as.numeric(difftime(d$dtarrol,d$dtnasc,units="days"))/365.25
d$female <- as.integer(d$sexo=="Feminino")
d$admtype <- factor(ifelse(d$proc_sepse=="Planejada","Elective",
                    ifelse(d$proc_sepse=="Não planejada, sem sepse","Non-sepsis","Sepsis")),
                    levels=c("Elective","Non-sepsis","Sepsis"))
d$death <- as.integer(d$ObitoPrimary_imput=="Óbito")
d$aki   <- as.integer(d$kdigo_obito3=="Sim")
d$rrt   <- as.integer(d$isubrhos=="Sim")
d <- add_postbaseline_exposures(d)
d$load1 <- d$study_load1
d$load2 <- d$study_load2
d$load3 <- d$study_load3
d$studyfl <- d$study_vol3
d$totalfl <- d$total_vol3
d$vol2 <- d$total_vol2
for (dd in 1:3) d[[paste0("dCl",dd)]] <- d[[paste0("clod",dd)]] - d$clobas

analytic <- !is.na(d$clobas) & (!is.na(d$clod1) | !is.na(d$clod2) | !is.na(d$clod3))
meas <- d[analytic,]                          # submitted analytic cohort
m2   <- d[!is.na(d$clobas)&!is.na(d$clod2),]  # day-2 paired
cat("Analytic cohort n =", nrow(meas), "| day-2 paired n =", nrow(m2), "\n")

# Construct-level validation: document why legacy Carga_cloro fields are not
# valid substitutes for the active post-randomization study-fluid exposure.
pre_enrollment_load <- 154*z(d$vosr24h)/1000 + 98*z(d$vosb24h)/1000
valid_legacy_day2 <- is.na(d$voflud2_total) | d$voflud2_total < 60000
exposure_validation <- data.frame(
  check=c(
    "Legacy day-0 load correlation with pre-enrollment 24-hour load",
    "Legacy day-1 residual correlation with non-study-fluid volume",
    "Legacy day-2 residual correlation with non-study-fluid volume",
    "Documented raw day-2 study-fluid outlier",
    "Corrected day-2 study-fluid value",
    "Maximum corrected daily study-fluid volume",
    "Maximum corrected daily total-fluid volume"
  ),
  value=c(
    cor(d$Carga_cloro_d0,pre_enrollment_load,use="complete.obs"),
    cor(d$Carga_cloro_d1-d$study_load_day1,d$nonstudy_volume_day1,use="complete.obs"),
    cor((d$Carga_cloro_d2-d$study_load_day2)[valid_legacy_day2],
        d$nonstudy_volume_day2[valid_legacy_day2],use="complete.obs"),
    max(d$voflud2_total,na.rm=TRUE),
    unique(d$study_volume_day2[which(d$voflud2_total==65000)]),
    max(unlist(d[paste0("study_volume_day",1:3)]),na.rm=TRUE),
    max(unlist(d[paste0("total_volume_day",1:3)]),na.rm=TRUE)
  ),
  stringsAsFactors=FALSE
)
write.csv(exposure_validation,file.path(P$reproduction,"exposure_validation.csv"),row.names=FALSE)
cat(sprintf("Exposure validation: legacy day-0/pre-enrollment correlation %.4f; day-1/day-2 legacy residual correlations with non-study volume %.4f/%.4f; raw/corrected day-2 maximum %.0f/%.0f mL.\n",
            exposure_validation$value[1],exposure_validation$value[2],exposure_validation$value[3],exposure_validation$value[4],exposure_validation$value[5]))

miq <- function(x){x<-x[!is.na(x)];sprintf("%.0f [%.0f-%.0f]",median(x),quantile(x,.25),quantile(x,.75))}
miq1<- function(x){x<-x[!is.na(x)];sprintf("%.1f [%.1f-%.1f]",median(x),quantile(x,.25),quantile(x,.75))}
pct <- function(x){sprintf("%d (%.0f%%)",sum(x,na.rm=TRUE),100*mean(x,na.rm=TRUE))}

# =====================================================================
cat("\n################ TABLE 1 — by arm (measured cohort) ################\n")
for (a in c("Saline","Plasma-Lyte")){ s<-meas[meas$arm==a,]
  cat(sprintf("\n--- %s (n=%d) ---\n",a,nrow(s)))
  cat("  Age:",miq(s$age)," | Female:",pct(s$female),"\n")
  cat("  SOFA:",miq1(s$sofa)," | APACHE II:",miq(s$apache),"\n")
  cat("  Admission — Elective/Non-sepsis/Sepsis:",
      paste(table(s$admtype),collapse=" / "),"\n")
  cat("  Baseline Cl:",miq1(s$clobas)," Creat:",miq1(s$creabas)," HCO3:",miq1(s$bsbas)," Na:",miq(s$sodbas),"\n")
  cat(sprintf("  Study fluid d1-3: %s  Total fluid d1-3: %s  Study-fluid chloride d1-3: %s\n",
              miq(s$studyfl), miq(s$totalfl), miq(s$load3)))
}
cat("\n  [selection] analytic vs rest: SOFA",round(mean(meas$sofa,na.rm=TRUE),1),"vs",
    round(mean(d$sofa[!analytic],na.rm=TRUE),1),
    "| %Saline",round(100*mean(meas$arm=="Saline"),1),"vs",round(100*mean(d$arm[!analytic]=="Saline"),1),"\n")

# =====================================================================
cat("\n################ KEY NUMBERS ################\n")
## Descriptive normalization of randomized between-group differences.
for (dd in 1:3){ ok<-!is.na(d$clobas)&!is.na(d[[paste0("clod",dd)]]); dcl<-d[[paste0("dCl",dd)]]; cum<-d[[paste0("load",dd)]]
  dd1<-mean(dcl[ok&d$arm=="Saline"])-mean(dcl[ok&d$arm=="Plasma-Lyte"]); ld<-mean(cum[ok&d$arm=="Saline"])-mean(cum[ok&d$arm=="Plasma-Lyte"])
  cat(sprintf("  Randomized day%d: ΔCl arm-diff=%+.2f, study-load diff=%.0f => descriptive ratio %.2f /100mmol between-group excess study-fluid chloride\n",dd,dd1,ld,100*dd1/ld)) }
ivslope <- { ok<-!is.na(d$clobas)&!is.na(d$clod2); dd1<-mean(d$dCl2[ok&d$arm=="Saline"])-mean(d$dCl2[ok&d$arm=="Plasma-Lyte"]); ld<-mean(d$load2[ok&d$arm=="Saline"])-mean(d$load2[ok&d$arm=="Plasma-Lyte"]); dd1/ld }
cat("  >>> Randomized normalization: ratio of randomized group differences; within-saline estimates: observational volume associations.\n")

dir.create(dirname(P$cache_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(d=d,meas=meas,m2=m2,ivslope=ivslope,col_s=col_s,col_p=col_p,col_acc=col_acc), P$cache_file)
cat("\nData cached. Figures built by companion script.\n")
