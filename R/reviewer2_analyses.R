# reviewer2_analyses.R — analyses added for the second review round (R2).
#   1. Descriptive 95% CIs for Results (reviewer 1, comment 13) and per-day
#      paired counts and load gaps supporting the day-2 rationale (comments 3, 11, 20).
#   2. Split of the measured-baseline block of the primary variance partition
#      (comments 5, 14).
#   3. Exploratory correlates of the patient-specific intercepts (comment 17).
#   4. Sensitivity refit with randomized arm added to the fluid block.
# Writes aggregate-only CSVs to reproduction/; the refit is saved to model_fits/
# (git-ignored). Run from revision_2/: Rscript R/reviewer2_analyses.R
suppressMessages(library(brms))
if (!exists("BASICS_PATHS", inherits = TRUE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grepl("^--file=", .args)][[1]])
  source(file.path(dirname(normalizePath(.file)), "project_paths.R"))
}
P <- get("BASICS_PATHS", inherits = TRUE)
source(file.path(P$root, "R", "exposure_construction.R"))
out <- function(x, name) write.csv(x, file.path(P$reproduction, name), row.names = FALSE)
q3 <- function(x) c(estimate = median(x), low = unname(quantile(x, .025)), high = unname(quantile(x, .975)))

## ---- 1. Descriptive intervals and day-2 rationale ----
C <- readRDS(P$cache_file); d <- C$d; m2 <- C$m2
for (k in 1:3) d[[paste0("dCl", k)]] <- d[[paste0("clod", k)]] - d$clobas
d <- add_postbaseline_exposures(d)
desc <- do.call(rbind, lapply(1:3, function(k) {
  x <- d[[paste0("dCl", k)]]; s <- x[d$arm == "Saline" & !is.na(x)]; p <- x[d$arm == "Plasma-Lyte" & !is.na(x)]
  tt <- t.test(s, p)
  has <- !is.na(d$clobas) & !is.na(d[[paste0("clod", k)]]); ld <- d[[paste0("study_load", k)]]
  lt <- t.test(ld[has & d$arm == "Saline"], ld[has & d$arm == "Plasma-Lyte"])
  data.frame(day = k, n_saline = length(s), n_plasmalyte = length(p), paired_n = sum(has),
             diff_change = mean(s) - mean(p), diff_low = tt$conf.int[1], diff_high = tt$conf.int[2],
             load_gap = unname(diff(rev(lt$estimate))), load_gap_low = lt$conf.int[1], load_gap_high = lt$conf.int[2])
}))
out(desc, "reviewer2_daily_differences.csv")

fdf <- m2[!is.na(m2$apache) & !is.na(m2$dCl2) & !is.na(m2$load2), ]
fdf$q <- cut(fdf$apache, quantile(fdf$apache, c(0, .25, .5, .75, 1)), include.lowest = TRUE, labels = FALSE)
quart <- do.call(rbind, lapply(c("Saline", "Plasma-Lyte"), function(a) do.call(rbind, lapply(1:4, function(q) {
  x <- fdf$dCl2[fdf$arm == a & fdf$q == q]; se <- sd(x) / sqrt(length(x))
  data.frame(arm = a, apache_quartile = q, n = length(x), mean_change = mean(x),
             low = mean(x) - 1.96 * se, high = mean(x) + 1.96 * se)
}))))
out(quart, "reviewer2_apache_quartile_changes.csv")

## ---- 2. Split of the measured-baseline block ----
fit <- readRDS(P$varpart_fit); long <- fit$data; post <- as_draws_df(fit)
V <- function(X, b) rowSums((b %*% cov(X)) * b)
Vb <- function(cols) V(as.matrix(long[, cols, drop = FALSE]), as.matrix(post[, paste0("b_", cols)]))
Xt <- model.matrix(~day, long)[, -1, drop = FALSE]; bt <- as.matrix(post[, c("b_day2", "b_day3")])
allc <- c("load", "vol", "clobas", "age", "creabas", "sofa")
Vfull <- V(cbind(as.matrix(long[, allc]), Xt), cbind(as.matrix(post[, paste0("b_", allc)]), bt))
Vtot <- Vfull + post$sd_id__Intercept^2 + post$sigma^2
Vp <- Vb(c("clobas", "age", "creabas", "sofa")); Vcl <- Vb("clobas"); Voth <- Vb(c("age", "creabas", "sofa"))
split <- rbind(
  c(component = "measured_baseline_block", 100 * q3(Vp / Vtot)),
  c(component = "baseline_chloride_alone", 100 * q3(Vcl / Vtot)),
  c(component = "age_creatinine_sofa", 100 * q3(Voth / Vtot)),
  c(component = "within_block_covariance", 100 * q3((Vp - Vcl - Voth) / Vtot)),
  c(component = "baseline_chloride_share_of_block", 100 * q3(Vcl / Vp)),
  c(component = "between_patient_sd_mmol", q3(post$sd_id__Intercept)),
  c(component = "residual_sd_mmol", q3(post$sigma))
)
out(as.data.frame(split), "reviewer2_baseline_block_split.csv")

## ---- 3. Correlates of the patient-specific intercepts (exploratory) ----
load(P$data_file)
dd <- dados; rm(dados)
dd$arm <- ifelse(dd$grupo == "Saline Solution", "Saline", "Plasma-Lyte")
dd$age <- as.numeric(difftime(dd$dtarrol, dd$dtnasc, units = "days")) / 365.25
dd <- add_postbaseline_exposures(dd)
lg <- build_varpart_long(dd); lg <- lg[complete.cases(lg), ]
stopifnot(nrow(lg) == nrow(long), isTRUE(all.equal(lg$clev, long$clev)),
          identical(as.integer(lg$id), as.integer(as.character(long$id))))
re <- ranef(fit)$id[, "Estimate", "Intercept"]; ids <- as.integer(names(re))
pt <- data.frame(u = unname(re), sodium = dd$sodbas[ids], ph = dd$phabas[ids], map = dd$pambas[ids],
                 weight = dd$peso[ids], preenrol_saline_L = zero_if_missing(dd$vosr24h[ids]) / 1000,
                 apache = dd$apache[ids], saline_arm = as.integer(dd$arm[ids] == "Saline"),
                 female = as.integer(dd$sexo[ids] == "Feminino"),
                 vasopressor = as.integer(as.character(dd$vasopressorbas[ids]) == "Sim"),
                 admission = factor(dd$proc_sepse[ids]), site = factor(dd$HospitalId[ids]))
m <- lm(u ~ scale(sodium) + scale(ph) + scale(map) + scale(weight) + scale(preenrol_saline_L) +
          scale(apache) + saline_arm + female + vasopressor + admission, data = pt)
ms <- update(m, . ~ . + site)
cf <- summary(m)$coefficients; ci <- confint(m)
assoc <- data.frame(term = rownames(cf), estimate = cf[, 1], low = ci[, 1], high = ci[, 2], p = cf[, 4])
assoc$term <- sub("admissionNão planejada, sem sepse", "admission: unplanned, non-sepsis (vs elective)", assoc$term)
assoc$term <- sub("admissionNão planejada, com sepse", "admission: unplanned, sepsis (vs elective)", assoc$term)
assoc <- rbind(assoc, data.frame(term = c("R2_baseline_variables", "R2_with_site", "R2_site_alone", "n_patients_complete", "n_patients_total", "n_sites"),
  estimate = c(summary(m)$r.squared, summary(ms)$r.squared, summary(lm(u ~ site, pt))$r.squared, nobs(m), nrow(pt), nlevels(pt$site)),
  low = NA, high = NA, p = NA))
out(assoc, "reviewer2_intercept_correlates.csv")

## ---- 4. Sensitivity: randomized arm in the fluid block ----
long$saline <- as.integer(dd$arm[as.integer(as.character(long$id))] == "Saline")
arm_path <- file.path(P$root, "model_fits", "reviewer2_arm_varpart_fit.rds")
fa <- brm(clev ~ saline + load + vol + day + clobas + age + creabas + sofa + (1 | id), data = long,
          chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 1, refresh = 0,
          prior = prior(normal(0, 5), class = "b"))
saveRDS(fa, arm_path)
pa <- as_draws_df(fa)
fl <- c("saline", "load", "vol"); pb <- c("clobas", "age", "creabas", "sofa")
Vf <- V(as.matrix(long[, fl]), as.matrix(pa[, paste0("b_", fl)]))
Vpp <- V(as.matrix(long[, pb]), as.matrix(pa[, paste0("b_", pb)]))
Vt <- V(Xt, as.matrix(pa[, c("b_day2", "b_day3")]))
Vf2 <- V(cbind(as.matrix(long[, c(fl, pb)]), Xt), cbind(as.matrix(pa[, paste0("b_", c(fl, pb))]), as.matrix(pa[, c("b_day2", "b_day3")])))
Vt2 <- Vf2 + pa$sd_id__Intercept^2 + pa$sigma^2
dg <- posterior::summarise_draws(posterior::as_draws_array(fa, variable = c(paste0("b_", c(fl, pb)), "b_day2", "b_day3", "sd_id__Intercept", "sigma")), "rhat", "ess_bulk", "ess_tail")
nuts <- nuts_params(fa)
arm <- rbind(
  data.frame(component = c("fluid_incl_arm", "study_day", "measured_baseline", "covariance", "between_patient", "residual"),
             rbind(100 * q3(Vf / Vt2), 100 * q3(Vt / Vt2), 100 * q3(Vpp / Vt2), 100 * q3((Vf2 - Vf - Vpp - Vt) / Vt2),
                   100 * q3(pa$sd_id__Intercept^2 / Vt2), 100 * q3(pa$sigma^2 / Vt2))),
  data.frame(component = "b_saline_mmol", t(q3(pa$b_saline))),
  data.frame(component = c("max_rhat", "min_ess_bulk", "min_ess_tail", "divergences"),
             estimate = c(max(dg$rhat), min(dg$ess_bulk), min(dg$ess_tail), sum(nuts$Parameter == "divergent__" & nuts$Value == 1)),
             low = NA, high = NA))
out(arm, "reviewer2_arm_sensitivity.csv")
cat("Reviewer-2 analyses written to", P$reproduction, "\n")
