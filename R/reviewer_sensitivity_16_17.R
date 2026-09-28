# Reviewer sensitivity analyses for Comments 16 and 17.
#
# Comment 16: refit the longitudinal variance-partition model after excluding
# observations obtained during or after kidney-replacement therapy (KRT).
# Comment 17: repeat the partition using change from baseline serum chloride
# rather than the observed serum chloride level.

suppressMessages(library(brms))

if (!exists("BASICS_PATHS", inherits = TRUE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grepl("^--file=", .args)][[1]])
  source(file.path(dirname(normalizePath(.file)), "project_paths.R"))
}

P <- get("BASICS_PATHS", inherits = TRUE)
source(file.path(P$root, "R", "exposure_construction.R"))
C <- readRDS(P$cache_file)
d <- C$d
long <- build_varpart_long(d, include_krt = TRUE)
long <- long[complete.cases(long[, c(
  "clev", "dcl", "load", "vol", "clobas", "age", "creabas", "sofa"
)]), ]

scale_predictors <- function(dat) {
  scale_varpart_predictors(dat)
}

fit_or_load <- function(dat, outcome, path, seed) {
  refit <- identical(Sys.getenv("REFIT_REVIEWER_SENSITIVITY", unset = "0"), "1")
  if (!refit && file.exists(path)) return(readRDS(path))

  form <- as.formula(paste0(
    outcome,
    " ~ load + vol + day + clobas + age + creabas + sofa + (1 | id)"
  ))
  fit <- brm(
    form,
    data = dat,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = 4,
    seed = seed,
    refresh = 100,
    prior = prior(normal(0, 5), class = "b")
  )
  saveRDS(fit, path)
  fit
}

fit_diagnostics <- function(fit, analysis) {
  diag_vars <- c(
    "b_load", "b_vol", "b_day2", "b_day3", "b_clobas", "b_age", "b_creabas", "b_sofa",
    "sd_id__Intercept", "sigma"
  )
  diag <- posterior::summarise_draws(
    posterior::as_draws_array(fit, variable = diag_vars),
    "rhat", "ess_bulk", "ess_tail"
  )
  nuts <- brms::nuts_params(fit)
  data.frame(
    analysis = analysis,
    max_rhat = max(diag$rhat, na.rm = TRUE),
    min_ess_bulk = min(diag$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(diag$ess_tail, na.rm = TRUE),
    divergences = sum(nuts$Parameter == "divergent__" & nuts$Value == 1),
    stringsAsFactors = FALSE
  )
}

partition_draws <- function(fit, dat) {
  post <- as_draws_df(fit)
  Xf <- as.matrix(dat[, c("load", "vol")])
  Xp <- as.matrix(dat[, c("clobas", "age", "creabas", "sofa")])
  Xt <- model.matrix(~ day, dat)[, -1, drop = FALSE]
  bf <- as.matrix(post[, c("b_load", "b_vol")])
  bp <- as.matrix(post[, c("b_clobas", "b_age", "b_creabas", "b_sofa")])
  bt <- as.matrix(post[, c("b_day2", "b_day3")])
  Vf <- rowSums((bf %*% cov(Xf)) * bf)
  Vp <- rowSums((bp %*% cov(Xp)) * bp)
  Vt <- rowSums((bt %*% cov(Xt)) * bt)
  Xfull <- cbind(Xf, Xp, Xt)
  bfull <- cbind(bf, bp, bt)
  Vfull <- rowSums((bfull %*% cov(Xfull)) * bfull)
  Vcov <- Vfull - Vf - Vp - Vt
  Vbp <- post$sd_id__Intercept^2
  Vr <- post$sigma^2
  Vtot <- Vfull + Vbp + Vr
  data.frame(
    fluid = Vf / Vtot,
    study_day = Vt / Vtot,
    patient = Vp / Vtot,
    covariance = Vcov / Vtot,
    between_patient = Vbp / Vtot,
    residual = Vr / Vtot
  )
}

summarize_partition <- function(draws, analysis, outcome, dat) {
  components <- names(draws)
  do.call(rbind, lapply(components, function(component) {
    x <- 100 * draws[[component]]
    data.frame(
      analysis = analysis,
      outcome = outcome,
      component = component,
      observations = nrow(dat),
      patients = length(unique(dat$id)),
      estimate_percent = median(x),
      conf_low_percent = unname(quantile(x, 0.025)),
      conf_high_percent = unname(quantile(x, 0.975)),
      stringsAsFactors = FALSE
    )
  }))
}

primary_dat <- scale_predictors(long)
primary_fit <- readRDS(P$varpart_fit)
primary_result <- summarize_partition(
  partition_draws(primary_fit, primary_dat),
  "Primary level model",
  "Serum chloride level on days 1-3",
  primary_dat
)

# Restrict to observations with known absence of KRT through the measurement day.
no_krt <- long[!is.na(long$krt) & long$krt == "FALSE", ]
krt_exposed <- long[!is.na(long$krt) & long$krt == "TRUE", ]
krt_unknown <- long[is.na(long$krt), ]
no_krt <- scale_predictors(no_krt)
no_krt_path <- file.path(P$root, "model_fits", "reviewer16_no_krt_varpart_fit.rds")
no_krt_fit <- fit_or_load(no_krt, "clev", no_krt_path, 16017)
no_krt_result <- summarize_partition(
  partition_draws(no_krt_fit, no_krt),
  "Sensitivity excluding KRT-exposed observations",
  "Serum chloride level on days 1-3",
  no_krt
)

delta_dat <- scale_predictors(long)
delta_path <- file.path(P$root, "model_fits", "reviewer17_delta_varpart_fit.rds")
delta_fit <- fit_or_load(delta_dat, "dcl", delta_path, 17017)
delta_result <- summarize_partition(
  partition_draws(delta_fit, delta_dat),
  "Sensitivity using change from baseline",
  "Change in serum chloride from baseline",
  delta_dat
)

result <- rbind(primary_result, no_krt_result, delta_result)
out <- file.path(P$reproduction, "reviewer16_17_sensitivity_results.csv")
write.csv(result, out, row.names = FALSE)

diagnostics <- rbind(
  fit_diagnostics(primary_fit, "Primary level model"),
  fit_diagnostics(no_krt_fit, "Sensitivity excluding KRT-exposed observations"),
  fit_diagnostics(delta_fit, "Sensitivity using change from baseline")
)
diagnostic_out <- file.path(
  P$reproduction, "reviewer16_17_sensitivity_diagnostics.csv"
)
write.csv(diagnostics, diagnostic_out, row.names = FALSE)

krt_counts <- data.frame(
  total_observations = nrow(long),
  total_patients = length(unique(long$id)),
  excluded_total_observations = nrow(long) - nrow(no_krt),
  excluded_krt_observations = nrow(krt_exposed),
  patients_with_excluded_krt_observations = length(unique(krt_exposed$id)),
  excluded_unknown_krt_observations = nrow(krt_unknown),
  patients_with_unknown_krt_observations = length(unique(krt_unknown$id)),
  retained_observations = nrow(no_krt),
  retained_patients = length(unique(no_krt$id))
)
krt_count_out <- file.path(P$reproduction, "reviewer16_krt_exclusion_counts.csv")
write.csv(krt_counts, krt_count_out, row.names = FALSE)

print(result, row.names = FALSE, digits = 3)
print(diagnostics, row.names = FALSE, digits = 3)
cat("\nKRT-exposed observations excluded:", nrow(krt_exposed), "\n")
cat("Observations with unknown KRT status excluded:", nrow(krt_unknown), "\n")
cat("Results written to", out, "\n")
cat("Diagnostics written to", diagnostic_out, "\n")
cat("KRT exclusion counts written to", krt_count_out, "\n")
