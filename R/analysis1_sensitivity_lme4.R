#!/usr/bin/env Rscript

# Independent frequentist corroboration of the two qualitative sensitivity
# statements for the Bayesian longitudinal effect-modification analysis.
# Participant-level values are never written; output contains aggregate model
# dimensions and day-specific interaction estimates only.
suppressPackageStartupMessages(library(lme4))

args <- commandArgs(trailingOnly = FALSE)
this_file <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][[1]]), mustWork = TRUE)
revision_root <- dirname(dirname(this_file))
Sys.setenv(BASICS_REVISION_ROOT = revision_root)
source(file.path(revision_root, "R", "project_paths.R"))

C <- readRDS(BASICS_PATHS$cache_file)
d <- C$d
d$id <- seq_len(nrow(d))

egfr <- function(scr, age, female) {
  k <- ifelse(female == 1, 0.7, 0.9)
  a <- ifelse(female == 1, -0.241, -0.302)
  142 * pmin(scr / k, 1)^a * pmax(scr / k, 1)^(-1.2) * 0.9938^age * ifelse(female == 1, 1.012, 1)
}

d$eGFR <- egfr(d$creabas, d$age, d$female)
d$eGFR[!is.finite(d$eGFR) | d$eGFR < 5 | d$eGFR > 200] <- NA
d$eGFRc <- (d$eGFR - 90) / 15
d$creatc <- as.numeric(scale(d$creabas))

make_long <- function(data, kidney_variable) {
  pieces <- lapply(1:3, function(study_day) {
    chloride <- data[[paste0("clod", study_day)]]
    keep <- !is.na(chloride) & !is.na(data$clobas) & !is.na(data[[kidney_variable]])
    data.frame(
      id = data$id[keep],
      chloride = chloride[keep],
      day = study_day,
      dayf = factor(study_day, levels = 1:3),
      arm = factor(data$arm[keep], levels = c("Plasma-Lyte", "Saline")),
      baseline_chloride = data$clobas[keep] - 106,
      kidney = data[[kidney_variable]][keep]
    )
  })
  do.call(rbind, pieces)
}

fit_model <- function(data) {
  lmer(
    chloride ~ arm * dayf * kidney + baseline_chloride + (1 + day | id),
    data = data,
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
  )
}

day_interactions <- function(fit, label, data, scale_label) {
  coefficients <- fixef(fit)
  covariance <- vcov(fit)
  base_name <- "armSaline:kidney"
  rows <- lapply(2:3, function(study_day) {
    increment_name <- paste0("armSaline:dayf", study_day, ":kidney")
    contrast <- setNames(rep(0, length(coefficients)), names(coefficients))
    contrast[c(base_name, increment_name)] <- 1
    estimate <- sum(contrast * coefficients)
    standard_error <- sqrt(as.numeric(t(contrast) %*% covariance %*% contrast))
    data.frame(
      analysis = label,
      observations = nrow(data),
      participants = length(unique(data$id)),
      day = study_day,
      scale = scale_label,
      estimate = estimate,
      standard_error = standard_error,
      conf_low = estimate - 1.96 * standard_error,
      conf_high = estimate + 1.96 * standard_error,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

primary <- make_long(d, "eGFRc")
all_days_ids <- Reduce(intersect, lapply(1:3, function(study_day) {
  which(!is.na(d[[paste0("clod", study_day)]]) & !is.na(d$clobas) & !is.na(d$eGFRc))
}))
all_days <- primary[primary$id %in% all_days_ids, ]
creatinine <- make_long(d, "creatc")

results <- rbind(
  day_interactions(fit_model(primary), "Primary eGFR model", primary, "per 15 mL/min/1.73 m2 higher eGFR"),
  day_interactions(fit_model(all_days), "All-three-days eGFR sensitivity", all_days, "per 15 mL/min/1.73 m2 higher eGFR"),
  day_interactions(fit_model(creatinine), "Baseline-creatinine sensitivity", creatinine, "per 1 SD higher baseline creatinine")
)

output <- file.path(BASICS_PATHS$reproduction, "analysis1_sensitivity_lme4.csv")
write.csv(results, output, row.names = FALSE)
print(results, row.names = FALSE)
