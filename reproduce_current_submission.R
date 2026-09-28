#!/usr/bin/env Rscript

# Reproduce the current revision's aggregate cohort, Table 1, and Figures 1-4
# using the duplicated raw dataset and the saved Bayesian model objects.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (!length(file_arg)) stop("Run this file with Rscript.")
this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
revision_root <- dirname(this_file)
Sys.setenv(BASICS_REVISION_ROOT = revision_root)
source(file.path(revision_root, "R", "project_paths.R"))

dir.create(BASICS_PATHS$reproduction, recursive = TRUE, showWarnings = FALSE)
generated <- file.path(BASICS_PATHS$reproduction, "generated_figures")
dir.create(generated, recursive = TRUE, showWarnings = FALSE)

run_logged <- function(label, script_name, log_name) {
  cat(sprintf("[%s] starting\n", label))
  log_path <- file.path(BASICS_PATHS$reproduction, log_name)
  con <- file(log_path, open = "wt")
  sink(con)
  sink(con, type = "message")
  error <- NULL
  tryCatch(
    source(file.path(revision_root, "R", script_name), local = .GlobalEnv),
    error = function(e) error <<- e
  )
  sink(type = "message")
  sink()
  close(con)
  if (!is.null(error)) stop(sprintf("%s failed: %s", label, conditionMessage(error)))
  cat(sprintf("[%s] complete: %s\n", label, log_path))
}

run_logged("cohort and key numbers", "basics_chloride_analysis.R", "analysis_output.txt")
run_logged("Table 1", "table1_final.R", "table1_output.txt")

old_fig_dir <- Sys.getenv("BASICS_FIGURE_DIR", unset = NA_character_)
Sys.setenv(BASICS_FIGURE_DIR = generated)
run_logged("Figures 1-4 from saved fits", "basics_chloride_figures.R", "figures_output.txt")
if (is.na(old_fig_dir)) Sys.unsetenv("BASICS_FIGURE_DIR") else Sys.setenv(BASICS_FIGURE_DIR = old_fig_dir)

analysis_text <- paste(readLines(file.path(BASICS_PATHS$reproduction, "analysis_output.txt"), warn = FALSE), collapse = "\n")
table_text <- paste(readLines(file.path(BASICS_PATHS$reproduction, "table1_output.txt"), warn = FALSE), collapse = "\n")
if (!grepl("Analytic cohort n = 3482 \\| day-2 paired n = 2721", analysis_text)) {
  stop("Aggregate analysis did not reproduce the submitted cohort counts.")
}
if (!grepl("ANALYTIC COHORT n = 3482", table_text) ||
    !grepl("Saline 1736  Plasma-Lyte 1746", table_text) ||
    !grepl("PAIRED DAY-2 n = 2721", table_text)) {
  stop("Table 1 did not reproduce the submitted cohort and arm counts.")
}

volume_results_path <- file.path(BASICS_PATHS$reproduction, "analysis2_volume_results.csv")
if (!file.exists(volume_results_path)) stop("Analysis 2 volume-results file was not generated.")
volume_results <- read.csv(volume_results_path, check.names = FALSE)
expected_volume_results <- rbind(
  unweighted = c(0.3059921977, 0.2077739384, 0.4229671080),
  ipw_adjusted = c(0.3061768536, 0.2170426174, 0.4385250661)
)
observed_volume_results <- as.matrix(volume_results[, c("estimate", "conf_low", "conf_high")])
if (!identical(dim(observed_volume_results), dim(expected_volume_results)) ||
    max(abs(observed_volume_results - expected_volume_results)) > 1e-8) {
  stop("Analysis 2 saline-only volume estimates or confidence intervals did not reproduce.")
}

concentration_results_path <- file.path(BASICS_PATHS$reproduction, "analysis2_concentration_results.csv")
if (!file.exists(concentration_results_path)) stop("Analysis 2 concentration-results file was not generated.")
concentration_results <- read.csv(concentration_results_path, check.names = FALSE)
if (nrow(concentration_results) != 1L ||
    concentration_results$analysis_n != 2721L ||
    concentration_results$exposure_scale != "Per 100 mmol between-group excess post-randomization study-fluid chloride associated with randomized fluid assignment" ||
    abs(concentration_results$chloride_change_difference - 2.1454173753) > 1e-8 ||
    abs(concentration_results$chloride_load_difference - 145.7030329939) > 1e-8 ||
    abs(concentration_results$estimate - 1.4724589676) > 1e-8 ||
    abs(concentration_results$conf_low - 1.1725081520) > 1e-8 ||
    abs(concentration_results$conf_high - 1.8102120216) > 1e-8) {
  stop("Analysis 2 randomized between-group differences, normalized ratio, or bootstrap interval did not reproduce.")
}

writeLines(capture.output(sessionInfo()), file.path(BASICS_PATHS$reproduction, "session_info.txt"))
cat("Current corrected exposure analysis and participant-bootstrap Analysis 2 intervals reproduced.\n")
