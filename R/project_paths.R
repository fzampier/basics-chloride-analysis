# Central path resolver for the self-contained revision workspace.

basics_revision_root <- function() {
  configured <- Sys.getenv("BASICS_REVISION_ROOT", unset = "")
  if (nzchar(configured)) {
    return(normalizePath(configured, mustWork = TRUE))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
    candidate <- dirname(script)
    if (basename(candidate) == "R") candidate <- dirname(candidate)
    if (dir.exists(file.path(candidate, "raw_data"))) return(candidate)
  }

  candidate <- normalizePath(getwd(), mustWork = TRUE)
  if (basename(candidate) == "R") candidate <- dirname(candidate)
  if (dir.exists(file.path(candidate, "raw_data"))) return(candidate)

  stop("Cannot locate revision root. Run from revision/ or set BASICS_REVISION_ROOT.")
}

basics_paths <- function(root = basics_revision_root()) {
  list(
    root = root,
    data_file = file.path(root, "raw_data", "Dados_BASICS_MainPaper.RData"),
    cache_file = file.path(root, "derived_data", ".chl_cache.rds"),
    brms_fit = file.path(root, "model_fits", "chl_brms_fit.rds"),
    varpart_fit = file.path(root, "model_fits", "chl_varpart_fit.rds"),
    figures = file.path(root, "figures"),
    reproduction = file.path(root, "reproduction")
  )
}

BASICS_PATHS <- basics_paths()
