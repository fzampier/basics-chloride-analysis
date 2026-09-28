# Shared construction of post-randomization study-fluid exposures.
#
# The legacy Carga_cloro_d* fields combine study and non-study chloride, and
# Carga_cloro_d0 describes the 24 hours before enrollment. They are therefore
# not used as the randomized study-fluid exposure in the revision analyses.

zero_if_missing <- function(x) ifelse(is.na(x), 0, x)

correct_fluid_volume <- function(x) {
  out <- as.numeric(x)
  idx <- !is.na(out) & out > 60000
  out[idx] <- out[idx] - 60000
  out
}

add_postbaseline_exposures <- function(d) {
  if (!"arm" %in% names(d)) {
    d$arm <- ifelse(d$grupo == "Saline Solution", "Saline", "Plasma-Lyte")
  }
  chloride_concentration <- ifelse(d$arm == "Saline", 154, 98)

  for (day in 1:3) {
    study_day <- correct_fluid_volume(d[[paste0("voflud", day, "_total")]])
    total_day <- correct_fluid_volume(d[[paste0("VolumeTotal_d", day)]])
    nonstudy_day <- zero_if_missing(d[[paste0("vofludnd", day, "_total")]])

    d[[paste0("study_volume_day", day)]] <- zero_if_missing(study_day)
    d[[paste0("total_volume_day", day)]] <- zero_if_missing(total_day)
    d[[paste0("nonstudy_volume_day", day)]] <- nonstudy_day
    d[[paste0("study_load_day", day)]] <-
      d[[paste0("study_volume_day", day)]] * chloride_concentration / 1000

    days <- seq_len(day)
    d[[paste0("study_vol", day)]] <- Reduce(
      `+`, lapply(days, function(j) d[[paste0("study_volume_day", j)]])
    )
    d[[paste0("total_vol", day)]] <- Reduce(
      `+`, lapply(days, function(j) d[[paste0("total_volume_day", j)]])
    )
    d[[paste0("nonstudy_vol", day)]] <- Reduce(
      `+`, lapply(days, function(j) d[[paste0("nonstudy_volume_day", j)]])
    )
    d[[paste0("study_load", day)]] <- Reduce(
      `+`, lapply(days, function(j) d[[paste0("study_load_day", j)]])
    )
  }

  # Dataset-specific safeguards for the documented 60,000 mL transcription
  # artifact. These fail loudly if a future data extract changes unexpectedly.
  raw_study_day2 <- d$voflud2_total
  raw_total_day2 <- d$VolumeTotal_d2
  study_outliers <- which(!is.na(raw_study_day2) & raw_study_day2 > 60000)
  total_outliers <- which(!is.na(raw_total_day2) & raw_total_day2 > 60000)
  if (!identical(study_outliers, total_outliers) || length(study_outliers) != 1L ||
      raw_study_day2[study_outliers] != 65000 ||
      raw_total_day2[total_outliers] != 65000 ||
      d$study_volume_day2[study_outliers] != 5000 ||
      d$total_volume_day2[total_outliers] != 5000) {
    stop("Unexpected day-2 fluid-volume correction pattern.")
  }
  if (any(vapply(1:3, function(day) {
    any(d[[paste0("study_volume_day", day)]] > 60000, na.rm = TRUE) ||
      any(d[[paste0("total_volume_day", day)]] > 60000, na.rm = TRUE)
  }, logical(1)))) {
    stop("Uncorrected >60,000 mL daily fluid volume remains.")
  }

  d
}

build_varpart_long <- function(d, include_krt = FALSE) {
  pieces <- lapply(1:3, function(day) {
    out <- data.frame(
      id = seq_len(nrow(d)),
      day = factor(day, levels = 1:3),
      clev = d[[paste0("clod", day)]],
      dcl = d[[paste0("clod", day)]] - d$clobas,
      load = d[[paste0("study_load", day)]],
      vol = d[[paste0("total_vol", day)]],
      clobas = d$clobas,
      age = d$age,
      creabas = d$creabas,
      sofa = d$sofa
    )
    if (include_krt) out$krt <- d[[paste0("isubrd", day, "_acum")]]
    out
  })
  do.call(rbind, pieces)
}

scale_varpart_predictors <- function(dat) {
  for (v in c("load", "vol", "clobas", "age", "creabas", "sofa")) {
    dat[[v]] <- as.numeric(scale(dat[[v]]))
  }
  dat
}
