# _verify_clinical.R -- clinical-field derivation for the VERIFICATION suites only.
#
# WHY THIS FILE EXISTS, AND WHY IT MUST NOT BE "TIDIED AWAY"
#
# R/21 and R/22 check the manuscript's sample sizes and headline models. Until
# round 9 they did that by calling load_se() and .clinical_df() out of
# R/_helpers.R -- the module 31 analysis scripts also use. That made the check
# circular in the one direction that matters. The manuscript's numbers are
# produced BY the pipeline, so if .clinical_df had mis-derived a survival time
# or an age, the pipeline would have written a wrong number into the paper and
# the audit would have re-derived the SAME wrong number through the SAME
# function and reported a match. Two of this project's three original data bugs
# were in exactly this kind of clinical-field parsing.
#
# This file re-derives those fields from colData without touching _helpers.R,
# so a defect in the pipeline's parsing now shows up as a mismatch against the
# published values instead of being reproduced by the checker.
#
# RULES:
#   * Never source R/_helpers.R from here, and never make this file call it.
#   * Never make an analysis script (00-20, 23-40) source THIS file.
#   * If you change a definition here, change it because TCGA's schema changed,
#     not to make a failing check agree with _helpers.R. A disagreement between
#     the two IS the signal this file exists to produce.
#
# Sourced by R/21_assertions.R and R/22_stat_recompute.R. Not run on its own.

suppressPackageStartupMessages({ library(SummarizedExperiment) })

# Read the cached SummarizedExperiment straight off disk.
vload_se <- function(project) {
  f <- file.path("data", paste0(project, "_expr_se.rds"))
  if (!file.exists(f))
    stop("Missing ", f, " -- run R/01_download.R first.")
  readRDS(f)
}

# Overall survival and demographics, one row per column of `se`, in the same
# order as colnames(se). Derived here from first principles:
#
#   event = 1 iff vital_status is "Dead"
#   time  = days_to_death for a death, otherwise days_to_last_follow_up
#   age   = age_at_index in years, or age_at_diagnosis (days) / 365.25
#
# Deliberately does NOT derive grade: R/21 and R/22 each build grade from
# paper_Grade themselves, so the study-definition GBM rule in _helpers.R is
# not on this path at all.
vclinical <- function(se) {
  cd <- as.data.frame(colData(se))
  col <- function(...) {
    for (nm in c(...)) if (nm %in% names(cd)) return(cd[[nm]])
    rep(NA, nrow(cd))
  }
  num <- function(x) suppressWarnings(as.numeric(x))

  vital    <- as.character(col("vital_status"))
  is_dead  <- !is.na(vital) & vital == "Dead"
  d_death  <- num(col("days_to_death"))
  d_follow <- num(col("days_to_last_follow_up", "days_to_last_followup"))

  age_yr  <- num(col("age_at_index"))
  age_day <- num(col("age_at_diagnosis"))
  age     <- ifelse(!is.na(age_yr), age_yr, age_day / 365.25)

  bc <- colnames(se)
  data.frame(
    barcode     = bc,
    patient     = substr(bc, 1, 12),
    sample_type = substr(bc, 14, 15),
    time        = ifelse(is_dead, d_death, d_follow),
    event       = as.integer(is_dead),
    age         = age,
    sex         = as.character(col("gender")),
    stringsAsFactors = FALSE
  )
}
