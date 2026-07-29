# _helpers.R -- shared functions to turn a prepared SummarizedExperiment into a
# tidy per-patient table with expression, overall survival, and covariates.
# Sourced by 02/03/04. Not run on its own.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
})

load_se <- function(project) {
  f <- file.path("data", paste0(project, "_expr_se.rds"))
  if (!file.exists(f))
    stop("Missing ", f, " -- run R/01_download.R first.")
  readRDS(f)
}

# log2(TPM+1) expression for the requested gene symbols.
# Returns a data.frame: one row per sample barcode, one column per gene found.
.expr_for_genes <- function(se, genes) {
  tpm <- assay(se, "tpm_unstrand")
  sym <- rowData(se)$gene_name
  keep <- match(genes, sym)               # first match per symbol
  found <- genes[!is.na(keep)]
  missing <- genes[is.na(keep)]
  if (length(missing))
    message("  (genes not found, skipped: ", paste(missing, collapse = ", "), ")")
  mat <- tpm[keep[!is.na(keep)], , drop = FALSE]
  mat <- log2(mat + 1)
  out <- as.data.frame(t(mat))
  colnames(out) <- found
  out$barcode <- colnames(se)
  out
}

# Pull overall-survival + covariates out of colData, coping with the columns
# actually present (they vary a little between projects).
.clinical_df <- function(se, project) {
  cd <- as.data.frame(colData(se))
  pick <- function(cands) {
    hit <- cands[cands %in% names(cd)]
    if (length(hit)) cd[[hit[1]]] else rep(NA, nrow(cd))
  }
  vital <- as.character(pick(c("vital_status")))
  d_death <- suppressWarnings(as.numeric(pick(c("days_to_death"))))
  d_follow <- suppressWarnings(as.numeric(pick(c("days_to_last_follow_up",
                                                 "days_to_last_followup"))))
  age <- suppressWarnings(as.numeric(pick(c("age_at_index"))))
  age_dx <- suppressWarnings(as.numeric(pick(c("age_at_diagnosis"))))  # in days
  age <- ifelse(is.na(age) & !is.na(age_dx), age_dx / 365.25, age)
  sex <- as.character(pick(c("gender")))

  # Grade: try known columns; GBM is uniformly grade IV if nothing is recorded.
  grade_raw <- pick(c("paper_Grade", "paper_Histology", "tumor_grade",
                      "neoplasm_histologic_grade", "primary_diagnosis"))
  grade <- as.character(grade_raw)
  grade[grade %in% c("", "not reported", "NA", "Not Reported")] <- NA
  if (grepl("GBM", project)) grade[is.na(grade)] <- "G4"

  time <- ifelse(vital == "Dead", d_death, d_follow)
  data.frame(
    barcode = rownames(cd),
    patient = substr(rownames(cd), 1, 12),
    sample_type = substr(rownames(cd), 14, 15),
    time  = as.numeric(time),
    event = as.integer(vital == "Dead"),
    age   = age,
    sex   = sex,
    grade = grade,
    stringsAsFactors = FALSE
  )
}

# Full analysis table for one cohort: primary tumors only, one row per patient,
# valid positive survival time, with the requested genes attached.
build_analysis <- function(project, genes) {
  se   <- load_se(project)
  expr <- .expr_for_genes(se, genes)
  clin <- .clinical_df(se, project)
  df <- clin %>%
    left_join(expr, by = "barcode") %>%
    filter(sample_type == "01") %>%              # primary solid tumor only
    filter(!is.na(time) & time > 0 & !is.na(event)) %>%
    arrange(patient) %>%
    distinct(patient, .keep_all = TRUE)          # one sample per patient
  df$cohort <- project
  df
}

# Add a per-cohort median split for one gene: "High"/"Low" factor (Low = ref).
add_median_split <- function(df, gene) {
  stopifnot(gene %in% names(df))
  med <- median(df[[gene]], na.rm = TRUE)
  df[[paste0(gene, "_grp")]] <-
    factor(ifelse(df[[gene]] > med, "High", "Low"), levels = c("Low", "High"))
  df
}
