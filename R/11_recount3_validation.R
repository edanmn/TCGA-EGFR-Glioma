# 11_recount3_validation.R -- pipeline-robustness check via recount3.
#
# recount3 provides TCGA RNA-seq uniformly reprocessed FROM RAW READS by the
# Monorail pipeline (hg38/Gencode), independent of the GDC STAR pipeline used
# elsewhere in this study. We test whether, under this independent raw-read
# reprocessing: (a) EGFR passes the positive control (rises with grade, higher
# in IDH-wildtype), and (b) the IDH-confounding result replicates.
#
# This is a PROCESSING-robustness test on the SAME patients, not new-patient
# external validation (CGGA raw reads are controlled-access; see paper).
#
# Run from project root:  Rscript R/11_recount3_validation.R
# Writes: results/recount3_validation.csv

suppressPackageStartupMessages({
  library(recount3); library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")

# ---- clinical (IDH, grade, survival) from the cached GDC objects, by patient ----
clin_for <- function(project) {
  se <- load_se(project); cd <- as.data.frame(colData(se))
  base <- .clinical_df(se, project)
  base$idh <- as.character(cd$paper_IDH.status)
  base$grade_raw <- as.character(cd$paper_Grade)
  base %>% filter(sample_type == "01", !is.na(time), time > 0, !is.na(event)) %>%
    transmute(patient, time, event, age,
              idh_f = factor(ifelse(idh %in% c("WT","Mutant"), idh, NA), levels=c("WT","Mutant")),
              grade_f = factor(dplyr::case_when(grepl("G2|II", grade_raw) & !grepl("III|IV", grade_raw) ~ "II",
                                                grepl("G3|III", grade_raw) ~ "III",
                                                grepl("G4|IV", grade_raw) ~ "IV", TRUE ~ NA_character_),
                               levels=c("II","III","IV"))) %>%
    distinct(patient, .keep_all = TRUE)
}

# ---- recount3 expression: EGFR log2(CPM+1) by patient ----
egfr_recount3 <- function(tcga_project) {
  ap <- available_projects()
  pr <- subset(ap, project == tcga_project & file_source == "tcga")
  rse <- create_rse(pr)
  cts <- assay(rse, "raw_counts")
  cpm <- log2(sweep(cts, 2, colSums(cts), "/") * 1e6 + 1)
  sym <- rowData(rse)$gene_name
  egfr <- cpm[which(sym == "EGFR")[1], ]
  bc_col <- grep("barcode", tolower(colnames(colData(rse))), value = TRUE)[1]
  bc <- if (!is.na(bc_col)) as.character(colData(rse)[[bc_col]]) else colData(rse)$external_id
  data.frame(barcode = bc, stringsAsFactors = FALSE) %>%
    mutate(patient = substr(barcode, 1, 12), sample_type = substr(barcode, 14, 15),
           EGFR = egfr) %>%
    filter(sample_type == "01") %>% distinct(patient, .keep_all = TRUE)
}

run <- function(project, tcga_project) {
  clin <- clin_for(project)
  ex <- egfr_recount3(tcga_project)
  d <- inner_join(clin, ex[, c("patient","EGFR")], by = "patient")
  d$EGFR_z <- as.numeric(scale(d$EGFR))
  cat(sprintf("\n==== %s (recount3): n=%d events=%d matched ====\n", project, nrow(d), sum(d$event)))
  cat("EGFR by grade (positive control, expect rising):\n"); print(round(tapply(d$EGFR, d$grade_f, mean, na.rm=TRUE), 2))
  cat("EGFR by IDH (expect WT>Mutant):\n"); print(round(tapply(d$EGFR, d$idh_f, mean, na.rm=TRUE), 2))
  gr_ok <- nlevels(droplevels(d$grade_f[!is.na(d$grade_f)])) >= 2
  fits <- list(
    c("EGFR_z", "1. unadjusted"),
    c(if (gr_ok) "EGFR_z + age + grade_f" else "EGFR_z + age", if (gr_ok) "2. +age+grade" else "2. +age"),
    c(if (gr_ok) "EGFR_z + age + grade_f + idh_f" else "EGFR_z + age + idh_f", "3. +age+grade+IDH"))
  do.call(rbind, lapply(fits, function(f) {
    dd <- d; if (grepl("idh", f[1])) dd <- d[!is.na(d$idh_f), ]
    m <- tryCatch(coxph(as.formula(paste("Surv(time,event)~", f[1])), data = dd), error=function(e) NULL)
    if (is.null(m) || !"EGFR_z" %in% rownames(summary(m)$coefficients)) return(NULL)
    s <- summary(m); ci <- s$conf.int["EGFR_z", ]
    data.frame(cohort=project, model=f[2], HR=signif(ci[1],3), CI_low=signif(ci[3],3),
               CI_high=signif(ci[4],3), p=signif(s$coefficients["EGFR_z","Pr(>|z|)"],3),
               n=m$n, events=m$nevent, stringsAsFactors=FALSE)
  }))
}

res <- rbind(run("TCGA-LGG", "LGG"), run("TCGA-GBM", "GBM"))
write.csv(res, "results/recount3_validation.csv", row.names = FALSE)
cat("\n=== recount3 (Monorail raw-read reprocessing) EGFR Cox, per SD ===\n")
print(res, row.names = FALSE)
