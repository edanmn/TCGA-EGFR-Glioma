# 01_download.R -- pull LGG and GBM RNA-seq + clinical from the GDC via TCGAbiolinks.
#
# Downloads STAR gene-expression counts (which also carry TPM/FPKM) for the two
# TCGA brain-cancer studies and saves each as a prepared SummarizedExperiment
# (.rds) under data/. This is the slow, network-heavy step -- run it once, then
# the analysis scripts read the cached .rds files.
#
# Run from the project root:  Rscript R/01_download.R
#
# Cohorts:
#   TCGA-LGG  Lower Grade Glioma   (grade II/III, ~500 samples)
#   TCGA-GBM  Glioblastoma         (grade IV,     ~170 samples with STAR counts)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

dir.create("data", showWarnings = FALSE)

download_cohort <- function(project) {
  out_rds <- file.path("data", paste0(project, "_expr_se.rds"))
  if (file.exists(out_rds)) {
    message(project, ": cached ", out_rds, " already exists -- skipping.")
    return(invisible(out_rds))
  }
  message("==== ", project, ": querying GDC ====")
  query <- GDCquery(
    project       = project,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query, method = "api", files.per.chunk = 20)
  se <- GDCprepare(query)          # SummarizedExperiment
  saveRDS(se, out_rds)
  message(project, ": saved ", out_rds,
          " (", nrow(se), " genes x ", ncol(se), " samples)")
  invisible(out_rds)
}

for (p in c("TCGA-LGG", "TCGA-GBM")) download_cohort(p)

message("Download step complete. Prepared objects are in data/.")
