# 01b_download_cnv.R -- pull gene-level copy number for LGG and GBM from the GDC.
#
# WHY THIS FILE EXISTS
# --------------------
# R/14_egfr_amplification.R and R/21_assertions.R both read data/<PROJECT>_cnv.rds,
# but no script in this repository ever created it: the original download was run
# ad hoc and only its logs were committed (results/_cnv_download.log,
# results/_cnv_prepare.log). The repo could therefore reproduce its expression
# results but not its copy-number results -- three of the 28 assertions in
# R/21 depend on a file the pipeline could not rebuild.
#
# The query below is reconstructed from those two logs, which record the genome
# of reference (hg38), the data type (Gene Level Copy Number), the assay names
# (copy_number, min_copy_number, max_copy_number) and the shapes: TCGA-GBM
# 60623 x 511. R/21 reads assay(x) with no name, i.e. the first assay, which is
# copy_number -- matching what the log lists first.
#
# Run from the project root, AFTER R/01_download.R:  Rscript R/01b_download_cnv.R

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

dir.create("data", showWarnings = FALSE)

download_cnv <- function(project) {
  out_rds <- file.path("data", paste0(project, "_cnv.rds"))
  if (file.exists(out_rds)) {
    message(project, ": cached ", out_rds, " already exists -- skipping.")
    return(invisible(out_rds))
  }
  message("==== ", project, ": querying GDC for gene-level copy number ====")
  query <- GDCquery(
    project       = project,
    data.category = "Copy Number Variation",
    data.type     = "Gene Level Copy Number",
    workflow.type = "ASCAT3"
  )
  GDCdownload(query, method = "api", files.per.chunk = 20)
  se <- GDCprepare(query)
  stopifnot("copy_number" %in% assayNames(se))
  saveRDS(se, out_rds)
  message(project, ": saved ", out_rds, " (", nrow(se), " x ", ncol(se),
          "; assays: ", paste(assayNames(se), collapse = ", "), ")")
  invisible(out_rds)
}

for (p in c("TCGA-LGG", "TCGA-GBM")) download_cnv(p)
