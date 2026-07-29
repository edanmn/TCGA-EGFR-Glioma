# 04_pathway_screen.R -- screen the EGFR / RTK-PI3K pathway for survival
# association across cohorts, with Benjamini-Hochberg (FDR) correction for
# testing multiple genes.
#
# Run from the project root:  Rscript R/04_pathway_screen.R
# Produces: results/pathway_screen.csv
#
# For each gene we fit an age-adjusted Cox model (continuous log2 TPM) and
# collect the gene's hazard ratio and p-value, then adjust p across genes
# WITHIN each cohort.

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})
source("R/_helpers.R")

# EGFR and its core pathway partners in glioma.
GENES <- c("EGFR", "PDGFRA", "PTEN", "NF1", "PIK3CA", "PIK3R1",
           "RB1", "CDKN2A", "TP53", "IDH1")

dir.create("results", showWarnings = FALSE)

lgg <- build_analysis("TCGA-LGG", GENES)
gbm <- build_analysis("TCGA-GBM", GENES)
combined <- bind_rows(lgg, gbm); combined$cohort <- "LGG+GBM"
cohorts <- list("TCGA-LGG" = lgg, "TCGA-GBM" = gbm, "LGG+GBM" = combined)

screen_one <- function(df, gene) {
  if (!gene %in% names(df)) return(NULL)
  m <- tryCatch(
    coxph(Surv(time, event) ~ df[[gene]] + age, data = df),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)$coefficients
  ci <- summary(m)$conf.int
  data.frame(
    gene = gene,
    HR = signif(unname(ci[1, "exp(coef)"]), 3),
    p  = s[1, "Pr(>|z|)"],
    n = m$n, events = m$nevent,
    stringsAsFactors = FALSE
  )
}

all_rows <- list()
for (nm in names(cohorts)) {
  df <- cohorts[[nm]]
  res <- do.call(rbind, lapply(GENES, function(g) screen_one(df, g)))
  res$cohort <- nm
  res$p_BH <- p.adjust(res$p, method = "BH")     # FDR within cohort
  res <- res[order(res$p_BH), ]
  res$p <- signif(res$p, 3); res$p_BH <- signif(res$p_BH, 3)
  all_rows <- c(all_rows, list(res))
  cat("\n==== ", nm, " (age-adjusted, BH-corrected) ====\n")
  print(res[, c("gene", "HR", "p", "p_BH", "n", "events")], row.names = FALSE)
}

screen_tbl <- do.call(rbind, all_rows)[, c("cohort", "gene", "HR", "p", "p_BH", "n", "events")]
write.csv(screen_tbl, "results/pathway_screen.csv", row.names = FALSE)
cat("\nFull screen written to results/pathway_screen.csv\n")
cat("Significant after FDR: p_BH < 0.05.\n")
