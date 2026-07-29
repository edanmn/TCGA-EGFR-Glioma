# 02_survival_km.R -- Kaplan-Meier + log-rank for EGFR (median split) in each
# brain-cancer cohort, plus the combined LGG+GBM cohort.
#
# Run from the project root:  Rscript R/02_survival_km.R
# Produces: figures/KM_EGFR_<cohort>.pdf and results/km_logrank_summary.csv

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})
source("R/_helpers.R")

# survminer makes nicer KM plots (curve + risk table) but drags in a heavy
# dependency chain (nloptr/lme4/car). It's optional -- if it isn't installed we
# fall back to base-R plot.survfit. All statistics below use base `survival`.
HAVE_SURVMINER <- requireNamespace("survminer", quietly = TRUE)
if (!HAVE_SURVMINER)
  message("survminer not installed -- using base-R KM plots. ",
          "Install it (needs CMake) for curve + risk-table figures.")

GENE <- "EGFR"
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Build per-cohort tables and a combined one.
lgg <- build_analysis("TCGA-LGG", GENE)
gbm <- build_analysis("TCGA-GBM", GENE)
combined <- bind_rows(lgg, gbm)
combined$cohort <- "LGG+GBM"

cohorts <- list("TCGA-LGG" = lgg, "TCGA-GBM" = gbm, "LGG+GBM" = combined)

run_km <- function(df, label) {
  df <- add_median_split(df, GENE)
  # Use a real, fixed-name grouping column ("grp") and a plain formula.
  # ggsurvplot re-evaluates the fit's formula in `data`, so a `df[[var]]`
  # term (var being a local variable) breaks with 'object not found'.
  df$grp <- df[[paste0(GENE, "_grp")]]
  fit <- survfit(Surv(time, event) ~ grp, data = df)
  sd  <- survdiff(Surv(time, event) ~ grp, data = df)
  p   <- 1 - pchisq(sd$chisq, length(sd$n) - 1)

  # median survival (days) per group
  med <- summary(fit)$table
  med_low  <- med[grep("Low",  rownames(med)), "median"]
  med_high <- med[grep("High", rownames(med)), "median"]

  fig <- file.path("figures", paste0("KM_", GENE, "_",
                                     gsub("[^A-Za-z0-9]", "", label), ".pdf"))
  ttl <- paste0(label, ": overall survival by ", GENE,
                " expression (median split)")
  if (HAVE_SURVMINER) {
    g <- survminer::ggsurvplot(
      fit, data = df, pval = TRUE, risk.table = TRUE,
      legend.labs = c(paste0(GENE, " Low"), paste0(GENE, " High")),
      palette = c("#2c7fb8", "#d95f0e"),
      xlab = "Days", ylab = "Overall survival", title = ttl
    )
    # A ggsurvplot carries a risk table, so render it through its own print
    # method on a PDF device rather than ggsave().
    pdf(fig, width = 7, height = 8); print(g); dev.off()
  } else {
    pdf(fig, width = 7, height = 6)
    plot(fit, col = c("#2c7fb8", "#d95f0e"), lwd = 2,
         xlab = "Days", ylab = "Overall survival", main = ttl)
    legend("topright", bty = "n", lwd = 2, col = c("#2c7fb8", "#d95f0e"),
           legend = c(paste0(GENE, " Low"), paste0(GENE, " High")))
    legend("bottomleft", bty = "n",
           legend = paste0("log-rank p = ", signif(p, 3)))
    dev.off()
  }

  data.frame(
    cohort = label,
    gene = GENE,
    n_low = sum(df$grp == "Low"),
    n_high = sum(df$grp == "High"),
    events = sum(df$event),
    median_surv_low_days = as.numeric(med_low),
    median_surv_high_days = as.numeric(med_high),
    logrank_p = signif(p, 3),
    stringsAsFactors = FALSE
  )
}

summary_tbl <- do.call(rbind, Map(run_km, cohorts, names(cohorts)))
write.csv(summary_tbl, "results/km_logrank_summary.csv", row.names = FALSE)

cat("\n=== Kaplan-Meier / log-rank summary ===\n")
print(summary_tbl, row.names = FALSE)
cat("\nFigures written to figures/, summary to results/km_logrank_summary.csv\n")
