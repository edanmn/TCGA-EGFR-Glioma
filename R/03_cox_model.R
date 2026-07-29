# 03_cox_model.R -- Cox proportional-hazards models for EGFR, controlling for
# age and tumor grade. This is the honest version of the question: does EGFR
# expression predict survival AFTER accounting for the fact that it tracks with
# grade?
#
# Run from the project root:  Rscript R/03_cox_model.R
# Produces: results/cox_EGFR_summary.csv
#
# Models per cohort:
#   univariate:   Surv(time,event) ~ EGFR (continuous log2 TPM)
#   adjusted:     Surv(time,event) ~ EGFR + age + grade
# GBM alone is uniformly grade IV, so grade is dropped there automatically.

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})
source("R/_helpers.R")

GENE <- "EGFR"
dir.create("results", showWarnings = FALSE)

lgg <- build_analysis("TCGA-LGG", GENE)
gbm <- build_analysis("TCGA-GBM", GENE)
combined <- bind_rows(lgg, gbm); combined$cohort <- "LGG+GBM"
cohorts <- list("TCGA-LGG" = lgg, "TCGA-GBM" = gbm, "LGG+GBM" = combined)

# Collapse messy grade strings to an ordered-ish factor when possible.
tidy_grade <- function(g) {
  g <- toupper(as.character(g))
  g[grepl("G2|GRADE II$|GRADE 2| II ", g)] <- "II"
  g[grepl("G3|GRADE III|GRADE 3| III ", g)] <- "III"
  g[grepl("G4|GRADE IV|GRADE 4| IV ", g)] <- "IV"
  factor(g)
}

extract <- function(model, term, cohort, kind) {
  s <- summary(model)$coefficients
  if (!term %in% rownames(s)) return(NULL)
  ci <- summary(model)$conf.int[term, ]
  data.frame(
    cohort = cohort, model = kind, term = term,
    HR = signif(unname(ci["exp(coef)"]), 3),
    CI_low = signif(unname(ci["lower .95"]), 3),
    CI_high = signif(unname(ci["upper .95"]), 3),
    p = signif(s[term, "Pr(>|z|)"], 3),
    n = model$n, events = model$nevent,
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (nm in names(cohorts)) {
  df <- cohorts[[nm]]
  df$grade_f <- tidy_grade(df$grade)

  # Univariate
  # Reference the gene as a bare column name so the model term is named `GENE`
  # (a df[[GENE]] term would be named "df[[GENE]]" and break extract()).
  gcol <- paste0("`", GENE, "`")
  m1 <- coxph(as.formula(paste0("Surv(time, event) ~ ", gcol)), data = df)
  r1 <- extract(m1, GENE, nm, "univariate")

  # Adjusted -- include grade only if it actually varies in this cohort.
  use_grade <- nlevels(droplevels(df$grade_f[!is.na(df$grade_f)])) >= 2
  fml <- if (use_grade)
    as.formula(paste0("Surv(time, event) ~ ", gcol, " + age + grade_f")) else
    as.formula(paste0("Surv(time, event) ~ ", gcol, " + age"))
  m2 <- coxph(fml, data = df)
  r2 <- extract(m2, GENE, nm, if (use_grade) "adjusted(age+grade)" else "adjusted(age)")

  rows <- c(rows, list(r1, r2))
  cat("\n==== ", nm, " ====\n"); print(summary(m2))
}

cox_tbl <- do.call(rbind, rows)
write.csv(cox_tbl, "results/cox_EGFR_summary.csv", row.names = FALSE)

cat("\n=== Cox HR summary (EGFR term) ===\n")
print(cox_tbl, row.names = FALSE)
cat("\nSummary written to results/cox_EGFR_summary.csv\n")
cat("Read HR as hazard per 1-unit increase in log2(TPM+1) EGFR; HR>1 = worse survival.\n")
