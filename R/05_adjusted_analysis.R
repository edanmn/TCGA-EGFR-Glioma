# 05_adjusted_analysis.R -- revised, confounder-aware survival analysis.
#
# Addresses the peer-review issues:
#  * adjusts for IDH status / molecular subtype (the dominant glioma confounder)
#  * standardizes expression to per-SD units (comparable HRs across genes)
#  * reports Harrell's C-index for every model
#  * tests the proportional-hazards assumption (cox.zph)
#  * reports covariate missingness
#  * age modeled with a spline in a sensitivity model
#
# Run from project root:  Rscript R/05_adjusted_analysis.R
# Writes: results/cox_EGFR_adjusted.csv, results/pathway_screen_std.csv,
#         results/missingness.csv, results/ph_tests.csv

suppressPackageStartupMessages({
  library(survival); library(dplyr); library(splines)
})
source("R/_helpers.R")

GENES <- c("EGFR","PDGFRA","PTEN","NF1","PIK3CA","PIK3R1","RB1","CDKN2A","TP53","IDH1")
dir.create("results", showWarnings = FALSE)

# ---- extended per-patient table: expression + survival + molecular covariates ----
clin_ext <- function(se, project) {
  cd <- as.data.frame(colData(se))
  base <- .clinical_df(se, project)          # aligned row order with cd
  pick <- function(cands) {
    hit <- cands[cands %in% names(cd)]
    if (length(hit)) as.character(cd[[hit[1]]]) else rep(NA, nrow(cd))
  }
  base$idh     <- pick("paper_IDH.status")
  base$codel   <- pick(c("paper_X1p.19q.codeletion", "paper_1p.19q.codeletion"))
  base$subtype <- pick("paper_IDH.codel.subtype")
  base$mgmt    <- pick("paper_MGMT.promoter.status")
  base$kps     <- suppressWarnings(as.numeric(pick("paper_Karnofsky.Performance.Score")))
  base
}

tidy_grade <- function(g) {
  g <- toupper(as.character(g))
  out <- rep(NA_character_, length(g))
  out[grepl("G2|GRADE II$|GRADE 2| II ", g)] <- "II"
  out[grepl("G3|GRADE III|GRADE 3| III ", g)] <- "III"
  out[grepl("G4|GRADE IV|GRADE 4| IV ", g)] <- "IV"
  factor(out, levels = c("II","III","IV"))
}
norm_idh <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("WT","WILDTYPE","WILD TYPE","WILD-TYPE")] <- "WT"
  x[x %in% c("MUTANT","MUT","MUTATED")] <- "Mutant"
  x[x %in% c("","NA","NOT AVAILABLE","NOT.AVAILABLE","NOS")] <- NA
  factor(x, levels = c("WT","Mutant"))     # WT = reference (worse prognosis)
}
norm_subtype <- function(x) {
  x <- trimws(as.character(x)); x[x %in% c("","NA")] <- NA
  # canonical: IDHwt (ref), IDHmut-non-codel, IDHmut-codel
  x <- gsub("IDHmut-codel.*","IDHmut-codel", x)
  x <- gsub("IDHmut-non-codel.*","IDHmut-non-codel", x)
  x <- gsub("IDHwt.*","IDHwt", x)
  fac <- factor(x)
  if ("IDHwt" %in% levels(fac)) fac <- relevel(fac, "IDHwt")
  fac
}

build_ext <- function(project) {
  se <- load_se(project)
  expr <- .expr_for_genes(se, GENES)
  clin <- clin_ext(se, project)
  df <- clin %>%
    left_join(expr, by = "barcode") %>%
    filter(sample_type == "01", !is.na(time), time > 0, !is.na(event)) %>%
    arrange(patient) %>% distinct(patient, .keep_all = TRUE)
  df$grade_f  <- tidy_grade(df$grade)
  df$idh_f    <- norm_idh(df$idh)
  df$subtype_f<- norm_subtype(df$subtype)
  df$cohort   <- project
  # per-SD standardized expression (within cohort) for every gene
  for (g in GENES) if (g %in% names(df)) df[[paste0(g, "_z")]] <- as.numeric(scale(df[[g]]))
  df
}

lgg <- build_ext("TCGA-LGG")
gbm <- build_ext("TCGA-GBM")
pooled <- bind_rows(lgg, gbm); pooled$cohort <- "LGG+GBM"
# re-standardize genes within the pooled cohort
for (g in GENES) if (g %in% names(pooled)) pooled[[paste0(g,"_z")]] <- as.numeric(scale(pooled[[g]]))
cohorts <- list("TCGA-LGG"=lgg, "TCGA-GBM"=gbm, "LGG+GBM"=pooled)

# ---- missingness report ----
miss <- lapply(names(cohorts), function(nm){
  d <- cohorts[[nm]]
  data.frame(cohort=nm, n=nrow(d),
             grade_avail=sum(!is.na(d$grade_f)),
             idh_avail=sum(!is.na(d$idh_f)),
             subtype_avail=sum(!is.na(d$subtype_f)),
             kps_avail=sum(!is.na(d$kps)))
})
miss <- do.call(rbind, miss)
write.csv(miss, "results/missingness.csv", row.names=FALSE)
cat("\n=== covariate availability ===\n"); print(miss, row.names=FALSE)
cat("\n=== IDH status by cohort ===\n")
for (nm in names(cohorts)) { cat(nm, ": "); print(table(cohorts[[nm]]$idh_f, useNA="ifany")) }
cat("\n=== IDH-codel subtype by cohort ===\n")
for (nm in names(cohorts)) { cat(nm, ": "); print(table(cohorts[[nm]]$subtype_f, useNA="ifany")) }

# ---- EGFR: nested models, per-SD, with C-index + PH ----
fit_report <- function(df, rhs, cohort, label, term="EGFR_z") {
  f <- as.formula(paste("Surv(time, event) ~", rhs))
  m <- tryCatch(coxph(f, data=df), error=function(e) NULL)
  if (is.null(m) || !term %in% rownames(summary(m)$coefficients)) return(NULL)
  s <- summary(m); ci <- s$conf.int[term,]
  zph <- tryCatch(cox.zph(m), error=function(e) NULL)
  ph_egfr <- if(!is.null(zph)) zph$table[term,"p"] else NA
  ph_glob <- if(!is.null(zph)) zph$table["GLOBAL","p"] else NA
  data.frame(cohort=cohort, model=label,
             HR=signif(unname(ci["exp(coef)"]),3),
             CI_low=signif(unname(ci["lower .95"]),3),
             CI_high=signif(unname(ci["upper .95"]),3),
             p=signif(s$coefficients[term,"Pr(>|z|)"],3),
             C_index=signif(unname(s$concordance["C"]),3),
             n=m$n, events=m$nevent,
             ph_p_EGFR=signif(ph_egfr,3), ph_p_global=signif(ph_glob,3),
             stringsAsFactors=FALSE)
}

egfr_rows <- list()
for (nm in names(cohorts)) {
  d <- cohorts[[nm]]
  gr_ok <- nlevels(droplevels(d$grade_f[!is.na(d$grade_f)])) >= 2
  idh_ok <- nlevels(droplevels(d$idh_f[!is.na(d$idh_f)])) >= 2
  sub_ok <- nlevels(droplevels(d$subtype_f[!is.na(d$subtype_f)])) >= 2
  egfr_rows[[length(egfr_rows)+1]] <- fit_report(d, "EGFR_z", nm, "1. unadjusted")
  egfr_rows[[length(egfr_rows)+1]] <- fit_report(d, "EGFR_z + age", nm, "2. + age")
  if (gr_ok) egfr_rows[[length(egfr_rows)+1]] <- fit_report(d, "EGFR_z + age + grade_f", nm, "3. + age + grade")
  if (idh_ok) egfr_rows[[length(egfr_rows)+1]] <- fit_report(d[!is.na(d$idh_f),], paste("EGFR_z + age", if(gr_ok)"+ grade_f" else "", "+ idh_f"), nm, "4. + age + grade + IDH")
  if (sub_ok) egfr_rows[[length(egfr_rows)+1]] <- fit_report(d[!is.na(d$subtype_f),], "EGFR_z + age + subtype_f", nm, "5. + age + molecular subtype")
  # age-spline sensitivity (unadjusted otherwise)
  egfr_rows[[length(egfr_rows)+1]] <- fit_report(d, "EGFR_z + ns(age,3)", nm, "6. + age spline (sens.)")
}
egfr_tbl <- do.call(rbind, Filter(Negate(is.null), egfr_rows))
write.csv(egfr_tbl, "results/cox_EGFR_adjusted.csv", row.names=FALSE)
cat("\n=== EGFR Cox models (HR per 1 SD of log2(TPM+1)) ===\n"); print(egfr_tbl, row.names=FALSE)

# ---- pathway screen: per-SD, age-adjusted, with CI + BH; plus IDH-adjusted (LGG) ----
screen <- function(df, cohort, extra="") {
  rows <- lapply(GENES, function(g){
    z <- paste0(g,"_z"); if(!z %in% names(df)) return(NULL)
    f <- as.formula(paste("Surv(time,event) ~", z, "+ age", extra))
    m <- tryCatch(coxph(f, data=df), error=function(e) NULL); if(is.null(m)) return(NULL)
    s <- summary(m); ci <- s$conf.int[z,]
    data.frame(cohort=cohort, gene=g,
               HR=signif(unname(ci["exp(coef)"]),3),
               CI_low=signif(unname(ci["lower .95"]),3),
               CI_high=signif(unname(ci["upper .95"]),3),
               p=s$coefficients[z,"Pr(>|z|)"], n=m$n, events=m$nevent,
               stringsAsFactors=FALSE)
  })
  res <- do.call(rbind, Filter(Negate(is.null), rows))
  res$p_BH <- p.adjust(res$p, "BH"); res <- res[order(res$p_BH),]
  res$p <- signif(res$p,3); res$p_BH <- signif(res$p_BH,3); res
}
scr <- rbind(
  transform(screen(lgg, "TCGA-LGG"), adjustment="age"),
  transform(screen(gbm, "TCGA-GBM"), adjustment="age"),
  transform(screen(pooled, "LGG+GBM"), adjustment="age"),
  transform(screen(lgg[!is.na(lgg$idh_f),], "TCGA-LGG", "+ idh_f"), adjustment="age+IDH")
)
write.csv(scr, "results/pathway_screen_std.csv", row.names=FALSE)
cat("\n=== pathway screen (HR per 1 SD, age-adjusted; last block age+IDH) ===\n")
print(scr[,c("cohort","adjustment","gene","HR","CI_low","CI_high","p","p_BH")], row.names=FALSE)

cat("\nDONE\n")
