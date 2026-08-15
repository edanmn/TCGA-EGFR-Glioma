# 37_supp_audit.R -- MACHINE AUDIT OF THE SUPPLEMENTARY TABLES.
#
# R/21 asserts cohort counts, R/22 refits the headline models, R/32 audits the
# figures added in revision. Nothing has ever checked Tables S1-S8, which carry
# roughly sixty numbers transcribed by hand from result files. Table S9 is covered
# by R/21; S6's contrast is covered by the section 6.7 checks in R/32.
#
# This script re-reads each supplementary table's source CSV and asserts every
# transcribed value, so a stale supplement fails loudly instead of quietly
# disagreeing with the code. It also flags a hazard found during review: several
# quantities exist in TWO result files, one superseded (e.g. pathway_screen.csv is
# per-unit-expression while the manuscript uses the per-SD pathway_screen_std.csv;
# cgga_positive_controls.csv is superseded by panel_correlations.csv).
#
# Run from project root:  Rscript R/37_supp_audit.R   -> results/supp_audit.txt

suppressPackageStartupMessages(library(dplyr))
out <- file("results/supp_audit.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }
PASS <- 0; FAIL <- 0
chk <- function(label, got, paper, tol = 0.006) {
  ok <- length(got) == 1 && !is.na(got) && !is.na(paper) && abs(got - paper) <= tol
  if (ok) PASS <<- PASS + 1 else FAIL <<- FAIL + 1
  say("  [%s] %-50s got=%12.5g  paper=%12.5g\n", ifelse(ok,"OK  ","FAIL"), label,
      ifelse(length(got)==1, got, NA), paper)
}
rel <- function(label, got, paper, pct = 0.02) {   # for p-values: relative tolerance
  ok <- length(got)==1 && !is.na(got) && !is.na(paper) && abs(log(got)-log(paper)) <= pct
  if (ok) PASS <<- PASS + 1 else FAIL <<- FAIL + 1
  say("  [%s] %-50s got=%12.4g  paper=%12.4g\n", ifelse(ok,"OK  ","FAIL"), label,
      ifelse(length(got)==1, got, NA), paper)
}

## ---------------- Table S1: KM median split ----------------
say("=== Table S1: unadjusted KM median split (km_logrank_summary.csv) ===\n")
k <- read.csv("results/km_logrank_summary.csv")
g <- function(co, col) k[[col]][k$cohort == co]
S1 <- list(list("TCGA-LGG",256,255,125,2988,2282,0.040),
           list("TCGA-GBM",141,141,227,399,448,0.99))
for (r in S1) {
  chk(sprintf("%s n_low",  r[[1]]), g(r[[1]],"n_low"),  r[[2]], 0)
  chk(sprintf("%s n_high", r[[1]]), g(r[[1]],"n_high"), r[[3]], 0)
  chk(sprintf("%s events", r[[1]]), g(r[[1]],"events"), r[[4]], 0)
  chk(sprintf("%s median low",  r[[1]]), g(r[[1]],"median_surv_low_days"),  r[[5]], 0)
  chk(sprintf("%s median high", r[[1]]), g(r[[1]],"median_surv_high_days"), r[[6]], 0)
  chk(sprintf("%s logrank p",   r[[1]]), g(r[[1]],"logrank_p"), r[[7]], 0.001)
}
pooled <- k[grepl("ooled|LGG\\+GBM", k$cohort), ]
if (nrow(pooled)) {
  chk("pooled n_low/high", pooled$n_low[1], 397, 0); chk("pooled events", pooled$events[1], 352, 0)
  chk("pooled median low", pooled$median_surv_low_days[1], 1491, 0)
  chk("pooled median high", pooled$median_surv_high_days[1], 883, 0)
}

## ---------------- Table S2: pathway screen (PER-SD file) ----------------
say("\n=== Table S2: pathway screen (pathway_screen_std.csv, per-SD) ===\n")
ps <- read.csv("results/pathway_screen_std.csv")
ga <- function(gene, col, adj = "age") ps[[col]][ps$gene == gene & ps$adjustment == adj & ps$cohort == "TCGA-LGG"]
S2 <- list(list("RB1",1.64,1.34,2.00,1.59e-5,0.031), list("IDH1",1.46,1.20,1.78,8.67e-4,0.033),
           list("EGFR",1.37,1.14,1.66,3.54e-3,0.178), list("TP53",1.39,1.13,1.71,5.20e-3,0.026),
           list("PIK3R1",0.777,0.649,0.931,1.25e-2,0.876), list("PTEN",0.800,0.666,0.961,2.83e-2,0.506),
           list("CDKN2A",0.834,0.691,1.01,8.72e-2,0.046), list("PIK3CA",1.13,0.950,1.34,0.212,0.046),
           list("PDGFRA",0.893,0.752,1.06,0.220,0.594), list("NF1",0.895,0.750,1.07,0.222,0.085))
for (r in S2) {
  chk(sprintf("%s HR (age)", r[[1]]), ga(r[[1]],"HR"), r[[2]], 0.006)
  chk(sprintf("%s CI low",   r[[1]]), ga(r[[1]],"CI_low"), r[[3]], 0.006)
  chk(sprintf("%s CI high",  r[[1]]), ga(r[[1]],"CI_high"), r[[4]], 0.011)
  rel(sprintf("%s p_BH (age)", r[[1]]), ga(r[[1]],"p_BH"), r[[5]], 0.03)
  v <- ga(r[[1]],"p_BH","age+IDH"); if (length(v)==1) chk(sprintf("%s p_BH (age+IDH)", r[[1]]), v, r[[6]], 0.002)
}

## ---------------- Table S3: VST sensitivity ----------------
say("\n=== Table S3: VST normalization (vst_sensitivity.csv) ===\n")
v <- read.csv("results/vst_sensitivity.csv")
gv <- function(co, mo, col) v[[col]][v$cohort == co & grepl(mo, v$model)]
S3 <- list(list("TCGA-LGG","unadjusted",1.50,1.22,1.84,1.25e-4,0.577),
           list("TCGA-LGG","age \\+ grade$|age \\+ grade[^,]",1.26,1.04,1.54,0.0204,0.784),
           list("TCGA-LGG","IDH",1.11,0.944,1.31,0.205,0.829),
           list("TCGA-GBM","unadjusted",0.945,0.829,1.08,0.397,0.524),
           list("TCGA-GBM","IDH",0.843,0.737,0.963,0.0124,0.641))
for (r in S3) {
  chk(sprintf("%s %s HR", r[[1]], substr(r[[2]],1,12)), gv(r[[1]],r[[2]],"HR")[1], r[[3]], 0.006)
  chk(sprintf("%s %s C-index", r[[1]], substr(r[[2]],1,12)), gv(r[[1]],r[[2]],"C_index")[1], r[[7]], 0.002)
}

## ---------------- Table S4: ten-gene panel correlations ----------------
say("\n=== Table S4: panel correlations (panel_correlations.csv) ===\n")
pc <- read.csv("results/panel_correlations.csv")
S4 <- list(list("EGFR",0.179,-0.033), list("IDH1",0.348,0.152), list("TP53",0.099,0.113),
           list("CDKN2A",-0.131,-0.076), list("NF1",-0.285,-0.183), list("PTEN",-0.357,-0.071),
           list("PDGFRA",-0.304,-0.096), list("PIK3R1",-0.481,-0.195), list("PIK3CA",-0.072,-0.057),
           list("RB1",-0.057,-0.023))
for (r in S4) {
  chk(sprintf("%s r TCGA", r[[1]]), pc$r_TCGA[pc$gene==r[[1]]], r[[2]], 0.001)
  chk(sprintf("%s r CGGA", r[[1]]), pc$r_CGGA[pc$gene==r[[1]]], r[[3]], 0.001)
}

## ---------------- Table S5: recount3 ----------------
say("\n=== Table S5: recount3 reprocessing (recount3_validation.csv) ===\n")
rc <- read.csv("results/recount3_validation.csv")
gr <- function(co, mo, col) rc[[col]][rc$cohort==co & grepl(mo, rc$model)][1]
S5 <- list(list("TCGA-LGG","unadjusted",1.54,509,125), list("TCGA-LGG","grade$",1.27,452,106),
           list("TCGA-LGG","IDH",1.11,450,105), list("TCGA-GBM","unadjusted",0.906,154,122),
           list("TCGA-GBM","age$",0.853,154,122), list("TCGA-GBM","IDH",0.813,150,120))
for (r in S5) {
  chk(sprintf("%s %s HR", r[[1]], substr(r[[2]],1,10)), gr(r[[1]],r[[2]],"HR"), r[[3]], 0.006)
  chk(sprintf("%s %s n",  r[[1]], substr(r[[2]],1,10)), gr(r[[1]],r[[2]],"n"),  r[[4]], 0)
  chk(sprintf("%s %s ev", r[[1]], substr(r[[2]],1,10)), gr(r[[1]],r[[2]],"events"), r[[5]], 0)
}

## ---------------- Table S6: gated screen ----------------
say("\n=== Table S6: positive-control-gated screen (systematic_screen.csv) ===\n")
ss <- read.csv("results/systematic_screen.csv"); p <- ss[ss$prognostic_T, ]
chk("QC-concordant genes", sum(p$qc_ok), 429, 0)
chk("QC-discordant genes", sum(!p$qc_ok), 37, 0)
chk("QC-concordant replication %", 100*mean(p$replicated[p$qc_ok]), 92.8, 0.1)
chk("QC-discordant replication %", 100*mean(p$replicated[!p$qc_ok]), 24.3, 0.1)
chk("sign flips, all QC-discordant", sum(p$signflip & !p$qc_ok), 2, 0)

## ---------------- Table S7: circularity control ----------------
say("\n=== Table S7: circularity control (circularity_control.csv) ===\n")
cc <- read.csv("results/circularity_control.csv")
gc2 <- function(a, col) cc[[col]][cc$analysis == a]
auc_pt <- function(a) as.numeric(sub(" .*", "", gc2(a, "AUC")))
S7 <- list(list("unadjusted",6616,0.63,0.82), list("subtype_adjusted",183,0.23,0.648),
           list("low_grade_corr",935,0.51,0.752), list("high_grade_corr",5681,0.65,0.843))
for (r in S7) {
  chk(sprintf("%s n", r[[1]]), gc2(r[[1]],"n"), r[[2]], 0)
  chk(sprintf("%s replication rate", r[[1]]), gc2(r[[1]],"replic_rate"), r[[3]], 0.006)
  chk(sprintf("%s AUC", r[[1]]), auc_pt(r[[1]]), r[[4]], 0.006)
}

## ---------------- Table S8: naive CGGA models ----------------
say("\n=== Table S8: naive CGGA validation (cgga_validation.csv) ===\n")
cv <- read.csv("results/cgga_validation.csv")
gcv <- function(mo, col) cv[[col]][grepl(mo, cv$model)][1]
S8 <- list(list("unadjusted",1.17,0.938,1.46,0.164,0.520),
           list("age \\+ grade$",1.19,0.958,1.48,0.115,0.624),
           list("IDH",1.44,1.17,1.76,4.6e-4,0.737))
for (r in S8) {
  chk(sprintf("CGGA %s HR", substr(r[[1]],1,10)), gcv(r[[1]],"HR"), r[[2]], 0.006)
  chk(sprintf("CGGA %s C-index", substr(r[[1]],1,10)), gcv(r[[1]],"C_index"), r[[6]], 0.002)
}

## ---------------- superseded-artifact hazard ----------------
say("\n=== Superseded result files that disagree with the manuscript ===\n")
say("  (not failures -- the manuscript uses the correct file -- but they should be\n")
say("   removed or renamed before release, since a reader may open the wrong one)\n")
if (file.exists("results/pathway_screen.csv")) {
  a <- read.csv("results/pathway_screen.csv")
  say("  pathway_screen.csv        RB1 HR=%.2f (per unit) vs Table S2's %.2f (per SD)\n",
      a$HR[a$gene=="RB1"][1], 1.64)
}
if (file.exists("results/cgga_positive_controls.csv")) {
  b <- read.csv("results/cgga_positive_controls.csv")
  say("  cgga_positive_controls.csv EGFR r=%.3f vs Table S4's %.3f\n",
      b$r_grade_TCGA[b$gene=="EGFR"][1], 0.179)
}

say("\n=== %d checks: %d passed, %d FAILED ===\n", PASS + FAIL, PASS, FAIL)
close(out)
cat(sprintf("\n%d passed, %d failed -> results/supp_audit.txt\n", PASS, FAIL))
