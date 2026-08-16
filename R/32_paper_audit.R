# 32_paper_audit.R -- machine audit of the numbers added in revision.
#
# R/21 and R/22 verify the ORIGINAL manuscript's cohort sizes and headline models.
# Nothing verifies the numbers added in §6.4, §6.9, Table 5 and Table 6, which were
# transcribed by hand from the outputs of R/23-R/31. This script re-reads those
# result files and asserts every transcribed figure against them, so a stale edit
# fails loudly instead of silently disagreeing with the code.
#
# Run from project root:  Rscript R/32_paper_audit.R   -> results/paper_audit.txt
# Exits non-zero-equivalent by reporting "N FAILED" for the run_all.R gate to catch.

suppressPackageStartupMessages(library(dplyr))
out <- file("results/paper_audit.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }
PASS <- 0; FAIL <- 0
chk <- function(label, got, paper, tol = 0.0015) {
  ok <- !is.na(got) && !is.na(paper) && abs(got - paper) <= tol
  if (ok) PASS <<- PASS + 1 else FAIL <<- FAIL + 1
  say("  [%s] %-52s got=%9.4f  paper=%9.4f\n", ifelse(ok, "OK  ", "FAIL"), label, got, paper)
}

## ---------------- Table 5: matched null on Table 4's gene universe ----------------
# Table 5 reports R/34 (Table 4's universe, gene-level paired bootstrap). R/33's
# 2,000-gene run is the robustness check and is audited loosely at the end.
say("=== Table 5: matched-null calibration (from results/null_full_universe.csv) ===\n")
G <- read.csv("results/null_full_universe.csv")
gg <- function(set, met, col) G[[col]][G$setting == set & G$metric == met]
T5 <- list(
  list("glioma: TCGA->CGGA-693","prod",0.900,0.824,-0.075,-0.093,-0.058,5495),
  list("glioma: TCGA->CGGA-693","absd",0.530,0.669,+0.138,+0.102,+0.176,5495),
  list("glioma: TCGA->CGGA-693","eff", 0.820,0.785,-0.035,-0.056,-0.013,5495),
  list("glioma: TCGA->CGGA-325","prod",0.893,0.937,+0.044,+0.032,+0.055,6037),
  list("glioma: TCGA->CGGA-325","absd",0.611,0.796,+0.185,+0.159,+0.210,6037),
  list("glioma: TCGA->CGGA-325","eff", 0.844,0.767,-0.077,-0.093,-0.060,6037),
  list("glioma: TCGA->array-301","prod",0.793,0.924,+0.132,+0.120,+0.144,6381),
  list("glioma: TCGA->array-301","absd",0.613,0.760,+0.148,+0.129,+0.167,6381),
  list("glioma: TCGA->array-301","eff", 0.768,0.733,-0.036,-0.051,-0.020,6381),
  list("breast: METABRIC->BRCA","prod",0.716,0.566,-0.150,-0.207,-0.094,1021),
  list("breast: METABRIC->BRCA","absd",0.544,0.487,-0.057,-0.118,+0.004,1021),
  list("breast: METABRIC->BRCA","eff", 0.634,0.566,-0.068,-0.125,-0.011,1021))
for (r in T5) {
  lab <- sprintf("%s %s", sub(".*: ","",r[[1]]), r[[2]])
  chk(paste(lab,"null"), gg(r[[1]],r[[2]],"auc_null"), r[[3]])
  chk(paste(lab,"real"), gg(r[[1]],r[[2]],"auc_real"), r[[4]])
  chk(paste(lab,"diff"), gg(r[[1]],r[[2]],"diff"),     r[[5]])
  chk(paste(lab,"CI lo"), gg(r[[1]],r[[2]],"lo"), r[[6]], 0.012)
  chk(paste(lab,"CI hi"), gg(r[[1]],r[[2]],"hi"), r[[7]], 0.012)
  chk(paste(lab,"n genes"), gg(r[[1]],r[[2]],"n_genes"), r[[8]], 0)
}
say("\n  derived claims:\n")
gl <- G[G$metric=="absd" & grepl("glioma", G$setting), ]
gp <- G[G$metric=="prod" & grepl("glioma", G$setting), ]
ge <- G[G$metric=="eff" & is.finite(G$diff), ]
chk("difference form: glioma pairs with CI excluding 0 (paper: 3 of 3)", sum(gl$lo > 0), 3, 0)
chk("product form: glioma pairs with CI excluding 0 (paper: 2 of 3)", sum(gp$lo > 0), 2, 0)
chk("product form inverted in CGGA-693 (diff<0)", as.numeric(gg("glioma: TCGA->CGGA-693","prod","diff") < 0), 1, 0)
chk("effect-size control negative in all evaluable (paper: 4)", sum(ge$diff < 0), 4, 0)
chk("product null AUC, glioma min (paper 0.79)", min(gp$auc_null), 0.793, 0.002)
chk("product null AUC, glioma max (paper 0.90)", max(gp$auc_null), 0.900, 0.002)
chk("difference null AUC, glioma min (paper 0.53)", min(gl$auc_null), 0.530, 0.002)
chk("difference null AUC, glioma max (paper 0.61)", max(gl$auc_null), 0.613, 0.002)
chk("BRCA->METABRIC prognostic genes at 8,000 universe (paper: none)",
    sum(G$n_genes[G$setting=="breast: BRCA->METABRIC"], na.rm=TRUE), 0, 0)
# robustness: the 2,000-gene run (R/33)
say("\n  robustness check vs the 2,000-gene run (R/33):\n")
S33 <- read.csv("results/gene_level_inference.csv")
g33 <- S33[S33$metric=="absd" & grepl("glioma", S33$setting), ]
chk("2,000-gene run: glioma difference-form CIs excluding 0 (paper: 2 of 3)", sum(g33$lo > 0), 2, 0)
chk("2,000-gene CGGA-693 difference diff (paper +0.051)",
    S33$diff[S33$setting=="glioma: TCGA->CGGA-693" & S33$metric=="absd"], 0.051, 0.002)
chk("2,000-gene CGGA-693 difference CI lo (paper -0.019)",
    S33$lo[S33$setting=="glioma: TCGA->CGGA-693" & S33$metric=="absd"], -0.019, 0.012)
chk("2,000-gene BRCA->METABRIC gene count (paper 26)",
    S33$n_genes[S33$setting=="breast: BRCA->METABRIC"][1], 26, 0)

## ---------------- Table 6: corruption modes ----------------
say("\n=== Table 6: corruption modes (from results/corruption_modes.csv) ===\n")
cm <- read.csv("results/corruption_modes.csv") %>% group_by(mode, dose) %>%
  summarise(auc = mean(auc_detect), rc = mean(replic_corrupt), rl = mean(replic_clean),
            vf = mean(var_frac_corrupted), vpre = mean(var_frac_uncorrupted), .groups = "drop")
gv <- function(m, ds, col) cm[[col]][cm$mode == m & abs(cm$dose - ds) < 1e-9]
T6 <- list(
  list("compress",0.25,0.848,0.43,0.93,0.046), list("compress",0.50,0.730,0.80,0.95,0.136),
  list("compress",0.75,0.639,0.91,0.96,0.236), list("permute",0.25,0.630,0.92,0.97,0.197),
  list("permute",0.50,0.728,0.68,0.95,0.098),  list("permute",0.75,0.858,0.24,0.96,0.035),
  list("noise",0.50,0.567,0.95,0.96,0.268),    list("noise",2.00,0.779,0.54,0.93,0.079),
  list("contaminate",0.25,0.543,0.89,0.93,0.308), list("contaminate",0.75,0.734,0.49,0.94,0.329),
  list("floor",0.25,0.498,0.94,0.95,0.341),    list("floor",0.75,0.599,0.82,0.95,0.184))
for (r in T6) {
  chk(sprintf("%s dose=%.2f AUC", r[[1]], r[[2]]), gv(r[[1]], r[[2]], "auc"), r[[3]])
  chk(sprintf("%s dose=%.2f repl corrupt", r[[1]], r[[2]]), gv(r[[1]], r[[2]], "rc"), r[[4]], 0.006)
  chk(sprintf("%s dose=%.2f repl clean",   r[[1]], r[[2]]), gv(r[[1]], r[[2]], "rl"), r[[5]], 0.006)
  chk(sprintf("%s dose=%.2f var share",    r[[1]], r[[2]]), gv(r[[1]], r[[2]], "vf"), r[[6]], 0.006)
}
chk("uncorrupted variance-share baseline (paper: 0.339)", mean(cm$vpre), 0.339, 0.002)
chk("compress lambda=0.5 variance RATIO (paper: 0.40)", gv("compress",0.50,"vf")/mean(cm$vpre), 0.401, 0.004)
chk("flooring max AUC across doses (paper: <=0.60)", max(cm$auc[cm$mode=="floor"]), 0.599, 0.002)

## ---------------- §6.4 composition / variance decomposition ----------------
say("\n=== §6.4: composition and variance decomposition ===\n")
cc <- read.csv("results/cgga_composition.csv")
idhf <- cc[cc$analysis == "idh_mutant_fraction", ]
for (k in c("TCGA","CGGA-693","CGGA-325","CGGA-301")) {
  r <- idhf[idhf$cohort == k, ]
  paper4 <- c(TCGA=0.08, `CGGA-693`=0.18, `CGGA-325`=0.13, `CGGA-301`=0.13)[[k]]
  chk(sprintf("%s grade-IV IDH-mutant fraction", k), r$grade4, paper4, 0.005)
}
wi <- cc[cc$analysis == "r_egfr_grade_within_idh", ]
for (k in c("CGGA-693","CGGA-325","CGGA-301")) {
  paperv <- c(`CGGA-693`=0.163, `CGGA-325`=0.065, `CGGA-301`=0.079)[[k]]
  chk(sprintf("%s within-IDH-wt r(EGFR,grade)", k), wi$r_idhwt[wi$cohort == k], paperv, 0.002)
}
vd <- read.csv("results/egfr_variance_decomposition.csv")
for (k in c("TCGA","CGGA-693","CGGA-325","CGGA-301")) {
  paperv <- c(TCGA=0.106, `CGGA-693`=0.045, `CGGA-325`=0.030, `CGGA-301`=0.024)[[k]]
  chk(sprintf("%s EGFR variance share", k), vd$frac_between[vd$cohort == k], paperv, 0.002)
  paperw <- c(TCGA=1.88, `CGGA-693`=1.84, `CGGA-325`=1.84, `CGGA-301`=1.84)[[k]]
  chk(sprintf("%s EGFR within-stratum SD", k), vd$sd_within[vd$cohort == k], paperw, 0.006)
}
td <- read.csv("results/transport_decomposition.csv")
chk("composition share of discordance, min (paper 8%)",  min(td$frac_composition), 0.08, 0.006)
chk("composition share of discordance, max (paper 25%)", max(td$frac_composition), 0.25, 0.006)

## ---------------- §6.9 effect-size baseline ----------------
say("\n=== §6.9: effect-size baseline (results/effectsize_baseline.csv) ===\n")
eb <- read.csv("results/effectsize_baseline.csv")
u <- eb[eb$criterion == "unadjusted", ]
chk("effect-size baseline, unadjusted min (paper 0.74)", min(u$auc_absbT), 0.738, 0.002)
chk("effect-size baseline, unadjusted max (paper 0.78)", max(u$auc_absbT), 0.783, 0.002)
ai <- read.csv("results/anchor_increment.csv")
chk("anchor increment over effect size, min (paper +0.035)", min(ai$gain), 0.035, 0.002)
chk("anchor increment over effect size, max (paper +0.199)", max(ai$gain), 0.199, 0.002)

## ---------------- controlled composition sweep (§6.9 text) ----------------
say("\n=== §6.9: composition sweep (results/controlled_shift.csv) ===\n")
cs <- read.csv("results/controlled_shift.csv") %>% filter(frac == 0) %>%
  group_by(delta) %>% summarise(auc = mean(auc_replication), rep = mean(replic_rate),
                                cpl = mean(abs(coupling_A - coupling_B)), .groups = "drop")
chk("delta=0 AUC (paper 0.943)",   cs$auc[cs$delta == 0],   0.943)
chk("delta=0.4 AUC (paper 0.609)", cs$auc[cs$delta == 0.4], 0.609)
chk("delta=0 replication (paper 0.957)",   cs$rep[cs$delta == 0],   0.957)
chk("delta=0.4 replication (paper 0.450)", cs$rep[cs$delta == 0.4], 0.450)
chk("delta=0 coupling diff (paper 0.045)",   cs$cpl[cs$delta == 0],   0.045)
chk("delta=0.4 coupling diff (paper 0.796)", cs$cpl[cs$delta == 0.4], 0.796)

## ---------------- §6.10 breast ----------------
say("\n=== §6.10: breast (results/breast_verify.csv) ===\n")
bv <- read.csv("results/breast_verify.csv")
chk("forward ER AUC (paper 0.60)",  bv$auc_er[1], 0.598, 0.002)
chk("reverse ER AUC (paper 0.56)",  bv$auc_er[2], 0.556, 0.002)
chk("max basal AUC (paper 0.54)",   max(bv$auc_basal), 0.542, 0.002)
chk("ER concordance r (paper 0.84)", bv$anchor_concordance_er[1], 0.839, 0.002)
chk("reverse effect-size baseline (paper 0.57)", bv$auc_effsize[2], 0.574, 0.002)


## ---------------- cluster bootstrap (R/35), reported in §6.9 ----------------
say("\n=== §6.9: co-expression cluster bootstrap (results/cluster_bootstrap.csv) ===\n")
CB <- read.csv("results/cluster_bootstrap.csv")
chk("median interval inflation (paper 2.2x)", median(CB$inflation), 2.23, 0.05)
chk("min inflation (paper 1.6x)", min(CB$inflation), 1.63, 0.05)
chk("max inflation (paper 2.7x)", max(CB$inflation), 2.69, 0.05)
chk("difference form: glioma pairs still conclusive (paper 3)",
    sum(CB$clust_excludes0[CB$metric=="absd"]), 3, 0)
chk("product form: glioma pairs still conclusive (paper 3)",
    sum(CB$clust_excludes0[CB$metric=="prod"]), 3, 0)
cb693 <- CB[CB$setting=="glioma: TCGA->CGGA-693", ]
chk("CGGA-693 difference cluster CI lo (paper +0.076)", cb693$clust_lo[cb693$metric=="absd"], 0.076, 0.012)
chk("CGGA-693 difference cluster CI hi (paper +0.200)", cb693$clust_hi[cb693$metric=="absd"], 0.200, 0.012)
chk("CGGA-693 product cluster CI lo (paper -0.115)", cb693$clust_lo[cb693$metric=="prod"], -0.115, 0.012)
chk("CGGA-693 product cluster CI hi (paper -0.033)", cb693$clust_hi[cb693$metric=="prod"], -0.033, 0.012)

## ---------------- threshold / anchor sensitivity (R/36), reported in §6.9 ----------------
say("\n=== §6.9: analyst-choice grid (results/threshold_sensitivity.csv) ===\n")
TS <- read.csv("results/threshold_sensitivity.csv")
TS <- TS[is.finite(TS$diff), ]
exp_sign <- function(set, met) if (met=="prod" && grepl("693", set)) -1 else 1
TS$agree <- mapply(function(s,m,d) sign(d)==exp_sign(s,m), TS$setting, TS$metric, TS$diff)
chk("grid cells evaluated (paper 108)", nrow(TS), 108, 0)
chk("cells agreeing with reported sign (paper 93)", sum(TS$agree), 93, 0)
chk("grade-anchor cells agreeing (paper 54 of 54)", sum(TS$agree[TS$anchor=="grade"]), 54, 0)
chk("grade-anchor cells total (paper 54)", sum(TS$anchor=="grade"), 54, 0)
chk("IDH-anchor cells agreeing (paper 39 of 54)", sum(TS$agree[TS$anchor=="idh"]), 39, 0)
chk("CGGA-325 product/IDH cells agreeing (paper 0 of 9)",
    sum(TS$agree[TS$anchor=="idh" & TS$metric=="prod" & grepl("325", TS$setting)]), 0, 0)

## ---------------- split-variance combination (R/38), reported in §6.9 ----------------
say("\n=== §6.9: split + gene variance combination (results/split_variance.csv) ===\n")
SV <- read.csv("results/split_variance.csv")
gs <- function(set, met, col) SV[[col]][SV$setting == set & SV$metric == met]
SVP <- list(list("TCGA->CGGA-693","prod",-0.097,-0.183,-0.011), list("TCGA->CGGA-693","absd",+0.063,-0.036,+0.162),
            list("TCGA->CGGA-325","prod",+0.077,-0.015,+0.169), list("TCGA->CGGA-325","absd",+0.198,+0.091,+0.306),
            list("TCGA->array-301","prod",+0.068,-0.019,+0.155), list("TCGA->array-301","absd",+0.170,-0.008,+0.348))
for (r in SVP) {
  lab <- sprintf("%s %s combined", sub("TCGA->","",r[[1]]), r[[2]])
  chk(paste(lab,"est"), gs(r[[1]],r[[2]],"est"), r[[3]], 0.002)
  chk(paste(lab,"lo"),  gs(r[[1]],r[[2]],"lo"),  r[[4]], 0.003)
  chk(paste(lab,"hi"),  gs(r[[1]],r[[2]],"hi"),  r[[5]], 0.003)
}
sh <- SV$sd_between^2 / SV$sd_total^2
chk("between-split share of variance, min (paper 77%)", min(sh), 0.765, 0.01)
chk("between-split share of variance, max (paper 84%)", max(sh), 0.844, 0.01)
chk("interval inflation min (paper 3.1x)", min(SV$inflation_vs_within), 3.06, 0.06)
chk("interval inflation max (paper 8.2x)", max(SV$inflation_vs_within), 8.18, 0.06)
ab <- SV[SV$metric=="absd", ]
chk("difference form: glioma pairs excluding 0 after combining (paper 1 of 3)", sum(ab$excludes0), 1, 0)
chk("difference form: all point estimates positive (paper 3 of 3)", sum(ab$est > 0), 3, 0)
chk("split-level difference-form estimate range, min (paper +0.030)", min(ab$est_min), 0.030, 0.004)
chk("split-level difference-form estimate range, max (paper +0.283)", max(ab$est_max), 0.283, 0.004)

say("\n=== %d checks: %d passed, %d FAILED ===\n", PASS + FAIL, PASS, FAIL)
close(out)
cat(sprintf("\n%d passed, %d failed -> results/paper_audit.txt\n", PASS, FAIL))
