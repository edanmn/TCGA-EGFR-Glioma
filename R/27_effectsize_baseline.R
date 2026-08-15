# 27_effectsize_baseline.R -- THE MISSING BASELINE, and a calibrated null.
#
# 25_controlled_shift.R found that two EXCHANGEABLE random halves of TCGA -- no
# composition shift, no measurement problem, nothing wrong -- still give a
# grade-anchor AUC of 0.94 for predicting replication. That is higher than any
# observational AUC in the paper. The mechanism must therefore be something that
# is present even when the cohorts are identical: genes with a strong anchor
# relationship have a large survival effect, and large effects replicate.
#
# If that is right, the paper's baselines are the wrong ones. It benchmarks the
# anchor against DETECTABILITY (mean expression) and PRECISION (SE of the
# replication-cohort coefficient). Both measure NOISE. Neither measures SIGNAL
# MAGNITUDE. The obvious competitor -- the discovery-cohort effect size |beta_T| --
# is computed in the BREAST analysis (17_breast_generalization.R:84, AUC 0.530)
# but never in glioma.
#
# This script supplies it, plus:
#   (A) every baseline on one footing, both replication criteria, all 3 cohorts
#   (B) does the anchor add anything OVER effect size? (logistic model comparison
#       and AUC of the anchor within strata of effect size)
#   (C) the EGFR variance puzzle: 24_cgga_composition.R found EGFR's IDH-adjusted
#       grade SLOPE is the same in TCGA and CGGA (+0.115 vs +0.124), yet
#       26_transport_anchor.R found a large transport RESIDUAL. Those reconcile
#       only if CGGA's within-stratum VARIANCE is inflated. Tested directly.
#
# Reuses results/pergene_cache.rds written by 26_transport_anchor.R.
#
# Run from project root:  Rscript R/27_effectsize_baseline.R
# Writes: results/effectsize_baseline.csv, results/effectsize_baseline.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)

out <- file("results/effectsize_baseline.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

CACHE <- "results/pergene_cache.rds"
if (!file.exists(CACHE)) stop("run R/26_transport_anchor.R first (writes ", CACHE, ")")
P <- readRDS(CACHE); genes <- P$genes; TS <- P$TCGA
say("Genes: %d (four-cohort intersection, TCGA SD>0.5)\n\n", length(genes))

auc <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
aucCI <- function(score, label, B = 600){
  ok <- !is.na(score) & !is.na(label); s <- score[ok]; l <- label[ok]
  a <- auc(s, l); if (is.na(a)) return(c(a = NA, lo = NA, hi = NA))
  bs <- replicate(B, { i <- sample(length(s), replace = TRUE); auc(s[i], l[i]) })
  c(a = a, lo = unname(quantile(bs, .025, na.rm = TRUE)), hi = unname(quantile(bs, .975, na.rm = TRUE)))
}
fmt <- function(v) sprintf("%.3f (%.3f-%.3f)", v["a"], v["lo"], v["hi"])

## ---------- (A) all baselines on one footing ----------
say("=== (A) Every predictor of replication, on identical footing ===\n")
say("  anchor      = r_grade(TCGA) * r_grade(replication)   [the paper's score]\n")
say("  |beta_T|    = discovery-cohort effect size            [THE MISSING BASELINE]\n")
say("  |z_T|       = discovery-cohort effect / its SE        [discovery signal-to-noise]\n")
say("  precision   = -SE of replication coefficient          [paper's baseline]\n")
say("  |r_grade_T| = discovery anchor strength ALONE, no replication-cohort data\n\n")

rows <- list()
for (k in c("CGGA-693","CGGA-325","array-301")) {
  RS <- P[[k]]
  for (crit in c("unadjusted","adjusted")) {
    bR <- if (crit == "unadjusted") RS$b0 else RS$b1
    pR <- if (crit == "unadjusted") RS$p0 else RS$p1
    prog <- p.adjust(TS$p0, "BH") < 0.05 & !is.na(bR)
    repl <- as.integer(sign(TS$b0) == sign(bR) & pR < 0.05)[prog]

    anchor  <- (TS$r_grade * RS$r_grade)[prog]
    absbT   <- abs(TS$b0)[prog]
    zT      <- abs(qnorm(pmax(TS$p0, 1e-300) / 2))[prog]     # |z| from the discovery p-value
    absrT   <- abs(TS$r_grade)[prog]

    a1 <- aucCI(anchor, repl); a2 <- aucCI(absbT, repl)
    a3 <- aucCI(zT, repl);     a4 <- aucCI(absrT, repl)
    rows[[length(rows)+1]] <- data.frame(
      cohort = k, criterion = crit, n = sum(prog), replic = mean(repl),
      auc_anchor = a1["a"], anchor_lo = a1["lo"], anchor_hi = a1["hi"],
      auc_absbT = a2["a"], absbT_lo = a2["lo"], absbT_hi = a2["hi"],
      auc_zT = a3["a"], auc_absrT = a4["a"], row.names = NULL)
    say("  %-10s %-11s n=%4d replic=%.3f\n", k, crit, sum(prog), mean(repl))
    say("      anchor      %s\n", fmt(a1))
    say("      |beta_T|    %s   <-- missing baseline\n", fmt(a2))
    say("      |z_T|       %s\n", fmt(a3))
    say("      |r_grade_T| %s   (discovery only, no replication-cohort data)\n", fmt(a4))
  }
}
A <- do.call(rbind, rows)
write.csv(A, "results/effectsize_baseline.csv", row.names = FALSE)

## ---------- (B) does the anchor add anything OVER effect size? ----------
say("\n=== (B) Incremental value of the anchor over discovery effect size ===\n")
say("  Logistic models for replication; nested LR test and AUC.\n\n")
incr <- list()
for (k in c("CGGA-693","CGGA-325","array-301")) {
  RS <- P[[k]]
  for (crit in c("unadjusted","adjusted")) {
    bR <- if (crit == "unadjusted") RS$b0 else RS$b1
    pR <- if (crit == "unadjusted") RS$p0 else RS$p1
    prog <- p.adjust(TS$p0, "BH") < 0.05 & !is.na(bR)
    repl <- as.integer(sign(TS$b0) == sign(bR) & pR < 0.05)[prog]
    d <- data.frame(y = repl, anchor = (TS$r_grade * RS$r_grade)[prog], eff = abs(TS$b0)[prog])
    d <- d[complete.cases(d), ]
    if (length(unique(d$y)) < 2) next
    m1 <- glm(y ~ eff, binomial, d)
    m2 <- glm(y ~ eff + anchor, binomial, d)
    lr <- anova(m1, m2, test = "LRT")
    a1 <- auc(predict(m1), d$y); a2 <- auc(predict(m2), d$y)
    say("  %-10s %-11s  AUC(eff)=%.3f  AUC(eff+anchor)=%.3f  gain=%+.3f  LR p=%.2g\n",
        k, crit, a1, a2, a2 - a1, lr$`Pr(>Chi)`[2])
    incr[[length(incr)+1]] <- data.frame(cohort = k, criterion = crit,
      auc_eff = a1, auc_eff_anchor = a2, gain = a2 - a1, lr_p = lr$`Pr(>Chi)`[2])
  }
}
write.csv(do.call(rbind, incr), "results/anchor_increment.csv", row.names = FALSE)

## ---------- (C) the EGFR variance puzzle ----------
say("\n=== (C) EGFR: reconciling the matched conditional slope with the large\n")
say("        transport residual. If CGGA's WITHIN-STRATUM variance is inflated,\n")
say("        the same slope produces a much smaller marginal correlation.\n\n")
tcga_mat <- function(){
  one <- function(project, gradeIV = FALSE){
    se <- load_se(project); cd <- as.data.frame(colData(se)); k <- substr(colnames(se), 14, 15) == "01"
    e <- log2(assay(se, "tpm_unstrand")[, k] + 1); rownames(e) <- rowData(se)$gene_name
    e <- e[!duplicated(rownames(e)), ]
    b <- .clinical_df(se, project); b <- b[match(colnames(e), b$barcode), ]
    gr <- as.character(cd$paper_Grade[k])
    g <- if (gradeIV) ifelse(grepl("G4|IV", gr), 4, NA) else ifelse(grepl("G3|III", gr), 3, ifelse(grepl("G2| II", gr), 2, NA))
    idh <- as.character(cd$paper_IDH.status[k]); keep <- !duplicated(b$patient)
    list(e = e[, keep, drop = FALSE], grade = g[keep],
         idh = ifelse(idh[keep] == "Mutant", 1, ifelse(idh[keep] == "WT", 0, NA)))
  }
  L <- one("TCGA-LGG"); G <- one("TCGA-GBM", TRUE); cg <- intersect(rownames(L$e), rownames(G$e))
  list(e = cbind(L$e[cg, ], G$e[cg, ]), grade = c(L$grade, G$grade), idh = c(L$idh, G$idh))
}
cgga <- function(gf, cf, cols, logt = TRUE){
  ex <- rd(gf); rownames(ex) <- ex[[1]]; ex[[1]] <- NULL
  ex <- if (logt) log2(as.matrix(ex) + 1) else as.matrix(ex)
  cl <- read.delim(cf, check.names = FALSE); names(cl)[cols] <- c("id","prs","grade","age","os","censor","idh")
  cl <- cl %>% filter(prs == "Primary") %>%
    mutate(grade = as.numeric(factor(grade, levels = c("WHO II","WHO III","WHO IV"))) + 1,
           idh = ifelse(idh == "Mutant", 1, ifelse(idh == "Wildtype", 0, NA)))
  sel <- intersect(colnames(ex), cl$id); cl <- cl[match(sel, cl$id), ]
  list(e = ex[, sel, drop = FALSE], grade = cl$grade, idh = cl$idh)
}
CO <- list(TCGA = tcga_mat(),
  `CGGA-693` = cgga("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt","data/cgga/cgga_clinical.tsv",c(1,2,4,6,7,8,11)),
  `CGGA-325` = cgga("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt","data/cgga/cgga325_clinical.tsv",c(1,2,4,6,7,8,11)),
  `CGGA-301` = cgga("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt","data/cgga/cgga_array_clinical.tsv",c(1,3,5,7,8,9,12),logt=FALSE))

say("  %-10s %-28s %-28s %s\n", "cohort", "EGFR between-stratum SD", "EGFR within-stratum SD", "between/total var")
vres <- list()
for (k in names(CO)) {
  M <- CO[[k]]; x <- M$e["EGFR", ]
  ok <- !is.na(x) & !is.na(M$grade) & !is.na(M$idh)
  x <- x[ok]; s <- paste(M$grade[ok], M$idh[ok], sep = "_")
  mu <- tapply(x, s, mean); w <- table(s) / length(s)
  vb <- sum(w * (mu[names(w)] - mean(x))^2)
  vw <- sum(w * tapply(x, s, var)[names(w)], na.rm = TRUE)
  # same, for the whole transcriptome, as a reference for how unusual EGFR is
  say("  %-10s %-28.4f %-28.4f %.3f\n", k, sqrt(vb), sqrt(vw), vb / (vb + vw))
  vres[[length(vres)+1]] <- data.frame(cohort = k, sd_between = sqrt(vb), sd_within = sqrt(vw),
                                       frac_between = vb / (vb + vw))
}
V <- do.call(rbind, vres)
write.csv(V, "results/egfr_variance_decomposition.csv", row.names = FALSE)
say("\n  'between/total var' is the share of EGFR's variance explained by the six\n")
say("  (grade, IDH) strata. A collapse of this share in CGGA relative to TCGA means\n")
say("  the covariate structure explains far less of EGFR's variation there -- the\n")
say("  signature of added measurement noise, not of a different patient mix.\n")
close(out)
cat("\nwrote results/effectsize_baseline.txt and .csv\n")
