# 26_transport_anchor.R -- TRANSPORT-ADJUSTED ANCHOR CONCORDANCE.
#
# 24_cgga_composition.R showed that for EGFR the raw anchor discordance between
# TCGA and CGGA is almost entirely explained by a difference in the cohorts'
# joint (grade, IDH) composition, not by anything about EGFR itself. That is a
# statement about one gene. This script turns it into a method and applies it to
# the whole transcriptome.
#
# The raw anchor score r_T * r_R conflates two very different things:
#   (1) COMPOSITION SHIFT -- the cohorts have different (grade, IDH) mixes, so a
#       gene measured identically in both still shows different MARGINAL anchor
#       correlations. Nothing is wrong with the gene.
#   (2) CONDITIONAL DISCORDANCE -- the gene's relationship to the anchor differs
#       WITHIN strata. This is the part that implicates the gene.
#
# Transport step. Grade is constant within a (grade, IDH) stratum, so the
# correlation implied by TCGA's conditional structure under cohort R's stratum
# weights w_R has a closed form (law of total covariance; the within-stratum
# term vanishes because grade does not vary inside a stratum):
#
#   E[X]   = sum_s w_R(s) mu_T(s)
#   Cov    = sum_s w_R(s) (mu_T(s) - E[X]) (g(s) - E[G])
#   Var[X] = sum_s w_R(s) [ sigma^2_T(s) + (mu_T(s) - E[X])^2 ]
#   r_pred = Cov / sqrt(Var[X] Var[G])
#
# r_pred answers: "if cohort R measured this gene exactly as TCGA does, but with
# R's patient mix, what anchor correlation would we see?" Correlation is scale
# invariant, so the differing expression units across platforms cancel.
#
# Two new scores follow, and both are tested against the raw score:
#   RESIDUAL  = -|r_actual_R - r_pred|   (discordance beyond what composition explains)
#   EXPECTED-SHIFT = -|r_pred - r_T|     (how much composition ALONE distorts this
#                                         gene's cross-cohort anchor comparison)
#
# EXPECTED-SHIFT is the practically interesting one: it uses TCGA expression and
# the replication cohort's CLINICAL TABLE ONLY. No replication-cohort expression
# is required, so it can be computed BEFORE the validation data exist.
#
# Run from project root:  Rscript R/26_transport_anchor.R   (~5 min, then cached)
# Writes: results/transport_anchor.csv, results/transport_anchor.txt,
#         results/pergene_cache.rds

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
set.seed(20260815)   # bootstrap CIs below are resampled; seed for reproducibility
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)

out <- file("results/transport_anchor.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

## ---------- cohorts (grade coded 2/3/4 everywhere in this script) ----------
tcga_mat <- function(){
  one <- function(project, gradeIV = FALSE){
    se <- load_se(project); cd <- as.data.frame(colData(se)); k <- substr(colnames(se), 14, 15) == "01"
    e <- log2(assay(se, "tpm_unstrand")[, k] + 1); rownames(e) <- rowData(se)$gene_name
    e <- e[!duplicated(rownames(e)), ]
    b <- .clinical_df(se, project); b <- b[match(colnames(e), b$barcode), ]
    gr <- as.character(cd$paper_Grade[k])
    g <- if (gradeIV) ifelse(grepl("G4|IV", gr), 4, NA) else ifelse(grepl("G3|III", gr), 3, ifelse(grepl("G2| II", gr), 2, NA))
    idh <- as.character(cd$paper_IDH.status[k])
    keep <- !is.na(b$time) & b$time > 0 & !is.na(b$event) & !duplicated(b$patient)
    list(e = e[, keep], time = b$time[keep], event = b$event[keep], age = b$age[keep],
         grade = g[keep], idh = ifelse(idh[keep] == "Mutant", 1, ifelse(idh[keep] == "WT", 0, NA)))
  }
  L <- one("TCGA-LGG"); G <- one("TCGA-GBM", TRUE)
  cg <- intersect(rownames(L$e), rownames(G$e))
  list(e = cbind(L$e[cg, ], G$e[cg, ]),
       clin = data.frame(time = c(L$time, G$time), event = c(L$event, G$event),
                         age = c(L$age, G$age), grade = c(L$grade, G$grade), idh = c(L$idh, G$idh)))
}
cgga_mat <- function(gfile, cfile, cols, logt = TRUE){
  ex <- rd(gfile); rownames(ex) <- ex[[1]]; ex[[1]] <- NULL
  ex <- if (logt) log2(as.matrix(ex) + 1) else as.matrix(ex)
  cl <- read.delim(cfile, check.names = FALSE); names(cl)[cols] <- c("id","prs","grade","age","os","censor","idh")
  cl <- cl %>% filter(prs == "Primary") %>%
    mutate(time = as.numeric(os), event = as.integer(censor), age = as.numeric(age),
           grade = as.numeric(factor(grade, levels = c("WHO II","WHO III","WHO IV"))) + 1,
           idh = ifelse(idh == "Mutant", 1, ifelse(idh == "Wildtype", 0, NA)))
  sel <- intersect(colnames(ex), cl$id); cl <- cl[match(sel, cl$id), ]
  list(e = ex[, sel], clin = cl[, c("time","event","age","grade","idh")])
}

CACHE <- "results/pergene_cache.rds"
if (file.exists(CACHE)) {
  say("loading cached per-gene statistics from %s\n", CACHE); P <- readRDS(CACHE)
} else {
  TC   <- tcga_mat()
  C693 <- cgga_mat("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", "data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
  C325 <- cgga_mat("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", "data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
  CARR <- cgga_mat("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt", "data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt = FALSE)

  genes <- Reduce(intersect, list(rownames(TC$e), rownames(C693$e), rownames(C325$e), rownames(CARR$e)))
  sdT <- apply(TC$e[genes, ], 1, sd); genes <- genes[sdT > 0.5]
  if (length(genes) > 8000) genes <- names(sort(sdT[genes], decreasing = TRUE))[1:8000]

  # per-gene Cox (unadjusted + subtype-adjusted), Pearson anchor r, and the
  # per-stratum mean/variance the transport formula needs
  cohort_stats <- function(M, gs, lab){
    cl <- M$clin; e <- M$e[gs, , drop = FALSE]
    cc <- !is.na(cl$idh) & !is.na(cl$grade)
    S0 <- Surv(cl$time, cl$event); a0 <- cl$age
    S1 <- Surv(cl$time[cc], cl$event[cc]); a1 <- cl$age[cc]
    g1 <- factor(cl$grade[cc]); i1 <- factor(cl$idh[cc])
    cat(sprintf("  %s: n=%d (adjusted n=%d)\n", lab, nrow(cl), sum(cc)))
    b0 <- p0 <- b1 <- p1 <- numeric(length(gs))
    for (i in seq_along(gs)) {
      z <- as.numeric(scale(e[i, ]))
      m0 <- tryCatch(coxph(S0 ~ z + a0), error = function(e) NULL)
      zc <- z[cc]
      m1 <- tryCatch(coxph(S1 ~ zc + a1 + g1 + i1), error = function(e) NULL)
      s0 <- if (is.null(m0)) c(NA,NA) else summary(m0)$coefficients["z",  c("coef","Pr(>|z|)")]
      s1 <- if (is.null(m1)) c(NA,NA) else summary(m1)$coefficients["zc", c("coef","Pr(>|z|)")]
      b0[i] <- s0[1]; p0[i] <- s0[2]; b1[i] <- s1[1]; p1[i] <- s1[2]
    }
    okg <- !is.na(cl$grade)
    r_grade <- suppressWarnings(apply(e, 1, function(x) cor(x[okg], cl$grade[okg])))
    # stratum means / variances, strata = (grade, idh)
    key <- paste(cl$grade, cl$idh, sep = "_")
    kk <- key[cc]
    lev <- sort(unique(kk))
    mu <- sapply(lev, function(L) rowMeans(e[, cc, drop = FALSE][, kk == L, drop = FALSE]))
    v  <- sapply(lev, function(L) apply(e[, cc, drop = FALSE][, kk == L, drop = FALSE], 1, var))
    w  <- sapply(lev, function(L) mean(kk == L))
    list(gene = gs, b0 = b0, p0 = p0, b1 = b1, p1 = p1, r_grade = r_grade,
         mu = mu, v = v, w = w, lev = lev, n = nrow(cl), n_adj = sum(cc))
  }
  cat("computing per-gene statistics (4 cohorts)...\n")
  P <- list(genes = genes,
            TCGA = cohort_stats(TC, genes, "TCGA"),
            `CGGA-693` = cohort_stats(C693, genes, "CGGA-693"),
            `CGGA-325` = cohort_stats(C325, genes, "CGGA-325"),
            `array-301` = cohort_stats(CARR, genes, "array-301"))
  saveRDS(P, CACHE)
  say("computed and cached per-gene statistics for %d genes\n", length(genes))
}
genes <- P$genes
say("Genes: %d (four-cohort intersection, TCGA SD>0.5)\n\n", length(genes))

## ---------- transport: predict r_R from TCGA conditional structure + R weights ----------
transport_r <- function(Tst, Rst){
  lev <- intersect(Tst$lev, Rst$lev)
  w <- Rst$w[match(lev, Rst$lev)]; w <- w / sum(w)
  mu <- Tst$mu[, match(lev, Tst$lev), drop = FALSE]
  v  <- Tst$v [, match(lev, Tst$lev), drop = FALSE]
  g  <- as.numeric(sub("_.*", "", lev))
  Eg <- sum(w * g); Vg <- sum(w * (g - Eg)^2)
  Ex <- as.numeric(mu %*% w)
  cov <- as.numeric(sweep(mu, 1, Ex, "-") %*% (w * (g - Eg)))
  vx  <- as.numeric((v + sweep(mu, 1, Ex, "-")^2) %*% w)
  list(r = cov / sqrt(vx * Vg), weight_kept = sum(Rst$w[match(lev, Rst$lev)]), lev = lev)
}

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

TS <- P$TCGA
rows <- list(); dec <- list()
for (k in c("CGGA-693","CGGA-325","array-301")) {
  RS <- P[[k]]
  tr <- transport_r(TS, RS)
  r_pred <- tr$r; r_act <- RS$r_grade; r_T <- TS$r_grade

  raw       <- r_T * r_act              # current paper score
  residual  <- -abs(r_act - r_pred)     # discordance beyond composition
  expected  <- -abs(r_pred - r_T)       # composition-only distortion (no R expression)

  for (crit in c("unadjusted","adjusted")) {
    bT <- if (crit == "unadjusted") TS$b0 else TS$b0   # discovery always the unadjusted set
    pT <- TS$p0
    bR <- if (crit == "unadjusted") RS$b0 else RS$b1
    pR <- if (crit == "unadjusted") RS$p0 else RS$p1
    prog <- p.adjust(pT, "BH") < 0.05 & !is.na(bR)
    repl <- as.integer(sign(bT) == sign(bR) & pR < 0.05)[prog]
    A <- aucCI(raw[prog], repl); B2 <- aucCI(residual[prog], repl); C <- aucCI(expected[prog], repl)
    rows[[length(rows)+1]] <- data.frame(
      cohort = k, criterion = crit, n = sum(prog), replic_rate = mean(repl),
      auc_raw = A["a"], raw_lo = A["lo"], raw_hi = A["hi"],
      auc_residual = B2["a"], res_lo = B2["lo"], res_hi = B2["hi"],
      auc_expected = C["a"], exp_lo = C["lo"], exp_hi = C["hi"], row.names = NULL)
  }
  # how much of the total anchor discordance is composition?
  tot <- abs(r_act - r_T); comp <- abs(r_pred - r_T); res <- abs(r_act - r_pred)
  dec[[length(dec)+1]] <- data.frame(
    cohort = k, weight_kept = tr$weight_kept, strata = length(tr$lev),
    mean_total_discord = mean(tot, na.rm = TRUE),
    mean_composition   = mean(comp, na.rm = TRUE),
    mean_residual      = mean(res, na.rm = TRUE),
    frac_composition   = mean(comp, na.rm = TRUE) / mean(tot, na.rm = TRUE),
    cor_pred_actual    = cor(r_pred, r_act, use = "complete.obs"),
    cor_T_actual       = cor(r_T,    r_act, use = "complete.obs"),
    egfr_r_T = r_T[genes == "EGFR"], egfr_r_pred = r_pred[genes == "EGFR"],
    egfr_r_act = r_act[genes == "EGFR"],
    egfr_resid = (r_act - r_pred)[genes == "EGFR"],
    egfr_resid_pctile = mean(abs(r_act - r_pred) <= abs((r_act - r_pred)[genes == "EGFR"]), na.rm = TRUE))
}
R <- do.call(rbind, rows); D <- do.call(rbind, dec)
write.csv(R, "results/transport_anchor.csv", row.names = FALSE)
write.csv(D, "results/transport_decomposition.csv", row.names = FALSE)

say("=== Does the transport prediction actually work? ===\n")
say("  correlation between PREDICTED and ACTUAL anchor r across %d genes:\n", length(genes))
for (i in seq_len(nrow(D)))
  say("    %-10s  cor(pred, actual) = %.3f    vs cor(TCGA r, actual) = %.3f  [%d strata, %.0f%% of R's patients]\n",
      D$cohort[i], D$cor_pred_actual[i], D$cor_T_actual[i], D$strata[i], 100*D$weight_kept[i])
say("\n  If cor(pred, actual) > cor(TCGA r, actual), transporting TCGA's conditional\n")
say("  structure to R's composition predicts R's anchor correlations better than\n")
say("  TCGA's raw correlations do -- i.e. composition shift is real and estimable.\n")

say("\n=== Decomposition of cross-cohort anchor discordance ===\n")
say("%-10s %12s %14s %12s %14s\n", "cohort", "|r_R - r_T|", "composition", "residual", "%% composition")
for (i in seq_len(nrow(D)))
  say("%-10s %12.4f %14.4f %12.4f %13.0f%%\n", D$cohort[i], D$mean_total_discord[i],
      D$mean_composition[i], D$mean_residual[i], 100*D$frac_composition[i])

say("\n=== EGFR specifically ===\n")
for (i in seq_len(nrow(D)))
  say("  %-10s r_TCGA=%+0.3f  r_predicted=%+0.3f  r_actual=%+0.3f  residual=%+0.3f (%.0f%% of genes have a smaller residual)\n",
      D$cohort[i], D$egfr_r_T[i], D$egfr_r_pred[i], D$egfr_r_act[i], D$egfr_resid[i],
      100*D$egfr_resid_pctile[i])

say("\n=== Predicting replication: three scores ===\n")
say("%-10s %-11s %6s %8s  %-22s %-22s %-22s\n", "cohort","criterion","n","replic",
    "raw r_T*r_R","transport residual","expected shift (no R expr)")
for (i in seq_len(nrow(R)))
  say("%-10s %-11s %6d %8.3f  %5.3f (%.3f-%.3f)      %5.3f (%.3f-%.3f)      %5.3f (%.3f-%.3f)\n",
      R$cohort[i], R$criterion[i], R$n[i], R$replic_rate[i],
      R$auc_raw[i], R$raw_lo[i], R$raw_hi[i],
      R$auc_residual[i], R$res_lo[i], R$res_hi[i],
      R$auc_expected[i], R$exp_lo[i], R$exp_hi[i])
say("\n  'expected shift' uses TCGA expression and the replication cohort's CLINICAL\n")
say("  TABLE ONLY -- no replication-cohort expression. Any AUC above chance there is\n")
say("  a prospective capability the raw score does not have.\n")
close(out)
cat("\nwrote results/transport_anchor.csv, transport_decomposition.csv, .txt\n")
