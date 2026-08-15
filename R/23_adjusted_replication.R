# 23_adjusted_replication.R -- decompose the circularity concern into its two
# separate components, across all three replication cohorts.
#
# The concern: survival in pooled glioma is grade-dominated, so sign(Cox coef) is
# close to sign(r(expression, grade)) in BOTH cohorts. The grade-anchor QC score
# (r_T * r_R) is then near-algebraically tied to the replication criterion
# (sign agreement + p<0.05), and the headline AUC would be partly tautological.
#
# 18_circularity_control.R tests this by restricting to the subtype-adjusted
# prognostic set, but that changes TWO things at once -- which genes are in the
# evaluation set AND how replication is defined -- on only 183 genes. This script
# separates them into a 2x2, so each factor is isolated on a well-powered set:
#
#                             replication criterion
#                             unadjusted        grade+IDH-adjusted
#   discovery  unadjusted     (A) headline      (B) <-- KEY missing cell
#   set        adjusted       (C)               (D) = Table S7 row 2
#
# Cell B keeps the full ~6,500-gene evaluation set and changes ONLY the criterion:
# if the AUC survives there, the QC is not merely re-reading the grade signal.
#
# Uses the SAME four-cohort gene universe as 16_method_framework.R (Table 4), so
# cell A reproduces the headline AUCs exactly and the cells are comparable.
#
# Run from project root:  Rscript R/23_adjusted_replication.R
# Writes: results/adjusted_replication.csv, results/adjusted_replication.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)

out <- file("results/adjusted_replication.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

## ---------- cohort assembly (identical to 16_method_framework.R) ----------
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
# mRNA-array_301 is log-ratio scaled with negative values; log2(x+1) would drop
# every value <= -1 (see paper section 4.6).
cgga_mat <- function(gfile, cfile, cols, logt = TRUE){
  ex <- rd(gfile); rownames(ex) <- ex[[1]]; ex[[1]] <- NULL
  ex <- if (logt) log2(as.matrix(ex) + 1) else as.matrix(ex)
  cl <- read.delim(cfile, check.names = FALSE); names(cl)[cols] <- c("id","prs","grade","age","os","censor","idh")
  cl <- cl %>% filter(prs == "Primary") %>%
    mutate(time = as.numeric(os), event = as.integer(censor), age = as.numeric(age),
           grade = as.numeric(factor(grade, levels = c("WHO II","WHO III","WHO IV"))),
           idh = ifelse(idh == "Mutant", 1, ifelse(idh == "Wildtype", 0, NA)))
  sel <- intersect(colnames(ex), cl$id); cl <- cl[match(sel, cl$id), ]
  list(e = ex[, sel], clin = cl[, c("time","event","age","grade","idh")])
}

TC   <- tcga_mat()
C693 <- cgga_mat("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", "data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
C325 <- cgga_mat("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", "data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
CARR <- cgga_mat("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt", "data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt = FALSE)

## ---------- gene universe: four-cohort intersection, as in Table 4 ----------
genes <- Reduce(intersect, list(rownames(TC$e), rownames(C693$e), rownames(C325$e), rownames(CARR$e)))
sdT <- apply(TC$e[genes, ], 1, sd); genes <- genes[sdT > 0.5]
if (length(genes) > 8000) genes <- names(sort(sdT[genes], decreasing = TRUE))[1:8000]
say("Gene universe (four-cohort intersection, TCGA SD>0.5): %d\n", length(genes))

## ---------- per-gene: unadjusted AND subtype-adjusted Cox, in one pass ----------
# Complete-case index for the adjusted model is computed once per cohort, not per gene.
stats_cohort <- function(M, gs, lab){
  cl <- M$clin; e <- M$e
  cc <- !is.na(cl$idh) & !is.na(cl$grade)          # adjusted model complete cases
  S0 <- Surv(cl$time, cl$event); age0 <- cl$age
  S1 <- Surv(cl$time[cc], cl$event[cc])
  age1 <- cl$age[cc]; gr1 <- factor(cl$grade[cc]); idh1 <- factor(cl$idh[cc])
  say("  %-10s n=%d (adjusted-model complete cases n=%d, events=%d)\n",
      lab, nrow(cl), sum(cc), sum(cl$event[cc], na.rm = TRUE))
  do.call(rbind, lapply(gs, function(g){
    x <- e[g, ]; z <- as.numeric(scale(x))
    m0 <- tryCatch(coxph(S0 ~ z + age0), error = function(e) NULL)
    zc <- z[cc]
    m1 <- tryCatch(coxph(S1 ~ zc + age1 + gr1 + idh1), error = function(e) NULL)
    b0 <- if (is.null(m0)) c(NA,NA,NA) else summary(m0)$coefficients["z",  c("coef","se(coef)","Pr(>|z|)")]
    b1 <- if (is.null(m1)) c(NA,NA,NA) else summary(m1)$coefficients["zc", c("coef","se(coef)","Pr(>|z|)")]
    data.frame(gene = g,
               b0 = b0[1], se0 = b0[2], p0 = b0[3],
               b1 = b1[1], se1 = b1[2], p1 = b1[3],
               rg = suppressWarnings(cor(x, cl$grade, method = "spearman", use = "complete.obs")),
               ri = suppressWarnings(cor(x, cl$idh,   method = "spearman", use = "complete.obs")),
               meanx = mean(x, na.rm = TRUE))
  }))
}
say("\nFitting per-gene models (2 Cox fits x %d genes x 4 cohorts)...\n", length(genes))
sT <- stats_cohort(TC,   genes, "TCGA")
s6 <- stats_cohort(C693, genes, "CGGA-693")
s3 <- stats_cohort(C325, genes, "CGGA-325")
sA <- stats_cohort(CARR, genes, "array-301")

## ---------- AUC with bootstrap CI ----------
auc <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
aucCI <- function(score, label, B = 1000){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  a <- auc(score, label)
  if (is.na(a)) return(c(a = NA, lo = NA, hi = NA))
  bs <- replicate(B, { i <- sample(length(score), replace = TRUE); auc(score[i], label[i]) })
  c(a = a, lo = unname(quantile(bs, .025, na.rm = TRUE)), hi = unname(quantile(bs, .975, na.rm = TRUE)))
}

## ---------- the 2x2 ----------
# disc: which TCGA column defines "prognostic" (BH-adjusted within the universe)
# repl: which replication-cohort column defines "replicated"
cell <- function(sR, cohort, disc, repl){
  D <- merge(sT, sR, by = "gene", suffixes = c("_T","_R"))
  bT <- D[[paste0(disc, "_T")]]; pT <- D[[if (disc == "b0") "p0_T" else "p1_T"]]
  bR <- D[[paste0(repl, "_R")]]; pR <- D[[if (repl == "b0") "p0_R" else "p1_R"]]
  seR <- D[[if (repl == "b0") "se0_R" else "se1_R"]]
  prog <- p.adjust(pT, "BH") < 0.05 & !is.na(bR)
  D <- D[prog, ]; bT <- bT[prog]; bR <- bR[prog]; pR <- pR[prog]; seR <- seR[prog]
  replicated <- as.integer(sign(bT) == sign(bR) & pR < 0.05)
  qc_grade <- D$rg_T * D$rg_R
  qc_idh   <- D$ri_T * D$ri_R
  detect   <- pmin(D$meanx_T, D$meanx_R)
  precis   <- -seR
  g <- aucCI(qc_grade, replicated); i <- aucCI(qc_idh, replicated)
  data.frame(cohort = cohort,
             discovery   = if (disc == "b0") "unadjusted" else "grade+IDH-adjusted",
             replication = if (repl == "b0") "unadjusted" else "grade+IDH-adjusted",
             n = nrow(D), replic_rate = round(mean(replicated), 3),
             AUC_grade = round(g["a"], 3), grade_lo = round(g["lo"], 3), grade_hi = round(g["hi"], 3),
             AUC_idh   = round(i["a"], 3),
             AUC_detect = round(auc(detect, replicated), 3),
             AUC_precision = round(auc(precis, replicated), 3),
             row.names = NULL)
}

cohorts <- list(`CGGA-693` = s6, `CGGA-325` = s3, `array-301` = sA)
res <- do.call(rbind, lapply(names(cohorts), function(k)
  do.call(rbind, list(
    cell(cohorts[[k]], k, "b0", "b0"),   # A: headline
    cell(cohorts[[k]], k, "b0", "b1"),   # B: KEY -- criterion changed only
    cell(cohorts[[k]], k, "b1", "b0"),   # C: set changed only
    cell(cohorts[[k]], k, "b1", "b1")    # D: both (= Table S7 row 2 regime)
  ))))

write.csv(res, "results/adjusted_replication.csv", row.names = FALSE)

say("\n=== 2x2 circularity decomposition (grade-anchor QC) ===\n")
say("%-10s %-19s %-19s %6s %7s  %-22s %7s %7s %7s\n",
    "cohort","discovery set","replication crit.","n","rep.rate","AUC grade (95% CI)","AUC IDH","detect","precis")
for (i in seq_len(nrow(res))) {
  say("%-10s %-19s %-19s %6d %7.3f  %5.3f (%.3f-%.3f)      %7.3f %7.3f %7.3f\n",
      res$cohort[i], res$discovery[i], res$replication[i], res$n[i], res$replic_rate[i],
      res$AUC_grade[i], res$grade_lo[i], res$grade_hi[i],
      res$AUC_idh[i], res$AUC_detect[i], res$AUC_precision[i])
}

say("\nRead cell B (unadjusted discovery set, adjusted replication criterion) against\n")
say("cell A: it holds the evaluation set fixed and removes the shared grade signal\n")
say("from the replication definition only. Cell D reproduces the Table S7 regime.\n")
close(out)
cat("\nwrote results/adjusted_replication.csv and .txt\n")
