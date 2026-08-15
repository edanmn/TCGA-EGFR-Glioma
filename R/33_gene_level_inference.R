# 33_gene_level_inference.R -- CORRECT THE INFERENCE IN TABLE 5.
#
# 30_multisetting_null.R reported Welch t-tests comparing the null and real arms
# across 8 replicates. That inference is invalid, and the review that found it can
# be reproduced in one line: re-running the same comparison with 3, 4, ... 8 of the
# replicates moves p from 0.32 to 0.010 while the effect size barely moves
# (+0.030 -> +0.037). The replicates are RESAMPLES OF ONE FIXED COHORT, so their
# spread is Monte-Carlo error of the estimate, which shrinks as 1/sqrt(reps).
# The p-value therefore measures how many replicates we chose to run, not how
# strong the evidence is, and would approach 0 for any nonzero difference.
#
# The correct sampling unit is the GENE: the AUC is a statistic over genes, and the
# scientific question ("does this metric rank genes by replication better against a
# real cohort than against a same-population null?") is a question about the gene
# population. Both arms are scored on the SAME genes, so the comparison is paired.
#
# This script therefore:
#   * runs ONE clean replicate per setting and retains the per-gene vectors,
#   * computes the null-vs-real AUC difference with a PAIRED BOOTSTRAP OVER GENES
#     (2,000 resamples), giving a CI whose width reflects gene-level sampling
#     variability rather than an arbitrary replicate count,
#   * applies Benjamini-Hochberg across the settings tested, which 30 did not do,
#   * fixes the co-expression metric's self-correlation leakage (a gene's own
#     correlation with itself, r=1, was included when it fell in the hub set).
#
# Run from project root:  Rscript R/33_gene_level_inference.R   (~5 min)
# Writes: results/gene_level_inference.csv, results/gene_level_inference.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)
set.seed(20260815)

out <- file("results/gene_level_inference.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }
NGENE <- 2000; NHUB <- 100; NBOOT <- 2000

## ---------------- cohort loaders (identical to R/30) ----------------
tcga_glioma <- function(){
  one <- function(project, gradeIV = FALSE){
    se <- load_se(project); cd <- as.data.frame(colData(se)); k <- substr(colnames(se), 14, 15) == "01"
    e <- log2(assay(se, "tpm_unstrand")[, k] + 1); rownames(e) <- rowData(se)$gene_name
    e <- e[!duplicated(rownames(e)), ]
    b <- .clinical_df(se, project); b <- b[match(colnames(e), b$barcode), ]
    gr <- as.character(cd$paper_Grade[k])
    g <- if (gradeIV) ifelse(grepl("G4|IV", gr), 4, NA) else ifelse(grepl("G3|III", gr), 3, ifelse(grepl("G2| II", gr), 2, NA))
    keep <- !is.na(b$time) & b$time > 0 & !is.na(b$event) & !duplicated(b$patient)
    list(e = e[, keep], cl = data.frame(time = b$time[keep], event = b$event[keep],
                                        age = b$age[keep], anchor = g[keep]))
  }
  L <- one("TCGA-LGG"); G <- one("TCGA-GBM", TRUE)
  cg <- intersect(rownames(L$e), rownames(G$e))
  list(e = cbind(L$e[cg, ], G$e[cg, ]), cl = rbind(L$cl, G$cl))
}
cgga_load <- function(gfile, cfile, cols, logt = TRUE){
  ex <- rd(gfile); rownames(ex) <- ex[[1]]; ex[[1]] <- NULL
  ex <- if (logt) log2(as.matrix(ex) + 1) else as.matrix(ex)
  cl <- read.delim(cfile, check.names = FALSE); names(cl)[cols] <- c("id","prs","grade","age","os","censor","idh")
  cl <- cl %>% filter(prs == "Primary") %>%
    mutate(time = as.numeric(os), event = as.integer(censor), age = as.numeric(age),
           anchor = as.numeric(factor(grade, levels = c("WHO II","WHO III","WHO IV"))) + 1)
  sel <- intersect(colnames(ex), cl$id); cl <- cl[match(sel, cl$id), ]
  ok <- !is.na(cl$time) & cl$time > 0 & !is.na(cl$event) & !is.na(cl$anchor) & !is.na(cl$age)
  list(e = ex[, sel][, ok, drop = FALSE], cl = cl[ok, c("time","event","age","anchor")])
}
breast_load <- function(){
  dd <- "data/brca/moanna_data"
  le <- function(f){ m <- rd(f); rownames(m) <- m[[1]]; m[[1]] <- NULL; m <- as.matrix(m)
    ex <- m[, grep("_EXPR$", colnames(m)), drop = FALSE]; colnames(ex) <- sub("_EXPR$","",colnames(ex)); ex }
  lb <- function(f){ l <- rd(f); rownames(l) <- l[[1]]; l[[1]] <- NULL; l }
  mbE <- rbind(le(file.path(dd,"training/moanna_training_data.tsv")), le(file.path(dd,"training/moanna_validation_data.tsv")))
  tcE <- le(file.path(dd,"testing/moanna_testing_data.tsv"))
  mbL <- rbind(lb(file.path(dd,"training/moanna_training_label.tsv")), lb(file.path(dd,"training/moanna_validation_label.tsv")))
  tcL <- lb(file.path(dd,"testing/moanna_testing_label.tsv"))
  mbc <- read.delim("data/brca/mb_clin.txt", comment.char="#", check.names=FALSE)
  mbc <- data.frame(id=mbc$PATIENT_ID, time=as.numeric(mbc$OS_MONTHS),
                    event=as.integer(grepl("DECEASED", mbc$OS_STATUS)), age=suppressWarnings(as.numeric(mbc$AGE_AT_DIAGNOSIS)))
  tcc <- read.delim("data/brca/brca_tcga_clin.txt", comment.char="#", check.names=FALSE)
  tcc <- data.frame(id=tcc$PATIENT_ID, time=as.numeric(tcc$OS_MONTHS),
                    event=as.integer(grepl("DECEASED", tcc$OS_STATUS)), age=suppressWarnings(as.numeric(tcc$AGE)))
  asm <- function(E, L, clin, idmap){
    ids <- rownames(E); cl <- clin[match(idmap(ids), clin$id), ]
    ok <- !is.na(cl$time) & cl$time > 0 & !is.na(cl$event) & !is.na(cl$age)
    list(e = t(E[ok, , drop = FALSE]),
         cl = data.frame(time=cl$time[ok], event=cl$event[ok], age=cl$age[ok],
                         anchor=as.numeric(L[ids,"ERStatus"])[ok]))
  }
  list(MB = asm(mbE, mbL, mbc, function(x) x), TC = asm(tcE, tcL, tcc, function(x) substr(x,1,12)))
}

say("loading cohorts...\n")
GL <- tcga_glioma()
C693 <- cgga_load("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt","data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
C325 <- cgga_load("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt","data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
C301 <- cgga_load("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt","data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt = FALSE)
BR <- breast_load()
okg <- !is.na(GL$cl$anchor) & !is.na(GL$cl$age); GL$e <- GL$e[, okg]; GL$cl <- GL$cl[okg, ]

SETTINGS <- list(
  list(lab="glioma: TCGA->CGGA-693",  D=GL, R=C693),
  list(lab="glioma: TCGA->CGGA-325",  D=GL, R=C325),
  list(lab="glioma: TCGA->array-301", D=GL, R=C301),
  list(lab="breast: BRCA->METABRIC",  D=BR$TC, R=BR$MB),
  list(lab="breast: METABRIC->BRCA",  D=BR$MB, R=BR$TC))

cox_stats <- function(e, cl){
  S <- Surv(cl$time, cl$event); age <- cl$age
  b <- p <- numeric(nrow(e))
  for (i in seq_len(nrow(e))) {
    z <- as.numeric(scale(e[i, ]))
    f <- tryCatch(coxph(S ~ z + age), error = function(err) NULL)
    s <- if (is.null(f)) c(NA,NA) else summary(f)$coefficients["z", c("coef","Pr(>|z|)")]
    b[i] <- s[1]; p[i] <- s[2]
  }
  list(b = b, p = p)
}
# LEAKAGE FIX: drop a gene's own column from its hub correlation profile
coexpr <- function(eD, eR, hub){
  cD <- suppressWarnings(cor(t(eD), t(eD[hub, , drop = FALSE])))
  cR <- suppressWarnings(cor(t(eR), t(eR[hub, , drop = FALSE])))
  vapply(seq_len(nrow(eD)), function(i){
    a <- cD[i, ]; b <- cR[i, ]
    j <- which(hub == i); if (length(j)) { a <- a[-j]; b <- b[-j] }   # <-- remove self-correlation
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 10) NA_real_ else suppressWarnings(cor(a[ok], b[ok]))
  }, numeric(1))
}
auc_fast <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1*(n1+1)/2) / (n1*n0)
}

rows <- list()
for (S in SETTINGS) {
  D <- S$D; R <- S$R
  genes <- intersect(rownames(D$e), rownames(R$e))
  sdD <- apply(D$e[genes, , drop = FALSE], 1, sd)
  genes <- names(sort(sdD[is.finite(sdD) & sdD > 0], decreasing = TRUE))
  genes <- genes[seq_len(min(NGENE, length(genes)))]
  eD <- D$e[genes, , drop = FALSE]; eR <- R$e[genes, , drop = FALSE]
  hub <- seq_len(min(NHUB, length(genes)))
  nR <- min(nrow(R$cl), floor(nrow(D$cl)/2))

  nullidx <- sample.int(nrow(D$cl), nR); discidx <- setdiff(seq_len(nrow(D$cl)), nullidx)
  realidx <- sample.int(nrow(R$cl), nR)
  say("\n=== %s ===\n  discovery n=%d, both arms n=%d, genes=%d\n",
      S$lab, length(discidx), nR, length(genes))

  sD <- cox_stats(eD[, discidx, drop = FALSE], D$cl[discidx, ])
  prog <- p.adjust(sD$p, "BH") < 0.05
  arm <- function(eA, clA){
    sA <- cox_stats(eA, clA)
    rD <- suppressWarnings(apply(eD[, discidx, drop=FALSE], 1, function(x) cor(x, D$cl$anchor[discidx], method="spearman", use="complete.obs")))
    rR <- suppressWarnings(apply(eA, 1, function(x) cor(x, clA$anchor, method="spearman", use="complete.obs")))
    list(repl = as.integer(sign(sD$b) == sign(sA$b) & sA$p < 0.05),
         prod = rD*rR, absd = -abs(rD-rR), coex = coexpr(eD[, discidx, drop=FALSE], eA, hub),
         eff = abs(sD$b), ok = !is.na(sA$b))
  }
  AN <- arm(eD[, nullidx, drop=FALSE], D$cl[nullidx, ])
  AR <- arm(eR[, realidx, drop=FALSE], R$cl[realidx, ])
  keep <- which(prog & AN$ok & AR$ok)
  say("  prognostic genes scored in both arms: %d\n", length(keep))

  for (met in c("prod","absd","coex","eff")) {
    sN <- AN[[met]][keep]; sR <- AR[[met]][keep]
    lN <- AN$repl[keep];   lR <- AR$repl[keep]
    obs <- auc_fast(sR, lR) - auc_fast(sN, lN)
    bs <- replicate(NBOOT, { i <- sample(length(keep), replace = TRUE)
      auc_fast(sR[i], lR[i]) - auc_fast(sN[i], lN[i]) })
    ci <- quantile(bs, c(.025,.975), na.rm = TRUE)
    pb <- 2*min(mean(bs <= 0, na.rm=TRUE), mean(bs >= 0, na.rm=TRUE))   # bootstrap two-sided p
    rows[[length(rows)+1]] <- data.frame(setting=S$lab, metric=met,
      auc_null=auc_fast(sN,lN), auc_real=auc_fast(sR,lR), diff=obs,
      lo=unname(ci[1]), hi=unname(ci[2]), p_boot=max(pb, 1/NBOOT), n_genes=length(keep))
    say("  %-6s null=%.3f real=%.3f diff=%+.3f (95%% CI %+.3f to %+.3f) p=%.4g\n",
        met, auc_fast(sN,lN), auc_fast(sR,lR), obs, ci[1], ci[2], max(pb,1/NBOOT))
  }
}
G <- do.call(rbind, rows)
G$fdr <- NA
for (met in unique(G$metric)) G$fdr[G$metric==met] <- p.adjust(G$p_boot[G$metric==met], "BH")
write.csv(G, "results/gene_level_inference.csv", row.names = FALSE)

say("\n\n=== Gene-level paired bootstrap, BH-adjusted within metric ===\n")
say("%-26s %-6s %8s %8s %9s %22s %9s %8s\n","setting","metric","null","real","diff","95%% CI","p","FDR")
for (i in seq_len(nrow(G)))
  say("%-26s %-6s %8.3f %8.3f %+9.3f   (%+.3f, %+.3f) %9.4g %8.4g\n",
      G$setting[i], G$metric[i], G$auc_null[i], G$auc_real[i], G$diff[i], G$lo[i], G$hi[i], G$p_boot[i], G$fdr[i])

say("\n=== What changes relative to the replicate-based t-tests in R/30 ===\n")
gl <- G[G$metric=="absd" & grepl("glioma", G$setting), ]
say("  difference form, glioma: %d of 3 settings with CI excluding zero\n", sum(gl$lo > 0))
gp <- G[G$metric=="prod" & grepl("glioma", G$setting), ]
say("  product form,    glioma: %d of 3 settings with CI excluding zero\n", sum(gp$lo > 0))
say("  These CIs are the inference the paper should report. The replicate-based\n")
say("  p-values are Monte-Carlo artifacts and have been removed from the manuscript.\n")
close(out)
cat("\nwrote results/gene_level_inference.csv and .txt\n")
