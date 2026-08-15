# 35_cluster_bootstrap.R -- IS THE GENE BOOTSTRAP VALID?
#
# R/33 and R/34 put confidence intervals on the null-vs-real AUC difference by
# resampling GENES independently with replacement. That assumes genes are
# exchangeable and independent. They are not: transcriptomes are strongly
# co-expressed, so genes come in correlated blocks and an independent gene
# bootstrap treats correlated observations as if they carried independent
# information. The standard consequence is intervals that are too NARROW, which
# would make §6.9's conclusions look firmer than the data support.
#
# This is the same class of error as the one already corrected once in this
# project (replicate t-tests, R/33), so it is worth testing rather than assuming.
#
# Test: cluster the genes by expression profile in the discovery cohort, then
# resample whole CLUSTERS with replacement instead of individual genes. A cluster
# bootstrap respects the correlation blocks and is conservative. If the interval
# widths are similar, the naive bootstrap was adequate; if they widen enough to
# cross zero, §6.9 must be re-stated.
#
# Reported for the three glioma settings and both anchor forms, on Table 4's
# universe (the manuscript's primary analysis).
#
# Run from project root:  Rscript R/35_cluster_bootstrap.R   (~12 min)
# Writes: results/cluster_bootstrap.csv, results/cluster_bootstrap.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)
set.seed(20260815)

out <- file("results/cluster_bootstrap.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }
NCAP <- 8000; NBOOT <- 2000; NCLUST <- 200

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

say("loading cohorts...\n")
GL <- tcga_glioma()
C693 <- cgga_load("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt","data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
C325 <- cgga_load("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt","data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
C301 <- cgga_load("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt","data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt = FALSE)
okg <- !is.na(GL$cl$anchor) & !is.na(GL$cl$age); GL$e <- GL$e[, okg]; GL$cl <- GL$cl[okg, ]
uni <- local({
  g <- Reduce(intersect, list(rownames(GL$e), rownames(C693$e), rownames(C325$e), rownames(C301$e)))
  sdT <- apply(GL$e[g, , drop = FALSE], 1, sd); g <- g[is.finite(sdT) & sdT > 0.5]
  if (length(g) > NCAP) g <- names(sort(sdT[g], decreasing = TRUE))[1:NCAP]
  g })

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
auc_fast <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1*(n1+1)/2) / (n1*n0)
}

SET <- list(list("glioma: TCGA->CGGA-693", C693), list("glioma: TCGA->CGGA-325", C325),
            list("glioma: TCGA->array-301", C301))
rows <- list()
for (S in SET) {
  lab <- S[[1]]; R <- S[[2]]
  genes <- intersect(uni, rownames(R$e))
  eD <- GL$e[genes, , drop = FALSE]; eR <- R$e[genes, , drop = FALSE]
  nR <- min(nrow(R$cl), floor(nrow(GL$cl)/2))
  nullidx <- sample.int(nrow(GL$cl), nR); discidx <- setdiff(seq_len(nrow(GL$cl)), nullidx)
  realidx <- sample.int(nrow(R$cl), nR)
  say("\n=== %s ===\n  genes=%d discovery n=%d arms n=%d\n", lab, length(genes), length(discidx), nR)

  sD <- cox_stats(eD[, discidx, drop = FALSE], GL$cl[discidx, ])
  prog <- p.adjust(sD$p, "BH") < 0.05
  arm <- function(eA, clA){
    sA <- cox_stats(eA, clA)
    rD <- suppressWarnings(apply(eD[, discidx, drop=FALSE], 1, function(x) cor(x, GL$cl$anchor[discidx], method="spearman", use="complete.obs")))
    rR <- suppressWarnings(apply(eA, 1, function(x) cor(x, clA$anchor, method="spearman", use="complete.obs")))
    list(repl = as.integer(sign(sD$b) == sign(sA$b) & sA$p < 0.05),
         prod = rD*rR, absd = -abs(rD-rR), ok = !is.na(sA$b))
  }
  AN <- arm(eD[, nullidx, drop=FALSE], GL$cl[nullidx, ])
  AR <- arm(eR[, realidx, drop=FALSE], R$cl[realidx, ])
  keep <- which(prog & AN$ok & AR$ok)

  ## cluster the KEPT genes by discovery-cohort expression profile
  Z <- t(scale(t(eD[keep, discidx, drop = FALSE])))
  Z[!is.finite(Z)] <- 0
  km <- kmeans(Z, centers = min(NCLUST, nrow(Z) - 1), iter.max = 50, nstart = 3)
  cl_id <- km$cluster
  idx_by_cluster <- split(seq_along(keep), cl_id)
  # how much correlation structure is there? mean within-cluster |r| vs overall
  say("  %d genes in %d clusters (median size %.0f)\n",
      length(keep), length(idx_by_cluster), median(lengths(idx_by_cluster)))

  for (met in c("prod","absd")) {
    sN <- AN[[met]][keep]; sR2 <- AR[[met]][keep]
    lN <- AN$repl[keep];   lR <- AR$repl[keep]
    obs <- auc_fast(sR2, lR) - auc_fast(sN, lN)
    # (a) naive gene bootstrap, as in R/33-34
    bs_g <- replicate(NBOOT, { i <- sample(length(keep), replace = TRUE)
      auc_fast(sR2[i], lR[i]) - auc_fast(sN[i], lN[i]) })
    ci_g <- quantile(bs_g, c(.025,.975), na.rm = TRUE)
    # (b) cluster bootstrap: resample whole co-expression clusters
    K <- length(idx_by_cluster)
    bs_c <- replicate(NBOOT, {
      cs <- sample.int(K, K, replace = TRUE)
      i <- unlist(idx_by_cluster[cs], use.names = FALSE)
      auc_fast(sR2[i], lR[i]) - auc_fast(sN[i], lN[i]) })
    ci_c <- quantile(bs_c, c(.025,.975), na.rm = TRUE)
    wg <- diff(ci_g); wc <- diff(ci_c)
    say("  %-5s diff=%+.3f | gene CI (%+.3f,%+.3f) w=%.3f | cluster CI (%+.3f,%+.3f) w=%.3f | inflation %.2fx | %s\n",
        met, obs, ci_g[1], ci_g[2], wg, ci_c[1], ci_c[2], wc, wc/wg,
        ifelse(sign(ci_c[1]) == sign(ci_c[2]), "still excludes 0", "NOW SPANS 0"))
    rows[[length(rows)+1]] <- data.frame(setting=lab, metric=met, diff=obs,
      gene_lo=unname(ci_g[1]), gene_hi=unname(ci_g[2]), gene_width=unname(wg),
      clust_lo=unname(ci_c[1]), clust_hi=unname(ci_c[2]), clust_width=unname(wc),
      inflation=unname(wc/wg), clust_excludes0=sign(ci_c[1])==sign(ci_c[2]),
      n_genes=length(keep), n_clusters=K)
  }
}
G <- do.call(rbind, rows)
write.csv(G, "results/cluster_bootstrap.csv", row.names = FALSE)

say("\n\n=== VERDICT ===\n")
say("  median interval inflation from respecting co-expression blocks: %.2fx\n", median(G$inflation))
say("  range: %.2fx to %.2fx\n", min(G$inflation), max(G$inflation))
ab <- G[G$metric == "absd", ]
say("  difference form: %d of %d glioma pairs still exclude zero under the cluster bootstrap\n",
    sum(ab$clust_excludes0), nrow(ab))
pr <- G[G$metric == "prod", ]
say("  product form:    %d of %d still exclude zero\n", sum(pr$clust_excludes0), nrow(pr))
say("\n  An inflation factor near 1 means the naive gene bootstrap was adequate despite\n")
say("  co-expression. A large factor, or any conclusion crossing zero, means §6.9's\n")
say("  intervals are optimistic and the manuscript must report the cluster version.\n")
close(out)
cat("\nwrote results/cluster_bootstrap.csv and .txt\n")
