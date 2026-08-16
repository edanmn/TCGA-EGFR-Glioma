# 38_split_variance.R -- COMBINE THE TWO SOURCES OF VARIABILITY.
#
# Limitation twelve currently concedes something the manuscript never fixed: the
# intervals in Table 5 come from ONE split per setting, bootstrapped over genes,
# so they express gene-level sampling only. Split-to-split variability is real --
# two independent splits gave array-301 difference-form estimates of +0.148 and
# +0.227 -- and is not in those intervals.
#
# This script closes it. K independent splits per setting; within each split the
# gene-level paired bootstrap gives a within-split variance; across splits we get
# a between-split variance. They combine by the standard Rubin rule for repeated
# sampling:
#
#     Var_total = mean(Var_within) + (1 + 1/K) * Var_between
#
# The combined interval is the honest one: it is wider than either component and
# is what the manuscript should report if the conclusions survive it.
#
# Run on Table 4's universe so it speaks directly to the primary analysis.
#
# Run from project root:  Rscript R/38_split_variance.R   (~35 min)
# Writes: results/split_variance.csv, results/split_variance.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)
set.seed(20260815)

out <- file("results/split_variance.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }
NCAP <- 8000; NBOOT <- 400; K <- 6

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
  if (length(g) > NCAP) g <- names(sort(sdT[g], decreasing = TRUE))[1:NCAP]; g })
say("universe: %d genes; K=%d splits per setting; B=%d bootstrap per split\n\n", length(uni), K, NBOOT)

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
anc <- function(e, a) suppressWarnings(apply(e, 1, function(x) cor(x, a, method="spearman", use="complete.obs")))
auc_fast <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label==1); n0 <- sum(label==0); if (n1==0||n0==0) return(NA_real_)
  r <- rank(score); (sum(r[label==1]) - n1*(n1+1)/2)/(n1*n0)
}

SET <- list(list("TCGA->CGGA-693", C693), list("TCGA->CGGA-325", C325), list("TCGA->array-301", C301))
rows <- list()
for (S in SET) {
  lab <- S[[1]]; R <- S[[2]]
  genes <- intersect(uni, rownames(R$e))
  eD <- GL$e[genes, , drop=FALSE]; eR <- R$e[genes, , drop=FALSE]
  nR <- min(nrow(R$cl), floor(nrow(GL$cl)/2))
  per <- list(prod = list(est=numeric(0), var=numeric(0)), absd = list(est=numeric(0), var=numeric(0)))
  for (k in seq_len(K)) {
    cat(sprintf("  %s split %d/%d\n", lab, k, K))
    nullidx <- sample.int(nrow(GL$cl), nR); discidx <- setdiff(seq_len(nrow(GL$cl)), nullidx)
    realidx <- sample.int(nrow(R$cl), nR)
    sD <- cox_stats(eD[, discidx, drop=FALSE], GL$cl[discidx, ])
    sN <- cox_stats(eD[, nullidx, drop=FALSE], GL$cl[nullidx, ])
    sR <- cox_stats(eR[, realidx, drop=FALSE], R$cl[realidx, ])
    aD <- anc(eD[, discidx, drop=FALSE], GL$cl$anchor[discidx])
    aN <- anc(eD[, nullidx, drop=FALSE], GL$cl$anchor[nullidx])
    aR <- anc(eR[, realidx, drop=FALSE], R$cl$anchor[realidx])
    keep <- which(p.adjust(sD$p,"BH") < 0.05 & !is.na(sN$b) & !is.na(sR$b))
    lN <- as.integer(sign(sD$b)==sign(sN$b) & sN$p<0.05)[keep]
    lR <- as.integer(sign(sD$b)==sign(sR$b) & sR$p<0.05)[keep]
    for (met in c("prod","absd")) {
      f <- function(x,y) if (met=="prod") x*y else -abs(x-y)
      dN <- f(aD,aN)[keep]; dR <- f(aD,aR)[keep]
      est <- auc_fast(dR,lR) - auc_fast(dN,lN)
      bs <- replicate(NBOOT, { i <- sample(length(keep), replace=TRUE)
        auc_fast(dR[i],lR[i]) - auc_fast(dN[i],lN[i]) })
      per[[met]]$est <- c(per[[met]]$est, est)
      per[[met]]$var <- c(per[[met]]$var, var(bs, na.rm=TRUE))
    }
  }
  for (met in c("prod","absd")) {
    e <- per[[met]]$est; w <- per[[met]]$var
    Vw <- mean(w, na.rm=TRUE); Vb <- var(e, na.rm=TRUE)
    Vt <- Vw + (1 + 1/K) * Vb
    m <- mean(e, na.rm=TRUE); se <- sqrt(Vt)
    lo <- m - 1.96*se; hi <- m + 1.96*se
    say("%-18s %-5s  est=%+.3f  within SD=%.4f  between SD=%.4f  combined SD=%.4f  95%% CI (%+.3f, %+.3f)  %s\n",
        lab, met, m, sqrt(Vw), sqrt(Vb), se, lo, hi,
        ifelse(sign(lo)==sign(hi), "excludes 0", "SPANS 0"))
    rows[[length(rows)+1]] <- data.frame(setting=lab, metric=met, K=K, est=m,
      sd_within=sqrt(Vw), sd_between=sqrt(Vb), sd_total=se, lo=lo, hi=hi,
      excludes0=sign(lo)==sign(hi), inflation_vs_within=se/sqrt(Vw),
      est_min=min(e), est_max=max(e))
  }
}
G <- do.call(rbind, rows)
write.csv(G, "results/split_variance.csv", row.names = FALSE)

say("\n\n=== VERDICT ===\n")
say("  Between-split share of total VARIANCE: %.0f%%-%.0f%%  (of total SD: %.0f%%-%.0f%%)\n",
    100*min(G$sd_between^2/G$sd_total^2), 100*max(G$sd_between^2/G$sd_total^2),
    100*min(G$sd_between/G$sd_total), 100*max(G$sd_between/G$sd_total))
say("  Interval inflation vs gene-bootstrap-only: %.2fx to %.2fx (median %.2fx)\n",
    min(G$inflation_vs_within), max(G$inflation_vs_within), median(G$inflation_vs_within))
ab <- G[G$metric=="absd", ]; pr <- G[G$metric=="prod", ]
say("  difference form: %d of %d glioma pairs still exclude zero\n", sum(ab$excludes0), nrow(ab))
say("  product form:    %d of %d still conclusive\n", sum(pr$excludes0), nrow(pr))
say("  point-estimate range across splits (difference form): %+.3f to %+.3f\n",
    min(ab$est_min), max(ab$est_max))
say("\n  This is the interval limitation twelve says the manuscript does not report.\n")
say("  If the conclusions survive it, limitation twelve can be downgraded from an\n")
say("  unquantified caveat to a quantified one.\n")
close(out)
cat("\nwrote results/split_variance.csv and .txt\n")
