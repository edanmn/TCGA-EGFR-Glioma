# 36_threshold_sensitivity.R -- DO THE CONCLUSIONS SURVIVE MY ARBITRARY CHOICES?
#
# Section 6.9 rests on a chain of analyst decisions that were made once, by us,
# without justification from the data:
#
#   * "prognostic" = BH-adjusted p < 0.05 in the discovery half
#   * "replicated" = same sign AND nominal p < 0.05 in the replication cohort
#   * anchor       = tumour grade (IDH was reported for Table 4 but never for the
#                    matched null)
#   * gene universe cap, cluster count, corruption fraction (elsewhere)
#
# Each was defensible; none was pre-registered. A conclusion that holds only at
# our particular settings is a garden-of-forking-paths result, and this script
# is the check. Because the thresholds are applied AFTER the Cox models, the
# whole grid costs one model pass per arm.
#
# Grid: FDR in {0.01, 0.05, 0.10} x replication p in {0.01, 0.05, 0.10}
#       x anchor in {grade, IDH} x metric in {product, difference}
#       for the three glioma settings, on Table 4's universe.
#
# The pre-specified summary is the fraction of the grid in which each metric's
# null-vs-real difference keeps the sign the manuscript reports.
#
# Run from project root:  Rscript R/36_threshold_sensitivity.R   (~12 min)
# Writes: results/threshold_sensitivity.csv, results/threshold_sensitivity.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)
set.seed(20260815)

out <- file("results/threshold_sensitivity.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }
NCAP <- 8000

tcga_glioma <- function(){
  one <- function(project, gradeIV = FALSE){
    se <- load_se(project); cd <- as.data.frame(colData(se)); k <- substr(colnames(se), 14, 15) == "01"
    e <- log2(assay(se, "tpm_unstrand")[, k] + 1); rownames(e) <- rowData(se)$gene_name
    e <- e[!duplicated(rownames(e)), ]
    b <- .clinical_df(se, project); b <- b[match(colnames(e), b$barcode), ]
    gr <- as.character(cd$paper_Grade[k])
    g <- if (gradeIV) ifelse(grepl("G4|IV", gr), 4, NA) else ifelse(grepl("G3|III", gr), 3, ifelse(grepl("G2| II", gr), 2, NA))
    idh <- as.character(cd$paper_IDH.status[k])
    keep <- !is.na(b$time) & b$time > 0 & !is.na(b$event) & !duplicated(b$patient)
    list(e = e[, keep], cl = data.frame(time = b$time[keep], event = b$event[keep], age = b$age[keep],
         grade = g[keep], idh = ifelse(idh[keep]=="Mutant",1,ifelse(idh[keep]=="WT",0,NA))))
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
           grade = as.numeric(factor(grade, levels = c("WHO II","WHO III","WHO IV"))) + 1,
           idh = ifelse(idh == "Mutant", 1, ifelse(idh == "Wildtype", 0, NA)))
  sel <- intersect(colnames(ex), cl$id); cl <- cl[match(sel, cl$id), ]
  ok <- !is.na(cl$time) & cl$time > 0 & !is.na(cl$event) & !is.na(cl$grade) & !is.na(cl$age)
  list(e = ex[, sel][, ok, drop = FALSE], cl = cl[ok, c("time","event","age","grade","idh")])
}
say("loading cohorts...\n")
GL <- tcga_glioma()
C693 <- cgga_load("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt","data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
C325 <- cgga_load("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt","data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
C301 <- cgga_load("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt","data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt = FALSE)
okg <- !is.na(GL$cl$grade) & !is.na(GL$cl$age); GL$e <- GL$e[, okg]; GL$cl <- GL$cl[okg, ]
uni <- local({
  g <- Reduce(intersect, list(rownames(GL$e), rownames(C693$e), rownames(C325$e), rownames(C301$e)))
  sdT <- apply(GL$e[g, , drop = FALSE], 1, sd); g <- g[is.finite(sdT) & sdT > 0.5]
  if (length(g) > NCAP) g <- names(sort(sdT[g], decreasing = TRUE))[1:NCAP]; g })

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
anchors <- function(e, cl, which){
  a <- if (which == "grade") cl$grade else cl$idh
  suppressWarnings(apply(e, 1, function(x) cor(x, a, method = "spearman", use = "complete.obs")))
}
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
  nullidx <- sample.int(nrow(GL$cl), nR); discidx <- setdiff(seq_len(nrow(GL$cl)), nullidx)
  realidx <- sample.int(nrow(R$cl), nR)
  say("\n=== %s ===  genes=%d discovery=%d arms=%d\n", lab, length(genes), length(discidx), nR)

  sD <- cox_stats(eD[, discidx, drop=FALSE], GL$cl[discidx, ])
  sN <- cox_stats(eD[, nullidx, drop=FALSE], GL$cl[nullidx, ])
  sR <- cox_stats(eR[, realidx, drop=FALSE], R$cl[realidx, ])
  A <- list()
  for (an in c("grade","idh")) {
    A[[an]] <- list(D = anchors(eD[, discidx, drop=FALSE], GL$cl[discidx, ], an),
                    N = anchors(eD[, nullidx, drop=FALSE], GL$cl[nullidx, ], an),
                    R = anchors(eR[, realidx, drop=FALSE], R$cl[realidx, ], an))
  }
  fdrBH <- p.adjust(sD$p, "BH")
  for (fdr in c(0.01, 0.05, 0.10)) for (pr in c(0.01, 0.05, 0.10)) for (an in c("grade","idh")) {
    prog <- fdrBH < fdr & !is.na(sN$b) & !is.na(sR$b)
    if (sum(prog) < 100) next
    lN <- as.integer(sign(sD$b)==sign(sN$b) & sN$p < pr)[prog]
    lR <- as.integer(sign(sD$b)==sign(sR$b) & sR$p < pr)[prog]
    for (met in c("prod","absd")) {
      sc <- function(x, y) if (met=="prod") x*y else -abs(x-y)
      dN <- sc(A[[an]]$D, A[[an]]$N)[prog]; dR <- sc(A[[an]]$D, A[[an]]$R)[prog]
      d <- auc_fast(dR, lR) - auc_fast(dN, lN)
      rows[[length(rows)+1]] <- data.frame(setting=lab, fdr=fdr, repl_p=pr, anchor=an,
        metric=met, n_prog=sum(prog), auc_null=auc_fast(dN,lN), auc_real=auc_fast(dR,lR), diff=d)
    }
  }
}
G <- do.call(rbind, rows)
write.csv(G, "results/threshold_sensitivity.csv", row.names = FALSE)

say("\n\n=== Sign stability across the analyst-choice grid ===\n")
say("  manuscript reports: difference form POSITIVE in all 3 glioma pairs;\n")
say("  product form positive in CGGA-325 and array-301, NEGATIVE (inverted) in CGGA-693.\n\n")
say("%-18s %-6s %-7s %7s %9s %9s %s\n","setting","metric","anchor","cells","sign +","sign -","matches manuscript")
expect <- list(`TCGA->CGGA-693`=list(prod=-1, absd=+1),
               `TCGA->CGGA-325`=list(prod=+1, absd=+1),
               `TCGA->array-301`=list(prod=+1, absd=+1))
tot <- 0; agree <- 0
for (S in unique(G$setting)) for (met in c("prod","absd")) for (an in c("grade","idh")) {
  g <- G[G$setting==S & G$metric==met & G$anchor==an & is.finite(G$diff), ]
  if (!nrow(g)) next
  e <- expect[[S]][[met]]
  ok <- sum(sign(g$diff) == e); tot <- tot + nrow(g); agree <- agree + ok
  say("%-18s %-6s %-7s %7d %9d %9d %5.0f%%\n", S, met, an, nrow(g),
      sum(g$diff>0), sum(g$diff<0), 100*ok/nrow(g))
}
say("\n  overall: %d of %d grid cells agree with the manuscript's reported sign (%.0f%%)\n",
    agree, tot, 100*agree/tot)
gr <- G[G$anchor=="grade" & is.finite(G$diff), ]; id <- G[G$anchor=="idh" & is.finite(G$diff), ]
say("\n  by anchor: grade %.0f%% agree, IDH %.0f%% agree\n",
    100*mean(mapply(function(s,m,d) sign(d)==expect[[s]][[m]], gr$setting, gr$metric, gr$diff)),
    100*mean(mapply(function(s,m,d) sign(d)==expect[[s]][[m]], id$setting, id$metric, id$diff)))
say("\n  A high agreement rate means the conclusions are not an artifact of our\n")
say("  particular thresholds. Cells that disagree are reported, not dropped.\n")
close(out)
cat("\nwrote results/threshold_sensitivity.csv and .txt\n")
