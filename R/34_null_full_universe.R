# 34_null_full_universe.R -- ROBUSTNESS OF THE MATCHED NULL TO THE GENE UNIVERSE.
#
# R/30 and R/33 ran the matched-null calibration on the 2,000 most-variable shared
# genes, for runtime. Table 4 -- the table the calibration is meant to reinterpret --
# uses a different universe: the four-cohort intersection with TCGA SD > 0.5, capped
# at the 8,000 most variable, of which 6,512 were prognostic at full discovery n.
# The manuscript therefore reinterprets one gene set using a calibration measured on
# a smaller, more variable one, and flags that as untested (limitation N12).
#
# This script closes that gap. The design is identical to R/33 -- one discovery split
# per setting, a disjoint same-population null arm and a size-matched real arm, paired
# bootstrap over genes -- but run on TABLE 4's UNIVERSE:
#
#   glioma settings : four-cohort intersection (TCGA n CGGA-693 n CGGA-325 n array-301),
#                     TCGA SD > 0.5, capped at 8,000 most variable  <- exactly Table 4's rule
#   breast settings : pairwise intersection under the same SD/cap rule
#
# The prognostic subset is whatever the split yields; because the discovery half is
# smaller than the full cohort, it will not be 6,512, and that count is reported.
#
# The question is narrow and pre-specified: do the SIGNS and the INTERVALS of the
# null-vs-real differences hold when the universe changes? Point estimates are
# expected to move, since a larger universe admits lower-variance, weaker-signal
# genes. Conclusions flipping would invalidate §6.9.
#
# Run from project root:  Rscript R/34_null_full_universe.R   (~25 min)
# Writes: results/null_full_universe.csv, results/null_full_universe.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)
set.seed(20260815)

out <- file("results/null_full_universe.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }
NCAP <- 8000; NHUB <- 100; NBOOT <- 2000

## ---------------- cohort loaders (identical to R/33) ----------------
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

## Table 4's universe: FOUR-cohort intersection, TCGA SD>0.5, cap 8,000 most variable
glioma_universe <- local({
  g <- Reduce(intersect, list(rownames(GL$e), rownames(C693$e), rownames(C325$e), rownames(C301$e)))
  sdT <- apply(GL$e[g, , drop = FALSE], 1, sd); g <- g[is.finite(sdT) & sdT > 0.5]
  if (length(g) > NCAP) g <- names(sort(sdT[g], decreasing = TRUE))[1:NCAP]
  g
})
say("Glioma universe (Table 4 rule: 4-cohort intersection, TCGA SD>0.5, cap %d): %d genes\n",
    NCAP, length(glioma_universe))

SETTINGS <- list(
  list(lab="glioma: TCGA->CGGA-693",  D=GL, R=C693, uni=glioma_universe),
  list(lab="glioma: TCGA->CGGA-325",  D=GL, R=C325, uni=glioma_universe),
  list(lab="glioma: TCGA->array-301", D=GL, R=C301, uni=glioma_universe),
  list(lab="breast: METABRIC->BRCA",  D=BR$MB, R=BR$TC, uni=NULL),
  # appended LAST so the RNG stream for the settings above is unchanged; this
  # direction was unevaluable at the 2,000-gene universe (26 prognostic genes)
  # and is re-tested here rather than assumed.
  list(lab="breast: BRCA->METABRIC",  D=BR$TC, R=BR$MB, uni=NULL))

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
coexpr <- function(eD, eR, hub){
  cD <- suppressWarnings(cor(t(eD), t(eD[hub, , drop = FALSE])))
  cR <- suppressWarnings(cor(t(eR), t(eR[hub, , drop = FALSE])))
  vapply(seq_len(nrow(eD)), function(i){
    a <- cD[i, ]; b <- cR[i, ]
    j <- which(hub == i); if (length(j)) { a <- a[-j]; b <- b[-j] }   # no self-correlation
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
  genes <- if (is.null(S$uni)) {
    g <- intersect(rownames(D$e), rownames(R$e))
    sdD <- apply(D$e[g, , drop = FALSE], 1, sd); g <- g[is.finite(sdD) & sdD > 0.5]
    if (length(g) > NCAP) g <- names(sort(sdD[g], decreasing = TRUE))[1:NCAP]
    g
  } else intersect(S$uni, rownames(R$e))
  eD <- D$e[genes, , drop = FALSE]; eR <- R$e[genes, , drop = FALSE]
  hub <- seq_len(min(NHUB, length(genes)))
  nR <- min(nrow(R$cl), floor(nrow(D$cl)/2))

  nullidx <- sample.int(nrow(D$cl), nR); discidx <- setdiff(seq_len(nrow(D$cl)), nullidx)
  realidx <- sample.int(nrow(R$cl), nR)
  say("\n=== %s ===\n  universe=%d genes, discovery n=%d, both arms n=%d\n",
      S$lab, length(genes), length(discidx), nR)

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
  say("  prognostic in the discovery half and scored in both arms: %d of %d\n", length(keep), length(genes))

  for (met in c("prod","absd","coex","eff")) {
    sN <- AN[[met]][keep]; sR2 <- AR[[met]][keep]
    lN <- AN$repl[keep];   lR <- AR$repl[keep]
    obs <- auc_fast(sR2, lR) - auc_fast(sN, lN)
    bs <- replicate(NBOOT, { i <- sample(length(keep), replace = TRUE)
      auc_fast(sR2[i], lR[i]) - auc_fast(sN[i], lN[i]) })
    ci <- quantile(bs, c(.025,.975), na.rm = TRUE)
    pb <- 2*min(mean(bs <= 0, na.rm=TRUE), mean(bs >= 0, na.rm=TRUE))
    rows[[length(rows)+1]] <- data.frame(setting=S$lab, metric=met, universe=length(genes),
      n_genes=length(keep), auc_null=auc_fast(sN,lN), auc_real=auc_fast(sR2,lR), diff=obs,
      lo=unname(ci[1]), hi=unname(ci[2]), p_boot=max(pb, 1/NBOOT))
    say("  %-6s null=%.3f real=%.3f diff=%+.3f (95%% CI %+.3f to %+.3f)\n",
        met, auc_fast(sN,lN), auc_fast(sR2,lR), obs, ci[1], ci[2])
  }
}
G <- do.call(rbind, rows)
for (met in unique(G$metric)) G$fdr[G$metric==met] <- p.adjust(G$p_boot[G$metric==met], "BH")
write.csv(G, "results/null_full_universe.csv", row.names = FALSE)

## ---------------- side-by-side against the 2,000-gene run ----------------
say("\n\n=== Robustness: 2,000-gene run (R/33) vs Table 4's universe (this script) ===\n")
S33 <- read.csv("results/gene_level_inference.csv")
say("%-26s %-6s %20s %20s %10s\n", "setting","metric","2,000 genes","Table 4 universe","sign held")
flip <- 0; tested <- 0
for (i in seq_len(nrow(G))) {
  a <- S33[S33$setting==G$setting[i] & S33$metric==G$metric[i], ]
  if (!nrow(a) || !is.finite(a$diff[1])) next
  tested <- tested + 1
  same <- sign(a$diff[1]) == sign(G$diff[i])
  concl_a <- if (a$lo[1] > 0) "+" else if (a$hi[1] < 0) "-" else "0"
  concl_b <- if (G$lo[i] > 0) "+" else if (G$hi[i] < 0) "-" else "0"
  if (concl_a != concl_b) flip <- flip + 1
  say("%-26s %-6s %+9.3f (%+.3f,%+.3f) %+9.3f (%+.3f,%+.3f) %10s\n",
      G$setting[i], G$metric[i], a$diff[1], a$lo[1], a$hi[1], G$diff[i], G$lo[i], G$hi[i],
      ifelse(same, "yes", "NO"))
}
say("\n  comparisons: %d ; conclusion (sign of interval) changed in %d\n", tested, flip)
gl <- G[G$metric=="absd" & grepl("glioma", G$setting), ]
gp <- G[G$metric=="prod" & grepl("glioma", G$setting), ]
say("  difference form, glioma pairs with CI excluding 0: %d of %d (2,000-gene run: 2 of 3)\n",
    sum(gl$lo > 0), nrow(gl))
say("  product form,    glioma pairs with CI excluding 0: %d of %d (2,000-gene run: 2 of 3)\n",
    sum(gp$lo > 0), nrow(gp))
ge <- G[G$metric=="eff", ]
say("  effect-size control negative in %d of %d settings (expected: all)\n", sum(ge$diff < 0), nrow(ge))
say("  product-form null AUC range, glioma: %.3f-%.3f (2,000-gene run: 0.830-0.931)\n",
    min(gp$auc_null), max(gp$auc_null))
say("  difference-form null AUC range, glioma: %.3f-%.3f (2,000-gene run: 0.521-0.579)\n",
    min(gl$auc_null), max(gl$auc_null))
close(out)
cat("\nwrote results/null_full_universe.csv and .txt\n")
