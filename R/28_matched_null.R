# 28_matched_null.R -- MATCHED-SIZE NULL CALIBRATION.
#
# SUPERSEDED. Retained for provenance only. This script compares the two arms
# with a t-test across resampling replicates; R/33 shows that inference is
# invalid (replicates are resamples of one fixed cohort, so the p-value tracks
# the replicate count, not the evidence). The manuscript reports R/30's five-
# setting design with R/33's gene-level bootstrap intervals instead. Nothing in
# the paper depends on this script's p-values.
#
# 25_controlled_shift.R showed that splitting TCGA in half gives a grade-anchor
# AUC of 0.94 -- higher than any observational AUC the paper reports. But that
# comparison is not quite fair: the split halves have ~358 patients each, while
# the paper's discovery cohort has 716. Discovery power differs, so the two AUCs
# are not directly comparable and the floor could be an artefact of the design.
#
# This script removes that objection. ONE discovery set is used for everything,
# and it is compared against two replication cohorts OF THE SAME SIZE:
#
#   * NULL replication      -- a held-out, disjoint subsample of TCGA itself.
#                              Same population, same platform, same processing.
#                              No artefact of any kind can be present.
#   * REAL replication      -- CGGA-693, subsampled to exactly the same n.
#
# Discovery is the remaining TCGA patients, disjoint from the null replication
# set, so nothing is evaluated on its own data.
#
# If AUC(null) is as high as AUC(real), then the anchor score's apparent ability
# to "predict replication" carries no information about cross-cohort artefacts:
# it is a signal-strength statistic that behaves the same whether or not anything
# is wrong. That is the claim this script is built to test, and it is the single
# number that decides whether Table 4 supports its interpretation.
#
# Run from project root:  Rscript R/28_matched_null.R   (~5 min)
# Writes: results/matched_null.csv, results/matched_null.txt, figures/matched_null.png

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)

set.seed(20260815)
out <- file("results/matched_null.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

## ---------- cohorts ----------
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
  ok <- !is.na(cl$time) & cl$time > 0 & !is.na(cl$event)
  list(e = ex[, sel][, ok, drop = FALSE], clin = cl[ok, c("time","event","age","grade","idh")])
}

TC <- tcga_mat()
CG <- cgga_mat("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", "data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
okT <- !is.na(TC$clin$grade) & !is.na(TC$clin$age)
TE <- TC$e[, okT]; TCL <- TC$clin[okT, ]
okC <- !is.na(CG$clin$grade) & !is.na(CG$clin$age)
CE <- CG$e[, okC]; CCL <- CG$clin[okC, ]

genes <- intersect(rownames(TE), rownames(CE))
sdT <- apply(TE[genes, ], 1, sd)
genes <- names(sort(sdT[sdT > 0.5], decreasing = TRUE))[1:2000]
TE <- TE[genes, ]; CE <- CE[genes, ]
say("TCGA (graded, OS-evaluable): n=%d, events=%d\n", nrow(TCL), sum(TCL$event))
say("CGGA-693 (graded, OS-evaluable): n=%d, events=%d\n", nrow(CCL), sum(CCL$event))
say("Genes: %d\n\n", length(genes))

NREP <- nrow(CCL)                       # matched replication-cohort size
say("Matched replication-cohort size for BOTH arms: n=%d\n", NREP)
say("Discovery = the %d TCGA patients not in the null replication set.\n\n", nrow(TCL) - NREP)

stats <- function(e, cl){
  S <- Surv(cl$time, cl$event); age <- cl$age
  b <- p <- numeric(nrow(e))
  for (i in seq_len(nrow(e))) {
    z <- as.numeric(scale(e[i, ]))
    f <- tryCatch(coxph(S ~ z + age), error = function(err) NULL)
    s <- if (is.null(f)) c(NA, NA) else summary(f)$coefficients["z", c("coef","Pr(>|z|)")]
    b[i] <- s[1]; p[i] <- s[2]
  }
  rg <- suppressWarnings(apply(e, 1, function(x) cor(x, cl$grade, method = "spearman")))
  list(b = b, p = p, rg = rg)
}
auc <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

REPS <- 10
res <- list()
for (r in seq_len(REPS)) {
  cat(sprintf("  rep %d/%d\n", r, REPS))
  nullidx <- sample.int(nrow(TCL), NREP)          # held-out TCGA = null replication
  discidx <- setdiff(seq_len(nrow(TCL)), nullidx) # disjoint discovery
  cggidx  <- sample.int(nrow(CCL), NREP)          # CGGA subsampled to the same n

  sD <- stats(TE[, discidx, drop = FALSE], TCL[discidx, ])
  sN <- stats(TE[, nullidx, drop = FALSE], TCL[nullidx, ])
  sR <- stats(CE[, cggidx,  drop = FALSE], CCL[cggidx,  ])

  prog <- p.adjust(sD$p, "BH") < 0.05
  for (arm in c("null (held-out TCGA)", "real (CGGA-693)")) {
    s <- if (grepl("null", arm)) sN else sR
    keep <- prog & !is.na(s$b)
    if (sum(keep) < 50) next
    repl <- as.integer(sign(sD$b) == sign(s$b) & s$p < 0.05)[keep]
    res[[length(res)+1]] <- data.frame(
      rep = r, arm = arm, n_disc = length(discidx), n_repl = NREP,
      n_prog = sum(keep), replic_rate = mean(repl),
      auc_anchor = auc((sD$rg * s$rg)[keep], repl),
      auc_effsize = auc(abs(sD$b)[keep], repl))
  }
}
R <- do.call(rbind, res)
write.csv(R, "results/matched_null.csv", row.names = FALSE)

sm <- R %>% group_by(arm) %>%
  summarise(reps = dplyr::n(), n_prog = mean(n_prog),
            replic = mean(replic_rate), replic_sd = sd(replic_rate),
            anchor = mean(auc_anchor), anchor_sd = sd(auc_anchor),
            effsize = mean(auc_effsize), effsize_sd = sd(auc_effsize),
            .groups = "drop") %>% as.data.frame()

say("=== Matched-size comparison (%d replicates) ===\n\n", REPS)
say("%-22s %8s %16s %20s %20s\n", "replication cohort", "n prog", "replication rate",
    "AUC anchor", "AUC |beta_disc|")
for (i in seq_len(nrow(sm)))
  say("%-22s %8.0f %10.3f (%.3f) %12.3f (%.3f) %12.3f (%.3f)\n",
      sm$arm[i], sm$n_prog[i], sm$replic[i], sm$replic_sd[i],
      sm$anchor[i], sm$anchor_sd[i], sm$effsize[i], sm$effsize_sd[i])

nl <- sm[grepl("null", sm$arm), ]; rl <- sm[grepl("real", sm$arm), ]
d <- rl$anchor - nl$anchor
tt <- t.test(R$auc_anchor[grepl("real", R$arm)], R$auc_anchor[grepl("null", R$arm)])
say("\n  anchor AUC, real minus null: %+.3f (Welch t-test p=%.3g)\n", d, tt$p.value)
say("\n=== INTERPRETATION ===\n")
if (d >= 0.02) {
  say("  The anchor performs BETTER against a genuinely independent cohort than\n")
  say("  against a same-population null of identical size. The extra discrimination\n")
  say("  is attributable to cross-cohort discordance, which is what the paper claims.\n")
} else if (d <= -0.02) {
  say("  The anchor performs WORSE against the real cohort than against a\n")
  say("  same-population null of identical size. Its apparent ability to predict\n")
  say("  replication is therefore NOT evidence that it detects cross-cohort\n")
  say("  artefacts: an AUC of this magnitude arises with no artefact present at all.\n")
  say("  Table 4's AUCs must be reported against this null, not against 0.5.\n")
} else {
  say("  The anchor performs the SAME against a real independent cohort as against\n")
  say("  a same-population null of identical size. Its apparent ability to predict\n")
  say("  replication therefore carries no information about cross-cohort artefacts.\n")
  say("  Table 4's AUCs must be reported against this null, not against 0.5.\n")
}
say("\n  Note the replication RATE does separate the arms (%.3f null vs %.3f real):\n",
    nl$replic, rl$replic)
say("  the cohorts genuinely differ. It is the ANCHOR'S AUC that fails to register it.\n")

suppressPackageStartupMessages(library(ggplot2))
p <- ggplot(R, aes(arm, auc_anchor, fill = arm)) +
  geom_boxplot(width = 0.5, show.legend = FALSE) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.7, show.legend = FALSE) +
  labs(x = NULL, y = "AUC for predicting replication",
       title = "The anchor score cannot tell a real cohort from a same-population null",
       subtitle = sprintf("identical discovery set, identical replication-cohort size (n=%d), %d replicates", NREP, REPS)) +
  theme_bw(base_size = 11)
ggsave("figures/matched_null.png", p, width = 6.4, height = 4.2, dpi = 300)
say("\nwrote figures/matched_null.png\n")
close(out)
cat("\nwrote results/matched_null.csv and .txt\n")
