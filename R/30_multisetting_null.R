# 30_multisetting_null.R -- EARN (or refute) THE GENERALIZATION CLAIM.
#
# 28_matched_null.R showed that a same-population null cohort produces a HIGHER
# anchor AUC than a real replication cohort -- in ONE setting (TCGA -> CGGA-693,
# grade anchor, glioma). On that basis the manuscript asserted that the problem
# "applies to any per-gene cross-cohort quality metric validated this way". That
# assertion is not supported by one cohort pair and one metric. This script tests
# it directly, on two axes:
#
#   AXIS 1 -- SETTINGS. The matched-null calibration is repeated in five cohort
#     pairs spanning two tumour types and three platforms:
#       TCGA -> CGGA-693, TCGA -> CGGA-325, TCGA -> array-301 (glioma),
#       TCGA-BRCA -> METABRIC and METABRIC -> TCGA-BRCA (breast).
#
#   AXIS 2 -- METRICS. Five per-gene cross-cohort quality metrics, not just the
#     paper's anchor, each scored the same way:
#       anchor_product   r_D(x,anchor) * r_R(x,anchor)      [the paper's gate]
#       anchor_absdiff   -|r_D - r_R|                        [same anchor, different form]
#       level_concord    -|mean_D(z) - mean_R(z)|            [expression-level agreement]
#       range_concord    -|log(sd_D / sd_R)|                 [dynamic-range agreement]
#       coexpr_preserv   cor over 100 hub genes of x's       [co-expression preservation,
#                        correlation profile, D vs R          a standard cross-study QC idea]
#
#   AXIS 3 (control) -- RATE MATCHING. The null arm replicates at a much higher
#     rate than the real arm, and although AUC is base-rate invariant in
#     expectation, the two arms are not otherwise identical. We therefore also
#     report a rate-matched variant in which the null arm's replication threshold
#     is tightened until its replication rate equals the real arm's, so the two
#     AUCs are computed at the same base rate.
#
# If null >= real for EVERY metric in EVERY setting, the generalization is earned
# empirically. If it holds only for the anchor in glioma, the manuscript's claim
# must be narrowed to exactly that.
#
# Run from project root:  Rscript R/30_multisetting_null.R   (~25 min)
# Writes: results/multisetting_null.csv, results/multisetting_null.txt,
#         figures/multisetting_null.png

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)

set.seed(20260815)
out <- file("results/multisetting_null.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

NGENE <- 2000; NHUB <- 100; REPS <- 8

## ---------------- cohort loaders: all return list(e = genes x samples, cl) ----------------
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
  load_expr <- function(f){ m <- rd(f); rownames(m) <- m[[1]]; m[[1]] <- NULL; m <- as.matrix(m)
    ex <- m[, grep("_EXPR$", colnames(m)), drop = FALSE]; colnames(ex) <- sub("_EXPR$","",colnames(ex)); ex }
  lab <- function(f){ l <- rd(f); rownames(l) <- l[[1]]; l[[1]] <- NULL; l }
  mbE <- rbind(load_expr(file.path(dd,"training/moanna_training_data.tsv")),
               load_expr(file.path(dd,"training/moanna_validation_data.tsv")))
  tcE <- load_expr(file.path(dd,"testing/moanna_testing_data.tsv"))
  mbL <- rbind(lab(file.path(dd,"training/moanna_training_label.tsv")),
               lab(file.path(dd,"training/moanna_validation_label.tsv")))
  tcL <- lab(file.path(dd,"testing/moanna_testing_label.tsv"))
  mbc <- read.delim("data/brca/mb_clin.txt", comment.char="#", check.names=FALSE)
  mbc <- data.frame(id=mbc$PATIENT_ID, time=as.numeric(mbc$OS_MONTHS),
                    event=as.integer(grepl("DECEASED", mbc$OS_STATUS)),
                    age=suppressWarnings(as.numeric(mbc$AGE_AT_DIAGNOSIS)))
  tcc <- read.delim("data/brca/brca_tcga_clin.txt", comment.char="#", check.names=FALSE)
  tcc <- data.frame(id=tcc$PATIENT_ID, time=as.numeric(tcc$OS_MONTHS),
                    event=as.integer(grepl("DECEASED", tcc$OS_STATUS)),
                    age=suppressWarnings(as.numeric(tcc$AGE)))
  asm <- function(E, L, clin, idmap){
    ids <- rownames(E); cl <- clin[match(idmap(ids), clin$id), ]
    ok <- !is.na(cl$time) & cl$time > 0 & !is.na(cl$event) & !is.na(cl$age)
    list(e = t(E[ok, , drop = FALSE]),
         cl = data.frame(time = cl$time[ok], event = cl$event[ok], age = cl$age[ok],
                         anchor = as.numeric(L[ids, "ERStatus"])[ok]))
  }
  list(MB = asm(mbE, mbL, mbc, function(x) x),
       TC = asm(tcE, tcL, tcc, function(x) substr(x, 1, 12)))
}

say("loading cohorts...\n")
GL <- tcga_glioma()
C693 <- cgga_load("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt","data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
C325 <- cgga_load("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt","data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
C301 <- cgga_load("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt","data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt = FALSE)
BR <- breast_load()
okg <- !is.na(GL$cl$anchor) & !is.na(GL$cl$age)
GL$e <- GL$e[, okg]; GL$cl <- GL$cl[okg, ]

SETTINGS <- list(
  list(lab = "glioma: TCGA->CGGA-693",  D = GL,    R = C693, anchor = "grade"),
  list(lab = "glioma: TCGA->CGGA-325",  D = GL,    R = C325, anchor = "grade"),
  list(lab = "glioma: TCGA->array-301", D = GL,    R = C301, anchor = "grade"),
  list(lab = "breast: BRCA->METABRIC",  D = BR$TC, R = BR$MB, anchor = "ER"),
  list(lab = "breast: METABRIC->BRCA",  D = BR$MB, R = BR$TC, anchor = "ER")
)

## ---------------- per-gene statistics ----------------
cox_stats <- function(e, cl){
  S <- Surv(cl$time, cl$event); age <- cl$age
  b <- p <- numeric(nrow(e))
  for (i in seq_len(nrow(e))) {
    z <- as.numeric(scale(e[i, ]))
    f <- tryCatch(coxph(S ~ z + age), error = function(err) NULL)
    s <- if (is.null(f)) c(NA, NA) else summary(f)$coefficients["z", c("coef","Pr(>|z|)")]
    b[i] <- s[1]; p[i] <- s[2]
  }
  list(b = b, p = p)
}
qc_metrics <- function(eD, clD, eR, clR, hub){
  rD <- suppressWarnings(apply(eD, 1, function(x) cor(x, clD$anchor, method = "spearman", use = "complete.obs")))
  rR <- suppressWarnings(apply(eR, 1, function(x) cor(x, clR$anchor, method = "spearman", use = "complete.obs")))
  zD <- t(scale(t(eD))); zR <- t(scale(t(eR)))
  mD <- rowMeans(eD); mR <- rowMeans(eR)
  sD <- apply(eD, 1, sd);  sR <- apply(eR, 1, sd)
  # co-expression preservation over a fixed hub set
  cD <- suppressWarnings(cor(t(eD), t(eD[hub, , drop = FALSE])))
  cR <- suppressWarnings(cor(t(eR), t(eR[hub, , drop = FALSE])))
  coex <- vapply(seq_len(nrow(eD)), function(i) {
    a <- cD[i, ]; b <- cR[i, ]; ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 10) NA_real_ else suppressWarnings(cor(a[ok], b[ok]))
  }, numeric(1))
  list(anchor_product = rD * rR,
       anchor_absdiff = -abs(rD - rR),
       level_concord  = -abs(scale(mD)[, 1] - scale(mR)[, 1]),
       range_concord  = -abs(log(pmax(sD, 1e-9) / pmax(sR, 1e-9))),
       coexpr_preserv = coex)
}
auc <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

res <- list()
for (S in SETTINGS) {
  D <- S$D; R <- S$R
  genes <- intersect(rownames(D$e), rownames(R$e))
  sdD <- apply(D$e[genes, , drop = FALSE], 1, sd)
  genes <- names(sort(sdD[is.finite(sdD) & sdD > 0], decreasing = TRUE))
  genes <- genes[seq_len(min(NGENE, length(genes)))]
  eD_all <- D$e[genes, , drop = FALSE]; eR_all <- R$e[genes, , drop = FALSE]
  hub <- seq_len(min(NHUB, length(genes)))

  nR <- min(nrow(R$cl), floor(nrow(D$cl) / 2))
  say("\n=== %s ===\n", S$lab)
  say("  discovery cohort n=%d (events %d) | replication cohort n=%d (events %d)\n",
      nrow(D$cl), sum(D$cl$event), nrow(R$cl), sum(R$cl$event))
  say("  matched replication size for both arms: n=%d; discovery = %d\n", nR, nrow(D$cl) - nR)

  for (rep in seq_len(REPS)) {
    cat(sprintf("  %s rep %d/%d\n", S$lab, rep, REPS))
    nullidx <- sample.int(nrow(D$cl), nR)
    discidx <- setdiff(seq_len(nrow(D$cl)), nullidx)
    realidx <- sample.int(nrow(R$cl), nR)

    sD <- cox_stats(eD_all[, discidx, drop = FALSE], D$cl[discidx, ])
    prog <- p.adjust(sD$p, "BH") < 0.05
    if (sum(prog, na.rm = TRUE) < 40) next

    # both arms computed together so the null can be rate-matched to the real arm
    arms <- lapply(c("null","real"), function(arm){
      eA  <- if (arm == "null") eD_all[, nullidx, drop = FALSE] else eR_all[, realidx, drop = FALSE]
      clA <- if (arm == "null") D$cl[nullidx, ] else R$cl[realidx, ]
      sA <- cox_stats(eA, clA)
      M  <- qc_metrics(eD_all[, discidx, drop = FALSE], D$cl[discidx, ], eA, clA, hub)
      keep <- prog & !is.na(sA$b)
      list(arm = arm, sA = sA, M = M, keep = keep,
           agree = (sign(sD$b) == sign(sA$b))[keep], pR = sA$p[keep])
    })
    names(arms) <- c("null","real")
    rate_real <- mean(as.integer(arms$real$agree & arms$real$pR < 0.05))

    for (arm in c("null","real","null_ratematched")) {
      A <- if (arm == "real") arms$real else arms$null
      if (arm == "null_ratematched") {
        # tighten the null arm's p threshold until its replication rate matches the real arm
        cand <- sort(A$pR[A$agree])
        k <- max(1, min(length(cand), floor(rate_real * length(A$pR))))
        thr <- if (length(cand)) cand[k] else 0.05
        repl <- as.integer(A$agree & A$pR <= thr)
      } else {
        repl <- as.integer(A$agree & A$pR < 0.05)
      }
      row <- data.frame(setting = S$lab, anchor = S$anchor, rep = rep, arm = arm,
                        n_disc = length(discidx), n_repl = nR, n_prog = sum(A$keep),
                        replic_rate = mean(repl))
      for (nm in names(A$M)) row[[paste0("auc_", nm)]] <- auc(A$M[[nm]][A$keep], repl)
      row$auc_effsize <- auc(abs(sD$b)[A$keep], repl)
      res[[length(res)+1]] <- row
    }
  }
}
R0 <- do.call(rbind, res)
write.csv(R0, "results/multisetting_null.csv", row.names = FALSE)

METRICS <- c("anchor_product","anchor_absdiff","level_concord","range_concord","coexpr_preserv","effsize")
say("\n\n=== NULL vs REAL, five settings x six metrics (%d replicates each) ===\n", REPS)
say("  positive difference = the metric discriminates BETTER against a real cohort\n")
say("  than against a same-population null, i.e. it registers genuine discordance.\n\n")
say("%-26s %-16s %8s %8s %9s %9s\n", "setting", "metric", "null", "real", "real-null", "p")
summ <- list()
for (S in unique(R0$setting)) {
  for (m in METRICS) {
    col <- paste0("auc_", m)
    a <- R0[[col]][R0$setting == S & R0$arm == "null"]
    b <- R0[[col]][R0$setting == S & R0$arm == "real"]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 3 || length(b) < 3) next
    tt <- tryCatch(t.test(b, a), error = function(e) list(p.value = NA))
    say("%-26s %-16s %8.3f %8.3f %+9.3f %9.3g\n", S, m, mean(a), mean(b), mean(b) - mean(a), tt$p.value)
    summ[[length(summ)+1]] <- data.frame(setting = S, metric = m, null = mean(a), real = mean(b),
                                         diff = mean(b) - mean(a), p = tt$p.value)
  }
  say("\n")
}
SM <- do.call(rbind, summ)
write.csv(SM, "results/multisetting_null_summary.csv", row.names = FALSE)

say("=== VERDICT ===\n")
n_tot <- nrow(SM); n_neg <- sum(SM$diff < 0, na.rm = TRUE)
say("  metric x setting combinations tested: %d\n", n_tot)
say("  combinations where the null scores AT LEAST AS HIGH as the real cohort: %d (%.0f%%)\n",
    n_neg, 100 * n_neg / n_tot)
ap <- SM[SM$metric == "anchor_product", ]
say("  anchor_product specifically: null >= real in %d of %d settings\n", sum(ap$diff < 0), nrow(ap))
say("\n  A high proportion across metrics AND tumour types earns the general claim.\n")
say("  A result confined to the anchor in glioma requires the claim to be narrowed\n")
say("  to exactly that, which is how it will be written if this is what we see.\n")

suppressPackageStartupMessages({ library(ggplot2); library(tidyr) })
pl <- SM %>% filter(is.finite(diff))
p <- ggplot(pl, aes(metric, diff, fill = diff < 0)) +
  geom_col() + facet_wrap(~setting, ncol = 2) +
  geom_hline(yintercept = 0, linewidth = 0.4) + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "#1b9e77"), guide = "none") +
  labs(x = NULL, y = "AUC(real cohort) - AUC(same-population null)",
       title = "Does any cross-cohort quality metric register a real cohort as different?",
       subtitle = "orange = the metric scores no better on a real cohort than on a null where nothing is wrong") +
  theme_bw(base_size = 10)
ggsave("figures/multisetting_null.png", p, width = 9, height = 6, dpi = 300)
say("\nwrote figures/multisetting_null.png\n")
close(out)
cat("\nwrote results/multisetting_null.csv and .txt\n")
