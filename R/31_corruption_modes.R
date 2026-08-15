# 31_corruption_modes.R -- WHICH failure modes does the gate actually catch?
#
# 25_controlled_shift.R injected exactly one kind of corruption: permutation of a
# fraction of a gene's values. That is a single failure mode, and the manuscript
# concedes as much. It is also NOT the failure mode EGFR exhibits in CGGA.
# 27_effectsize_baseline.R found EGFR's signature there to be SIGNAL COMPRESSION:
# the share of variance explained by the (grade, IDH) strata collapses from 0.106
# to 0.024-0.045 while the WITHIN-stratum SD is unchanged (1.88 vs ~1.84).
# Permutation does not reproduce that -- it degrades everything at once.
#
# So the controlled experiment, as run, cannot tell us whether the gate would have
# caught the thing we claim it caught. This script fixes that by injecting five
# distinct, mechanistically different failure modes at matched doses:
#
#   permute      shuffle a fraction of values across patients   (the original)
#   noise        x + N(0, k*sd)                                 (inflates within-stratum var)
#   compress     shrink the BETWEEN-stratum component by lambda (the EGFR signature:
#                x' = mu + lambda*(mu_s - mu) + (x - mu_s)       between-var falls,
#                                                                within-var preserved)
#   floor        censor values below a quantile to that quantile (detection limit / saturation)
#   contaminate  (1-w)*x + w*(an unrelated gene)                 (paralog / mixed feature)
#
# For each we record the AUC for identifying the genes we actually corrupted, the
# replication rate of corrupted vs clean genes, and the resulting variance
# signature -- so the mode that reproduces EGFR's observed signature can be
# identified and its detectability read off directly.
#
# No composition shift is applied (delta = 0) so the only thing varying is the
# measurement failure.
#
# Run from project root:  Rscript R/31_corruption_modes.R   (~15 min)
# Writes: results/corruption_modes.csv, results/corruption_modes.txt,
#         figures/corruption_modes.png

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
set.seed(20260815)

out <- file("results/corruption_modes.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

## ---------- TCGA pooled with complete grade + IDH (strata needed for `compress`) ----------
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
TC <- tcga_mat()
ok <- !is.na(TC$clin$grade) & !is.na(TC$clin$idh) & !is.na(TC$clin$age)
E <- TC$e[, ok]; CL <- TC$clin[ok, ]
sdT <- apply(E, 1, sd); genes <- names(sort(sdT[sdT > 0.5], decreasing = TRUE))[1:2000]
E <- E[genes, ]
STRAT <- paste(CL$grade, CL$idh, sep = "_")
say("Cohort: n=%d patients, %d events, %d genes, %d (grade,IDH) strata\n\n",
    nrow(CL), sum(CL$event), length(genes), length(unique(STRAT)))

## ---------- variance signature: share of variance explained by the strata ----------
var_signature <- function(x, s){
  mu <- tapply(x, s, mean); w <- table(s) / length(s); lev <- names(w)
  vb <- sum(w * (mu[lev] - mean(x))^2)
  vw <- sum(w * tapply(x, s, var)[lev], na.rm = TRUE)
  c(between = vb, within = vw, frac = vb / (vb + vw))
}

## ---------- the five corruption modes, applied to half B only ----------
apply_mode <- function(mat, gidx, mode, dose, s){
  if (!length(gidx)) return(mat)
  n <- ncol(mat)
  for (g in gidx) {
    x <- mat[g, ]
    mat[g, ] <- switch(mode,
      permute = { k <- max(2, round(dose * n)); i <- sample.int(n, k); y <- x; y[i] <- x[sample(i)]; y },
      noise   = x + rnorm(n, 0, dose * sd(x)),
      compress = {                              # EGFR-like: shrink between-stratum signal
        mu <- mean(x); ms <- tapply(x, s, mean)[s]
        mu + dose * (ms - mu) + (x - ms)
      },
      floor   = { q <- quantile(x, dose); pmax(x, q) },
      contaminate = {                            # mix in an unrelated gene, scale-matched
        d <- sample(setdiff(seq_len(nrow(mat)), g), 1)
        y <- mat[d, ]; y <- (y - mean(y)) / sd(y) * sd(x) + mean(x)
        (1 - dose) * x + dose * y
      },
      x)
  }
  mat
}

half_stats <- function(idx, mat){
  cl <- CL[idx, ]; m <- mat[, idx, drop = FALSE]
  S <- Surv(cl$time, cl$event); age <- cl$age
  b <- p <- numeric(nrow(m))
  for (i in seq_len(nrow(m))) {
    z <- as.numeric(scale(m[i, ]))
    f <- tryCatch(coxph(S ~ z + age), error = function(e) NULL)
    s <- if (is.null(f)) c(NA, NA) else summary(f)$coefficients["z", c("coef","Pr(>|z|)")]
    b[i] <- s[1]; p[i] <- s[2]
  }
  rg <- suppressWarnings(apply(m, 1, function(x) cor(x, cl$grade, method = "spearman")))
  list(b = b, p = p, rg = rg)
}
auc <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## dose grids chosen so each mode spans mild -> severe
DESIGN <- rbind(
  expand.grid(mode = "permute",     dose = c(0.25, 0.50, 0.75), stringsAsFactors = FALSE),
  expand.grid(mode = "noise",       dose = c(0.50, 1.00, 2.00), stringsAsFactors = FALSE),
  expand.grid(mode = "compress",    dose = c(0.75, 0.50, 0.25), stringsAsFactors = FALSE),
  expand.grid(mode = "floor",       dose = c(0.25, 0.50, 0.75), stringsAsFactors = FALSE),
  expand.grid(mode = "contaminate", dose = c(0.25, 0.50, 0.75), stringsAsFactors = FALSE)
)
REPS <- 4
say("Design: %d mode x dose conditions, %d replicates, 20%% of genes corrupted in half B\n\n",
    nrow(DESIGN), REPS)

res <- list(); ctr <- 0
for (i in seq_len(nrow(DESIGN))) {
  for (r in seq_len(REPS)) {
    ctr <- ctr + 1
    cat(sprintf("  [%d/%d] %s dose=%.2f rep=%d\n", ctr, nrow(DESIGN)*REPS,
                DESIGN$mode[i], DESIGN$dose[i], r))
    a <- runif(nrow(CL)) < 0.5
    A <- which(a); B <- which(!a)
    if (length(A) < 120 || length(B) < 120) next
    gidx <- sample.int(length(genes), round(0.20 * length(genes)))
    EB <- apply_mode(E, gidx, DESIGN$mode[i], DESIGN$dose[i], STRAT)

    sA <- half_stats(A, E); sB <- half_stats(B, EB)
    prog <- p.adjust(sA$p, "BH") < 0.05 & !is.na(sB$b)
    if (sum(prog) < 50) next
    repl <- as.integer(sign(sA$b) == sign(sB$b) & sB$p < 0.05)
    qc <- sA$rg * sB$rg
    iscorr <- as.integer(seq_along(genes) %in% gidx)

    # variance signature of corrupted genes in half B (vs their uncorrupted selves in A)
    sB_str <- STRAT[B]
    fr_corrupt <- mean(vapply(gidx[1:min(200, length(gidx))],
      function(g) var_signature(EB[g, B], sB_str)["frac"], numeric(1)), na.rm = TRUE)
    fr_clean_A <- mean(vapply(gidx[1:min(200, length(gidx))],
      function(g) var_signature(E[g, B], sB_str)["frac"], numeric(1)), na.rm = TRUE)

    res[[length(res)+1]] <- data.frame(
      mode = DESIGN$mode[i], dose = DESIGN$dose[i], rep = r,
      n_prog = sum(prog),
      replic_all = mean(repl[prog]),
      replic_corrupt = mean(repl[prog & iscorr == 1]),
      replic_clean   = mean(repl[prog & iscorr == 0]),
      auc_detect = auc(-qc, iscorr),
      var_frac_corrupted = fr_corrupt, var_frac_uncorrupted = fr_clean_A)
  }
}
R <- do.call(rbind, res)
write.csv(R, "results/corruption_modes.csv", row.names = FALSE)

sm <- R %>% group_by(mode, dose) %>%
  summarise(reps = dplyr::n(), replic_corrupt = mean(replic_corrupt),
            replic_clean = mean(replic_clean),
            auc_sd = sd(auc_detect),                 # computed BEFORE auc_detect is overwritten
            auc_detect = mean(auc_detect), var_frac = mean(var_frac_corrupted),
            var_frac_pre = mean(var_frac_uncorrupted), .groups = "drop") %>% as.data.frame()

say("=== Detectability by failure mode (ground truth = genes actually corrupted) ===\n\n")
say("%-13s %6s %11s %11s %20s %14s\n", "mode", "dose", "repl corrupt", "repl clean",
    "AUC detect (sd)", "var frac after")
for (i in seq_len(nrow(sm)))
  say("%-13s %6.2f %11.3f %11.3f %12.3f (%.3f) %14.3f\n",
      sm$mode[i], sm$dose[i], sm$replic_corrupt[i], sm$replic_clean[i],
      sm$auc_detect[i], sm$auc_sd[i], sm$var_frac[i])
say("\n  'var frac after' = share of a corrupted gene's variance explained by the\n")
say("  (grade, IDH) strata after corruption; uncorrupted baseline = %.3f.\n", mean(sm$var_frac_pre))

## ---------- which mode reproduces EGFR's observed signature? ----------
say("\n=== Which mode reproduces EGFR's signature in CGGA? ===\n")
say("  Observed for EGFR: variance share 0.106 (TCGA) -> 0.024-0.045 (CGGA),\n")
say("  a ratio of 0.23-0.42, with WITHIN-stratum SD essentially unchanged.\n\n")
sm$ratio <- sm$var_frac / sm$var_frac_pre
say("%-13s %6s %10s %14s\n", "mode", "dose", "var ratio", "AUC detect")
for (i in order(abs(sm$ratio - 0.32)))
  say("%-13s %6.2f %10.3f %14.3f%s\n", sm$mode[i], sm$dose[i], sm$ratio[i], sm$auc_detect[i],
      ifelse(sm$ratio[i] >= 0.23 & sm$ratio[i] <= 0.42, "   <-- matches EGFR", ""))

best <- sm[which.min(abs(sm$ratio - 0.32)), ]
say("\n  Closest match: %s at dose %.2f (variance ratio %.3f), detected with AUC %.3f.\n",
    best$mode, best$dose, best$ratio, best$auc_detect)
say("  This is the number that says whether the gate could have caught the EGFR\n")
say("  failure -- as opposed to catching permutation, which is not what EGFR shows.\n")

suppressPackageStartupMessages(library(ggplot2))
p <- ggplot(R, aes(factor(dose), auc_detect, fill = mode)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  facet_wrap(~mode, scales = "free_x", nrow = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  labs(x = "dose (mode-specific)", y = "AUC for identifying corrupted genes",
       title = "The gate detects some measurement failure modes and not others",
       subtitle = "ground truth = the genes actually corrupted; dashed line = chance") +
  theme_bw(base_size = 10)
ggsave("figures/corruption_modes.png", p, width = 10, height = 3.6, dpi = 300)
say("\nwrote figures/corruption_modes.png\n")
close(out)
cat("\nwrote results/corruption_modes.csv and .txt\n")
