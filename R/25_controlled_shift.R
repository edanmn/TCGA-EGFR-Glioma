# 25_controlled_shift.R -- CONTROLLED GROUND-TRUTH EXPERIMENT.
#
# Everything the paper says about what the positive-control gate detects is
# inferred from observational cohorts, where the truth is unknown: when TCGA and
# CGGA disagree we cannot see whether the cause was measurement or population
# structure. 24_cgga_composition.R showed the two are genuinely confusable for
# EGFR. This script removes the ambiguity by SPLITTING TCGA INTO TWO HALVES AND
# INJECTING EACH MECHANISM AT A KNOWN DOSE, so the gate can be scored against a
# ground truth we control.
#
# Two orthogonal factors, applied to the same patients and the same genes:
#
#   FACTOR 1 -- confounding-structure shift (delta).
#     Patients are assigned to half A or half B with a bias that makes grade and
#     IDH more tightly coupled in A than in B. delta=0 is a plain random split
#     (halves exchangeable); delta=0.4 is a strong decoupling. This reproduces,
#     under control, the exact difference measured between TCGA and CGGA.
#
#   FACTOR 2 -- measurement corruption (frac).
#     A randomly chosen 20% of genes have a fraction `frac` of their values
#     permuted across patients in half B ONLY. This mimics mis-quantification:
#     the gene's relationship to every covariate is degraded in one cohort while
#     the clinical data stay perfect. Ground truth = which genes were corrupted.
#
# The three numbers that matter:
#   (a) AUC at delta=0, frac=0 -- two exchangeable cohorts. NOT 0.5: genes with
#       strong anchor relationships have strong signal and replicate more often
#       regardless of any artifact. This is the floor the observational AUCs must
#       be read against, and the paper currently has no such floor.
#   (b) AUC as delta rises with frac=0 -- pure composition shift, zero measurement
#       problem. If AUC climbs, the gate provably fires on population structure.
#   (c) detection AUC for corrupted genes at delta=0 -- pure measurement problem,
#       zero composition shift. If it is high, the gate provably fires on
#       measurement too.
#
# Together (b) and (c) settle what the gate is: a detector of anchor-relationship
# discordance from EITHER cause, not a measurement-fidelity test.
#
# Run from project root:  Rscript R/25_controlled_shift.R   (~10 min)
# Writes: results/controlled_shift.csv, results/controlled_shift.txt,
#         figures/controlled_shift.png

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")

set.seed(20260815)   # this script is the one randomised analysis in the pipeline

out <- file("results/controlled_shift.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

## ---------- TCGA pooled, complete grade + IDH (the split needs both) ----------
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
say("Controlled-experiment cohort: n=%d patients, %d events\n", nrow(CL), sum(CL$event))

sdT <- apply(E, 1, sd)
NG <- 2000
genes <- names(sort(sdT[sdT > 0.5], decreasing = TRUE))[1:NG]
E <- E[genes, ]
say("Genes: %d (most variable, TCGA SD>0.5)\n", length(genes))
say("Baseline grade-IDH coupling in the full cohort: r=%.3f\n\n",
    cor(CL$grade, CL$idh, method = "spearman"))

## ---------- biased split: makes grade-IDH coupling stronger in A than in B ----------
# "typical" = the configuration that carries the canonical coupling
# (grade IV & IDH-wildtype, or lower grade & IDH-mutant).
typical <- (CL$grade == 4 & CL$idh == 0) | (CL$grade < 4 & CL$idh == 1)
split_halves <- function(delta){
  pA <- ifelse(typical, 0.5 + delta, 0.5 - delta)
  a <- runif(nrow(CL)) < pA
  list(A = which(a), B = which(!a))
}

## ---------- corruption: permute a fraction of each chosen gene's values ----------
corrupt <- function(mat, gidx, frac){
  if (frac <= 0 || !length(gidx)) return(mat)
  n <- ncol(mat)
  k <- max(2, round(frac * n))
  for (g in gidx) {
    i <- sample.int(n, k)
    mat[g, i] <- mat[g, sample(i)]
  }
  mat
}

## ---------- per-gene stats within one half ----------
half_stats <- function(idx, mat){
  cl <- CL[idx, ]; m <- mat[, idx, drop = FALSE]
  S <- Surv(cl$time, cl$event); age <- cl$age
  b <- numeric(nrow(m)); p <- numeric(nrow(m))
  for (i in seq_len(nrow(m))) {
    z <- as.numeric(scale(m[i, ]))
    fit <- tryCatch(coxph(S ~ z + age), error = function(e) NULL)
    if (is.null(fit)) { b[i] <- NA; p[i] <- NA } else {
      s <- summary(fit)$coefficients["z", c("coef", "Pr(>|z|)")]; b[i] <- s[1]; p[i] <- s[2]
    }
  }
  rg <- suppressWarnings(apply(m, 1, function(x) cor(x, cl$grade, method = "spearman")))
  data.frame(b = b, p = p, rg = rg)
}

auc <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---------- one replicate of one condition ----------
run_one <- function(delta, frac, rep){
  h <- split_halves(delta)
  A <- h$A; B <- h$B
  if (length(A) < 120 || length(B) < 120) return(NULL)
  # ground truth for factor 2
  ncorr <- round(0.20 * length(genes))
  gidx <- if (frac > 0) sample.int(length(genes), ncorr) else integer(0)
  EB <- corrupt(E, gidx, frac)

  sA <- half_stats(A, E)          # discovery half: never corrupted
  sB <- half_stats(B, EB)         # replication half: corrupted genes only
  prog <- p.adjust(sA$p, "BH") < 0.05 & !is.na(sB$b)
  if (sum(prog) < 50) return(NULL)
  replicated <- as.integer(sign(sA$b) == sign(sB$b) & sB$p < 0.05)
  qc <- sA$rg * sB$rg
  is_corrupt <- as.integer(seq_along(genes) %in% gidx)

  data.frame(
    delta = delta, frac = frac, rep = rep,
    nA = length(A), nB = length(B),
    coupling_A = cor(CL$grade[A], CL$idh[A], method = "spearman"),
    coupling_B = cor(CL$grade[B], CL$idh[B], method = "spearman"),
    n_prog = sum(prog),
    replic_rate = mean(replicated[prog]),
    # (b) does the gate predict replication?
    auc_replication = auc(qc[prog], replicated[prog]),
    # (c) does the gate identify the genes we actually corrupted?
    auc_detect_corrupt = if (frac > 0) auc(-qc, is_corrupt) else NA_real_,
    replic_corrupt = if (frac > 0) mean(replicated[prog & is_corrupt == 1]) else NA_real_,
    replic_clean   = mean(replicated[prog & is_corrupt == 0])
  )
}

## ---------- design ----------
REPS <- 6
design <- rbind(
  expand.grid(delta = c(0, 0.1, 0.2, 0.3, 0.4), frac = 0),        # composition sweep, no corruption
  expand.grid(delta = 0, frac = c(0.25, 0.50, 0.75))              # corruption sweep, no shift
)
say("Design: %d conditions x %d replicates\n", nrow(design), REPS)
say("  composition sweep: delta=0,0.1,0.2,0.3,0.4 at frac=0\n")
say("  corruption sweep:  frac=0.25,0.50,0.75 at delta=0\n\n")

res <- list(); ctr <- 0
for (i in seq_len(nrow(design))) {
  for (r in seq_len(REPS)) {
    ctr <- ctr + 1
    cat(sprintf("  [%d/%d] delta=%.1f frac=%.2f rep=%d\n", ctr, nrow(design) * REPS,
                design$delta[i], design$frac[i], r))
    res[[length(res) + 1]] <- run_one(design$delta[i], design$frac[i], r)
  }
}
R <- do.call(rbind, res)
write.csv(R, "results/controlled_shift.csv", row.names = FALSE)

## ---------- summarise ----------
sm <- R %>% group_by(delta, frac) %>%
  summarise(reps = dplyr::n(),
            d_coupling = mean(abs(coupling_A - coupling_B)),
            n_prog = mean(n_prog),
            replic = mean(replic_rate),
            auc_rep = mean(auc_replication), auc_rep_sd = sd(auc_replication),
            auc_det = mean(auc_detect_corrupt),
            rep_corrupt = mean(replic_corrupt), rep_clean = mean(replic_clean),
            .groups = "drop") %>% as.data.frame()

say("=== FACTOR 1: confounding-structure shift, no measurement problem ===\n")
say("%7s %12s %8s %8s %18s\n", "delta", "|coupling A-B|", "n prog", "replic", "AUC replication")
f1 <- sm[sm$frac == 0, ]
for (i in seq_len(nrow(f1)))
  say("%7.1f %12.3f %8.0f %8.3f %10.3f (sd %.3f)\n", f1$delta[i], f1$d_coupling[i],
      f1$n_prog[i], f1$replic[i], f1$auc_rep[i], f1$auc_rep_sd[i])
say("\n  delta=0 is the exchangeable-cohorts FLOOR: the AUC the anchor earns from\n")
say("  signal strength alone, with no artifact of any kind present.\n")

say("\n=== FACTOR 2: measurement corruption, no composition shift (delta=0) ===\n")
say("%7s %10s %10s %12s %18s\n", "frac", "replic", "AUC repl.", "corrupt/clean", "AUC detect corrupt")
f2 <- sm[sm$frac > 0, ]
for (i in seq_len(nrow(f2)))
  say("%7.2f %10.3f %10.3f %5.2f/%.2f %14.3f\n", f2$frac[i], f2$replic[i], f2$auc_rep[i],
      f2$rep_corrupt[i], f2$rep_clean[i], f2$auc_det[i])
say("\n  AUC-detect is scored against the genes we actually corrupted (ground truth).\n")

base <- f1$auc_rep[f1$delta == 0]
say("\n=== VERDICT ===\n")
say("  floor (exchangeable cohorts)      AUC = %.3f\n", base)
say("  max composition shift, no noise   AUC = %.3f  (%+.3f over floor)\n",
    max(f1$auc_rep), max(f1$auc_rep) - base)
say("  max corruption, no shift          AUC = %.3f  (%+.3f over floor)\n",
    max(f2$auc_rep), max(f2$auc_rep) - base)
say("  corrupted-gene detection          AUC = %.3f\n", max(f2$auc_det, na.rm = TRUE))
say("\n  If both factors raise the AUC above the floor, the gate is a detector of\n")
say("  anchor-relationship discordance from EITHER cause -- it cannot be described\n")
say("  as a measurement-fidelity test, and the observational AUCs in Table 4 must\n")
say("  be read against the floor rather than against 0.5.\n")

## ---------- figure ----------
suppressPackageStartupMessages(library(ggplot2))
p1 <- ggplot(R[R$frac == 0, ], aes(factor(delta), auc_replication)) +
  geom_boxplot(fill = "#9ecae1", outlier.size = 0.7) +
  geom_hline(yintercept = base, linetype = 2, colour = "grey40") +
  labs(x = expression(paste("confounding-structure shift ", delta)),
       y = "AUC for predicting replication",
       title = "Composition shift alone raises the gate's AUC",
       subtitle = "no measurement problem present; dashed line = exchangeable-cohort floor") +
  theme_bw(base_size = 11)
p2 <- ggplot(R[R$frac > 0, ], aes(factor(frac), auc_detect_corrupt)) +
  geom_boxplot(fill = "#fc9272", outlier.size = 0.7) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  labs(x = "fraction of values permuted in the replication half",
       y = "AUC for identifying corrupted genes",
       title = "Measurement corruption alone is also detected",
       subtitle = "no composition shift present; ground truth = injected corruption") +
  theme_bw(base_size = 11)
gt <- tryCatch({
  if (requireNamespace("patchwork", quietly = TRUE)) { library(patchwork); p1 | p2 } else NULL
}, error = function(e) NULL)
ggsave("figures/controlled_shift.png", if (is.null(gt)) p1 else gt,
       width = if (is.null(gt)) 5.5 else 10, height = 4.2, dpi = 300)
say("\nwrote figures/controlled_shift.png\n")
close(out)
cat("\nwrote results/controlled_shift.csv and .txt\n")
