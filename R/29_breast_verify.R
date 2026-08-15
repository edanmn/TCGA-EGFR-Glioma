# 29_breast_verify.R -- close the reproducibility gap in section 6.9.
#
# Two numbers in section 6.9 are produced by NO committed script:
#   (i)  "The basal anchor reached at most 0.54"
#   (ii) "genome-wide agreement of expression-ER correlations r=0.84"
# 17_breast_generalization.R computes the FORWARD direction only, and its basal
# AUC is 0.474. The reverse-direction ER AUC (0.556) survives only in a stray log
# (results/_breast_ci.log); the reverse basal AUC and the r=0.84 concordance
# figure appear nowhere. Both are cited in the paper.
#
# This script recomputes every section 6.9 number in both directions from the
# source data, so all of them are reproducible from the pipeline.
#
# Run from project root:  Rscript R/29_breast_verify.R
# Writes: results/breast_verify.csv, results/breast_verify.txt

suppressPackageStartupMessages({ library(dplyr); library(survival) })
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)
dd <- "data/brca/moanna_data"

out <- file("results/breast_verify.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

load_expr <- function(f){
  m <- rd(f); rownames(m) <- m[[1]]; m[[1]] <- NULL; m <- as.matrix(m)
  ex <- m[, grep("_EXPR$", colnames(m)), drop = FALSE]
  colnames(ex) <- sub("_EXPR$", "", colnames(ex)); ex
}
lab <- function(f){ l <- rd(f); rownames(l) <- l[[1]]; l[[1]] <- NULL; l }

mbE <- rbind(load_expr(file.path(dd, "training/moanna_training_data.tsv")),
             load_expr(file.path(dd, "training/moanna_validation_data.tsv")))
tcE <- load_expr(file.path(dd, "testing/moanna_testing_data.tsv"))
mbL <- rbind(lab(file.path(dd, "training/moanna_training_label.tsv")),
             lab(file.path(dd, "training/moanna_validation_label.tsv")))
tcL <- lab(file.path(dd, "testing/moanna_testing_label.tsv"))

mbc <- read.delim("data/brca/mb_clin.txt", comment.char = "#", check.names = FALSE)
mbc <- data.frame(id = mbc$PATIENT_ID, time = as.numeric(mbc$OS_MONTHS),
                  event = as.integer(grepl("DECEASED", mbc$OS_STATUS)),
                  age = suppressWarnings(as.numeric(mbc$AGE_AT_DIAGNOSIS)))
tcc <- read.delim("data/brca/brca_tcga_clin.txt", comment.char = "#", check.names = FALSE)
tcc <- data.frame(id = tcc$PATIENT_ID, time = as.numeric(tcc$OS_MONTHS),
                  event = as.integer(grepl("DECEASED", tcc$OS_STATUS)),
                  age = suppressWarnings(as.numeric(tcc$AGE)))

assemble <- function(E, L, clin, idmap){
  ids <- rownames(E); pid <- idmap(ids); cl <- clin[match(pid, clin$id), ]
  d <- list(e = E, time = cl$time, event = cl$event, age = cl$age,
            er = as.numeric(L[ids, "ERStatus"]), basal = as.numeric(L[ids, "BasalNonBasal"]))
  ok <- !is.na(d$time) & d$time > 0 & !is.na(d$event)
  for (k in c("time","event","age","er","basal")) d[[k]] <- d[[k]][ok]
  d$e <- d$e[ok, , drop = FALSE]; d
}
MB <- assemble(mbE, mbL, mbc, function(x) x)
TC <- assemble(tcE, tcL, tcc, function(x) substr(x, 1, 12))
say("TCGA-BRCA: n=%d events=%d | METABRIC: n=%d events=%d\n",
    length(TC$time), sum(TC$event), length(MB$time), sum(MB$event))

genes <- intersect(colnames(TC$e), colnames(MB$e))
sdT <- apply(TC$e[, genes], 2, sd); genes <- genes[is.finite(sdT) & sdT > 0]
say("shared genes: %d\n\n", length(genes))

stats <- function(D, gs){
  S <- Surv(D$time, D$event); age <- D$age
  b <- p <- numeric(length(gs))
  for (i in seq_along(gs)) {
    z <- as.numeric(scale(D$e[, gs[i]]))
    f <- tryCatch(coxph(S ~ z + age), error = function(e) NULL)
    s <- if (is.null(f)) c(NA, NA) else summary(f)$coefficients["z", c("coef","Pr(>|z|)")]
    b[i] <- s[1]; p[i] <- s[2]
  }
  r_er    <- suppressWarnings(apply(D$e[, gs, drop = FALSE], 2, function(x) cor(x, D$er,    method = "spearman", use = "complete.obs")))
  r_basal <- suppressWarnings(apply(D$e[, gs, drop = FALSE], 2, function(x) cor(x, D$basal, method = "spearman", use = "complete.obs")))
  data.frame(gene = gs, b = b, p = p, r_er = r_er, r_basal = r_basal)
}
say("computing per-gene statistics (both cohorts)...\n")
sT <- stats(TC, genes); sM <- stats(MB, genes)

auc <- function(score, label){
  ok <- !is.na(score) & !is.na(label); score <- score[ok]; label <- label[ok]
  n1 <- sum(label == 1); n0 <- sum(label == 0); if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
aucCI <- function(score, label, B = 600){
  ok <- !is.na(score) & !is.na(label); s <- score[ok]; l <- label[ok]
  a <- auc(s, l); if (is.na(a)) return(c(a = NA, lo = NA, hi = NA))
  bs <- replicate(B, { i <- sample(length(s), replace = TRUE); auc(s[i], l[i]) })
  c(a = a, lo = unname(quantile(bs, .025, na.rm = TRUE)), hi = unname(quantile(bs, .975, na.rm = TRUE)))
}

## ---- (i) the r=0.84 concordance claim ----
say("\n=== Genome-wide agreement of expression-anchor correlations ===\n")
r_er_cor    <- cor(sT$r_er,    sM$r_er,    use = "complete.obs")
r_basal_cor <- cor(sT$r_basal, sM$r_basal, use = "complete.obs")
say("  cor of expression-ER correlations,    TCGA-BRCA vs METABRIC: r=%.3f   [paper: 0.84]\n", r_er_cor)
say("  cor of expression-basal correlations, TCGA-BRCA vs METABRIC: r=%.3f\n", r_basal_cor)

## ---- (ii) both directions, both anchors ----
say("\n=== Both directions, both anchors, with effect-size baseline ===\n")
rows <- list()
run <- function(sD, sR, lab){
  prog <- p.adjust(sD$p, "BH") < 0.05 & !is.na(sR$b)
  repl <- as.integer(sign(sD$b) == sign(sR$b) & sR$p < 0.05)[prog]
  er <- aucCI((sD$r_er * sR$r_er)[prog], repl)
  ba <- aucCI((sD$r_basal * sR$r_basal)[prog], repl)
  ef <- aucCI(abs(sD$b)[prog], repl)
  say("  %-26s n=%4d  replic=%.3f\n", lab, sum(prog), mean(repl))
  say("      ER anchor      %.3f (%.3f-%.3f)\n", er["a"], er["lo"], er["hi"])
  say("      basal anchor   %.3f (%.3f-%.3f)\n", ba["a"], ba["lo"], ba["hi"])
  say("      |beta_disc|    %.3f (%.3f-%.3f)\n", ef["a"], ef["lo"], ef["hi"])
  data.frame(direction = lab, n = sum(prog), replic = mean(repl),
             auc_er = er["a"], er_lo = er["lo"], er_hi = er["hi"],
             auc_basal = ba["a"], basal_lo = ba["lo"], basal_hi = ba["hi"],
             auc_effsize = ef["a"], row.names = NULL)
}
rows[[1]] <- run(sT, sM, "TCGA-BRCA -> METABRIC")
rows[[2]] <- run(sM, sT, "METABRIC -> TCGA-BRCA")
B <- do.call(rbind, rows)
B$anchor_concordance_er <- r_er_cor
write.csv(B, "results/breast_verify.csv", row.names = FALSE)

say("\n=== Verification against the paper's section 6.9 ===\n")
chk <- function(lab, got, paper, tol) say("  [%s] %-46s got=%.3f  paper=%s\n",
  ifelse(abs(got - paper) <= tol, "OK  ", "FAIL"), lab, got, format(paper))
chk("forward ER AUC",  B$auc_er[1],    0.60, 0.01)
chk("reverse ER AUC",  B$auc_er[2],    0.56, 0.01)
chk("max basal AUC",   max(B$auc_basal), 0.54, 0.01)
chk("ER concordance r", r_er_cor,      0.84, 0.01)
close(out)
cat("\nwrote results/breast_verify.csv and .txt\n")
