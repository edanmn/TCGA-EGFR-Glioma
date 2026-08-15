# 24_cgga_composition.R -- test the leading BIOLOGICAL alternative to the
# "measurement artifact" interpretation of CGGA's EGFR positive-control failure.
#
# The objection: CGGA is enriched for younger patients and secondary glioblastoma,
# so its grade IV stratum carries a far higher IDH-mutant fraction than TCGA's.
# EGFR is high in IDH-wildtype tumours; if CGGA grade IV is diluted with IDH-mutant
# (secondary) GBM, EGFR would fail to rise with grade for reasons of COHORT
# COMPOSITION -- genuine biology of a different patient population, not a
# measurement problem. Under that explanation the flat r(EGFR, grade) is real and
# the paper's "artifact" framing is wrong.
#
# The composition explanation makes a sharp, falsifiable prediction: it operates
# ENTIRELY THROUGH IDH mixing, so WITHIN a single IDH stratum it has nothing left
# to explain, and EGFR should rise with grade in CGGA just as it does in TCGA.
#
#   (A) IDH-mutant fraction by grade, every cohort -- how large is the imbalance?
#   (B) r(EGFR, grade) computed WITHIN IDH-wildtype and WITHIN IDH-mutant tumours,
#       with Fisher r-to-z vs the matching TCGA stratum. <-- the decisive test
#   (C) EGFR by IDH status within grade IV only -- the secondary-GBM stratum where
#       the composition account is most concentrated.
#   (D) Panel-wide check: is the TCGA-vs-CGGA attenuation in Table S4 uniform
#       (points to a cohort-level/range-restriction effect) or is EGFR an outlier
#       (points to a gene-specific measurement problem)?
#
# Run from project root:  Rscript R/24_cgga_composition.R
# Writes: results/cgga_composition.csv, results/cgga_composition.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)

out <- file("results/cgga_composition.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }

## ---------- cohort assembly (matches 20_control_table.R: grade coded 2/3/4) ----------
tcga <- function(){
  one <- function(project, gradeIV = FALSE){
    se <- load_se(project); cd <- as.data.frame(colData(se)); k <- substr(colnames(se), 14, 15) == "01"
    e <- log2(assay(se, "tpm_unstrand")[, k] + 1); rownames(e) <- rowData(se)$gene_name
    e <- e[!duplicated(rownames(e)), ]
    b <- .clinical_df(se, project); b <- b[match(colnames(e), b$barcode), ]
    gr <- as.character(cd$paper_Grade[k])
    g <- if (gradeIV) ifelse(grepl("G4|IV", gr), 4, NA) else ifelse(grepl("G3|III", gr), 3, ifelse(grepl("G2| II", gr), 2, NA))
    idh <- as.character(cd$paper_IDH.status[k])
    keep <- !duplicated(b$patient)   # all graded primary tumours, as in Table 2
    list(e = e[, keep, drop = FALSE], grade = g[keep],
         idh = ifelse(idh[keep] == "Mutant", 1, ifelse(idh[keep] == "WT", 0, NA)))
  }
  L <- one("TCGA-LGG"); G <- one("TCGA-GBM", TRUE)
  cg <- intersect(rownames(L$e), rownames(G$e))
  list(e = cbind(L$e[cg, ], G$e[cg, ]), grade = c(L$grade, G$grade), idh = c(L$idh, G$idh))
}
cgga <- function(gfile, cfile, cols, logt = TRUE){
  ex <- rd(gfile); rownames(ex) <- ex[[1]]; ex[[1]] <- NULL
  ex <- if (logt) log2(as.matrix(ex) + 1) else as.matrix(ex)
  cl <- read.delim(cfile, check.names = FALSE); names(cl)[cols] <- c("id","prs","grade","age","os","censor","idh")
  cl <- cl %>% filter(prs == "Primary") %>%
    mutate(grade = as.numeric(factor(grade, levels = c("WHO II","WHO III","WHO IV"))) + 1,
           idh = ifelse(idh == "Mutant", 1, ifelse(idh == "Wildtype", 0, NA)))
  sel <- intersect(colnames(ex), cl$id); cl <- cl[match(sel, cl$id), ]
  list(e = ex[, sel, drop = FALSE], grade = cl$grade, idh = cl$idh)
}

TC   <- tcga()
D693 <- cgga("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", "data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
D325 <- cgga("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", "data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
D301 <- cgga("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt", "data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt = FALSE)
COH  <- list(`TCGA` = TC, `CGGA-693` = D693, `CGGA-325` = D325, `CGGA-301` = D301)

fisher_z <- function(r1, n1, r2, n2){
  if (any(is.na(c(r1, n1, r2, n2))) || n1 < 4 || n2 < 4) return(NA_real_)
  2 * pnorm(-abs((atanh(r1) - atanh(r2)) / sqrt(1/(n1-3) + 1/(n2-3))))
}

## ---------- (A) IDH-mutant fraction by grade ----------
say("=== (A) IDH-mutant fraction by WHO grade (primary tumours with grade + IDH) ===\n")
say("    If CGGA grade IV is diluted with IDH-mutant (secondary) GBM, the grade-IV\n")
say("    column is where it must show.\n\n")
say("  %-10s %-16s %-16s %-16s\n", "cohort", "grade II", "grade III", "grade IV")
compA <- do.call(rbind, lapply(names(COH), function(k){
  M <- COH[[k]]; ok <- !is.na(M$grade) & !is.na(M$idh)
  cells <- sapply(2:4, function(g){
    s <- M$idh[ok & M$grade == g]
    if (!length(s)) return(c(NA, 0))
    c(mean(s == 1), length(s))
  })
  say("  %-10s %-16s %-16s %-16s\n", k,
      sprintf("%.0f%% (n=%d)", 100*cells[1,1], cells[2,1]),
      sprintf("%.0f%% (n=%d)", 100*cells[1,2], cells[2,2]),
      sprintf("%.0f%% (n=%d)", 100*cells[1,3], cells[2,3]))
  data.frame(cohort = k, analysis = "idh_mutant_fraction",
             grade2 = cells[1,1], grade3 = cells[1,2], grade4 = cells[1,3],
             n2 = cells[2,1], n3 = cells[2,2], n4 = cells[2,3])
}))

## ---------- (B) r(EGFR, grade) WITHIN each IDH stratum -- the decisive test ----------
say("\n=== (B) Spearman r(EGFR, grade) WITHIN IDH strata ===\n")
say("    Composition prediction: within a single IDH stratum the mixing explanation\n")
say("    is exhausted, so CGGA should recover TCGA's positive EGFR-grade slope.\n")
say("    Measurement-artifact prediction: CGGA stays flat/negative even here.\n\n")
strat_r <- function(M, stratum){
  x <- M$e["EGFR", ]
  ok <- !is.na(x) & !is.na(M$grade) & !is.na(M$idh) & M$idh == stratum
  if (sum(ok) < 10) return(c(r = NA, n = sum(ok)))
  c(r = suppressWarnings(cor(x[ok], M$grade[ok], method = "spearman")), n = sum(ok))
}
say("  %-10s %-24s %-24s\n", "cohort", "IDH-wildtype", "IDH-mutant")
refs <- list(wt = strat_r(TC, 0), mut = strat_r(TC, 1))
compB <- do.call(rbind, lapply(names(COH), function(k){
  M <- COH[[k]]; w <- strat_r(M, 0); m <- strat_r(M, 1)
  pw <- if (k == "TCGA") NA else fisher_z(w["r"], w["n"], refs$wt["r"], refs$wt["n"])
  pm <- if (k == "TCGA") NA else fisher_z(m["r"], m["n"], refs$mut["r"], refs$mut["n"])
  say("  %-10s %-24s %-24s\n", k,
      sprintf("r=%+0.3f n=%3d %s", w["r"], w["n"], if (is.na(pw)) "(ref)   " else sprintf("p=%.1e", pw)),
      sprintf("r=%+0.3f n=%3d %s", m["r"], m["n"], if (is.na(pm)) "(ref)   " else sprintf("p=%.1e", pm)))
  data.frame(cohort = k, analysis = "r_egfr_grade_within_idh",
             r_idhwt = unname(w["r"]), n_idhwt = unname(w["n"]), p_vs_TCGA_idhwt = pw,
             r_idhmut = unname(m["r"]), n_idhmut = unname(m["n"]), p_vs_TCGA_idhmut = pm)
}))

## ---------- (C) EGFR by IDH status within grade IV only ----------
say("\n=== (C) mean EGFR by IDH status, grade IV only (the secondary-GBM stratum) ===\n")
say("    Expected biology: EGFR higher in IDH-wildtype. A reversal here cannot be\n")
say("    produced by grade/IDH mixing, since grade is held fixed.\n\n")
compC <- do.call(rbind, lapply(names(COH), function(k){
  M <- COH[[k]]; x <- M$e["EGFR", ]
  ok <- !is.na(x) & !is.na(M$idh) & !is.na(M$grade) & M$grade == 4
  mw <- mean(x[ok & M$idh == 0]); mm <- mean(x[ok & M$idh == 1])
  nw <- sum(ok & M$idh == 0);     nm <- sum(ok & M$idh == 1)
  p <- if (nw >= 3 && nm >= 3) suppressWarnings(wilcox.test(x[ok & M$idh == 0], x[ok & M$idh == 1])$p.value) else NA
  say("  %-10s IDH-wt %6.2f (n=%3d)   IDH-mut %6.2f (n=%3d)   diff %+0.2f  %s\n",
      k, mw, nw, mm, nm, mw - mm, if (is.na(p)) "(too few)" else sprintf("Wilcoxon p=%.2g", p))
  data.frame(cohort = k, analysis = "egfr_by_idh_grade4",
             mean_idhwt = mw, n_idhwt = nw, mean_idhmut = mm, n_idhmut = nm,
             diff = mw - mm, p_wilcox = p)
}))

## ---------- (D) is EGFR an outlier, or is the whole panel attenuated? ----------
say("\n=== (D) EGFR vs the genome-wide attenuation background (CGGA-693 vs TCGA) ===\n")
say("    Uniform attenuation across all genes would favour a cohort-level effect;\n")
say("    EGFR sitting in the tail of the residual distribution favours a\n")
say("    gene-specific measurement problem.\n\n")
gs <- intersect(rownames(TC$e), rownames(D693$e))
sdT <- apply(TC$e[gs, ], 1, sd); gs <- gs[sdT > 0.5]
rT <- suppressWarnings(apply(TC$e[gs, ], 1, function(x) cor(x, TC$grade, method = "spearman", use = "complete.obs")))
rC <- suppressWarnings(apply(D693$e[gs, ], 1, function(x) cor(x, D693$grade, method = "spearman", use = "complete.obs")))
keep <- !is.na(rT) & !is.na(rC)
rT <- rT[keep]; rC <- rC[keep]
fit <- lm(rC ~ rT)                     # attenuation slope: 1 = no attenuation
resid_z <- as.numeric(scale(residuals(fit)))
egfr_z <- if ("EGFR" %in% names(rT)) resid_z[which(names(rT) == "EGFR")] else NA
egfr_pct <- if (is.na(egfr_z)) NA else mean(resid_z <= egfr_z)
say("  genes analysed: %d\n", length(rT))
say("  overall attenuation slope (CGGA ~ TCGA): %.3f  (1.0 = no attenuation), R^2=%.3f\n",
    coef(fit)[2], summary(fit)$r.squared)
say("  correlation of the two correlation vectors: r=%.3f\n", cor(rT, rC))
say("  EGFR: r_TCGA=%+0.3f  r_CGGA=%+0.3f  residual z=%+0.2f  (%.1f%% of genes lie below it)\n",
    rT[["EGFR"]], rC[["EGFR"]], egfr_z, 100 * egfr_pct)
compD <- data.frame(cohort = "CGGA-693", analysis = "attenuation_background",
                    n_genes = length(rT), slope = unname(coef(fit)[2]),
                    r_of_r = cor(rT, rC), egfr_resid_z = egfr_z, egfr_percentile = egfr_pct)

res <- dplyr::bind_rows(compA, compB, compC, compD)
write.csv(res, "results/cgga_composition.csv", row.names = FALSE)
close(out)
cat("\nwrote results/cgga_composition.csv and .txt\n")
