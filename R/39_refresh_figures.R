# 39_refresh_figures.R -- BRING THE THREE MAIN FIGURES INTO LINE WITH THE TEXT.
#
# A visual review of the published figures found each asserting something the
# manuscript no longer claims, or referring to an element that is not drawn:
#
#   Figure 1 (positive_control_EGFR_grade.png)
#     Subtitle says CGGA is "flat". It is not flat: it dips from grade II to III
#     (-0.015 -> -0.055) and rises from III to IV (-> +0.070). The Spearman r is
#     near zero because the pattern is NON-MONOTONIC, not because the line is
#     level, and the grade-IV error bars overlap grade II's. Redrawn with an
#     accurate subtitle.
#
#   Figure 2 (systematic_method.png)
#     The manuscript's caption directs the reader to "the diagonal" -- which is
#     not drawn on the plot. The subtitle asserts flagged genes are "measured
#     inconsistently", the strong causal reading section 6.4 has since qualified.
#     Redrawn with the identity line actually present, an attenuation slope, and
#     a subtitle describing the pattern (flagged genes sit near zero on the CGGA
#     axis) rather than its presumed cause.
#
#   Figure 3 (method_auc.png)
#     Draws a dashed reference line at 0.5 labelled as chance. Section 6.9 shows
#     0.5 is the WRONG null: a same-population cohort reaches 0.79-0.90 with the
#     same product-form score. As drawn, the figure visually asserts exactly the
#     interpretation the paper retracts. Redrawn with the matched-null level
#     shown per cohort, so the bar-versus-null comparison is the visible one.
#
# Run from project root:  Rscript R/39_refresh_figures.R

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(ggplot2); library(tidyr)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep = "\t")) else read.delim(f, check.names = FALSE)

## ---------------- Figure 1 ----------------
tcga <- function(){
  one <- function(project, gradeIV = FALSE){
    se <- load_se(project); cd <- as.data.frame(colData(se)); k <- substr(colnames(se), 14, 15) == "01"
    e <- log2(assay(se, "tpm_unstrand")[, k] + 1); rownames(e) <- rowData(se)$gene_name
    e <- e[!duplicated(rownames(e)), ]
    b <- .clinical_df(se, project); b <- b[match(colnames(e), b$barcode), ]
    gr <- as.character(cd$paper_Grade[k])
    g <- if (gradeIV) ifelse(grepl("G4|IV", gr), 4, NA) else ifelse(grepl("G3|III", gr), 3, ifelse(grepl("G2| II", gr), 2, NA))
    keep <- !duplicated(b$patient)
    list(egfr = e["EGFR", keep], grade = g[keep])
  }
  L <- one("TCGA-LGG"); G <- one("TCGA-GBM", TRUE)
  data.frame(egfr = c(L$egfr, G$egfr), grade = c(L$grade, G$grade))
}
cgga <- function(){
  ex <- rd("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt"); rownames(ex) <- ex[[1]]; ex[[1]] <- NULL
  ex <- log2(as.matrix(ex) + 1)
  cl <- read.delim("data/cgga/cgga_clinical.tsv", check.names = FALSE)
  names(cl)[c(1,2,4)] <- c("id","prs","grade")
  cl <- cl %>% filter(prs == "Primary") %>%
    mutate(grade = as.numeric(factor(grade, levels = c("WHO II","WHO III","WHO IV"))) + 1)
  sel <- intersect(colnames(ex), cl$id); cl <- cl[match(sel, cl$id), ]
  data.frame(egfr = ex["EGFR", sel], grade = cl$grade)
}
mk <- function(d, lab){ d <- d[!is.na(d$egfr) & !is.na(d$grade), ]
  data.frame(cohort = lab, grade = factor(d$grade, levels = 2:4, labels = c("II","III","IV")),
             z = as.numeric(scale(d$egfr))) }
TC <- tcga(); CG <- cgga()
P <- rbind(mk(TC, "TCGA (TPM)"), mk(CG, "CGGA (RSEM)"))
S <- P %>% group_by(cohort, grade) %>% summarise(m = mean(z), se = sd(z)/sqrt(dplyr::n()), .groups = "drop")
rT <- cor(TC$egfr, TC$grade, method="spearman", use="complete.obs")
rC <- cor(CG$egfr, CG$grade, method="spearman", use="complete.obs")
g1 <- ggplot(S, aes(grade, m, colour = cohort, group = cohort)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  geom_errorbar(aes(ymin = m-se, ymax = m+se), width = 0.12) +
  labs(x = "WHO grade", y = "EGFR expression (per-SD units)", colour = "Cohort",
       title = "Positive control: EGFR expression vs grade",
       subtitle = sprintf("monotonic in TCGA (Spearman r = %+.2f);\nnon-monotonic and attenuated in CGGA (r = %+.2f)", rT, rC)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 10))
ggsave("figures/positive_control_EGFR_grade.png", g1, width = 7.0, height = 4.6, dpi = 300)
cat(sprintf("Figure 1 redrawn. TCGA r=%+.3f, CGGA r=%+.3f; CGGA means by grade: %s\n",
    rT, rC, paste(sprintf("%.3f", S$m[S$cohort=="CGGA (RSEM)"]), collapse=", ")))

## ---------------- Figure 2 ----------------
ss <- read.csv("results/systematic_screen.csv")
ss$QC <- ifelse(ss$qc_ok, "QC-concordant", "QC-discordant (flagged)")
fit <- lm(r_grade_C ~ r_grade_T, ss)
g2 <- ggplot(ss, aes(r_grade_T, r_grade_C, colour = QC)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey35") +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], colour = "grey20", linewidth = 0.5) +
  geom_point(alpha = 0.65, size = 1.5) +
  annotate("text", x = -0.72, y = 0.55, hjust = 0, size = 3.1, colour = "grey25",
           label = sprintf("dashed: identity\nsolid: fitted slope %.2f (attenuation)", coef(fit)[2])) +
  scale_colour_manual(values = c("QC-concordant" = "#3182bd", "QC-discordant (flagged)" = "#e6550d")) +
  labs(x = "Expression–grade correlation (TCGA)", y = "Expression–grade correlation (CGGA)",
       colour = "Positive control",
       title = "Cross-cohort expression QC for TCGA-prognostic genes",
       subtitle = "flagged genes sit near zero on the CGGA axis;\nCGGA correlations are attenuated overall") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 10))
ggsave("figures/systematic_method.png", g2, width = 7.4, height = 5.3, dpi = 300)
cat(sprintf("Figure 2 redrawn. Attenuation slope = %.3f; flagged genes: %d of %d\n",
    coef(fit)[2], sum(!ss$qc_ok), nrow(ss)))

## ---------------- Figure 3 ----------------
res <- read.csv("results/method_auc.csv")
nullv <- read.csv("results/null_full_universe.csv")
nl <- nullv[nullv$metric == "prod" & grepl("glioma", nullv$setting),
            c("setting","auc_null")]
nl$cohort <- sub(".*-> *", "", nl$setting)
nl$cohort <- ifelse(nl$cohort == "array-301", "array-301", nl$cohort)
pl <- res %>% select(cohort, AUC_grade, AUC_idh, AUC_detect, AUC_precision) %>%
  pivot_longer(-cohort, names_to = "metric", values_to = "AUC") %>%
  mutate(metric = recode(metric, AUC_grade = "Grade anchor", AUC_idh = "IDH anchor",
                         AUC_detect = "Detectability", AUC_precision = "Precision (SE of coefficient)"))
pl$cohort <- factor(pl$cohort, levels = c("CGGA-693","CGGA-325","array-301"))
nl$cohort <- factor(nl$cohort, levels = levels(pl$cohort))
g3 <- ggplot(pl, aes(cohort, AUC, fill = metric)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey60") +
  geom_segment(data = nl, aes(x = as.numeric(cohort) - 0.45, xend = as.numeric(cohort) + 0.45,
                              y = auc_null, yend = auc_null), inherit.aes = FALSE,
               colour = "#c0392b", linewidth = 0.9) +
  annotate("text", x = 0.55, y = 0.985, hjust = 0, size = 3.1, colour = "#c0392b",
           label = "red bar = same-population null (§6.9): the correct comparison") +
  annotate("text", x = 0.55, y = 0.44, hjust = 0, size = 3.0, colour = "grey45",
           label = "dotted 0.5 = chance, which is NOT the appropriate null here") +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(x = NULL, y = "AUC for predicting cross-cohort replication", fill = "Positive-control metric",
       title = "Anchor scores against chance and against a same-population null",
       subtitle = "bars exceed chance; the product-form anchor does not exceed\nits matched null in CGGA-693") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", legend.text = element_text(size = 8),
        legend.title = element_text(size = 9), plot.subtitle = element_text(size = 9)) +
  guides(fill = guide_legend(nrow = 2))
ggsave("figures/method_auc.png", g3, width = 7.6, height = 5.2, dpi = 300)
cat(sprintf("Figure 3 redrawn with matched-null levels: %s\n",
    paste(sprintf("%s=%.3f", nl$cohort, nl$auc_null), collapse = ", ")))
