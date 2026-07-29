# 13_systematic_method.R -- Systematic QC-gated cross-cohort prognostic screen.
#
# Generalizes the EGFR cautionary tale into a method: for a panel of genes,
# test the prognostic association in TCGA and its replication in CGGA, and gate
# on a per-gene expression-level positive control (does the gene's expression
# track WHO grade consistently across cohorts?). Quantifies how positive-control
# gating changes apparent cross-cohort replication.
#
# Run from project root:  Rscript R/13_systematic_method.R
# Writes: results/systematic_screen.csv, figures/systematic_method.png

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival); library(ggplot2)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly = TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep="\t")) else read.delim(f, check.names=FALSE)

## ---------- TCGA (LGG+GBM): expression matrix + clinical ----------
tcga_one <- function(project, gradeIV=FALSE) {
  se <- load_se(project); cd <- as.data.frame(colData(se))
  k <- substr(colnames(se),14,15)=="01"
  expr <- log2(assay(se,"tpm_unstrand")[,k]+1); rownames(expr) <- rowData(se)$gene_name
  b <- .clinical_df(se, project); b <- b[match(colnames(expr), b$barcode),]
  g <- if (gradeIV) rep(4, ncol(expr)) else ifelse(grepl("G3|III", cd$paper_Grade[k]), 3, 2)
  keep <- !is.na(b$time) & b$time>0 & !is.na(b$event) & !duplicated(b$patient)
  list(expr=expr[,keep], time=b$time[keep], event=b$event[keep], age=b$age[keep], grade=g[keep])
}
L <- tcga_one("TCGA-LGG"); G <- tcga_one("TCGA-GBM", gradeIV=TRUE)
common_sym <- intersect(rownames(L$expr), rownames(G$expr))
Texpr <- cbind(L$expr[common_sym,], G$expr[common_sym,])
Tclin <- data.frame(time=c(L$time,G$time), event=c(L$event,G$event),
                    age=c(L$age,G$age), grade=c(L$grade,G$grade))
Texpr <- Texpr[!duplicated(rownames(Texpr)), ]

## ---------- CGGA-693: full matrix + clinical ----------
cg <- rd("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
rownames(cg) <- cg[[1]]; cg[[1]] <- NULL; cg <- log2(as.matrix(cg)+1)
ccl <- read.delim("data/cgga/cgga_clinical.tsv", check.names=FALSE)
names(ccl)[c(1,2,4,6,7,8,11)] <- c("id","prs","grade","age","os","censor","idh")
ccl <- ccl %>% filter(prs=="Primary") %>%
  mutate(time=as.numeric(os), event=as.integer(censor), age=as.numeric(age),
         gnum=as.numeric(factor(grade, levels=c("WHO II","WHO III","WHO IV"))))
csel <- intersect(colnames(cg), ccl$id)
cg <- cg[, csel]; ccl <- ccl[match(csel, ccl$id),]

## ---------- panel: top-variable TCGA genes present in CGGA ----------
gv <- apply(Texpr, 1, var)
panel <- names(sort(gv, decreasing=TRUE))
panel <- intersect(panel, rownames(cg))
panel <- head(panel[!is.na(panel) & panel!=""], 500)
cat("Panel size:", length(panel), " (top-variable TCGA genes present in CGGA)\n")

## ---------- per-gene: prognostic (age-adj Cox) + grade positive control ----------
cox_beta <- function(expr_vec, clin) {
  z <- as.numeric(scale(expr_vec))
  m <- tryCatch(coxph(Surv(clin$time, clin$event) ~ z + clin$age), error=function(e) NULL)
  if (is.null(m)) return(c(NA,NA))
  s <- summary(m)$coefficients; c(s["z","coef"], s["z","Pr(>|z|)"])
}
res <- lapply(panel, function(g) {
  bT <- cox_beta(Texpr[g,], Tclin); bC <- cox_beta(cg[g,], ccl)
  rT <- suppressWarnings(cor(Texpr[g,], Tclin$grade, method="spearman", use="complete.obs"))
  rC <- suppressWarnings(cor(cg[g,], ccl$gnum, method="spearman", use="complete.obs"))
  data.frame(gene=g, beta_T=bT[1], p_T=bT[2], beta_C=bC[1], p_C=bC[2],
             r_grade_T=rT, r_grade_C=rC, stringsAsFactors=FALSE)
})
res <- do.call(rbind, res)
res$p_T_BH <- p.adjust(res$p_T, "BH")
# QC: gene measured consistently vs grade across cohorts
res$qc_ok <- with(res, sign(r_grade_T)==sign(r_grade_C) & abs(r_grade_T)>0.1 & abs(r_grade_C)>0.1)
# prognostic in TCGA (FDR) and replication in CGGA
res$prognostic_T <- res$p_T_BH < 0.05
res$replicated   <- res$prognostic_T & (sign(res$beta_T)==sign(res$beta_C)) & (res$p_C < 0.05)
res$signflip     <- res$prognostic_T & (sign(res$beta_T)!=sign(res$beta_C)) & (res$p_C < 0.05)
write.csv(res, "results/systematic_screen.csv", row.names=FALSE)

## ---------- summary ----------
prog <- res[res$prognostic_T & !is.na(res$beta_C), ]
tab <- prog %>% group_by(qc_ok) %>%
  summarise(n=n(),
            replicated_pct = round(100*mean(replicated, na.rm=TRUE),1),
            signflip_pct   = round(100*mean(signflip, na.rm=TRUE),1), .groups="drop")
cat("\n=== Of genes prognostic in TCGA (FDR<0.05), replication in CGGA by QC status ===\n")
print(as.data.frame(tab), row.names=FALSE)
cat(sprintf("\nEGFR: prognostic_T=%s  QC_ok=%s  beta_T=%.2f beta_C=%.2f  r_grade_T=%.2f r_grade_C=%.2f\n",
    res$prognostic_T[res$gene=="EGFR"], res$qc_ok[res$gene=="EGFR"],
    res$beta_T[res$gene=="EGFR"], res$beta_C[res$gene=="EGFR"],
    res$r_grade_T[res$gene=="EGFR"], res$r_grade_C[res$gene=="EGFR"]))
cat(sprintf("Overall QC pass rate in panel: %.1f%% (%d/%d)\n",
    100*mean(res$qc_ok,na.rm=TRUE), sum(res$qc_ok,na.rm=TRUE), sum(!is.na(res$qc_ok))))

## ---------- figure: TCGA vs CGGA grade-correlation, prognostic genes ----------
p <- ggplot(prog, aes(r_grade_T, r_grade_C, colour=qc_ok)) +
  geom_hline(yintercept=0, linewidth=.3, colour="grey70") + geom_vline(xintercept=0, linewidth=.3, colour="grey70") +
  geom_point(alpha=.6, size=1.6) +
  scale_colour_manual(values=c("FALSE"="#e31a1c","TRUE"="#1f78b4"),
                      labels=c("QC-discordant","QC-concordant"), name="Positive control") +
  labs(x="Expression–grade correlation (TCGA)", y="Expression–grade correlation (CGGA)",
       title="Cross-cohort expression QC for TCGA-prognostic genes\n(points off the diagonal are measured inconsistently)") +
  theme_bw(base_size=12) + theme(legend.position="top", plot.title=element_text(size=11,face="bold"))
ggsave("figures/systematic_method.png", p, width=6.5, height=5, dpi=300)
cat("\nwrote figures/systematic_method.png\n")
