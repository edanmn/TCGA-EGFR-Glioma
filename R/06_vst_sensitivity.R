# 06_vst_sensitivity.R -- sensitivity analysis re-running the key EGFR models on
# DESeq2 variance-stabilized (VST) expression instead of log2(TPM+1), to check
# that the conclusion (EGFR not independent of IDH) is not an artifact of TPM.
#
# Run from project root:  Rscript R/06_vst_sensitivity.R
# Writes: results/vst_sensitivity.csv

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(DESeq2); library(survival); library(dplyr)
})
source("R/_helpers.R")

GENES <- c("EGFR","PDGFRA","PTEN","NF1","PIK3CA","PIK3R1","RB1","CDKN2A","TP53","IDH1")

# VST-transform raw STAR counts and return a per-sample gene matrix (log2-scale VST).
vst_expr <- function(se, genes) {
  cts <- assay(se, "unstranded")
  mode(cts) <- "integer"
  cts[is.na(cts)] <- 0L
  keep <- rowSums(cts) >= 10                      # drop essentially-unexpressed genes
  cts <- cts[keep, ]
  vsd <- suppressMessages(vst(cts, blind = TRUE)) # fast, design-blind
  sym <- rowData(se)$gene_name[keep]
  idx <- match(genes, sym)
  found <- genes[!is.na(idx)]
  mat <- vsd[idx[!is.na(idx)], , drop = FALSE]
  out <- as.data.frame(t(mat)); colnames(out) <- found; out$barcode <- colnames(se)
  out
}

# clinical + molecular covariates (mirror 05)
clin_ext <- function(se, project) {
  cd <- as.data.frame(colData(se)); base <- .clinical_df(se, project)
  pick <- function(c) { h <- c[c %in% names(cd)]; if (length(h)) as.character(cd[[h[1]]]) else rep(NA, nrow(cd)) }
  base$idh <- pick("paper_IDH.status"); base$subtype <- pick("paper_IDH.codel.subtype"); base
}
tidy_grade <- function(g){g<-toupper(as.character(g));o<-rep(NA_character_,length(g));o[grepl("G2| II ",g)]<-"II";o[grepl("G3| III ",g)]<-"III";o[grepl("G4| IV ",g)]<-"IV";factor(o,levels=c("II","III","IV"))}
norm_idh <- function(x){x<-toupper(trimws(as.character(x)));x[x=="WILDTYPE"|x=="WILD-TYPE"]<-"WT";x[x=="MUTANT"|x=="MUT"]<-"Mutant";x[!x%in%c("WT","Mutant")]<-NA;factor(x,levels=c("WT","Mutant"))}

build_vst <- function(project) {
  se <- load_se(project)
  df <- clin_ext(se, project) %>%
    left_join(vst_expr(se, GENES), by = "barcode") %>%
    filter(sample_type == "01", !is.na(time), time > 0, !is.na(event)) %>%
    arrange(patient) %>% distinct(patient, .keep_all = TRUE)
  df$grade_f <- tidy_grade(df$grade); df$idh_f <- norm_idh(df$idh)
  for (g in GENES) if (g %in% names(df)) df[[paste0(g,"_z")]] <- as.numeric(scale(df[[g]]))
  df
}

fit_row <- function(df, rhs, cohort, label, term="EGFR_z") {
  m <- tryCatch(coxph(as.formula(paste("Surv(time,event)~",rhs)), data=df), error=function(e) NULL)
  if (is.null(m) || !term %in% rownames(summary(m)$coefficients)) return(NULL)
  s <- summary(m); ci <- s$conf.int[term,]
  data.frame(cohort=cohort, model=label, HR=signif(unname(ci["exp(coef)"]),3),
             CI_low=signif(unname(ci["lower .95"]),3), CI_high=signif(unname(ci["upper .95"]),3),
             p=signif(s$coefficients[term,"Pr(>|z|)"],3),
             C_index=signif(unname(s$concordance["C"]),3), n=m$n, events=m$nevent,
             stringsAsFactors=FALSE)
}

lgg <- build_vst("TCGA-LGG"); gbm <- build_vst("TCGA-GBM")
rows <- list(
  fit_row(lgg, "EGFR_z", "TCGA-LGG", "1. unadjusted"),
  fit_row(lgg, "EGFR_z + age + grade_f", "TCGA-LGG", "2. + age + grade"),
  fit_row(lgg[!is.na(lgg$idh_f),], "EGFR_z + age + grade_f + idh_f", "TCGA-LGG", "3. + age + grade + IDH"),
  fit_row(gbm, "EGFR_z", "TCGA-GBM", "1. unadjusted"),
  fit_row(gbm[!is.na(gbm$idh_f),], "EGFR_z + age + idh_f", "TCGA-GBM", "3. + age + IDH")
)
tbl <- do.call(rbind, Filter(Negate(is.null), rows))
tbl$normalization <- "VST"
write.csv(tbl, "results/vst_sensitivity.csv", row.names=FALSE)
cat("\n=== VST sensitivity: EGFR HR per 1 SD (compare to TPM Table 3) ===\n")
print(tbl, row.names=FALSE)
