# 14_egfr_amplification.R -- Does EGFR *amplification* (copy number) predict
# survival WITHIN IDH-wildtype glioma? Complements the expression result (12).
#
# Run from project root:  Rscript R/14_egfr_amplification.R
# Writes: results/egfr_amplification.csv

suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr); library(survival) })
source("R/_helpers.R")

# EGFR copy number per patient (primary tumors) from the gene-level CNV SE
egfr_cn <- function(project) {
  se <- readRDS(file.path("data", paste0(project, "_cnv.rds")))
  an <- assayNames(se); cn_assay <- an[grep("copy_number", an, ignore.case=TRUE)][1]
  if (is.na(cn_assay)) cn_assay <- an[1]
  gi <- which(rowData(se)$gene_name == "EGFR")[1]
  cn <- as.numeric(assay(se, cn_assay)[gi, ])
  data.frame(barcode = colnames(se), cn = cn, stringsAsFactors = FALSE) %>%
    mutate(patient = substr(barcode, 1, 12), sample_type = substr(barcode, 14, 15)) %>%
    filter(sample_type == "01", !is.na(cn)) %>%
    group_by(patient) %>% summarise(cn = max(cn), .groups = "drop")   # one per patient
}

# clinical (IDH, grade, survival) from cached expression objects
clin_for <- function(project) {
  se <- load_se(project); cd <- as.data.frame(colData(se))
  b <- .clinical_df(se, project)
  b$idh <- as.character(cd$paper_IDH.status)[match(b$barcode, rownames(cd))]
  b$grade_raw <- as.character(cd$paper_Grade)[match(b$barcode, rownames(cd))]
  b %>% filter(sample_type=="01", !is.na(time), time>0, !is.na(event)) %>%
    arrange(patient) %>% distinct(patient,.keep_all=TRUE) %>%
    transmute(patient, time, event, age, cohort=project, idh,
              grade_f=factor(dplyr::case_when(grepl("G2",grade_raw)~"II",grepl("G3",grade_raw)~"III",
                                              grepl("G4",grade_raw)~"IV",TRUE~NA_character_),levels=c("II","III","IV")))
}

cn  <- bind_rows(egfr_cn("TCGA-GBM"), egfr_cn("TCGA-LGG"))
clin<- bind_rows(clin_for("TCGA-GBM"), clin_for("TCGA-LGG"))
d <- inner_join(clin, cn, by="patient")
d$amp <- factor(ifelse(d$cn >= 6, "Amplified", "Not"), levels=c("Not","Amplified"))  # high-level amp
d$cn_log2 <- log2(d$cn + 1)

cat(sprintf("Matched CN+clinical: n=%d\n", nrow(d)))
cat("EGFR high-level amplification (CN>=6) by IDH:\n")
print(round(prop.table(table(d$idh, d$amp), 1), 3))

wt <- d %>% filter(idh=="WT"); wt$cn_z <- as.numeric(scale(wt$cn_log2))
cat(sprintf("\nIDH-wildtype: n=%d events=%d, amplified=%d (%.0f%%)\n",
    nrow(wt), sum(wt$event), sum(wt$amp=="Amplified"), 100*mean(wt$amp=="Amplified")))

row <- function(rhs, lab, term, df=wt) {
  m <- tryCatch(coxph(as.formula(paste("Surv(time,event)~",rhs)), data=df), error=function(e) NULL)
  if (is.null(m) || !term %in% rownames(summary(m)$coefficients)) return(NULL)
  s<-summary(m); ci<-s$conf.int[term,]
  data.frame(test=lab, HR=signif(ci[1],3), CI_low=signif(ci[3],3), CI_high=signif(ci[4],3),
             p=signif(s$coefficients[term,"Pr(>|z|)"],3), n=m$n, events=m$nevent, stringsAsFactors=FALSE)
}
res <- do.call(rbind, list(
  row("amp", "IDH-wt: amplified vs not (unadj)", "ampAmplified"),
  row("amp + age", "IDH-wt: amplified +age", "ampAmplified"),
  row("amp + age + grade_f", "IDH-wt: amplified +age+grade", "ampAmplified"),
  row("cn_z", "IDH-wt: CN (per SD, continuous)", "cn_z"),
  row("amp + age", "IDH-wt GBM: amplified +age", "ampAmplified", df=wt[wt$cohort=="TCGA-GBM",]) ))
write.csv(res, "results/egfr_amplification.csv", row.names=FALSE)
cat("\n=== EGFR amplification vs survival WITHIN IDH-wildtype ===\n"); print(res, row.names=FALSE)
