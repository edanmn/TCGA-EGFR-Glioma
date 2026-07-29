# 12_within_idhwt.R -- Does EGFR expression carry prognostic signal WITHIN
# IDH-wildtype glioma (the aggressive subtype), where the subtype-confounding
# explanation no longer applies? Uses TCGA (LGG+GBM), IDH-wildtype only.
#
# Run from project root:  Rscript R/12_within_idhwt.R
# Writes: results/within_idhwt.csv

suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr); library(survival) })
source("R/_helpers.R")

build <- function(project) {
  se <- load_se(project); cd <- as.data.frame(colData(se))
  tpm <- assay(se,"tpm_unstrand")
  egfr <- log2(tpm[match("EGFR",rowData(se)$gene_name),]+1)
  b <- .clinical_df(se, project)
  b$EGFR <- egfr[match(b$barcode, colnames(se))]
  b$idh <- as.character(cd$paper_IDH.status)[match(b$barcode, rownames(cd))]
  b$grade_raw <- as.character(cd$paper_Grade)[match(b$barcode, rownames(cd))]
  b %>% filter(sample_type=="01", !is.na(time), time>0, !is.na(event)) %>%
    arrange(patient) %>% distinct(patient,.keep_all=TRUE) %>%
    mutate(cohort=project,
           grade_f=factor(case_when(grepl("G2",grade_raw)~"II",grepl("G3",grade_raw)~"III",
                                    grepl("G4",grade_raw)~"IV",TRUE~NA_character_),levels=c("II","III","IV")))
}
d <- bind_rows(build("TCGA-LGG"), build("TCGA-GBM"))
wt <- d %>% filter(idh=="WT")            # IDH-wildtype only
wt$EGFR_z <- as.numeric(scale(wt$EGFR))
cat(sprintf("IDH-wildtype glioma: n=%d events=%d (LGG %d, GBM %d)\n",
    nrow(wt), sum(wt$event), sum(wt$cohort=="TCGA-LGG"), sum(wt$cohort=="TCGA-GBM")))
cat("Grade distribution (IDH-wt):\n"); print(table(wt$grade_f, useNA="ifany"))

cidx <- function(m){cc<-concordance(m); s<-sqrt(unname(cc$var)); sprintf("%.3f (%.3f-%.3f)",unname(cc$concordance),unname(cc$concordance)-1.96*s,unname(cc$concordance)+1.96*s)}
row <- function(rhs,lab,df=wt){
  m<-tryCatch(coxph(as.formula(paste("Surv(time,event)~",rhs)),data=df),error=function(e)NULL)
  if(is.null(m)||!"EGFR_z"%in%rownames(summary(m)$coefficients)) return(NULL)
  s<-summary(m); ci<-s$conf.int["EGFR_z",]
  data.frame(subset=lab, HR=signif(ci[1],3), CI_low=signif(ci[3],3), CI_high=signif(ci[4],3),
             p=signif(s$coefficients["EGFR_z","Pr(>|z|)"],3), C_index=cidx(m), n=m$n, events=m$nevent,
             stringsAsFactors=FALSE)
}
res <- do.call(rbind, list(
  row("EGFR_z","IDH-wt: unadjusted"),
  row("EGFR_z + age","IDH-wt: +age"),
  row("EGFR_z + age + grade_f","IDH-wt: +age+grade"),
  row("EGFR_z + age + cohort","IDH-wt: +age+cohort(LGG/GBM)"),
  row("EGFR_z + age", "IDH-wt GBM only: +age", df=within(wt[wt$cohort=="TCGA-GBM",],{EGFR_z<-as.numeric(scale(EGFR))})),
  row("EGFR_z + age", "IDH-wt LGG only: +age", df=within(wt[wt$cohort=="TCGA-LGG",],{EGFR_z<-as.numeric(scale(EGFR))}))
))
write.csv(res, "results/within_idhwt.csv", row.names=FALSE)
cat("\n=== EGFR (per SD) within IDH-wildtype glioma ===\n"); print(res, row.names=FALSE)
