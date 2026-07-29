# 10_review_fixes.R -- address round-3 review points:
#  (1) common-sample nested TCGA models (attribute attenuation to IDH, not sample change)
#  (2) C-index with 95% CI; likelihood-ratio test for EGFR's incremental value
#  (3) formal test of the TCGA-vs-CGGA EGFR-grade correlation difference (Fisher r-to-z)
#  (4) CGGA-325 independent-batch confirmation of the positive-control failure
# Run from project root:  Rscript R/10_review_fixes.R
# Writes: results/common_sample_cox.csv, results/cgga325_controls.csv

suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr); library(survival) })
source("R/_helpers.R")

## ---------- (1)+(2) common-sample nested models, C-index CI, LR test ----------
se <- load_se("TCGA-LGG"); cd <- as.data.frame(colData(se))
tpm <- assay(se,"tpm_unstrand"); egfr <- log2(tpm[match("EGFR",rowData(se)$gene_name),]+1)
base <- .clinical_df(se,"TCGA-LGG")
base$idh <- as.character(cd$paper_IDH.status); base$EGFR <- egfr[match(base$barcode, colnames(se))]
tg <- function(g){g<-toupper(as.character(g));o<-rep(NA,length(g));o[grepl("G2| II ",g)]<-"II";o[grepl("G3| III ",g)]<-"III";factor(o,levels=c("II","III"))}
d <- base %>% filter(sample_type=="01", !is.na(time), time>0, !is.na(event)) %>%
     arrange(patient) %>% distinct(patient,.keep_all=TRUE)
d$grade_f <- tg(d$grade)
d$idh_f <- factor(ifelse(d$idh %in% c("WT","Mutant"), d$idh, NA), levels=c("WT","Mutant"))
# COMMON complete-case sample: age, grade, IDH all present
dc <- d %>% filter(!is.na(age), !is.na(grade_f), !is.na(idh_f))
dc$EGFR_z <- as.numeric(scale(dc$EGFR))
cat(sprintf("Common-sample n=%d events=%d\n", nrow(dc), sum(dc$event)))

cidx <- function(m){ cc<-concordance(m); c<-unname(cc$concordance); se<-sqrt(unname(cc$var))
  sprintf("%.3f (%.3f-%.3f)", c, c-1.96*se, c+1.96*se) }
egfr_row <- function(m,lab){ s<-summary(m); ci<-s$conf.int["EGFR_z",]
  data.frame(model=lab, HR=sprintf("%.2f (%.2f-%.2f)",ci[1],ci[3],ci[4]),
             p=signif(s$coefficients["EGFR_z","Pr(>|z|)"],3), C_index=cidx(m), stringsAsFactors=FALSE)}
m1 <- coxph(Surv(time,event)~EGFR_z, data=dc)
m2 <- coxph(Surv(time,event)~EGFR_z+age+grade_f, data=dc)
m3 <- coxph(Surv(time,event)~EGFR_z+age+grade_f+idh_f, data=dc)
common <- rbind(egfr_row(m1,"1. unadjusted"), egfr_row(m2,"2. +age+grade"), egfr_row(m3,"3. +age+grade+IDH"))
write.csv(common, "results/common_sample_cox.csv", row.names=FALSE)
cat("\n=== (1) Common-sample nested EGFR models (identical rows) ===\n"); print(common, row.names=FALSE)
# LR test: does EGFR add anything over age+grade+IDH?
m_full_noE <- coxph(Surv(time,event)~age+grade_f+idh_f, data=dc)
lr <- anova(m_full_noE, m3)
cat(sprintf("\n(2) LR test EGFR incremental over age+grade+IDH: chisq=%.2f df=%d p=%.3f\n",
    lr$Chisq[2], lr$Df[2], lr$`Pr(>|Chi|)`[2]))
cat(sprintf("    C-index age+grade+IDH WITHOUT EGFR = %s ; WITH EGFR = %s\n", cidx(m_full_noE), cidx(m3)))

## ---------- (3) Fisher r-to-z: TCGA vs CGGA EGFR-grade correlation ----------
fisher_z <- function(r1,n1,r2,n2){ z<-(atanh(r1)-atanh(r2))/sqrt(1/(n1-3)+1/(n2-3)); 2*pnorm(-abs(z)) }
# TCGA EGFR-grade (LGG II/III + GBM IV)
gb <- load_se("TCGA-GBM"); egG <- log2(assay(gb,"tpm_unstrand")[match("EGFR",rowData(gb)$gene_name),]+1)
kG <- substr(colnames(gb),14,15)=="01"
tcg <- data.frame(egfr=c(d$EGFR, egG[kG]),
                  grade=c(ifelse(d$grade=="G3"|grepl("III",d$grade),3,2), rep(4,sum(kG))))
rT <- cor(tcg$egfr,tcg$grade,method="spearman",use="complete.obs"); nT<-sum(complete.cases(tcg))

## ---------- (4) CGGA-325 confirmation ----------
read_cgga <- function(gfile, cfile){
  ex<-read.delim(gfile,header=TRUE,row.names=1,check.names=FALSE)
  ep<-as.data.frame(t(log2(as.matrix(ex)+1))); ep$id<-rownames(ep)
  cl<-read.delim(cfile,header=TRUE,check.names=FALSE); names(cl)[c(1,2,4,7,8,11)]<-c("id","prs","grade","os","censor","idh")
  ep%>%inner_join(cl,by="id")%>%filter(prs=="Primary")%>%
    mutate(gnum=as.numeric(factor(grade,levels=c("WHO II","WHO III","WHO IV"))),
           time=as.numeric(os),event=as.integer(censor),
           idh_f=factor(ifelse(idh%in%c("Mutant","Wildtype"),idh,NA),levels=c("Wildtype","Mutant")))
}
c693 <- read_cgga("data/cgga/cgga_genes.tsv","data/cgga/cgga_clinical.tsv")
c325 <- read_cgga("data/cgga/cgga325_genes.tsv","data/cgga/cgga325_clinical.tsv")
ctrls <- function(x,lab){
  r<-cor(x$EGFR,x$gnum,method="spearman",use="complete.obs"); n<-sum(!is.na(x$EGFR)&!is.na(x$gnum))
  idhHR<-tryCatch(summary(coxph(Surv(time,event)~idh_f,data=x[!is.na(x$idh_f),]))$conf.int[1,1],error=function(e)NA)
  grHR<-tryCatch(summary(coxph(Surv(time,event)~factor(gnum),data=x))$conf.int[2,1],error=function(e)NA)
  data.frame(cohort=lab, EGFR_grade_r=round(r,3), n=n,
             p_vs_TCGA=signif(fisher_z(rT,nT,r,n),3),
             IDH_HR=round(idhHR,3), gradeIV_HR=round(grHR,2), stringsAsFactors=FALSE)
}
tab <- rbind(data.frame(cohort="TCGA (reference)",EGFR_grade_r=round(rT,3),n=nT,p_vs_TCGA=NA,IDH_HR=NA,gradeIV_HR=NA),
             ctrls(c693,"CGGA-693"), ctrls(c325,"CGGA-325"))
write.csv(tab, "results/cgga325_controls.csv", row.names=FALSE)
cat("\n=== (3)+(4) EGFR-grade positive control: TCGA vs CGGA batches ===\n"); print(tab, row.names=FALSE)
cat("\nCGGA-325 EGFR mean by grade:\n"); print(round(tapply(c325$EGFR,c325$gnum,mean),2))
cat("CGGA-325 EGFR mean by IDH:\n"); print(round(tapply(c325$EGFR,c325$idh_f,mean),2))
