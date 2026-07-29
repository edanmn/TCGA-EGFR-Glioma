# 16_method_framework.R -- develop the positive-control-gated cross-cohort
# validation into a proper framework and evaluate it at genome scale.
#
#   * Discovery cohort: TCGA (LGG+GBM). Replication cohorts: CGGA-693, -325, array-301.
#   * For every gene: age-adjusted Cox (prognostic signal) + biological "positive
#     control" anchors (expression vs grade; expression vs IDH) + detectability.
#   * Question: does a gene's positive-control concordance across cohorts PREDICT
#     whether its TCGA prognostic association replicates? Quantify with AUC, and
#     BENCHMARK the grade anchor against IDH-anchor and detectability baselines.
#
# Run from project root:  Rscript R/16_method_framework.R
# Writes: results/method_framework.csv (per-gene), results/method_auc.csv (summary),
#         figures/method_auc.png

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival); library(ggplot2)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly=TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep="\t")) else read.delim(f, check.names=FALSE)

## ---------- assemble a cohort: expression matrix (genes x samples), clinical ----------
tcga_mat <- function(){
  one <- function(project, gradeIV=FALSE){
    se<-load_se(project); cd<-as.data.frame(colData(se)); k<-substr(colnames(se),14,15)=="01"
    e<-log2(assay(se,"tpm_unstrand")[,k]+1); rownames(e)<-rowData(se)$gene_name
    e<-e[!duplicated(rownames(e)),]
    b<-.clinical_df(se,project); b<-b[match(colnames(e),b$barcode),]
    g<-if(gradeIV) rep(4,ncol(e)) else ifelse(grepl("G3|III",cd$paper_Grade[k]),3,2)
    idh<-as.character(cd$paper_IDH.status[k])
    keep<-!is.na(b$time)&b$time>0&!is.na(b$event)&!duplicated(b$patient)
    list(e=e[,keep], time=b$time[keep], event=b$event[keep], age=b$age[keep],
         grade=g[keep], idh=ifelse(idh[keep]=="Mutant",1,ifelse(idh[keep]=="WT",0,NA)))
  }
  L<-one("TCGA-LGG"); G<-one("TCGA-GBM",TRUE)
  cg<-intersect(rownames(L$e),rownames(G$e))
  list(e=cbind(L$e[cg,],G$e[cg,]),
       clin=data.frame(time=c(L$time,G$time),event=c(L$event,G$event),
                       age=c(L$age,G$age),grade=c(L$grade,G$grade),idh=c(L$idh,G$idh)))
}
cgga_mat <- function(gfile, cfile, cols){
  ex<-rd(gfile); rownames(ex)<-ex[[1]]; ex[[1]]<-NULL; ex<-log2(as.matrix(ex)+1)
  cl<-read.delim(cfile,check.names=FALSE); names(cl)[cols]<-c("id","prs","grade","age","os","censor","idh")
  cl<-cl%>%filter(prs=="Primary")%>%mutate(time=as.numeric(os),event=as.integer(censor),age=as.numeric(age),
    grade=as.numeric(factor(grade,levels=c("WHO II","WHO III","WHO IV"))),
    idh=ifelse(idh=="Mutant",1,ifelse(idh=="Wildtype",0,NA)))
  sel<-intersect(colnames(ex),cl$id); cl<-cl[match(sel,cl$id),]
  list(e=ex[,sel], clin=cl[,c("time","event","age","grade","idh")])
}

TC<-tcga_mat()
C693<-cgga_mat("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt","data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
C325<-cgga_mat("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt","data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
CARR<-cgga_mat("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt","data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12))

## ---------- gene universe: common, adequately expressed; cap for runtime ----------
genes<-Reduce(intersect,list(rownames(TC$e),rownames(C693$e),rownames(C325$e),rownames(CARR$e)))
sdT<-apply(TC$e[genes,],1,sd); genes<-genes[sdT>0.5]
if(length(genes)>8000) genes<-names(sort(sdT[genes],decreasing=TRUE))[1:8000]
cat("Gene universe:", length(genes), "\n")

## ---------- per-gene stats in a cohort ----------
stats_cohort <- function(M, gs){
  cl<-M$clin; e<-M$e
  do.call(rbind, lapply(gs, function(g){
    x<-e[g,]; z<-as.numeric(scale(x))
    m<-tryCatch(coxph(Surv(cl$time,cl$event)~z+cl$age), error=function(e) NULL)
    bp<-if(is.null(m)) c(NA,NA) else summary(m)$coefficients["z",c("coef","Pr(>|z|)")]
    data.frame(gene=g, beta=bp[1], p=bp[2],
               r_grade=suppressWarnings(cor(x,cl$grade,method="spearman",use="complete.obs")),
               r_idh=suppressWarnings(cor(x,cl$idh,method="spearman",use="complete.obs")),
               meanx=mean(x,na.rm=TRUE))
  }))
}
cat("computing TCGA...\n");  sT<-stats_cohort(TC,genes)
cat("computing CGGA-693...\n"); s6<-stats_cohort(C693,genes)
cat("computing CGGA-325...\n"); s3<-stats_cohort(C325,genes)
cat("computing array-301...\n"); sA<-stats_cohort(CARR,genes)
sT$p_BH<-p.adjust(sT$p,"BH")

## ---------- AUC (Mann-Whitney) ----------
auc<-function(score,label){ ok<-!is.na(score)&!is.na(label); score<-score[ok]; label<-label[ok]
  n1<-sum(label==1); n0<-sum(label==0); if(n1==0||n0==0) return(NA)
  r<-rank(score); (sum(r[label==1])-n1*(n1+1)/2)/(n1*n0) }

evaluate<-function(sR, lab){
  D<-merge(sT, sR, by="gene", suffixes=c("_T","_R"))
  prog<-D$p_BH<0.05 & !is.na(D$beta_R)
  D<-D[prog,]
  replicated<-as.integer(sign(D$beta_T)==sign(D$beta_R) & D$p_R<0.05)
  # QC scores (higher = more trustworthy)
  qc_grade<-D$r_grade_T*D$r_grade_R            # positive when grade-concordant
  qc_idh  <-D$r_idh_T*D$r_idh_R                # alternative biological anchor
  detect  <-pmin(D$meanx_T,D$meanx_R)          # naive detectability baseline
  data.frame(cohort=lab, n=nrow(D), replic_rate=round(mean(replicated),3),
             AUC_grade=round(auc(qc_grade,replicated),3),
             AUC_idh=round(auc(qc_idh,replicated),3),
             AUC_detect=round(auc(detect,replicated),3))
}
res<-rbind(evaluate(s6,"CGGA-693"), evaluate(s3,"CGGA-325"), evaluate(sA,"array-301"))
write.csv(res,"results/method_auc.csv",row.names=FALSE)
cat("\n=== QC score predicts cross-cohort replication (AUC) ===\n"); print(res,row.names=FALSE)

# precision/recall of a simple grade-QC flag (concordant sign & both|r|>0.1) in CGGA-693
D<-merge(sT,s6,by="gene",suffixes=c("_T","_R")); D<-D[D$p_BH<0.05 & !is.na(D$beta_R),]
rep6<-sign(D$beta_T)==sign(D$beta_R) & D$p_R<0.05
qcflag<-sign(D$r_grade_T)==sign(D$r_grade_R) & abs(D$r_grade_T)>0.1 & abs(D$r_grade_R)>0.1
cat(sprintf("\nGrade-QC flag (CGGA-693): of QC-fail genes, %.0f%% fail to replicate (precision of flag);\n  of non-replicating genes, %.0f%% are QC-flagged (recall).\n",
    100*mean(!rep6[!qcflag]), 100*mean(!qcflag[!rep6])))
write.csv(cbind(D[,c("gene","beta_T","p_BH","beta_R","p_R","r_grade_T","r_grade_R","r_idh_T","r_idh_R")],
                replicated=as.integer(rep6), qc_flag=as.integer(qcflag)),
          "results/method_framework.csv", row.names=FALSE)

## ---------- figure: AUC by QC metric across replication cohorts ----------
library(tidyr)
pl<-res%>%select(cohort,AUC_grade,AUC_idh,AUC_detect)%>%
  pivot_longer(-cohort,names_to="metric",values_to="AUC")%>%
  mutate(metric=recode(metric,AUC_grade="Grade anchor",AUC_idh="IDH anchor",AUC_detect="Detectability"))
p<-ggplot(pl,aes(cohort,AUC,fill=metric))+geom_col(position="dodge")+
  geom_hline(yintercept=0.5,linetype=2,colour="grey50")+coord_cartesian(ylim=c(0.4,1))+
  labs(x=NULL,y="AUC for predicting cross-cohort replication",fill="Positive-control metric",
       title="Biological positive controls predict cross-cohort replication")+
  theme_bw(base_size=12)+theme(legend.position="top")
ggsave("figures/method_auc.png",p,width=7,height=4.2,dpi=300)
cat("\nwrote figures/method_auc.png\n")
