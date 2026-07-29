# 18_circularity_control.R -- decisive control for the reviewer's circularity
# concern: does the grade-anchor QC predict cross-cohort replication because it
# captures MEASUREMENT FIDELITY, or merely because both "prognostic" status and
# "replication" are driven by the same grade signal the QC measures?
#
# Three tests (TCGA discovery -> CGGA-693 replication):
#   (0) sanity: unadjusted (age-only) prognostic set -> AUC (reproduce ~0.82)
#   (1) SUBTYPE-ADJUSTED: prognostic AFTER grade+IDH adjustment; replication of the
#       adjusted effect. If QC-AUC collapses to ~0.5, the effect was grade-driven.
#   (2) STRATIFIED: AUC within low- vs high-|grade-correlation| genes.
#   Bootstrap 95% CIs throughout.
#
# Run from project root:  Rscript R/18_circularity_control.R
# Writes: results/circularity_control.csv

suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr); library(survival) })
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly=TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep="\t")) else read.delim(f, check.names=FALSE)

## ---- TCGA (LGG+GBM) ----
tcga <- function(){
  one<-function(p,gIV=FALSE){se<-load_se(p);cd<-as.data.frame(colData(se));k<-substr(colnames(se),14,15)=="01"
    e<-log2(assay(se,"tpm_unstrand")[,k]+1);rownames(e)<-rowData(se)$gene_name;e<-e[!duplicated(rownames(e)),]
    b<-.clinical_df(se,p);b<-b[match(colnames(e),b$barcode),]
    g<-if(gIV) rep(4,ncol(e)) else ifelse(grepl("G3|III",cd$paper_Grade[k]),3,2)
    idh<-ifelse(cd$paper_IDH.status[k]=="Mutant",1,ifelse(cd$paper_IDH.status[k]=="WT",0,NA))
    keep<-!is.na(b$time)&b$time>0&!is.na(b$event)&!duplicated(b$patient)
    list(e=e[,keep],time=b$time[keep],event=b$event[keep],age=b$age[keep],grade=g[keep],idh=idh[keep])}
  L<-one("TCGA-LGG");G<-one("TCGA-GBM",TRUE);cg<-intersect(rownames(L$e),rownames(G$e))
  list(e=cbind(L$e[cg,],G$e[cg,]),clin=data.frame(time=c(L$time,G$time),event=c(L$event,G$event),
    age=c(L$age,G$age),grade=c(L$grade,G$grade),idh=c(L$idh,G$idh)))}
cgga <- function(){
  ex<-rd("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt");rownames(ex)<-ex[[1]];ex[[1]]<-NULL;ex<-log2(as.matrix(ex)+1)
  cl<-read.delim("data/cgga/cgga_clinical.tsv",check.names=FALSE);names(cl)[c(1,2,4,6,7,8,11)]<-c("id","prs","grade","age","os","censor","idh")
  cl<-cl%>%filter(prs=="Primary")%>%mutate(time=as.numeric(os),event=as.integer(censor),age=as.numeric(age),
    grade=as.numeric(factor(grade,levels=c("WHO II","WHO III","WHO IV"))),idh=ifelse(idh=="Mutant",1,ifelse(idh=="Wildtype",0,NA)))
  sel<-intersect(colnames(ex),cl$id);cl<-cl[match(sel,cl$id),];list(e=ex[,sel],clin=cl[,c("time","event","age","grade","idh")])}

TC<-tcga(); CG<-cgga()
genes<-intersect(rownames(TC$e),rownames(CG$e)); sdT<-apply(TC$e[genes,],1,sd)
genes<-genes[sdT>0.5]; if(length(genes)>8000) genes<-names(sort(sdT[genes],decreasing=TRUE))[1:8000]
cat("genes:",length(genes),"\n")

# per-gene: unadjusted (age) and subtype-adjusted (age+grade+idh) Cox coef+p, and grade corr
stat<-function(M,gs){cl<-M$clin;e<-M$e; cc<-cl[!is.na(cl$idh),]
  do.call(rbind,lapply(gs,function(g){z<-as.numeric(scale(e[g,]))
    m0<-tryCatch(coxph(Surv(cl$time,cl$event)~z+cl$age),error=function(e)NULL)
    zc<-z[!is.na(cl$idh)]
    m1<-tryCatch(coxph(Surv(cc$time,cc$event)~zc+cc$age+factor(cc$grade)+factor(cc$idh)),error=function(e)NULL)
    b0<-if(is.null(m0))c(NA,NA) else summary(m0)$coef["z",c("coef","Pr(>|z|)")]
    b1<-if(is.null(m1))c(NA,NA) else summary(m1)$coef["zc",c("coef","Pr(>|z|)")]
    data.frame(gene=g,b0=b0[1],p0=b0[2],b1=b1[1],p1=b1[2],
               rg=suppressWarnings(cor(e[g,],cl$grade,method="spearman",use="complete.obs")))}))}
cat("TCGA stats...\n");sT<-stat(TC,genes); cat("CGGA stats...\n");sC<-stat(CG,genes)
D<-merge(sT,sC,by="gene",suffixes=c("_T","_C"))
D$qc<-D$rg_T*D$rg_C                          # grade-anchor QC score
D$p0_BH<-p.adjust(D$p0_T,"BH"); D$p1_BH<-p.adjust(D$p1_T,"BH")

auc<-function(s,l){ok<-!is.na(s)&!is.na(l);s<-s[ok];l<-l[ok];n1<-sum(l==1);n0<-sum(l==0);if(n1==0||n0==0)return(NA)
  r<-rank(s);(sum(r[l==1])-n1*(n1+1)/2)/(n1*n0)}
aucCI<-function(s,l,B=1000){ok<-!is.na(s)&!is.na(l);s<-s[ok];l<-l[ok];a<-auc(s,l)
  bs<-replicate(B,{i<-sample(length(s),replace=TRUE);auc(s[i],l[i])});
  sprintf("%.3f (%.3f-%.3f)",a,quantile(bs,.025,na.rm=TRUE),quantile(bs,.975,na.rm=TRUE))}

# (0) unadjusted sanity
p0<-D$p0_BH<0.05 & !is.na(D$b0_C)
rep0<-as.integer(sign(D$b0_T)==sign(D$b0_C)&D$p0_C<0.05)
# (1) subtype-adjusted: prognostic independent of grade+IDH; replication of adjusted effect
p1<-D$p1_BH<0.05 & !is.na(D$b1_C)
rep1<-as.integer(sign(D$b1_T)==sign(D$b1_C)&D$p1_C<0.05)
# (2) stratify unadjusted-prognostic by |grade corr| in TCGA
lowg <- p0 & abs(D$rg_T)<0.2
higg <- p0 & abs(D$rg_T)>=0.2

cat("\n=== CIRCULARITY CONTROL (TCGA -> CGGA-693, grade-anchor QC) ===\n")
cat(sprintf("(0) unadjusted prognostic  n=%d  replic=%.2f  AUC=%s\n", sum(p0), mean(rep0[p0]), aucCI(D$qc[p0],rep0[p0])))
cat(sprintf("(1) subtype-ADJUSTED prog  n=%d  replic=%.2f  AUC=%s   <-- key\n", sum(p1), mean(rep1[p1]), aucCI(D$qc[p1],rep1[p1])))
cat(sprintf("(2) low grade-corr (|r|<0.2) n=%d replic=%.2f AUC=%s\n", sum(lowg), mean(rep0[lowg]), aucCI(D$qc[lowg],rep0[lowg])))
cat(sprintf("    high grade-corr(|r|>=0.2) n=%d replic=%.2f AUC=%s\n", sum(higg), mean(rep0[higg]), aucCI(D$qc[higg],rep0[higg])))
out<-data.frame(
  analysis=c("unadjusted","subtype_adjusted","low_grade_corr","high_grade_corr"),
  n=c(sum(p0),sum(p1),sum(lowg),sum(higg)),
  replic_rate=round(c(mean(rep0[p0]),mean(rep1[p1]),mean(rep0[lowg]),mean(rep0[higg])),3),
  AUC=c(aucCI(D$qc[p0],rep0[p0]),aucCI(D$qc[p1],rep1[p1]),aucCI(D$qc[lowg],rep0[lowg]),aucCI(D$qc[higg],rep0[higg])))
write.csv(out,"results/circularity_control.csv",row.names=FALSE)
cat("\nwrote results/circularity_control.csv\n")
