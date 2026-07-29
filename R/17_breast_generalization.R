# 17_breast_generalization.R -- does the positive-control framework generalize
# to a SECOND, unrelated tumor type? Discovery = TCGA-BRCA, replication = METABRIC
# (independent cohorts, RNA-seq vs microarray). Anchors: ER status, basal subtype.
#
# Expression + ER/subtype labels: Moanna processed set (Zenodo 4326602).
# Survival: METABRIC (cBioPortal datahub clinical), TCGA-BRCA (datahub clinical).
#
# Run from project root:  Rscript R/17_breast_generalization.R
# Writes: results/breast_method_auc.csv

suppressPackageStartupMessages({ library(dplyr); library(survival) })
fread_ok <- requireNamespace("data.table", quietly=TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep="\t")) else read.delim(f, check.names=FALSE)
dd <- "data/brca/moanna_data"

## ---- expression (EXPR features) + labels ----
load_expr <- function(f){
  m<-rd(f); rownames(m)<-m[[1]]; m[[1]]<-NULL; m<-as.matrix(m)
  ex<-m[, grep("_EXPR$", colnames(m)), drop=FALSE]
  colnames(ex)<-sub("_EXPR$","",colnames(ex)); ex
}
mbE <- rbind(load_expr(file.path(dd,"training/moanna_training_data.tsv")),
             load_expr(file.path(dd,"training/moanna_validation_data.tsv")))
tcE <- load_expr(file.path(dd,"testing/moanna_testing_data.tsv"))
lab<-function(f){l<-rd(f); rownames(l)<-l[[1]]; l[[1]]<-NULL; l}
mbL <- rbind(lab(file.path(dd,"training/moanna_training_label.tsv")),
             lab(file.path(dd,"training/moanna_validation_label.tsv")))
tcL <- lab(file.path(dd,"testing/moanna_testing_label.tsv"))

## ---- survival ----
mbc <- read.delim("data/brca/mb_clin.txt", comment.char="#", check.names=FALSE)
mbc <- data.frame(id=mbc$PATIENT_ID, time=as.numeric(mbc$OS_MONTHS),
                  event=as.integer(grepl("DECEASED", mbc$OS_STATUS)),
                  age=suppressWarnings(as.numeric(mbc$AGE_AT_DIAGNOSIS)))
tcc <- read.delim("data/brca/brca_tcga_clin.txt", comment.char="#", check.names=FALSE)
tcc <- data.frame(id=tcc$PATIENT_ID, time=as.numeric(tcc$OS_MONTHS),
                  event=as.integer(grepl("DECEASED", tcc$OS_STATUS)),
                  age=suppressWarnings(as.numeric(tcc$AGE)))

## ---- assemble a cohort: expression (samples x genes), clinical (time,event,age,ER,basal) ----
assemble <- function(E, L, clin, idmap){
  ids<-rownames(E); pid<-idmap(ids)
  cl<-clin[match(pid, clin$id),]
  data<-list(e=E, time=cl$time, event=cl$event, age=cl$age,
             er=as.numeric(L[ids,"ERStatus"]), basal=as.numeric(L[ids,"BasalNonBasal"]))
  ok<-!is.na(data$time)&data$time>0&!is.na(data$event)
  for(k in c("time","event","age","er","basal")) data[[k]]<-data[[k]][ok]
  data$e<-data$e[ok,,drop=FALSE]; data
}
MB <- assemble(mbE, mbL, mbc, function(x) x)                       # MB-#### match directly
TC <- assemble(tcE, tcL, tcc, function(x) substr(x,1,12))          # TCGA-..-....-01 -> patient
cat(sprintf("TCGA-BRCA (discovery): n=%d events=%d | METABRIC (replication): n=%d events=%d\n",
    length(TC$time),sum(TC$event),length(MB$time),sum(MB$event)))

genes <- intersect(colnames(TC$e), colnames(MB$e))
sdT <- apply(TC$e[,genes],2,sd); genes<-genes[is.finite(sdT)&sdT>0]
cat("shared genes:", length(genes), "\n")

## ---- per-gene stats in a cohort ----
stats<-function(D, gs){
  do.call(rbind, lapply(gs, function(g){
    x<-D$e[,g]; z<-as.numeric(scale(x))
    m<-tryCatch(coxph(Surv(D$time,D$event)~z+D$age), error=function(e) NULL)
    bp<-if(is.null(m)) c(NA,NA) else summary(m)$coefficients["z",c("coef","Pr(>|z|)")]
    data.frame(gene=g, beta=bp[1], p=bp[2],
               r_er=suppressWarnings(cor(x,D$er,method="spearman",use="complete.obs")),
               r_basal=suppressWarnings(cor(x,D$basal,method="spearman",use="complete.obs")))
  }))
}
cat("computing discovery (TCGA-BRCA)...\n"); sT<-stats(TC,genes)
cat("computing replication (METABRIC)...\n"); sM<-stats(MB,genes)
sT$p_BH<-p.adjust(sT$p,"BH")

auc<-function(score,label){ok<-!is.na(score)&!is.na(label);score<-score[ok];label<-label[ok]
  n1<-sum(label==1);n0<-sum(label==0); if(n1==0||n0==0) return(NA)
  r<-rank(score);(sum(r[label==1])-n1*(n1+1)/2)/(n1*n0)}

D<-merge(sT,sM,by="gene",suffixes=c("_T","_R"))
prog<-D$p_BH<0.05 & !is.na(D$beta_R)
Dp<-D[prog,]
replicated<-as.integer(sign(Dp$beta_T)==sign(Dp$beta_R) & Dp$p_R<0.05)
qc_er   <-Dp$r_er_T*Dp$r_er_R
qc_basal<-Dp$r_basal_T*Dp$r_basal_R
effsize <-abs(Dp$beta_T)                       # baseline predictor: discovery effect size
res<-data.frame(pair="TCGA-BRCA -> METABRIC", n_prognostic=nrow(Dp),
                replic_rate=round(mean(replicated),3),
                AUC_ER=round(auc(qc_er,replicated),3),
                AUC_basal=round(auc(qc_basal,replicated),3),
                AUC_effectsize=round(auc(effsize,replicated),3))
write.csv(res,"results/breast_method_auc.csv",row.names=FALSE)
cat("\n=== BREAST generalization: positive-control concordance predicts replication (AUC) ===\n")
print(res,row.names=FALSE)
cat(sprintf("\n(baseline: discovery effect-size |beta| AUC=%.3f)\n", res$AUC_effectsize))
