# 19_overlap_and_baselines.R -- analyses added in response to peer review:
#  (A) CGGA cross-dataset PATIENT OVERLAP, and a de-duplicated array-301 replication
#      (array-301 restricted to patients absent from CGGA-325) -- tests whether the
#      cross-platform result depends on shared patients.
#  (B) precision baseline: does replication-cohort ESTIMATION PRECISION (SE of the
#      Cox coefficient) predict replication as well as the biological anchor?
#      Guards against the "QC is just an attenuation proxy" objection.
#  (C) proportional-hazards sensitivity for the headline TCGA-LGG model (time-varying
#      EGFR effect; landmark analysis; period-specific HRs).
#  (D) TCGA IDH-status category counts for the cohort table.
#
# Run from project root:  Rscript R/19_overlap_and_baselines.R
# Writes: results/review_fixes.txt, results/dedup_array301.csv

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival)
})
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly=TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep="\t")) else read.delim(f, check.names=FALSE)
out <- file("results/review_fixes.txt", open="wt")
say <- function(...) { cat(sprintf(...), file=out); cat(sprintf(...)) }

## ---------------- cohort assembly (mirrors 16_method_framework.R, keeping IDs) ----------------
tcga_mat <- function(){
  one <- function(project, gradeIV=FALSE){
    se<-load_se(project); cd<-as.data.frame(colData(se)); k<-substr(colnames(se),14,15)=="01"
    e<-log2(assay(se,"tpm_unstrand")[,k]+1); rownames(e)<-rowData(se)$gene_name
    e<-e[!duplicated(rownames(e)),]
    b<-.clinical_df(se,project); b<-b[match(colnames(e),b$barcode),]
    gr<-as.character(cd$paper_Grade[k])
    g<-if(gradeIV) ifelse(grepl("G4|IV",gr),4,NA) else ifelse(grepl("G3|III",gr),3,ifelse(grepl("G2| II",gr),2,NA))
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
# NOTE: the mRNA-array_301 matrix is already on a log-ratio scale and contains
# negative values; applying log2(x+1) to it (as earlier scripts did) yields NaN for
# every value <= -1 and silently drops those samples. logt=FALSE for that cohort.
cgga_mat <- function(gfile, cfile, cols, logt=TRUE){
  ex<-rd(gfile); rownames(ex)<-ex[[1]]; ex[[1]]<-NULL
  ex<-if(logt) log2(as.matrix(ex)+1) else as.matrix(ex)
  cl<-read.delim(cfile,check.names=FALSE); names(cl)[cols]<-c("id","prs","grade","age","os","censor","idh")
  cl<-cl%>%filter(prs=="Primary")%>%mutate(time=as.numeric(os),event=as.integer(censor),age=as.numeric(age),
    grade=as.numeric(factor(grade,levels=c("WHO II","WHO III","WHO IV"))),
    idh=ifelse(idh=="Mutant",1,ifelse(idh=="Wildtype",0,NA)))
  sel<-intersect(colnames(ex),cl$id); cl<-cl[match(sel,cl$id),]
  list(e=ex[,sel], clin=cl[,c("time","event","age","grade","idh")], ids=sel)
}

TC  <- tcga_mat()
C693<- cgga_mat("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt","data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
C325<- cgga_mat("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt","data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
CARR<- cgga_mat("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt","data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt=FALSE)

## ================= (A) patient overlap between CGGA datasets =================
say("=== (A) CGGA cross-dataset patient overlap (analysed primary-tumour sets) ===\n")
say("  n: CGGA-693=%d  CGGA-325=%d  array-301=%d\n", length(C693$ids), length(C325$ids), length(CARR$ids))
say("  693 n 325 = %d\n", length(intersect(C693$ids, C325$ids)))
say("  693 n 301 = %d\n", length(intersect(C693$ids, CARR$ids)))
say("  325 n 301 = %d  (%.0f%% of the array-301 analysis set)\n",
    length(intersect(C325$ids, CARR$ids)),
    100*length(intersect(C325$ids, CARR$ids))/length(CARR$ids))

# de-duplicated array-301: drop any patient also present in CGGA-325 or CGGA-693
dup <- union(intersect(CARR$ids, C325$ids), intersect(CARR$ids, C693$ids))
keep <- !(CARR$ids %in% dup)
CARRD <- list(e=CARR$e[,keep], clin=CARR$clin[keep,], ids=CARR$ids[keep])
say("  de-duplicated array-301: %d -> %d patients (dropped %d)\n\n", length(CARR$ids), length(CARRD$ids), sum(!keep))

## ---------------- gene universe (identical rule to script 16) ----------------
genes<-Reduce(intersect,list(rownames(TC$e),rownames(C693$e),rownames(C325$e),rownames(CARR$e)))
sdT<-apply(TC$e[genes,],1,sd); genes<-genes[sdT>0.5]
if(length(genes)>8000) genes<-names(sort(sdT[genes],decreasing=TRUE))[1:8000]
say("Gene universe: %d\n\n", length(genes))

## ---------------- per-gene stats, now also returning the coefficient SE ----------------
stats_cohort <- function(M, gs){
  cl<-M$clin; e<-M$e
  do.call(rbind, lapply(gs, function(g){
    x<-e[g,]; z<-as.numeric(scale(x))
    m<-tryCatch(coxph(Surv(cl$time,cl$event)~z+cl$age), error=function(e) NULL)
    co<-if(is.null(m)) c(NA,NA,NA) else summary(m)$coefficients["z",c("coef","se(coef)","Pr(>|z|)")]
    data.frame(gene=g, beta=co[1], se=co[2], p=co[3],
               r_grade=suppressWarnings(cor(x,cl$grade,method="spearman",use="complete.obs")),
               r_idh=suppressWarnings(cor(x,cl$idh,method="spearman",use="complete.obs")),
               meanx=mean(x,na.rm=TRUE))
  }))
}
cat("computing TCGA...\n");            sT<-stats_cohort(TC,genes); sT$p_BH<-p.adjust(sT$p,"BH")
cat("computing CGGA-693...\n");        s6<-stats_cohort(C693,genes)
cat("computing array-301 (full)...\n");sA<-stats_cohort(CARR,genes)
cat("computing array-301 (dedup)...\n");sD<-stats_cohort(CARRD,genes)

auc<-function(score,label){ ok<-!is.na(score)&!is.na(label); score<-score[ok]; label<-label[ok]
  n1<-sum(label==1); n0<-sum(label==0); if(n1==0||n0==0) return(NA)
  r<-rank(score); (sum(r[label==1])-n1*(n1+1)/2)/(n1*n0) }
aucCI<-function(score,label,B=1000){ a<-auc(score,label)
  ok<-!is.na(score)&!is.na(label); s<-score[ok]; l<-label[ok]
  bs<-replicate(B,{i<-sample(length(s),replace=TRUE); auc(s[i],l[i])})
  sprintf("%.2f (%.2f-%.2f)", a, quantile(bs,.025,na.rm=TRUE), quantile(bs,.975,na.rm=TRUE)) }

evaluate<-function(sR, lab){
  D<-merge(sT, sR, by="gene", suffixes=c("_T","_R"))
  D<-D[D$p_BH<0.05 & !is.na(D$beta_R),]
  rep<-as.integer(sign(D$beta_T)==sign(D$beta_R) & D$p_R<0.05)
  data.frame(cohort=lab, n=nrow(D), replic_rate=round(mean(rep),3),
             AUC_grade=aucCI(D$r_grade_T*D$r_grade_R, rep),
             AUC_idh  =aucCI(D$r_idh_T*D$r_idh_R,     rep),
             AUC_detect=aucCI(pmin(D$meanx_T,D$meanx_R), rep),
             AUC_precision=aucCI(-D$se_R, rep))   # (B) higher precision = smaller SE
}
res<-rbind(evaluate(s6,"CGGA-693"), evaluate(sA,"array-301 (all)"), evaluate(sD,"array-301 (de-duplicated)"))
write.csv(res,"results/dedup_array301.csv",row.names=FALSE)
say("=== (A)+(B) replication prediction, with de-duplication and a precision baseline ===\n")
say("%s\n", paste(capture.output(print(res,row.names=FALSE)), collapse="\n"))

## ---------------- (A) EGFR positive control in the de-duplicated array-301 ----------------
egfr_ctrl <- function(M,lab){
  if(!"EGFR" %in% rownames(M$e)) { say("  %s: EGFR absent\n",lab); return(invisible()) }
  x<-M$e["EGFR",]; cl<-M$clin
  r<-suppressWarnings(cor(x,cl$grade,method="spearman",use="complete.obs"))
  idhHR<-tryCatch(exp(coef(coxph(Surv(cl$time,cl$event)~cl$idh))[1]), error=function(e) NA)
  k24<-cl$grade %in% c(1,3); g4<-as.integer(cl$grade[k24]==3)
  gHR<-tryCatch(exp(coef(coxph(Surv(cl$time[k24],cl$event[k24])~g4))[1]), error=function(e) NA)
  say("  %-28s r(EGFR,grade)=%+.3f  n=%d  IDH-mut HR=%.2f  gradeIV-vs-II HR=%.1f\n",
      lab, r, sum(!is.na(x)), idhHR, gHR)
  invisible(r)
}
say("\n=== (A) EGFR positive control, array-301 before vs after de-duplication ===\n")
rA<-egfr_ctrl(CARR,"array-301 (all)"); rD<-egfr_ctrl(CARRD,"array-301 (de-duplicated)")
# Fisher r-to-z vs TCGA reference (r=+0.185, n=826)
fz<-function(r1,n1,r2,n2){ z<-(atanh(r1)-atanh(r2))/sqrt(1/(n1-3)+1/(n2-3)); 2*pnorm(-abs(z)) }
if(!is.null(rD)) say("  Fisher r-to-z vs TCGA (r=+0.185,n=826): de-duplicated p=%.2g\n",
                     fz(0.185,826,rD,length(CARRD$ids)))

## ================= (C) proportional-hazards sensitivity, headline LGG model =================
say("\n=== (C) PH sensitivity for the headline TCGA-LGG model (EGFR + age + grade + IDH) ===\n")
se<-load_se("TCGA-LGG"); cd<-as.data.frame(colData(se))
tpm<-assay(se,"tpm_unstrand"); eg<-log2(tpm[match("EGFR",rowData(se)$gene_name),]+1)
b<-.clinical_df(se,"TCGA-LGG"); b$idh<-as.character(cd$paper_IDH.status); b$EGFR<-eg[match(b$barcode,colnames(se))]
tg<-function(g){g<-toupper(as.character(g));o<-rep(NA,length(g));o[grepl("G2| II ",g)]<-"II";o[grepl("G3| III ",g)]<-"III";factor(o,levels=c("II","III"))}
d<-b%>%filter(sample_type=="01",!is.na(time),time>0,!is.na(event))%>%arrange(patient)%>%distinct(patient,.keep_all=TRUE)
d$grade_f<-tg(d$grade); d$z<-as.numeric(scale(d$EGFR))
d<-d%>%filter(!is.na(grade_f), idh %in% c("Mutant","WT"), !is.na(age))
say("  analysis set: n=%d events=%d\n", nrow(d), sum(d$event))

m0<-coxph(Surv(time,event)~z+age+grade_f+idh, data=d)
s0<-summary(m0)$coefficients["z",]
say("  standard Cox            HR=%.2f (%.2f-%.2f) p=%.3f   [PH p=%.3f]\n",
    exp(s0[1]), exp(s0[1]-1.96*s0[3]), exp(s0[1]+1.96*s0[3]), s0[5],
    cox.zph(m0)$table["z","p"])

# time-varying EGFR effect
mtt<-coxph(Surv(time,event)~tt(z)+z+age+grade_f+idh, data=d,
           tt=function(x,t,...) x*log(t+1))
ctt<-summary(mtt)$coefficients
say("  + time-varying EGFR     EGFR main p=%.3f ; EGFR x log(t) p=%.3f (no time trend if n.s.)\n",
    ctt["z","Pr(>|z|)"], ctt["tt(z)","Pr(>|z|)"])

# landmark: restrict follow-up to first 3 years
d3<-d; d3$event<-ifelse(d3$time>1095,0,d3$event); d3$time<-pmin(d3$time,1095)
m3<-coxph(Surv(time,event)~z+age+grade_f+idh,data=d3); s3<-summary(m3)$coefficients["z",]
say("  landmark 0-3 years      HR=%.2f (%.2f-%.2f) p=%.3f  events=%d\n",
    exp(s3[1]), exp(s3[1]-1.96*s3[3]), exp(s3[1]+1.96*s3[3]), s3[5], sum(d3$event))

# period-specific HRs via survSplit
sp<-survSplit(Surv(time,event)~., data=d, cut=1095, episode="period")
for(p in sort(unique(sp$period))){
  dd<-sp[sp$period==p,]
  mm<-tryCatch(coxph(Surv(tstart,time,event)~z+age+grade_f+idh,data=dd),error=function(e)NULL)
  if(!is.null(mm)){ ss<-summary(mm)$coefficients["z",]
    say("  period %d (%s)        HR=%.2f (%.2f-%.2f) p=%.3f  events=%d\n", p,
        ifelse(p==1,"0-3y ","3y+  "), exp(ss[1]), exp(ss[1]-1.96*ss[3]), exp(ss[1]+1.96*ss[3]), ss[5], sum(dd$event)) }
}

## ================= (D) TCGA IDH-status category counts for the cohort table =================
say("\n=== (D) TCGA IDH-status categories (primary tumours, one per patient) ===\n")
for(pr in c("TCGA-LGG","TCGA-GBM")){
  se<-load_se(pr); cd<-as.data.frame(colData(se))
  bb<-.clinical_df(se,pr); bb$idh<-as.character(cd$paper_IDH.status)
  bb$sub<-as.character(cd$paper_IDH.codel.subtype)
  bb<-bb%>%filter(sample_type=="01",!is.na(time),time>0,!is.na(event))%>%arrange(patient)%>%distinct(patient,.keep_all=TRUE)
  t1<-table(bb$idh, useNA="ifany"); t2<-table(bb$sub, useNA="ifany")
  say("  %s  n=%d\n", pr, nrow(bb))
  say("    IDH status : %s\n", paste(sprintf("%s=%d", names(t1), as.integer(t1)), collapse="  "))
  say("    subtype    : %s\n", paste(sprintf("%s=%d", names(t2), as.integer(t2)), collapse="  "))
}
close(out)
cat("\nwrote results/review_fixes.txt and results/dedup_array301.csv\n")
