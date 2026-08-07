# 20_control_table.R -- reviewer fixes M8/M9:
#  (A) EGFR positive-control row for TCGA and all three CGGA datasets computed by a
#      SINGLE documented method (Spearman r with WHO grade; Fisher r-to-z vs TCGA;
#      IDH-mutant HR; grade IV vs grade II HR), replacing per-row ad-hoc computations.
#  (B) EGFR high-level amplification rate by IDH status (the "55% vs 2%" claim).
#
# Grade IV vs II HR uses only grade II and IV patients (not IV vs all others).
# Run from project root:  Rscript R/20_control_table.R   -> results/control_table.txt
suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr); library(survival) })
source("R/_helpers.R")
fread_ok <- requireNamespace("data.table", quietly=TRUE)
rd <- function(f) if (fread_ok) as.data.frame(data.table::fread(f, sep="\t")) else read.delim(f, check.names=FALSE)
out <- file("results/control_table.txt", open="wt")
say <- function(...) { cat(sprintf(...), file=out); cat(sprintf(...)) }

## ---------- TCGA reference (LGG+GBM primary tumours) ----------
tcga <- function(){
  one <- function(project, gradeIV=FALSE){
    se<-load_se(project); cd<-as.data.frame(colData(se)); k<-substr(colnames(se),14,15)=="01"
    e<-log2(assay(se,"tpm_unstrand")[,k]+1); rownames(e)<-rowData(se)$gene_name
    e<-e[!duplicated(rownames(e)),]
    b<-.clinical_df(se,project); b<-b[match(colnames(e),b$barcode),]
    gr<-as.character(cd$paper_Grade[k])
    g<-if(gradeIV) ifelse(grepl("G4|IV",gr),4,NA) else ifelse(grepl("G3|III",gr),3,ifelse(grepl("G2| II",gr),2,NA))
    idh<-as.character(cd$paper_IDH.status[k])
    keep<-!duplicated(b$patient)   # correlation uses all graded primary tumours; HRs drop OS-missing below
    list(e=e[,keep,drop=FALSE], grade=g[keep], time=b$time[keep], event=b$event[keep],
         idh=ifelse(idh[keep]=="Mutant",1,ifelse(idh[keep]=="WT",0,NA)), patient=b$patient[keep])
  }
  L<-one("TCGA-LGG"); G<-one("TCGA-GBM", TRUE)
  cg<-intersect(rownames(L$e),rownames(G$e))
  list(e=cbind(L$e[cg,],G$e[cg,]), grade=c(L$grade,G$grade), time=c(L$time,G$time),
       event=c(L$event,G$event), idh=c(L$idh,G$idh), patient=c(L$patient,G$patient))
}
as_df <- function(M) data.frame(egfr=M$e["EGFR",], grade=M$grade, time=M$time,
                                event=M$event, idh=M$idh, patient=M$patient)
## ---------- a CGGA dataset ----------
cgga <- function(gfile, cfile, cols, logt=TRUE){
  ex<-rd(gfile); rownames(ex)<-ex[[1]]; ex[[1]]<-NULL
  ex<-if(logt) log2(as.matrix(ex)+1) else as.matrix(ex)
  cl<-read.delim(cfile,check.names=FALSE); names(cl)[cols]<-c("id","prs","grade","age","os","censor","idh")
  cl<-cl%>%filter(prs=="Primary")%>%mutate(time=as.numeric(os),event=as.integer(censor),
      grade=as.numeric(factor(grade,levels=c("WHO II","WHO III","WHO IV")))+1,
      idh=ifelse(idh=="Mutant",1,ifelse(idh=="Wildtype",0,NA)))
  sel<-intersect(colnames(ex),cl$id); cl<-cl[match(sel,cl$id),]
  list(e=ex[,sel,drop=FALSE], grade=cl$grade, time=cl$time, event=cl$event, idh=cl$idh, patient=sel)
}

TCm <- tcga(); TC <- as_df(TCm)
D693m <- cgga("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt","data/cgga/cgga_clinical.tsv", c(1,2,4,6,7,8,11))
D325m <- cgga("data/cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt","data/cgga/cgga325_clinical.tsv", c(1,2,4,6,7,8,11))
D301m <- cgga("data/cgga/CGGA.mRNA_array_301_gene_level.20200506.txt","data/cgga/cgga_array_clinical.tsv", c(1,3,5,7,8,9,12), logt=FALSE)

D693 <- as_df(D693m); D325 <- as_df(D325m); D301 <- as_df(D301m)

row <- function(D, lab, ref=NULL){
  D$time[is.na(D$time) | D$time<=0] <- NA          # HRs use OS-evaluable patients only
  ok <- !is.na(D$egfr) & !is.na(D$grade)
  r  <- cor(D$egfr[ok], D$grade[ok], method="spearman"); n <- sum(ok)
  idhHR <- tryCatch(exp(coef(coxph(Surv(time,event)~idh, data=D[!is.na(D$idh)&!is.na(D$time)&!is.na(D$event),]))[1]), error=function(e) NA)
  k24 <- D$grade %in% c(2,4) & !is.na(D$time) & !is.na(D$event)
  gHR <- tryCatch(exp(coef(coxph(Surv(time,event)~I(grade==4), data=D[k24,]))[1]), error=function(e) NA)
  p <- if (is.null(ref)) NA else {
    z <- (atanh(r)-atanh(ref$r))/sqrt(1/(n-3)+1/(ref$n-3)); 2*pnorm(-abs(z)) }
  say("  %-22s r=%+0.3f  n=%3d  p_vs_TCGA=%-9s IDHmut HR=%s  gradeIV-vs-II HR=%s\n",
      lab, r, n, ifelse(is.na(p),"---",sprintf("%.1e",p)),
      ifelse(is.na(idhHR),"---",sprintf("%.2f",idhHR)), ifelse(is.na(gHR),"---",sprintf("%.1f",gHR)))
  list(r=r,n=n)
}
say("=== (A) EGFR positive control, one method for every cohort ===\n")
say("    Spearman r(EGFR expression, WHO grade); Fisher r-to-z vs TCGA;\n")
say("    IDH-mutant HR (univariable); grade IV vs grade II HR (grades II and IV only)\n\n")
ref <- row(TC,"TCGA (reference)")
invisible(row(D693,"CGGA-693 (RNA-seq)", ref))
invisible(row(D325,"CGGA-325 (RNA-seq)", ref))
invisible(row(D301,"CGGA-301 (microarray)", ref))

## ---------- (B) amplification rate by IDH status ----------
say("\n=== (B) EGFR high-level amplification (CN >= 6) by IDH status, TCGA LGG+GBM ===\n")
cn <- do.call(rbind, lapply(c("TCGA-LGG","TCGA-GBM"), function(p){
  f <- sprintf("data/%s_cnv.rds", p); if(!file.exists(f)) return(NULL)
  x <- readRDS(f); m <- assay(x); g <- rowData(x)$gene_name
  i <- which(g=="EGFR")[1]; if(is.na(i)) return(NULL)
  keep <- substr(colnames(x),14,15)=="01"          # primary tumours only (blood normals excluded)
  data.frame(barcode=colnames(x)[keep], cn=as.numeric(m[i,])[keep], project=p)
}))
if (!is.null(cn)) {
  cn$patient <- substr(cn$barcode,1,12)
  TC$pat <- substr(TC$patient,1,12)
  M <- merge(cn, TC[,c("pat","idh")], by.x="patient", by.y="pat")
  M <- M[!is.na(M$idh) & !is.na(M$cn),]
  M <- M %>% group_by(patient) %>% slice_max(cn, n=1, with_ties=FALSE) %>% ungroup() %>% as.data.frame()
  for (v in c(0,1)) {
    s <- M[M$idh==v,]
    say("  IDH-%-8s n=%3d   amplified (CN>=6) = %3d (%.1f%%)\n",
        ifelse(v==1,"mutant","wildtype"), nrow(s), sum(s$cn>=6), 100*mean(s$cn>=6))
  }
} else say("  copy-number object not found\n")

## ---------- (C) per-grade counts + regenerate Figure 1 from the SAME data ----------
say("\n=== (C) per-grade counts (same sets as the table above) ===\n")
cnt <- function(D,lab){ t<-table(D$grade[!is.na(D$egfr)&!is.na(D$grade)])
  say("  %-22s %s (total %d)\n", lab, paste(sprintf("grade %s: %d", names(t), as.integer(t)), collapse=", "), sum(t)) }
cnt(TC,"TCGA"); cnt(D693,"CGGA-693")

suppressPackageStartupMessages(library(ggplot2))
mk <- function(D,lab){ D<-D[!is.na(D$egfr)&!is.na(D$grade),]
  data.frame(cohort=lab, grade=factor(D$grade, levels=2:4, labels=c("II","III","IV")),
             z=as.numeric(scale(D$egfr))) }
P <- rbind(mk(TC,"TCGA (TPM)"), mk(D693,"CGGA (RSEM)"))
S <- P %>% group_by(cohort,grade) %>%
     summarise(m=mean(z), se=sd(z)/sqrt(dplyr::n()), .groups="drop")
g <- ggplot(S, aes(grade, m, colour=cohort, group=cohort)) +
  geom_line(linewidth=0.8) + geom_point(size=2.4) +
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=0.12) +
  labs(x="WHO grade", y="EGFR expression (per-SD units)", colour="Cohort",
       title="Positive control: EGFR expression vs grade",
       subtitle="rises with grade in TCGA but is flat in CGGA") +
  theme_bw(base_size=12) + theme(legend.position="top")
ggsave("figures/positive_control_EGFR_grade.png", g, width=6.4, height=4.2, dpi=300)
say("\n  regenerated figures/positive_control_EGFR_grade.png from these data\n")

## ---------- (D) ten-gene panel on the SAME sample sets as the EGFR row ----------
say("\n=== (D) ten-gene panel: expression-grade Spearman r, TCGA vs CGGA-693 ===\n")
panel <- c("EGFR","IDH1","TP53","CDKN2A","NF1","PTEN","PDGFRA","PIK3R1","PIK3CA","RB1")
rg <- function(M,g){ if(!g %in% rownames(M$e)) return(NA)
  x<-M$e[g,]; ok<-!is.na(x)&!is.na(M$grade); suppressWarnings(cor(x[ok],M$grade[ok],method="spearman")) }
tab <- data.frame(gene=panel,
                  r_TCGA=sapply(panel,function(g) rg(TCm,g)),
                  r_CGGA=sapply(panel,function(g) rg(D693m,g)))
write.csv(tab,"results/panel_correlations.csv",row.names=FALSE)
for(i in seq_len(nrow(tab))) say("  %-7s TCGA %+0.3f   CGGA-693 %+0.3f\n", tab$gene[i], tab$r_TCGA[i], tab$r_CGGA[i])
close(out); cat("\nwrote results/control_table.txt\n")
