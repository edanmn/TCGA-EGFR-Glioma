# 21_assertions.R -- independent consistency audit of every cohort size and
# category count the manuscript asserts. Rebuilds cohorts from the cached objects
# WITHOUT reusing any analysis script's helper logic, then checks each claim.
#
# Run from project root:  Rscript R/21_assertions.R  -> results/assertions.txt
suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr) })
source("R/_helpers.R")
out <- file("results/assertions.txt", open="wt")
say <- function(...) { cat(sprintf(...), file=out); cat(sprintf(...)) }
PASS <- 0; FAIL <- 0
chk <- function(label, got, want){
  ok <- isTRUE(all.equal(got, want))
  if (ok) PASS <<- PASS+1 else FAIL <<- FAIL+1
  say("  [%s] %-52s got=%s  paper=%s\n", ifelse(ok,"OK  ","FAIL"), label, paste(got,collapse="/"), paste(want,collapse="/"))
}

## ---- canonical rebuild: primary tumours, one per patient ----
build <- function(project, gradeIV=FALSE){
  se <- load_se(project); cd <- as.data.frame(colData(se))
  k  <- substr(colnames(se),14,15)=="01"
  b  <- .clinical_df(se, project); b <- b[match(colnames(se)[k], b$barcode),]
  gr <- as.character(cd$paper_Grade[k])
  g  <- if (gradeIV) ifelse(grepl("G4|IV",gr),4,NA) else ifelse(grepl("G3|III",gr),3,ifelse(grepl("G2| II",gr),2,NA))
  idh<- as.character(cd$paper_IDH.status[k])
  sub<- as.character(cd$paper_IDH.codel.subtype[k])
  tpm<- assay(se,"tpm_unstrand")[,k]; gn <- rowData(se)$gene_name
  d <- data.frame(patient=b$patient, time=b$time, event=b$event, age=b$age,
                  grade=g, idh=idh, subtype=sub, egfr=log2(tpm[which(gn=="EGFR")[1],]+1),
                  stringsAsFactors=FALSE)
  d[!duplicated(d$patient),]
}
L <- build("TCGA-LGG"); G <- build("TCGA-GBM", TRUE)
os <- function(d) d[!is.na(d$time) & d$time>0 & !is.na(d$event),]
Lo <- os(L); Go <- os(G); P <- rbind(Lo,Go)

say("=== TCGA cohort sizes (survival-evaluable, primary, one per patient) ===\n")
chk("TCGA-LGG patients / deaths", c(nrow(Lo), sum(Lo$event)), c(511,125))
chk("TCGA-GBM patients / deaths", c(nrow(Go), sum(Go$event)), c(282,227))
chk("pooled patients / deaths",   c(nrow(P),  sum(P$event)),  c(793,352))

say("\n=== Table S9 category counts ===\n")
chk("LGG IDH Mutant/WT/unknown", c(sum(Lo$idh=="Mutant",na.rm=TRUE), sum(Lo$idh=="WT",na.rm=TRUE),
                                   sum(!Lo$idh %in% c("Mutant","WT"))), c(414,94,3))
chk("GBM IDH Mutant/WT/unknown", c(sum(Go$idh=="Mutant",na.rm=TRUE), sum(Go$idh=="WT",na.rm=TRUE),
                                   sum(!Go$idh %in% c("Mutant","WT"))), c(21,244,17))
chk("LGG subtype codel/non-codel/wt",
    c(sum(Lo$subtype=="IDHmut-codel",na.rm=TRUE), sum(Lo$subtype=="IDHmut-non-codel",na.rm=TRUE),
      sum(Lo$subtype=="IDHwt",na.rm=TRUE)), c(166,248,94))
chk("GBM subtype non-codel/wt",
    c(sum(Go$subtype=="IDHmut-non-codel",na.rm=TRUE), sum(Go$subtype=="IDHwt",na.rm=TRUE)), c(19,238))
chk("LGG histologic grade available", sum(!is.na(Lo$grade)), 454)
chk("GBM histologic grade RECORDED", sum(!is.na(Go$grade)), 281)   # all GBM treated as IV by study definition

say("\n=== nested-model sample sizes (Table 1) ===\n")
chk("LGG unadjusted / +age n", nrow(Lo), 511)
chk("LGG +age+grade n / events", c(sum(!is.na(Lo$grade)), sum(Lo$event[!is.na(Lo$grade)])), c(454,106))
kg <- !is.na(Lo$grade) & Lo$idh %in% c("Mutant","WT")
chk("LGG +age+grade+IDH n / events", c(sum(kg), sum(Lo$event[kg])), c(452,105))
ks <- Lo$subtype %in% c("IDHmut-codel","IDHmut-non-codel","IDHwt")
chk("LGG +age+subtype n / events", c(sum(ks), sum(Lo$event[ks])), c(508,123))
chk("GBM +age+IDH n / events", c(sum(Go$idh %in% c("Mutant","WT")), sum(Go$event[Go$idh %in% c("Mutant","WT")])), c(265,214))
# pooled model uses the study-definition grade for GBM (all IV), per §4.1
Pg <- P; Pg$grade[is.na(Pg$grade) & Pg$patient %in% Go$patient] <- 4
kp <- !is.na(Pg$grade) & Pg$idh %in% c("Mutant","WT")
chk("pooled +age+grade+IDH n / events", c(sum(kp), sum(Pg$event[kp])), c(717,319))

say("\n=== within-IDH-wildtype (Table 3) ===\n")
W <- P[P$idh %in% "WT",]
chk("IDH-wt n / deaths", c(nrow(W), sum(W$event)), c(338,253))
chk("IDH-wt +grade n / events", c(sum(!is.na(W$grade)), sum(W$event[!is.na(W$grade)])), c(330,246))
Gw <- Go[Go$idh %in% "WT",]
chk("IDH-wt GBM n / events", c(nrow(Gw), sum(Gw$event)), c(244,203))

say("\n=== positive-control sets (Table 2, Figure 1) ===\n")
TCall <- rbind(L,G)
chk("TCGA graded primary tumours (Table 2 n)", sum(!is.na(TCall$grade)), 739)
chk("TCGA per-grade II/III/IV", as.integer(table(TCall$grade)), c(216,241,282))

## ---- CGGA-693 ----
cl <- read.delim("data/cgga/cgga_clinical.tsv", check.names=FALSE)
names(cl)[c(1,2,4,6,7,8,11)] <- c("id","prs","grade","age","os","censor","idh")
cl <- cl %>% filter(prs=="Primary") %>%
  mutate(gnum=as.numeric(factor(grade, levels=c("WHO II","WHO III","WHO IV")))+1,
         time=suppressWarnings(as.numeric(os)), event=suppressWarnings(as.integer(censor)))
ex <- data.table::fread("data/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", sep="\t")
ids <- intersect(names(ex)[-1], cl$id); cl <- cl[match(ids, cl$id),]
say("\n=== CGGA-693 (Table 2, §5) ===\n")
chk("primary tumours with grade", sum(!is.na(cl$gnum)), 422)
chk("per-grade II/III/IV", as.integer(table(cl$gnum)), c(138,144,140))
ok <- !is.na(cl$time) & cl$time>0 & !is.na(cl$event)
chk("OS-evaluable", sum(ok), 404)
chk("lower-grade (II/III) n / deaths", c(sum(ok & cl$gnum<4), sum(cl$event[ok & cl$gnum<4])), c(271,99))
chk("glioblastoma (IV) n / deaths",    c(sum(ok & cl$gnum==4), sum(cl$event[ok & cl$gnum==4])), c(133,110))

## ---- copy number ----
say("\n=== EGFR copy number (§6.6) ===\n")
cn <- do.call(rbind, lapply(c("TCGA-LGG","TCGA-GBM"), function(p){
  x <- readRDS(sprintf("data/%s_cnv.rds", p)); i <- which(rowData(x)$gene_name=="EGFR")[1]
  kk <- substr(colnames(x),14,15)=="01"
  data.frame(patient=substr(colnames(x)[kk],1,12), cn=as.numeric(assay(x)[i,])[kk])
}))
cn <- cn[!is.na(cn$cn),] %>% group_by(patient) %>% summarise(cn=max(cn), .groups="drop")
M  <- merge(rbind(L,G)[,c("patient","idh","time","event")], cn, by="patient")
wt <- M[M$idh %in% "WT",]; mu <- M[M$idh %in% "Mutant",]
chk("IDH-wt with CN / amplified %", c(nrow(wt), round(100*mean(wt$cn>=6),1)), c(167,55.1))
chk("IDH-mutant with CN / amplified %", c(nrow(mu), round(100*mean(mu$cn>=6),1)), c(222,2.3))
wto <- wt[!is.na(wt$time) & wt$time>0 & !is.na(wt$event),]
chk("IDH-wt CN + survival (Cox n / events)", c(nrow(wto), sum(wto$event)), c(166,127))

say("\n=== %d checks: %d passed, %d FAILED ===\n", PASS+FAIL, PASS, FAIL)
close(out); cat(sprintf("\n%d passed, %d failed -> results/assertions.txt\n", PASS, FAIL))
