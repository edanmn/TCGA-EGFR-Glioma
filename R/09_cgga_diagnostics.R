# 09_cgga_diagnostics.R -- positive-control QC of the CGGA external cohort.
#
# Motivation: the naive CGGA analysis (07) suggested an IDH-independent EGFR
# effect with a reversed EGFR-IDH direction. Before believing that, we test
# whether CGGA behaves canonically. Clinical controls (IDH, grade -> survival)
# and expression controls (EGFR should rise with grade and be higher in
# IDH-wildtype, as in TCGA) are checked.
#
# Run from project root:  Rscript R/09_cgga_diagnostics.R
# Writes: results/cgga_positive_controls.csv, figures/positive_control_EGFR_grade.png

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(survival); library(ggplot2)
})
source("R/_helpers.R")
GENES <- c("EGFR","PDGFRA","PTEN","NF1","PIK3CA","PIK3R1","RB1","CDKN2A","TP53","IDH1")

## ---- CGGA ----
ex <- read.delim("data/cgga/cgga_genes.tsv", header=TRUE, row.names=1, check.names=FALSE)
expr <- as.data.frame(t(log2(as.matrix(ex)+1))); expr$id <- rownames(expr)
cl <- read.delim("data/cgga/cgga_clinical.tsv", header=TRUE, check.names=FALSE)
names(cl)[c(1,2,4,6,7,8,11)] <- c("id","prs","grade","age","os","censor","idh")
cg <- expr %>% inner_join(cl, by="id") %>% filter(prs=="Primary")
cg$gnum <- as.numeric(factor(cg$grade, levels=c("WHO II","WHO III","WHO IV")))
cg$time <- as.numeric(cg$os); cg$event <- as.integer(cg$censor)
cg$idhf <- factor(ifelse(cg$idh %in% c("Mutant","Wildtype"), cg$idh, NA), levels=c("Wildtype","Mutant"))

cgi <- cg[!is.na(cg$idhf),]
idh_hr <- summary(coxph(Surv(time,event)~idhf, data=cgi))
grade_hr <- summary(coxph(Surv(time,event)~factor(gnum), data=cg))$conf.int[,1]

## ---- TCGA reference ----
lg <- load_se("TCGA-LGG"); gb <- load_se("TCGA-GBM")
gm <- function(se){tpm<-assay(se,"tpm_unstrand"); k<-substr(colnames(se),14,15)=="01"
  m<-as.data.frame(t(log2(tpm[match(GENES,rowData(se)$gene_name),k]+1))); colnames(m)<-GENES; m}
cdl <- as.data.frame(colData(lg)); kl <- substr(colnames(lg),14,15)=="01"
L <- gm(lg); L$gnum <- ifelse(grepl("G3|III", cdl$paper_Grade[kl]), 3, ifelse(grepl("G2| II", cdl$paper_Grade[kl]), 2, NA))
L$idh <- cdl$paper_IDH.status[kl]
G <- gm(gb); G$gnum <- 4; G$idh <- NA
TC <- bind_rows(L, G)

## ---- positive-control table: per-gene grade Spearman, EGFR-IDH means ----
ctrl <- data.frame(gene=GENES,
  r_grade_TCGA = sapply(GENES, function(g) cor(TC[[g]], TC$gnum, method="spearman", use="complete.obs")),
  r_grade_CGGA = sapply(GENES, function(g) cor(cg[[g]], cg$gnum, method="spearman", use="complete.obs")))
ctrl[,2:3] <- round(ctrl[,2:3],3)
write.csv(ctrl, "results/cgga_positive_controls.csv", row.names=FALSE)

cat("=== CGGA clinical positive controls (should be strong & correct) ===\n")
cat(sprintf("IDH-mutant HR = %.3f (p=%.1e)  [expect <<1]\n", idh_hr$conf.int[1,1], idh_hr$coefficients[1,5]))
cat(sprintf("Grade HR: III vs II = %.2f, IV vs II = %.2f  [expect rising >1]\n", grade_hr[1], grade_hr[2]))
cat("\n=== EGFR expression positive controls: TCGA (canonical) vs CGGA ===\n")
cat(sprintf("EGFR-grade Spearman r:  TCGA=%+.3f   CGGA=%+.3f  [expect positive]\n",
    ctrl$r_grade_TCGA[ctrl$gene=="EGFR"], ctrl$r_grade_CGGA[ctrl$gene=="EGFR"]))
cat("EGFR mean by grade -- TCGA:\n"); print(round(tapply(TC$EGFR, TC$gnum, mean),2))
cat("EGFR mean by grade -- CGGA:\n"); print(round(tapply(cg$EGFR, cg$gnum, mean),2))
cat("\n=== per-gene grade Spearman (TCGA vs CGGA) ===\n"); print(ctrl, row.names=FALSE)

## ---- figure: standardized EGFR by grade, TCGA vs CGGA ----
# Map each cohort's grade encoding to a common II/III/IV label (TCGA gnum=2/3/4,
# CGGA gnum=1/2/3) before combining, so the x-axis is aligned.
TC$glab <- factor(c("2"="II","3"="III","4"="IV")[as.character(TC$gnum)], levels=c("II","III","IV"))
cg$glab <- factor(c("1"="II","2"="III","3"="IV")[as.character(cg$gnum)], levels=c("II","III","IV"))
mk <- function(d, lab){ d$z <- as.numeric(scale(d$EGFR))
  d %>% filter(!is.na(glab)) %>% group_by(glab) %>%
    summarise(mean=mean(z), se=sd(z)/sqrt(n()), .groups="drop") %>% mutate(cohort=lab) }
pd <- bind_rows(mk(TC,"TCGA (TPM)"), mk(cg,"CGGA (RSEM)"))
pd$grade <- pd$glab
p <- ggplot(pd, aes(grade, mean, colour=cohort, group=cohort)) +
  geom_hline(yintercept=0, linetype=3, colour="grey60") +
  geom_line(linewidth=0.8) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.12) +
  geom_point(size=2.6) +
  scale_colour_manual(values=c("TCGA (TPM)"="#1f78b4","CGGA (RSEM)"="#e31a1c")) +
  labs(x="WHO grade", y="EGFR expression (per-cohort SD units)", colour="Cohort",
       title="Positive control: EGFR expression vs grade\nrises in TCGA but is flat in CGGA") +
  theme_bw(base_size=12) + theme(legend.position="top", plot.title=element_text(size=12,face="bold"))
ggsave("figures/positive_control_EGFR_grade.png", p, width=7, height=4.2, dpi=300)
cat("\nwrote figures/positive_control_EGFR_grade.png\n")
