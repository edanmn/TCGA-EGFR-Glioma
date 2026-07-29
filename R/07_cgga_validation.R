# 07_cgga_validation.R -- external validation in the CGGA mRNAseq_693 cohort.
#
# Tests whether the TCGA conclusion replicates: EGFR expression is associated
# with overall survival in lower-grade glioma before molecular adjustment, but
# the association does NOT survive adjustment for IDH status.
#
# Data: CGGA mRNAseq_693 (RSEM gene expression + clinical), primary tumors only.
#   Lower-grade analog = WHO II/III;  GBM analog = WHO IV.
# Run from project root:  Rscript R/07_cgga_validation.R
# Writes: results/cgga_validation.csv

suppressPackageStartupMessages({ library(survival); library(dplyr) })

GENES <- c("EGFR","PDGFRA","PTEN","NF1","PIK3CA","PIK3R1","RB1","CDKN2A","TP53","IDH1")
dd <- "data/cgga"

# ---- clinical (positional column names; header has spaces/parentheses) ----
clin <- read.delim(file.path(dd,"cgga_clinical.tsv"), header=TRUE,
                   stringsAsFactors=FALSE, check.names=FALSE)
names(clin)[1:13] <- c("id","prs","histology","grade","gender","age","os","censor",
                       "radio","chemo","idh","codel","mgmt")
clin <- clin %>%
  mutate(age=suppressWarnings(as.numeric(age)),
         time=suppressWarnings(as.numeric(os)),
         event=suppressWarnings(as.integer(censor)),
         idh_f=factor(ifelse(idh %in% c("Mutant","Wildtype"), idh, NA), levels=c("Wildtype","Mutant")),
         grade_f=factor(gsub("WHO ","",grade), levels=c("II","III","IV")),
         radio=suppressWarnings(as.numeric(radio)),
         chemo=suppressWarnings(as.numeric(chemo)))

# ---- expression: genes x samples -> samples x genes, log2(RSEM+1) ----
ex <- read.delim(file.path(dd,"cgga_genes.tsv"), header=TRUE, row.names=1,
                 stringsAsFactors=FALSE, check.names=FALSE)
expr <- as.data.frame(t(log2(as.matrix(ex)+1)))
expr$id <- rownames(expr)

df <- clin %>% inner_join(expr, by="id") %>%
  filter(prs=="Primary", !is.na(time), time>0, !is.na(event))

make_cohort <- function(d) { for (g in GENES) if (g %in% names(d)) d[[paste0(g,"_z")]] <- as.numeric(scale(d[[g]])); d }
lg  <- make_cohort(df %>% filter(grade_f %in% c("II","III")) %>% mutate(grade_f=droplevels(grade_f)))
gbm <- make_cohort(df %>% filter(grade_f=="IV"))

cat("CGGA primary cohorts: lower-grade(II/III) n=", nrow(lg),
    " events=", sum(lg$event), " | GBM(IV) n=", nrow(gbm), " events=", sum(gbm$event), "\n")
cat("Lower-grade IDH: "); print(table(lg$idh_f, useNA="ifany"))

fit_row <- function(d, rhs, cohort, label, term="EGFR_z") {
  m <- tryCatch(coxph(as.formula(paste("Surv(time,event)~",rhs)), data=d), error=function(e) NULL)
  if (is.null(m) || !term %in% rownames(summary(m)$coefficients)) return(NULL)
  s <- summary(m); ci <- s$conf.int[term,]
  data.frame(cohort=cohort, model=label, HR=signif(unname(ci["exp(coef)"]),3),
             CI_low=signif(unname(ci["lower .95"]),3), CI_high=signif(unname(ci["upper .95"]),3),
             p=signif(s$coefficients[term,"Pr(>|z|)"],3),
             C_index=signif(unname(s$concordance["C"]),3), n=m$n, events=m$nevent, stringsAsFactors=FALSE)
}

rows <- list(
  fit_row(lg,  "EGFR_z", "CGGA LGr(II/III)", "1. unadjusted"),
  fit_row(lg,  "EGFR_z + age + grade_f", "CGGA LGr(II/III)", "2. + age + grade"),
  fit_row(lg[!is.na(lg$idh_f),], "EGFR_z + age + grade_f + idh_f", "CGGA LGr(II/III)", "3. + age + grade + IDH"),
  fit_row(lg[!is.na(lg$idh_f) & !is.na(lg$radio) & !is.na(lg$chemo),],
          "EGFR_z + age + grade_f + idh_f + radio + chemo", "CGGA LGr(II/III)", "4. + treatment (sens.)"),
  fit_row(gbm, "EGFR_z", "CGGA GBM(IV)", "1. unadjusted"),
  fit_row(gbm[!is.na(gbm$idh_f),], "EGFR_z + age + idh_f", "CGGA GBM(IV)", "3. + age + IDH")
)
tbl <- do.call(rbind, Filter(Negate(is.null), rows))
write.csv(tbl, "results/cgga_validation.csv", row.names=FALSE)
cat("\n=== CGGA external validation: EGFR HR per 1 SD ===\n")
print(tbl, row.names=FALSE)
