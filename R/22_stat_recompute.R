# 22_stat_recompute.R -- independent recomputation of the manuscript's headline
# STATISTICS (not just sample sizes): every Table 1 row, the common-sample refit,
# the likelihood-ratio test, the C-index comparison, and every Table 3 row.
# Uses its own cohort construction and model fitting; shares no code with the
# analysis pipeline. That was NOT true until round 9, when this script still
# sourced R/_helpers.R for load_se()/.clinical_df() -- the same clinical parsing
# 31 pipeline scripts use, so a defect there would have been re-derived here and
# scored as a match. Now uses R/_verify_clinical.R. Do not repoint at _helpers.R.
#
# Run from project root:  Rscript R/22_stat_recompute.R -> results/stat_recompute.txt
suppressPackageStartupMessages({ library(SummarizedExperiment); library(dplyr); library(survival) })
source("R/_verify_clinical.R")
out <- file("results/stat_recompute.txt", open="wt")
say <- function(...) { cat(sprintf(...), file=out); cat(sprintf(...)) }
P <- 0; F <- 0
near <- function(a,b,tol) !is.na(a) && !is.na(b) && abs(a-b) <= tol
chk <- function(lbl, hr, lo, hi, pv, ci, paper){
  got <- sprintf("HR %.2f (%.2f-%.2f) p=%s C=%.3f", hr, lo, hi,
                 ifelse(pv<1e-3, sprintf("%.2e",pv), sprintf("%.3f",pv)), ci)
  ok <- near(hr,paper$hr,0.011) && near(pv,paper$p,max(1e-4,paper$p*0.05)) &&
        (is.na(paper$ci) || near(ci,paper$ci,0.0015))
  if (ok) P <<- P+1 else F <<- F+1
  say("  [%s] %-34s %s   | paper: HR %.2f p=%s C=%s\n", ifelse(ok,"OK  ","FAIL"), lbl, got,
      paper$hr, ifelse(paper$p<1e-3, sprintf("%.2e",paper$p), sprintf("%.3f",paper$p)),
      ifelse(is.na(paper$ci),"--",sprintf("%.3f",paper$ci)))
}

build <- function(project, gradeIV=FALSE){
  se <- vload_se(project); cd <- as.data.frame(colData(se)); k <- substr(colnames(se),14,15)=="01"
  b  <- vclinical(se); b <- b[match(colnames(se)[k], b$barcode),]
  gr <- as.character(cd$paper_Grade[k])
  g  <- if (gradeIV) rep("IV", sum(k)) else ifelse(grepl("G3|III",gr),"III",ifelse(grepl("G2| II",gr),"II",NA))
  tpm<- assay(se,"tpm_unstrand")[,k]; gn <- rowData(se)$gene_name
  d <- data.frame(patient=b$patient, time=b$time, event=b$event, age=b$age,
                  grade=g, idh=as.character(cd$paper_IDH.status[k]),
                  sub=as.character(cd$paper_IDH.codel.subtype[k]),
                  egfr=log2(tpm[which(gn=="EGFR")[1],]+1), stringsAsFactors=FALSE)
  d <- d[!duplicated(d$patient),]
  d[!is.na(d$time) & d$time>0 & !is.na(d$event),]
}
L <- build("TCGA-LGG"); G <- build("TCGA-GBM", TRUE); A <- rbind(L,G)
# per-SD standardisation WITHIN cohort, as the methods specify
L$z <- as.numeric(scale(L$egfr)); G$z <- as.numeric(scale(G$egfr)); A$z <- as.numeric(scale(A$egfr))

fit <- function(d, rhs, lbl, paper){
  f <- as.formula(paste("Surv(time,event) ~ z", rhs))
  m <- coxph(f, data=d); s <- summary(m); ci <- s$conf.int["z",]
  cidx <- unname(s$concordance[1])
  chk(lbl, ci[1], ci[3], ci[4], s$coefficients["z","Pr(>|z|)"], cidx, paper)
  invisible(m)
}
say("=== Table 1: nested Cox models, recomputed independently ===\n")
fit(L, "",                       "LGG unadjusted",      list(hr=1.59,p=6.85e-6,ci=0.592))
fit(L, "+ age",                  "LGG + age",           list(hr=1.37,p=1.06e-3,ci=0.744))
fit(L[!is.na(L$grade),], "+ age + factor(grade)", "LGG + age, grade", list(hr=1.31,p=8.22e-3,ci=0.786))
Lk <- L[!is.na(L$grade) & L$idh %in% c("Mutant","WT"),]
fit(Lk,"+ age + factor(grade) + factor(idh)", "LGG + age, grade, IDH", list(hr=1.13,p=0.164,ci=0.830))
Ls <- L[L$sub %in% c("IDHmut-codel","IDHmut-non-codel","IDHwt"),]
fit(Ls,"+ age + factor(sub)",    "LGG + age, subtype",  list(hr=1.14,p=0.096,ci=0.820))
fit(G, "",                       "GBM unadjusted",      list(hr=0.962,p=0.552,ci=0.523))
Gk <- G[G$idh %in% c("Mutant","WT"),]
fit(Gk,"+ age + factor(idh)",    "GBM + age, IDH",      list(hr=0.866,p=0.030,ci=0.639))
fit(A, "",                       "Pooled unadjusted",   list(hr=1.35,p=1.98e-7,ci=0.553))
Ag <- A; Ag$grade[is.na(Ag$grade) & Ag$patient %in% G$patient] <- "IV"
Ak <- Ag[!is.na(Ag$grade) & Ag$idh %in% c("Mutant","WT"),]
fit(Ak,"+ age + factor(grade) + factor(idh)", "Pooled + age, grade, IDH", list(hr=0.921,p=0.078,ci=0.837))

say("\n=== common-sample refit and incremental value (§6.2) ===\n")
# per §4.2 the common-sample refit standardises WITHIN the complete-case subset
Lc <- Lk; Lc$z <- as.numeric(scale(Lc$egfr))
fit(Lc,"+ age + factor(grade)",  "common-sample + age, grade", list(hr=1.308,p=0.00817,ci=0.786))
m1 <- coxph(Surv(time,event)~z+age+factor(grade)+factor(idh), data=Lk)
m0 <- coxph(Surv(time,event)~  age+factor(grade)+factor(idh), data=Lk)
lr <- 2*(m1$loglik[2]-m0$loglik[2]); pv <- pchisq(lr, 1, lower.tail=FALSE)
say("  LR test EGFR beyond age+grade+IDH: chisq=%.2f df=1 p=%.3f   | paper: 1.98, 0.16\n", lr, pv)
say("  C-index without EGFR=%.3f  with EGFR=%.3f   | paper: 0.834 vs 0.830\n",
    unname(summary(m0)$concordance[1]), unname(summary(m1)$concordance[1]))
say("  PH test (EGFR term, IDH-adjusted) p=%.3f   | paper: 0.004\n", cox.zph(m1)$table["z","p"])

say("\n=== Table 3: within IDH-wildtype ===\n")
W <- A[A$idh %in% "WT",]; W$z <- as.numeric(scale(W$egfr))
fit(W, "",                "IDH-wt unadjusted", list(hr=1.00,p=0.99,ci=NA_real_))
fit(W, "+ age",           "IDH-wt + age",      list(hr=0.93,p=0.22,ci=NA_real_))
Wg <- W[!is.na(W$grade),]
fit(Wg,"+ age + factor(grade)", "IDH-wt + age, grade", list(hr=0.90,p=0.10,ci=NA_real_))
Gw <- G[G$idh %in% "WT",]; Gw$z <- as.numeric(scale(Gw$egfr))
fit(Gw,"+ age",           "IDH-wt GBM + age",  list(hr=0.87,p=0.034,ci=NA_real_))

say("\n=== %d statistical checks: %d passed, %d FAILED ===\n", P+F, P, F)
close(out); cat(sprintf("\n%d passed, %d failed -> results/stat_recompute.txt\n", P, F))
