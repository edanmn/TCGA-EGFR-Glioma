---
title: "Supplementary Material for *Biological positive controls for cross-cohort validation of prognostic gene expression, with a case study of EGFR in glioma*"
author: "Rishik Kondadadi$^{1,2}$  \n$^1$University of Minnesota  \n$^2$Eastview High School"
date: "August 2026"
classoption: twocolumn
geometry: margin=0.85in
fontsize: 10pt
linkcolor: blue
urlcolor: blue
header-includes:
  - \usepackage{graphicx}
  - \usepackage{booktabs}
  - \usepackage{longtable}
  - \usepackage{float}
  - \usepackage{tikz}
  - \usetikzlibrary{arrows.meta,positioning}
  - \makeatletter\def\fps@figure{H}\makeatother
  - \raggedbottom
---

*This document is the supplementary material for the manuscript above. Every table and figure here is cited in the main text; numbering (S1, S2, ...) is continuous with those citations.*

### Supplementary Methods: Computational Environment

All analyses were performed in R version 4.6.0 (2026-04-24) on macOS Tahoe 26.5.1 (Darwin 25.5.1, aarch64-apple-darwin24.6.0). Key package versions: **TCGAbiolinks** 2.40.0 (Bioconductor), **SummarizedExperiment** 1.42.0 (Bioconductor), **DESeq2** 1.52.0 (Bioconductor; VST), **survival** 3.8-6 (CRAN; Cox models, Kaplan-Meier), **survminer** 0.5.2 (CRAN; survival plotting), **dplyr** 1.1.4, and **ggplot2** 3.5.1 (CRAN). Matrix operations used OpenBLAS 0.3.33 (BLAS) and R's internal LAPACK 3.12.1. Full session information is archived with the analysis scripts at https://github.com/edanmn/TCGA-EGFR-Glioma.

### Supplementary Figures

\begin{figure*}[t]
\centering
\resizebox{\textwidth}{!}{%
\begin{tikzpicture}[
  box/.style={draw, rounded corners, align=center, minimum height=11mm, inner sep=5pt, font=\small},
  arr/.style={-{Latex}, thick}]
\node[box] (a) {Data\\ TCGA-LGG/GBM (GDC),\\ CGGA mRNAseq\_693};
\node[box, right=7mm of a] (b) {Quality control\\ primary tumor,\\ 1/patient, OS $>0$};
\node[box, right=7mm of b] (c) {Expression\\ $\log_2(\cdot{+}1)$,\\ per-SD ($z$)};
\node[box, right=7mm of c] (d) {Cox models\\ unadj.\ $\to$ +age\\ $\to$ +grade\\ $\to$ +IDH/subtype};
\node[box, right=7mm of d] (e) {Robustness\\ VST norm.;\\ CGGA validation\\ + positive controls};
\draw[arr] (a)--(b); \draw[arr] (b)--(c); \draw[arr] (c)--(d); \draw[arr] (d)--(e);
\end{tikzpicture}}\\[8pt]
{\textbf{Figure S1.} Analysis pipeline. Standardized expression enters nested Cox models with progressive confounder adjustment, followed by VST-normalization sensitivity and a positive-control–gated CGGA external-validation attempt. Acquisition is cached separately from analysis.}
\end{figure*}

\begin{center}\includegraphics[width=\linewidth]{figures/KM_EGFR_TCGALGG.png}\end{center}

**Figure S2.** Unadjusted KM overall survival by EGFR group (median split), TCGA-LGG.

\begin{center}\includegraphics[width=\linewidth]{figures/KM_EGFR_TCGAGBM.png}\end{center}

**Figure S3.** Unadjusted KM overall survival by EGFR group, TCGA-GBM.

\begin{center}\includegraphics[width=\linewidth]{figures/KM_EGFR_LGGGBM.png}\end{center}

**Figure S4.** Unadjusted KM overall survival by EGFR group, pooled LGG+GBM.

\begin{center}\includegraphics[width=\linewidth]{figures/KM_LGG_subtype.png}\end{center}

**Figure S5.** TCGA-LGG overall survival by IDH/1p19q molecular subtype; subtype, not EGFR, drives the separation in Figure S2.

### Supplementary Tables

\begin{table*}[t]
\centering
{\textbf{Table S1.} Unadjusted Kaplan--Meier median-split results for EGFR in TCGA (descriptive).}\\[4pt]
\begin{tabular}{lllll}
\toprule
Cohort & $n$ low / high & Events & Median survival, low vs.\ high (days) & Log-rank $p$ \\
\midrule
TCGA-LGG & 256 / 255 & 125 & 2,988 vs.\ 2,282 & 0.040 \\
TCGA-GBM & 141 / 141 & 227 & 399 vs.\ 448 & 0.99 \\
Pooled (LGG+GBM) & 397 / 396 & 352 & 1,491 vs.\ 883 & $4.85\times10^{-3}$ \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table*}[t]
\centering
{\textbf{Table S2.} FDR-controlled pathway screen in TCGA-LGG (HR per 1 SD, age-adjusted), with FDR-adjusted $p$ before and after adding \emph{IDH} status.}\\[4pt]
\begin{tabular}{llll}
\toprule
Gene & HR (95\% CI), age-adjusted & $p_{\text{BH}}$ (age) & $p_{\text{BH}}$ (age $+$ IDH) \\
\midrule
RB1    & 1.64 (1.34--2.00) & $1.59\times10^{-5}$ & 0.031 \\
IDH1   & 1.46 (1.20--1.78) & $8.67\times10^{-4}$ & 0.033 \\
EGFR   & 1.37 (1.14--1.66) & $3.54\times10^{-3}$ & 0.178 \\
TP53   & 1.39 (1.13--1.71) & $5.20\times10^{-3}$ & 0.026 \\
PIK3R1 & 0.777 (0.649--0.931) & $1.25\times10^{-2}$ & 0.876 \\
PTEN   & 0.800 (0.666--0.961) & $2.83\times10^{-2}$ & 0.506 \\
CDKN2A & 0.834 (0.691--1.01) & $8.72\times10^{-2}$ & 0.046 \\
PIK3CA & 1.13 (0.950--1.34) & 0.212 & 0.046 \\
PDGFRA & 0.893 (0.752--1.06) & 0.220 & 0.594 \\
NF1    & 0.895 (0.750--1.07) & 0.222 & 0.085 \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table*}[t]
\centering
{\textbf{Table S3.} Normalization sensitivity: EGFR HR per 1 SD under DESeq2 VST (compare to the TPM results in Table 1).}\\[4pt]
\begin{tabular}{lllll}
\toprule
Cohort & Model & HR (95\% CI) & $p$ & C-index \\
\midrule
TCGA-LGG & Unadjusted          & 1.50 (1.22--1.84) & $1.25\times10^{-4}$ & 0.577 \\
TCGA-LGG & $+$ age, grade       & 1.26 (1.04--1.54) & 0.020               & 0.784 \\
TCGA-LGG & $+$ age, grade, IDH  & 1.11 (0.944--1.31) & 0.205              & 0.829 \\
TCGA-GBM & Unadjusted           & 0.945 (0.829--1.08) & 0.397             & 0.524 \\
TCGA-GBM & $+$ age, IDH         & 0.843 (0.737--0.963) & 0.012            & 0.641 \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table*}[t]
\centering
{\textbf{Table S4.} Positive controls: Spearman correlation of gene expression with WHO grade in TCGA vs.\ CGGA. EGFR (and the panel generally) tracks grade in TCGA but not CGGA, indicating attenuated expression–phenotype signal in the CGGA matrix. Correlations are computed on the same sample sets as Table 2, so the EGFR row matches that table exactly ($+0.179$ in TCGA, $n=739$; $-0.033$ in CGGA-693, $n=422$).}\\[4pt]
\begin{tabular}{lcc}
\toprule
Gene & $r$(expression, grade), TCGA & $r$(expression, grade), CGGA \\
\midrule
EGFR   & $+0.179$ & $-0.033$ \\
IDH1   & $+0.348$ & $+0.152$ \\
TP53   & $+0.099$ & $+0.113$ \\
CDKN2A & $-0.131$ & $-0.076$ \\
NF1    & $-0.285$ & $-0.183$ \\
PTEN   & $-0.357$ & $-0.071$ \\
PDGFRA & $-0.304$ & $-0.096$ \\
PIK3R1 & $-0.481$ & $-0.195$ \\
PIK3CA & $-0.072$ & $-0.057$ \\
RB1    & $-0.057$ & $-0.023$ \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table*}[t]
\centering
{\textbf{Table S5.} Pipeline robustness: nested EGFR Cox models (HR per 1 SD) on TCGA re-obtained from recount3 (independent Monorail raw-read reprocessing). Estimates match the GDC-pipeline results (Table 1).}\\[4pt]
\begin{tabular}{lllll}
\toprule
Cohort & Model & HR (95\% CI) & $p$ & $n$ / events \\
\midrule
TCGA-LGG & Unadjusted          & 1.54 (1.25--1.89) & $4.02\times10^{-5}$ & 509 / 125 \\
TCGA-LGG & $+$ age, grade       & 1.27 (1.04--1.55) & 0.021               & 452 / 106 \\
TCGA-LGG & $+$ age, grade, IDH  & 1.11 (0.94--1.31) & 0.220               & 450 / 105 \\
TCGA-GBM & Unadjusted           & 0.906 (0.755--1.09) & 0.289             & 154 / 122 \\
TCGA-GBM & $+$ age              & 0.853 (0.712--1.02) & 0.084             & 154 / 122 \\
TCGA-GBM & $+$ age, grade, IDH  & 0.813 (0.677--0.977) & 0.027            & 150 / 120 \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table*}[t]
\centering
{\textbf{Table S6.} Positive-control-gated cross-cohort replication for genes prognostic in TCGA (FDR $<0.05$), among 500 top-variable genes, in two cohort pairs that share no patients. Replication is far higher for QC-concordant (grade-concordant) than QC-discordant (grade-discordant) genes; Fisher exact $p<10^{-20}$ for both pairs. Both sign-flips observed in the panel fall among QC-discordant genes, though with only two such events this is suggestive rather than decisive.}\\[4pt]
\begin{tabular}{llcc}
\toprule
Cohort pair & Positive control & Genes & Replicated \% (95\% CI) \\
\midrule
TCGA $\to$ CGGA-693 & QC-concordant & 429 & 92.8 (89.9--95.0) \\
TCGA $\to$ CGGA-693 & QC-discordant & 37  & 24.3 (11.8--41.2) \\
TCGA $\to$ CGGA-325 & QC-concordant & 439 & 95.4 (93.1--97.2) \\
TCGA $\to$ CGGA-325 & QC-discordant & 26  & 19.2 (6.6--39.4) \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table*}[t]
\centering
{\textbf{Table S7.} Circularity control (TCGA $\to$ CGGA-693, grade-anchor QC; bootstrap 95\% CIs). The QC's predictive value is partly inflated by grade-correlated genes but persists for subtype-independent associations (grade+IDH-adjusted) and low-grade-correlation genes, excluding a purely tautological explanation. This control uses the TCGA$\cap$CGGA-693 gene universe (6,616 prognostic genes) rather than the four-cohort intersection used for the headline analysis (6,512 genes; Table 4), which is why $n$, the replication rate, and the AUC differ slightly from Table 4.}\\[4pt]
\begin{tabular}{lccc}
\toprule
Prognostic-gene set & Genes & Replication rate & AUC (95\% CI) \\
\midrule
Unadjusted (age only)                   & 6{,}616 & 0.63 & 0.82 (0.81--0.83) \\
Subtype-adjusted (grade + IDH)          & 183  & 0.23 & 0.65 (0.55--0.74) \\
Low grade-correlation ($|r|<0.2$)       & 935  & 0.51 & 0.75 (0.72--0.78) \\
High grade-correlation ($|r|\geq0.2$)   & 5{,}681 & 0.65 & 0.84 (0.83--0.85) \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table*}[t]
\centering
{\textbf{Table S8.} Naive CGGA external-validation models: EGFR HR per 1 SD in primary tumors. Retained for transparency only; these estimates are \emph{not} interpretable, because CGGA EGFR expression fails the positive control (Table 2, Table S4, Figure 1). The apparent IDH-adjusted and treatment-adjusted effects are artifacts of EGFR mis-measurement in CGGA, not evidence of a prognostic effect.}\\[4pt]
\begin{tabular}{lllll}
\toprule
Cohort & Model & HR (95\% CI) & $p$ & C-index \\
\midrule
CGGA II/III & Unadjusted             & 1.17 (0.938--1.46) & 0.164               & 0.520 \\
CGGA II/III & $+$ age, grade          & 1.19 (0.958--1.48) & 0.115               & 0.624 \\
CGGA II/III & $+$ age, grade, IDH     & 1.44 (1.17--1.76) & $4.6\times10^{-4}$   & 0.737 \\
CGGA II/III & $+$ treatment (sens.)   & 1.37 (1.09--1.71) & $6.8\times10^{-3}$   & 0.743 \\
CGGA IV     & $+$ age, IDH            & 1.12 (0.919--1.35) & 0.268               & 0.594 \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table}[t]
\centering
{\textbf{Table S9.} TCGA cohort characteristics and covariate availability after quality control. \emph{IDH} status is available for 508 (LGG) and 265 (GBM) patients and molecular subtype for 508 and 257; counts sum to the cohort total including the unknown category. The grade row counts patients with a \emph{recorded} histologic grade; all TCGA-GBM patients are treated as WHO grade IV by study definition, and one lacks a recorded grade field. Sex was not populated in the expression-object annotation.}\\[4pt]
\small
\begin{tabular}{lll}
\toprule
Characteristic & TCGA-LGG & TCGA-GBM \\
\midrule
Histologic grade & II--III & IV \\
Patients & 511 & 282 \\
Deaths (events) & 125 & 227 \\
Median age, years (IQR) & 41 (32--53) & 60 (51--69) \\
\emph{IDH} (Mutant / WT / unknown) & 414 / 94 / 3 & 21 / 244 / 17 \\
Subtype (codel / non-codel / wt) & 166 / 248 / 94 & --- / 19 / 238 \\
Histologic grade recorded & 454 & 281 \\
\bottomrule
\end{tabular}
\end{table}

