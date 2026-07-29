---
title: "Biological Positive Controls Detect Measurement-Artifact Discordances in Cross-Cohort Prognostic Analysis: A Framework Motivated by EGFR in Glioma"
author: "[TODO: Author name(s)]  \n[TODO: Affiliation]"
date: "July 2026"
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
---

## Abstract

Prognostic gene associations from one cohort frequently fail to replicate in another, and it is hard to know in advance which will. We propose a simple safeguard—a *biological positive control*—and show it predicts replicability at genome scale. The idea: a cross-cohort prognostic comparison for a gene is trustworthy only if both cohorts measure that gene faithfully, which can be checked by whether the gene's expression tracks a known biological covariate (here, tumor grade or *IDH* status) consistently across cohorts. Using The Cancer Genome Atlas (TCGA) as the discovery cohort and three independent Chinese Glioma Genome Atlas (CGGA) datasets (two RNA-seq batches, one microarray) for replication, we evaluated 6,512 TCGA-prognostic genes. A gene's positive-control concordance predicted whether its association replicated with **AUC $0.82$–$0.93$ (grade anchor) and $0.77$–$0.84$ (IDH anchor) across all three cohorts and both platforms**, whereas a naive detectability baseline was at or below chance (AUC $0.36$–$0.56$). Testing a second tumor type (breast; TCGA-BRCA$\leftrightarrow$METABRIC) delineated the method's scope: there the two cohorts measured genes concordantly (anchor agreement $r=0.84$) and the positive control did *not* predict replication (AUC $\approx0.55$–$0.60$), showing that the gate helps specifically when cross-cohort measurement fidelity varies—its intended use—rather than when non-replication is driven by limited power in well-measured cohorts. The framework was motivated by EGFR, a canonical case in which naive analysis produces a spurious cross-cohort discordance: in TCGA, EGFR expression is prognostic only until *IDH* is modeled (HR $1.31\to1.13$, $p=0.16$; no added discrimination; robust to normalization and to an independent raw-read reprocessing), and it carries no signal *within* IDH-wildtype glioma by expression (HR $\approx1.0$) or amplification (continuous copy number HR $1.00$, $0.84$–$1.20$; binary tests underpowered). EGFR expression fails the positive control reproducibly across all three CGGA datasets (Fisher $r$-to-$z$ $p<10^{-3}$ vs.\ TCGA) despite sound CGGA clinical data—so its apparent IDH-independent effect in CGGA is an artifact the positive control flags. We conclude that biological positive controls are an effective gate for detecting measurement-artifact-driven cross-cohort discordances—their designed purpose—with a scope bounded by measurement heterogeneity, and that EGFR is not an independent prognostic factor in glioma beyond marking the IDH-wildtype state.

---

## 1. Introduction

Diffuse gliomas vary widely in outcome, and their prognosis is defined primarily by molecular features—*IDH* mutation status and 1p/19q codeletion—rather than histology alone [5]. Against this backdrop, individual genes are frequently proposed as transcriptomic prognostic markers. EGFR is a natural candidate: a defining oncogene of the classical glioblastoma subtype and one of the most frequently altered genes in glioma [1, 6].

Any single-gene prognostic claim in glioma faces a specific threat: because molecular subtype governs both a tumor's transcriptional program and its outcome, a gene whose expression tracks subtype will appear prognostic even if it carries no information beyond subtype. A credible claim must therefore (i) survive adjustment for *IDH* status and molecular subtype, (ii) be robust to expression normalization, and (iii) replicate in an independent cohort—with the replication cohort itself verified to measure the gene faithfully.

Our primary contribution is **methodological**: a biological positive-control framework for cross-cohort prognostic validation, evaluated at genome scale (6,512 genes) across three replication cohorts and two platforms, in which a gene's expression–covariate concordance predicts whether its prognostic association replicates (AUC $0.82$–$0.93$), outperforming a detectability baseline and generalizing across two biological anchors (grade, *IDH*). We also delineate the method's **scope** in a second tumor type (breast, TCGA-BRCA$\leftrightarrow$METABRIC): where the two cohorts measure genes concordantly, the gate has little artifact to flag and does not predict replication—clarifying that the method targets measurement-fidelity–driven discordance specifically. The framework is motivated and stress-tested by a **case study of EGFR**, for which we additionally provide (i) nested Cox models with discrimination and proportional-hazards diagnostics showing the association vanishes after *IDH* adjustment, robust to normalization and to an independent raw-read reprocessing; (ii) a within-IDH-wildtype analysis of both expression and copy-number amplification showing no within-subtype signal; and (iii) a demonstration that EGFR's apparent cross-cohort discordance is a positive-control failure, not biology. The result pairs a reusable methodological safeguard with a definitive negative characterization of a widely studied marker.

The paper proceeds through background (§2), problem statement (§3), methods and pipeline (§4–5), setup (§6), results (§7), and discussion, limitations, future work, and conclusions (§8–11).

---

## 2. Background and Related Work

### 2.1 TCGA, CGGA, and the Genomic Data Commons

TCGA molecularly characterized thousands of tumors with matched survival data [1, 2, 3]; its glioma studies—glioblastoma [1] and lower-grade glioma [2]—include curated molecular classifications (*IDH* status, 1p/19q codeletion, subtype) [2]. CGGA is an independent resource providing RNA-sequencing and clinical data, including *IDH* status and treatment, for a large Chinese glioma cohort [17].

### 2.2 Molecular classification of glioma

Under the 2021 WHO classification, adult diffuse gliomas are defined molecularly, and *IDH* mutation is strongly favorable [5, 7]. The histologic label "lower-grade glioma" spans prognostically distinct molecular entities.

### 2.3 EGFR in glioma

EGFR amplification, mutation (EGFRvIII), and overexpression are recurrent in glioblastoma and mark the classical subtype [1, 6]. EGFR expression is characteristically higher in higher-grade and IDH-wildtype tumors—a fact we exploit as a positive control for whether a cohort measures EGFR faithfully.

### 2.4 Methodology

We use Kaplan–Meier estimation with the log-rank test [9, 10], Cox regression [11], Harrell's C-index [16], proportional-hazards (PH) testing via scaled Schoenfeld residuals [13], the Benjamini–Hochberg FDR procedure [12], and DESeq2 variance-stabilizing transformation [18].

### 2.5 Related work

Cross-study irreproducibility of transcriptomic prognostic signatures is well documented. Michiels et al. [24] re-analyzed seven microarray studies and found prognostic classifiers unstable under resampling; Venet et al. [25] showed that in ER-positive breast cancer almost any gene set—including biologically irrelevant ones—is significantly associated with outcome, because a dominant proliferation axis correlates with most genes. These results motivate two distinctions our framework makes explicit: *reproducibility is not validity* (a reproducible association may still be confounded, as with EGFR and *IDH*), and *non-replication has multiple causes*—limited power, population or biological differences, and measurement artifact. Prior cross-platform quality efforts target dataset-level concordance; by contrast we ask a per-gene, per-comparison question—does each cohort measure *this* gene faithfully, judged against known biology—and show the answer predicts replication specifically where measurement fidelity varies.

---

## 3. Problem Statement

For TCGA-LGG, TCGA-GBM, and (for validation) CGGA, let $z_{\text{EGFR}}$ be per-SD standardized EGFR expression. We ask:

1. **Association.** Is $z_{\text{EGFR}}$ associated with overall survival within each cohort?
2. **Independence.** Does any association survive adjustment for age, grade, and *IDH*/subtype?
3. **Discrimination.** Does EGFR improve the C-index beyond established covariates?
4. **Normalization robustness.** Are conclusions stable under VST rather than TPM?
5. **External reproducibility.** Do the TCGA conclusions replicate in CGGA—and does CGGA measure EGFR faithfully enough to test them?

---

## 4. Methodology

### 4.1 Data acquisition

TCGA RNA-sequencing and clinical data were retrieved from the Genomic Data Commons [3] with `TCGAbiolinks` [4] (harmonized data, GDC Data Release 45.0, 2025-12-04), using STAR-Counts quantifications [8]. Molecular annotations came from the TCGA marker-paper fields [1, 2]. The CGGA mRNAseq_693 RSEM matrix and clinical table were downloaded from the CGGA portal [17].

### 4.2 Expression and standardization

TCGA expression used the TPM assay, $\log_2(\text{TPM}+1)$; CGGA used $\log_2(\text{RSEM}+1)$. Within each cohort, expression was standardized to unit SD ($z$-scores). Genes were matched by HGNC symbol; in CGGA, EGFR was verified to be a unique feature.

### 4.3 Cohorts and endpoint

We retained primary solid tumors, one sample per patient, with positive overall-survival (OS) time (time to death as event, else censored at last follow-up). For CGGA, WHO grade II/III was the lower-grade analog and grade IV the glioblastoma analog.

### 4.4 Statistical analysis

Cox regression modeled $h(t\mid\mathbf{x}) = h_0(t)\exp(\boldsymbol{\beta}^{\top}\mathbf{x})$. For EGFR we fit nested models adding covariates stepwise: unadjusted; $+$age; $+$grade; $+$*IDH* status; $+$full IDH/1p19q subtype (TCGA). We report the EGFR HR per SD ($\text{HR}>1$ = worse survival), Harrell's C-index [16] with 95% CIs, and the PH test for the EGFR term [13]. To ensure that any attenuation across the nested models reflects the added covariate rather than the change in sample induced by complete-case fitting, the nested models were also refit on a single common *IDH*-complete subset. EGFR's incremental value beyond established covariates was assessed with a likelihood-ratio (LR) test and by the change in C-index. Because the primary nested tests span a small prespecified family across cohorts, their $p$-values are nominal. Kaplan–Meier curves (median split) are descriptive [9, 10]. We treat *IDH* status as a confounder of the EGFR–survival association; because *IDH* mutation is an upstream driver and EGFR expression is partly downstream of the IDH-wildtype state, the adjustment could also be removing a mediated path—either way, EGFR carrying no signal after conditioning on *IDH* means it adds no information beyond subtype.

### 4.5 Pathway screen

Ten EGFR/RTK–PI3K-axis genes were each entered per SD into age-adjusted and (in TCGA-LGG) age-plus-*IDH*-adjusted Cox models, with Benjamini–Hochberg FDR correction [12].

### 4.6 Normalization sensitivity and external validation

We recomputed the key TCGA EGFR models on DESeq2 VST-transformed STAR counts [18]. For CGGA we repeated the nested EGFR models (adding radiotherapy and chemotherapy status), and—before interpreting—assessed positive controls. The criterion was **specified a priori from established biology**: because EGFR expression is higher in higher-grade and IDH-wildtype glioma, a cohort that measures EGFR faithfully must reproduce a positive EGFR expression–grade correlation, as TCGA does. Clinical controls (IDH and grade must predict survival) verify the survival and annotation data independently of expression. The between-cohort difference in the EGFR–grade correlation was tested formally (Fisher $r$-to-$z$). To distinguish a genuine CGGA property from an analyst processing error, we (i) verified the expression–clinical sample join, (ii) reproduced the result by an independent computation, and (iii) repeated the positive controls in a second, independent CGGA batch (mRNAseq_325). Finally, to test robustness to the processing pipeline itself, we re-obtained TCGA RNA-seq from recount3 [19]—an independent reprocessing of the raw reads by the Monorail pipeline (hg38/Gencode)—and repeated the nested EGFR models. The CGGA microarray dataset (mRNA-array_301) was used as an additional cross-platform positive-control check.

### 4.6a Within-subtype and copy-number analysis

To test whether EGFR is prognostic within IDH-wildtype glioma, we restricted TCGA (LGG+GBM) to IDH-wildtype tumors and fit Cox models for EGFR expression (per SD) and, separately, for gene-level EGFR copy number obtained from the GDC (ASCAT3 gene-level calls). High-level amplification was pre-specified as copy number $\geq 6$, with sensitivity at CN $\geq 5$ and $\geq 7$ and a continuous (per-SD) analysis. Grade throughout is the original histologic grade recorded by TCGA; because under the WHO 2021 scheme EGFR amplification is itself a grade-4 criterion, grade-adjusted amplification analyses carry a partial circularity that the pre-2021 grading mitigates. For the amplification analysis we report a power/precision analysis (detectable HR at 80% power via the Schoenfeld relation) alongside the point estimates.

### 4.6b Systematic positive-control-gated screen

We generalized the positive control to the 500 most-variable genes shared by TCGA and CGGA. For each gene we recorded the age-adjusted Cox coefficient in each cohort and the expression–grade Spearman correlation in each cohort. A gene was "QC-concordant" if its grade correlation had the same sign and (pre-specified) magnitude $>0.1$ in both cohorts. Among genes prognostic in TCGA (Benjamini–Hochberg FDR $<0.05$), we computed cross-cohort replication (same sign, $p<0.05$ in CGGA) and sign-flip rates, stratified by QC status, with binomial 95% CIs and a Fisher exact test on the QC $\times$ replication table. We repeated the analysis in a second, independent cohort pair (TCGA$\to$CGGA-325), varied the QC magnitude threshold ($|r|>0.05$, $0.1$, $0.2$) as a sensitivity check, and compared against a no-gating baseline. We use the neutral term *grade-discordant* for QC-failing genes: this flags a candidate measurement artifact but does not by itself prove mis-measurement, as a gene's grade relationship could differ biologically across populations.

To evaluate the framework at genome scale, we extended it to all adequately expressed genes shared by TCGA and the CGGA cohorts (6,512 genes prognostic in TCGA at FDR $<0.05$). For each gene we scored positive-control concordance as the cross-cohort product of expression–grade Spearman correlations (grade anchor) and, as an alternative, of expression–*IDH* correlations (IDH anchor), and measured how well each score predicts replication (same-sign, $p<0.05$ in the replication cohort) by the area under the ROC curve (Mann–Whitney estimator). We benchmarked these biological anchors against a detectability baseline (minimum cross-cohort mean expression) and evaluated all three against three independent replication cohorts: CGGA mRNAseq_693, mRNAseq_325, and the mRNA-array_301 microarray. To probe scope in a second tumor type, we repeated the framework in breast cancer using TCGA-BRCA (discovery) and METABRIC (replication); standardized expression and ER/subtype labels came from a common processed resource [20], overall survival from the respective cohort clinical files [21], and ER status and basal/non-basal subtype served as the biological anchors (in place of grade).

### 4.7 Reproducibility

The pipeline is deterministic (random seeds not applicable). Software versions are in §6. Data provenance is pinned: TCGA via the GDC **Data Release 45.0 (2025-12-04)**, and CGGA datasets mRNAseq_693, mRNAseq_325, and mRNA-array_301 (all release 20200506). CGGA data were used under the CGGA data-use terms, which permit academic reuse with citation [17]. All positive-control and amplification thresholds ($|r|>0.1$; CN $\geq 6$) were pre-specified from biology, not tuned to outcomes. Analysis scripts (steps 00–15) accompany the manuscript and will be deposited in a public repository on publication https://github.com/edanmn/TCGA-EGFR-Glioma.

**Multiplicity.** The study reports one prespecified primary hypothesis (EGFR $\to$ overall survival, adjusted) evaluated across cohorts and robustness checks; these primary $p$-values are nominal and should be read as a small prespecified family. The pathway and systematic screens are explicitly exploratory and FDR-controlled within each analysis; the single nominally significant secondary observation (the inverse IDH-wt-GBM expression trend, $p=0.03$) is not corrected and is not interpreted as confirmatory.

---

## 5. Analysis Pipeline

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
\end{tikzpicture}}
\caption{Analysis pipeline. Standardized expression enters nested Cox models with progressive confounder adjustment, followed by VST-normalization sensitivity and a positive-control–gated CGGA external-validation attempt.}
\end{figure*}

Figure 1 summarizes the workflow; acquisition is cached separately from analysis.

---

## 6. Experimental Setup

### 6.1 Cohorts

After quality control, TCGA-LGG comprised 511 patients (125 deaths) and TCGA-GBM 282 (227 deaths); molecular composition differs sharply (Table 1). The CGGA validation used 271 primary lower-grade (WHO II/III) patients (99 deaths) and 133 primary glioblastoma (WHO IV) patients (110 deaths).

\begin{table*}[t]
\centering
\caption{TCGA cohort characteristics and covariate availability after quality control. Sex was not populated in the expression-object annotation.}
\begin{tabular}{lll}
\toprule
Characteristic & TCGA-LGG & TCGA-GBM \\
\midrule
Histologic grade & II--III & IV \\
Patients (KM / unadjusted Cox) & 511 & 282 \\
Deaths (events) & 125 & 227 \\
Median age, years (IQR) & 41 (32--53) & 60 (51--69) \\
\emph{IDH} status (Mutant / WT) & 408 / 93 & 21 / 243 \\
Molecular subtype (codel / non-codel / wt) & 166 / 242 / 93 & --- / 19 / 237 \\
\emph{IDH} status available & 508 & 265 \\
\bottomrule
\end{tabular}
\end{table*}

### 6.2 Software

R 4.6.0 [15] with `TCGAbiolinks` 2.40.0 [4], `SummarizedExperiment` 1.42.0, `DESeq2` [18], `survival` 3.8-6 [13], and `survminer` 0.5.2 [14].

---

## 7. Results

### 7.1 Unadjusted association in TCGA (descriptive)

Under a median split, high EGFR expression was associated with shorter survival in TCGA-LGG (log-rank $p=0.040$) and pooled ($p=4.85\times10^{-3}$), but not TCGA-GBM ($p=0.99$) (Table 2; Figures 2–4).

\begin{table*}[t]
\centering
\caption{Unadjusted Kaplan--Meier median-split results for EGFR in TCGA (descriptive).}
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

\begin{center}\includegraphics[width=\linewidth]{figures/KM_EGFR_TCGALGG.png}\end{center}

**Figure 2.** Unadjusted KM overall survival by EGFR group (median split), TCGA-LGG.

\begin{center}\includegraphics[width=\linewidth]{figures/KM_EGFR_TCGAGBM.png}\end{center}

**Figure 3.** Unadjusted KM overall survival by EGFR group, TCGA-GBM.

\begin{center}\includegraphics[width=\linewidth]{figures/KM_EGFR_LGGGBM.png}\end{center}

**Figure 4.** Unadjusted KM overall survival by EGFR group, pooled LGG+GBM.

### 7.2 In TCGA, EGFR is not independent of IDH status

The nested Cox models (Table 3) show that in TCGA-LGG the EGFR association—strong when unadjusted (HR $=1.59$ per SD) and after age and grade (HR $=1.31$; $p=8.2\times10^{-3}$)—was **abolished after adjustment for *IDH* status** (HR $=1.13$; $p=0.16$) and molecular subtype (HR $=1.14$; $p=0.10$). EGFR alone discriminated poorly (C-index $0.59$, 95% CI $0.51$–$0.67$) versus $0.83$ ($0.79$–$0.87$) for age, grade, and *IDH*.

To confirm that the attenuation reflects *IDH* rather than the sample change from complete-case fitting, we refit the models on the identical *IDH*-complete subset ($n=452$, 105 events): EGFR remained significant after age and grade (HR $=1.31$; $1.07$–$1.60$) and lost significance only when *IDH* was added (HR $=1.13$; $0.95$–$1.34$; $p=0.16$), so the attenuation is attributable to *IDH*. A likelihood-ratio test confirmed that EGFR adds nothing beyond age, grade, and *IDH* ($\chi^2=1.98$, df $=1$, $p=0.16$), and it did not improve discrimination: C-index $0.834$ ($0.79$–$0.87$) without EGFR versus $0.830$ ($0.79$–$0.87$) with it.

In TCGA-GBM, EGFR was not associated when unadjusted; a weak inverse association in the adjusted model (HR $\approx0.87$; $p=0.03$) was of borderline significance, in the opposite direction, and not robust across model choices, so we do not interpret it as a real effect. The pooled association (HR $=1.35$; $p=1.98\times10^{-7}$) vanished after age adjustment, demonstrating confounding. PH violations for the EGFR term appeared in the IDH-/subtype-adjusted LGG models, where EGFR is in any case non-significant.

\begin{table*}[t]
\centering
\caption{Nested Cox models for EGFR expression in TCGA (HR per 1 SD), with Harrell's C-index and the PH test $p$-value for the EGFR term. Grade is omitted for glioblastoma.}
\begin{tabular}{lllllll}
\toprule
Cohort & Model & HR (95\% CI) & $p$ & C-index & PH $p$ (EGFR) & $n$ / events \\
\midrule
TCGA-LGG & Unadjusted            & 1.59 (1.30--1.94) & $6.85\times10^{-6}$ & 0.592 & 0.98  & 511 / 125 \\
TCGA-LGG & $+$ age               & 1.37 (1.14--1.66) & $1.06\times10^{-3}$ & 0.744 & 0.88  & 511 / 125 \\
TCGA-LGG & $+$ age, grade        & 1.31 (1.07--1.59) & $8.22\times10^{-3}$ & 0.786 & 0.23  & 454 / 106 \\
TCGA-LGG & $+$ age, grade, IDH   & 1.13 (0.95--1.33) & 0.164               & 0.830 & 0.004 & 452 / 105 \\
TCGA-LGG & $+$ age, subtype      & 1.14 (0.98--1.34) & 0.096               & 0.820 & 0.012 & 508 / 123 \\
\midrule
TCGA-GBM & Unadjusted            & 0.962 (0.845--1.09) & 0.552             & 0.523 & 0.19  & 282 / 227 \\
TCGA-GBM & $+$ age, IDH          & 0.866 (0.761--0.986) & 0.030            & 0.639 & 0.40  & 265 / 214 \\
\midrule
Pooled   & Unadjusted            & 1.35 (1.20--1.51) & $1.98\times10^{-7}$ & 0.553 & 0.03  & 793 / 352 \\
Pooled   & $+$ age, grade, IDH   & 0.921 (0.840--1.01) & 0.078             & 0.837 & 0.56  & 717 / 319 \\
\bottomrule
\end{tabular}
\end{table*}

### 7.3 The confounder and pathway screen

In TCGA-LGG, survival separates strongly by molecular subtype (Figure 5), and EGFR expression is elevated in the IDH-wildtype group—so its unadjusted association reflects subtype. In the standardized age-adjusted pathway screen (Table 4), six genes were significant; after adding *IDH*, EGFR, PIK3R1, and PTEN lost significance while TP53, RB1, IDH1, CDKN2A, and PIK3CA were significant, consistent with subtype-driven signal. No gene was significant in TCGA-GBM. IDH1 *expression* is distinct from the *IDH* mutation.

\begin{center}\includegraphics[width=\linewidth]{figures/KM_LGG_subtype.png}\end{center}

**Figure 5.** TCGA-LGG overall survival by IDH/1p19q molecular subtype; subtype, not EGFR, drives the separation in Figure 2.

\begin{table*}[t]
\centering
\caption{FDR-controlled pathway screen in TCGA-LGG (HR per 1 SD, age-adjusted), with FDR-adjusted $p$ before and after adding \emph{IDH} status.}
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

### 7.4 The TCGA result is robust to normalization

Recomputing on DESeq2 VST-transformed counts reproduced the conclusion (Table 5): the lower-grade EGFR HR fell from $1.50$ (unadjusted) to $1.11$ ($p=0.21$) after age, grade, and *IDH*. The confounding pattern is not a TPM artifact.

\begin{table*}[t]
\centering
\caption{Normalization sensitivity: EGFR HR per 1 SD under DESeq2 VST (compare to the TPM results in Table 3).}
\begin{tabular}{lllll}
\toprule
Cohort & Model & HR (95\% CI) & $p$ & C-index \\
\midrule
TCGA-LGG & Unadjusted          & 1.50 (1.22--1.84) & $1.25\times10^{-4}$ & 0.577 \\
TCGA-LGG & $+$ age, grade       & 1.26 (1.04--1.54) & 0.020               & 0.784 \\
TCGA-LGG & $+$ age, grade, IDH  & 1.11 (0.944--1.31) & 0.205              & 0.829 \\
TCGA-GBM & $+$ age, IDH         & 0.843 (0.737--0.963) & 0.012            & 0.641 \\
\bottomrule
\end{tabular}
\end{table*}

### 7.5 External validation is inconclusive: CGGA fails EGFR positive controls

A naive CGGA analysis appeared to show the opposite of TCGA: EGFR was null unadjusted (HR $=1.17$; $p=0.16$) but *significant after* adjusting for age, grade, and *IDH* (HR $=1.44$; $p=4.6\times10^{-4}$), and after treatment (HR $=1.37$; $p=6.8\times10^{-3}$) (Table S1). Before interpreting this as a genuine discordance, we tested positive controls.

The CGGA **clinical** controls are strong and correct: IDH-mutant status is highly protective (HR $=0.23$; $p=2\times10^{-22}$) and grade predicts survival (grade IV vs.\ II HR $=8.5$), confirming that the survival, grade, and *IDH* annotations and the sample–clinical join are sound. However, the CGGA **EGFR-expression** controls fail. EGFR should rise with grade and be higher in IDH-wildtype tumors, as it does in TCGA; in CGGA it does neither (Table 6; Figure 6). EGFR expression is essentially flat across WHO grade (Spearman $r=-0.03$, vs.\ $+0.18$ in TCGA) and is *higher* in IDH-mutant tumors (mean $\log_2$ $4.17$ vs.\ $3.63$; the reverse of TCGA's $6.43$ vs.\ $7.44$). Moreover, expression–grade correlations are systematically attenuated across the whole gene panel in CGGA relative to TCGA (Table 6), indicating a general weakening of expression–phenotype signal in this matrix rather than an EGFR-specific mislabeling.

Two checks establish that this is a genuine property of CGGA's EGFR quantification rather than a processing error. First, all 693 expression samples matched clinical records (693/693), and the flat EGFR–grade relationship was reproduced by an independent computation. Second, the failure **reproduced across two further, independent CGGA datasets spanning a different platform**: a second RNA-seq batch (mRNAseq_325, $n=229$; EGFR flat across grade, means $5.00$/$5.22$/$5.11$) and a microarray dataset (mRNA-array_301, $n=264$). In all three datasets the EGFR–grade correlation differed significantly from TCGA (Fisher $r$-to-$z$ $p<10^{-3}$; Table 7), while clinical controls passed in every one (IDH-mutant HR $0.22$–$0.32$; high-grade HR $7$–$9$). The consistency across two RNA-seq batches and a microarray indicates the problem is not a single pipeline's quantification but a property of EGFR measurement in this resource. The positive-control criterion was specified a priori from biology, not chosen after seeing the disagreement; we nonetheless acknowledge the confirmation-bias risk and mitigate it through the a-priori criterion, the formal between-cohort test, and the independent-batch replication.

Because CGGA EGFR expression does not track the biology EGFR is known to follow—reproducibly, across two batches, with verified data handling—the apparent "IDH-independent EGFR effect" cannot be interpreted as genuine, and CGGA does not provide a valid external test of the TCGA finding for this gene. A naive reading—reporting an IDH-independent, treatment-robust EGFR effect—would have been a spurious positive that only the positive controls exposed. The mechanistic cause of CGGA's attenuated EGFR signal (normalization, gene model, or platform) remains unidentified and would require reprocessing from raw reads to resolve.

\begin{table*}[t]
\centering
\caption{Positive controls: Spearman correlation of gene expression with WHO grade in TCGA vs.\ CGGA. EGFR (and the panel generally) tracks grade in TCGA but not CGGA, indicating attenuated expression–phenotype signal in the CGGA matrix.}
\begin{tabular}{lcc}
\toprule
Gene & $r$(expression, grade), TCGA & $r$(expression, grade), CGGA \\
\midrule
EGFR   & $+0.180$ & $-0.033$ \\
IDH1   & $+0.364$ & $+0.152$ \\
TP53   & $+0.075$ & $+0.113$ \\
CDKN2A & $-0.158$ & $-0.076$ \\
NF1    & $-0.303$ & $-0.183$ \\
PTEN   & $-0.363$ & $-0.071$ \\
PDGFRA & $-0.319$ & $-0.096$ \\
PIK3R1 & $-0.504$ & $-0.195$ \\
PIK3CA & $-0.089$ & $-0.057$ \\
RB1    & $-0.046$ & $-0.023$ \\
\bottomrule
\end{tabular}
\end{table*}

\begin{table*}[t]
\centering
\caption{EGFR positive control across TCGA and three CGGA datasets spanning two platforms. The EGFR--grade Spearman correlation is positive in TCGA but null/negative in every CGGA dataset (Fisher $r$-to-$z$ vs.\ TCGA), whereas clinical controls (IDH, grade $\to$ survival) are strong and correct throughout. The expression-level failure is reproducible across batches and platforms; the clinical data are sound.}
\begin{tabular}{lccccc}
\toprule
Cohort & $r$(EGFR, grade) & $n$ & $p$ vs.\ TCGA & IDH-mutant HR & High-grade HR (vs.\ II) \\
\midrule
TCGA (reference)        & $+0.185$ & 826 & --- & --- & --- \\
CGGA-693 (RNA-seq)      & $-0.033$ & 422 & $2.4\times10^{-4}$ & 0.23 & 8.5 \\
CGGA-325 (RNA-seq)      & $-0.097$ & 229 & $1.6\times10^{-4}$ & 0.22 & 8.6 \\
CGGA-301 (microarray)   & $-0.086$ & 264 & $1.2\times10^{-4}$ & 0.32 & 7.1 \\
\bottomrule
\end{tabular}
\end{table*}

\begin{center}\includegraphics[width=\linewidth]{figures/positive_control_EGFR_grade.png}\end{center}

**Figure 6.** Positive control: standardized EGFR expression by WHO grade. EGFR rises with grade in TCGA (canonical) but is flat in CGGA, indicating that CGGA EGFR expression does not reflect expected EGFR biology.

### 7.6 The TCGA result is robust to the processing pipeline

To test whether the finding depends on the GDC pipeline, we re-obtained TCGA from recount3 [19], an independent Monorail reprocessing of the raw reads. Under this pipeline EGFR again passed the positive controls—expression rose with grade (lower-grade means $7.89$ at II, $8.41$ at III) and was higher in IDH-wildtype tumors ($8.92$ vs.\ $7.99$ in lower-grade; $9.20$ vs.\ $7.63$ in glioblastoma)—and the confounding result replicated: the lower-grade EGFR HR fell from $1.54$ (unadjusted) to $1.27$ (age$+$grade) to $1.11$ ($0.94$–$1.31$; $p=0.22$) after *IDH* (Table 8), matching the GDC-pipeline (Table 3) and VST (Table 5) estimates. The result is thus not a processing artifact, and a faithfully reprocessed cohort recovers the expected EGFR biology that CGGA lacks. Because CGGA raw reads are controlled-access, an equivalent raw-read reprocessing of CGGA—the only route to new-patient external validation here—was not possible and remains the key open step.

\begin{table*}[t]
\centering
\caption{Pipeline robustness: nested EGFR Cox models (HR per 1 SD) on TCGA re-obtained from recount3 (independent Monorail raw-read reprocessing). Estimates match the GDC-pipeline results (Table 3).}
\begin{tabular}{llllll}
\toprule
Cohort & Model & HR (95\% CI) & $p$ & $n$ / events \\
\midrule
TCGA-LGG & Unadjusted          & 1.54 (1.25--1.89) & $4.02\times10^{-5}$ & 509 / 125 \\
TCGA-LGG & $+$ age, grade       & 1.27 (1.04--1.55) & 0.021               & 452 / 106 \\
TCGA-LGG & $+$ age, grade, IDH  & 1.11 (0.94--1.31) & 0.220               & 450 / 105 \\
TCGA-GBM & Unadjusted           & 0.906 (0.755--1.09) & 0.289             & 154 / 122 \\
TCGA-GBM & $+$ age, grade, IDH  & 0.813 (0.677--0.977) & 0.027            & 150 / 120 \\
\bottomrule
\end{tabular}
\end{table*}

### 7.7 EGFR is not prognostic *within* IDH-wildtype glioma

The analyses above show EGFR expression is confounded by *IDH* status. A stronger test is whether EGFR carries any prognostic signal *within* the aggressive IDH-wildtype subtype, where subtype confounding no longer applies. Among IDH-wildtype gliomas (TCGA LGG+GBM; $n=338$, 253 deaths), EGFR expression was not prognostic (Table 9): unadjusted HR $=1.00$ ($0.89$–$1.13$; a well-powered null), age-adjusted $0.93$ ($p=0.22$), and $+$grade $0.90$ ($p=0.10$), with a weak inverse—not adverse—trend in IDH-wildtype glioblastoma (HR $=0.87$, $p=0.03$, not corrected for multiplicity). The proportional-hazards assumption held for the EGFR term ($p=0.08$). Grade here is the *original histologic* grade recorded by TCGA (pre-2021), not the molecular WHO 2021 grade; because under WHO 2021 EGFR amplification is itself a grade-4 criterion, grade-adjusted amplification analyses carry a partial circularity, mitigated here by the pre-2021 grading.

We obtained gene-level copy number for the same tumors. High-level EGFR amplification was strongly enriched in IDH-wildtype disease ($55\%$ vs.\ $2\%$ of IDH-mutant)—confirming the copy-number data are faithful, a positive control the CGGA expression data failed. Within IDH-wildtype glioma we found *no evidence* that amplification predicts survival: the continuous copy-number effect was null and reasonably powered (HR $=1.00$ per SD, $0.84$–$1.20$, age-adjusted), and binary high-level amplification was null across thresholds (CN $\geq 5$/$6$/$7$: HR $1.09$/$1.00$/$0.92$, all $p>0.6$; Table 9). The binary tests are, however, underpowered—with 127 events and $55\%$ amplified, 80% power extends only to HR $\geq 1.65$—so we can exclude a large but not a modest amplification effect. Taken together, EGFR carries no detectable within-subtype prognostic signal by expression, and no evidence of one by copy number: its apparent prognostic value is a between-subtype marker of the IDH-wildtype state.

\begin{table*}[t]
\centering
\caption{EGFR is not a within-subtype prognostic factor in IDH-wildtype glioma (TCGA). Expression (per 1 SD) shows a well-powered null; copy-number amplification shows no evidence of an effect (continuous, reasonably powered) and null across binary thresholds, though binary tests are underpowered (80\% power only for HR $\geq 1.65$). Grade is the pre-2021 histologic grade.}
\begin{tabular}{lllll}
\toprule
Marker & Model & HR (95\% CI) & $p$ & $n$ / events \\
\midrule
Expression   & Unadjusted                 & 1.00 (0.89--1.13) & 0.99 & 338 / 253 \\
Expression   & $+$ age                    & 0.93 (0.82--1.05) & 0.22 & 338 / 253 \\
Expression   & $+$ age, grade             & 0.90 (0.80--1.02) & 0.10 & 330 / 246 \\
Expression   & IDH-wt GBM, $+$ age        & 0.87 (0.76--0.99) & 0.034 & 244 / 203 \\
Amplification& Copy number (per SD), $+$age & 1.00 (0.84--1.20) & 0.98 & 166 / 127 \\
Amplification& CN $\geq 5$ vs.\ not, $+$age & 1.09 (0.76--1.57) & 0.63 & 166 / 127 \\
Amplification& CN $\geq 6$ vs.\ not, $+$age & 1.00 (0.70--1.42) & 0.98 & 166 / 127 \\
Amplification& CN $\geq 7$ vs.\ not, $+$age & 0.92 (0.65--1.32) & 0.66 & 166 / 127 \\
\bottomrule
\end{tabular}
\end{table*}

### 7.8 A positive-control-gated framework for cross-cohort validation

The CGGA episode motivates a general safeguard: before trusting a cross-cohort prognostic comparison for a gene, verify that both cohorts measure it faithfully via a biological positive control. We formalized this on the 500 most-variable genes shared by TCGA and CGGA. For each gene we tested the prognostic association (age-adjusted Cox) in TCGA and its replication in CGGA, and computed the expression–grade correlation in each cohort as a *pre-specified* positive control; a gene "passes" QC when its grade correlation has a consistent sign and magnitude $>0.1$ in both cohorts. Of the genes prognostic in TCGA (FDR $<0.05$), replication in CGGA-693 was $92.8\%$ (95% CI $89.9$–$95.0$) among QC-passing genes versus $24.3\%$ ($11.8$–$41.2$) among QC-failing genes (Fisher exact $p=5\times10^{-21}$, OR $=39$), with every sign-flip confined to QC-failing genes (Table 10; Figure 7). The contrast reproduced in an independent second cohort pair (TCGA$\to$CGGA-325: $95.4\%$ vs.\ $19.2\%$; Fisher $p=1\times10^{-20}$), was stable across QC thresholds ($|r|>0.05$, $0.1$, $0.2$), and exceeded the no-gating baseline ($87.3\%$ overall replication). Overall $88\%$ of the panel passed QC. We describe QC-failing genes as *grade-discordant*—a candidate signature of measurement artifact rather than proven mis-measurement, since a gene could in principle relate to grade differently across populations; for EGFR the artifact interpretation is independently supported by its canonical grade/IDH biology (§7.3, §7.5). Positive-control gating thus concentrates non-replicating and sign-flipping associations into a small, flaggable minority.

\begin{table*}[t]
\centering
\caption{Positive-control-gated cross-cohort replication for genes prognostic in TCGA (FDR $<0.05$), among 500 top-variable genes, in two independent cohort pairs. Replication is far higher for QC-passing (grade-concordant) than QC-failing (grade-discordant) genes; Fisher exact $p<10^{-20}$ for both pairs, and all sign-flip discordances fall in QC-failing genes.}
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

\begin{center}\includegraphics[width=0.86\linewidth]{figures/systematic_method.png}\end{center}

**Figure 7.** Cross-cohort expression QC for TCGA-prognostic genes. Each point is a gene; axes are its expression–grade correlation in TCGA (x) and CGGA (y). QC-concordant genes (blue) lie along the diagonal (measured consistently); QC-discordant genes (red), including EGFR, fall off it and are where cross-cohort prognostic replication breaks down.

### 7.9 The positive control predicts replication genome-wide

To test whether the positive control is a general safeguard rather than an EGFR-specific observation, we scaled it to **6,512 genes prognostic in TCGA** (FDR $<0.05$, age-adjusted) and asked, for each, whether its association replicated in three independent CGGA cohorts. For every gene we scored positive-control concordance as the cross-cohort product of its expression–grade correlations (and, as an alternative anchor, its expression–*IDH* correlations), and measured how well that score predicts replication by the area under the ROC curve (Table 11; Figure 8). The grade anchor predicted replication with **AUC $0.82$ (CGGA-693), $0.93$ (CGGA-325), and $0.90$ (microarray-301)**; the *IDH* anchor gave AUC $0.77$–$0.84$; and a naive detectability baseline (mean expression) was at or below chance (AUC $0.36$–$0.56$). A simple grade-concordance flag identified non-replicating associations with 70% precision and 60% recall in CGGA-693. Two conclusions follow. First, the safeguard generalizes—across three cohorts, two RNA-seq batches plus a microarray, and two independent biological anchors. Second, it works because it captures *faithful measurement*, not mere detectability: highly expressed genes are not more replicable (detectability is anti-predictive). EGFR falls in the flagged, grade-discordant minority, confirming that its cross-cohort behavior (§7.5) is a measurement artifact the framework catches, not a biological discordance.

\begin{table*}[t]
\centering
\caption{Genome-scale evaluation: a gene's biological positive-control concordance predicts whether its TCGA prognostic association replicates in an independent cohort (AUC). Grade and IDH anchors are strongly predictive across three cohorts and two platforms; a detectability baseline is at/below chance. $n=6,512$ TCGA-prognostic genes (FDR $<0.05$). Bootstrap 95\% CIs for the grade-anchor AUCs are tight given $n\approx6,500$ (CGGA-693: $0.82$, $0.81$–$0.83$); Table 12 gives the circularity control.}
\begin{tabular}{lcccc}
\toprule
Replication cohort & Replication rate & AUC (grade anchor) & AUC (IDH anchor) & AUC (detectability) \\
\midrule
CGGA-693 (RNA-seq)    & 0.62 & 0.82 & 0.77 & 0.36 \\
CGGA-325 (RNA-seq)    & 0.70 & 0.93 & 0.84 & 0.37 \\
CGGA-301 (microarray) & 0.59 & 0.90 & 0.83 & 0.56 \\
\bottomrule
\end{tabular}
\end{table*}

\begin{center}\includegraphics[width=0.9\linewidth]{figures/method_auc.png}\end{center}

**Figure 8.** Biological positive controls predict cross-cohort replication of prognostic associations. AUC for predicting replication of TCGA-prognostic genes in each CGGA cohort, by positive-control metric; grade and IDH anchors far exceed the detectability baseline (dashed line = chance).

**Controlling for circularity.** Because survival in pooled glioma is grade-driven, the prognostic set is enriched for grade-correlated genes, and the grade-anchor QC could predict replication partly *because* both prognostic status and replication share the grade signal it measures. We tested this directly (Table 12). Restricting to associations that remain prognostic *after adjustment for grade and IDH*—i.e., subtype-independent signal—the QC still predicted replication above chance (AUC $=0.67$, 95% CI $0.61$–$0.73$), as it did for genes with low grade-correlation ($|r|<0.2$; AUC $=0.76$, $0.73$–$0.79$). The AUC is higher for strongly grade-correlated genes ($0.85$) than for low-grade-correlation genes ($0.76$), so part of the headline value does reflect shared grade signal; but the effect clearly persists for subtype-independent associations, indicating the positive control captures measurement fidelity *beyond* the anchor covariate itself rather than acting tautologically.

\begin{table*}[t]
\centering
\caption{Circularity control (TCGA $\to$ CGGA-693, grade-anchor QC; bootstrap 95\% CIs). The QC's predictive value is partly inflated by grade-correlated genes but persists for subtype-independent associations (grade+IDH-adjusted) and low-grade-correlation genes, excluding a purely tautological explanation.}
\begin{tabular}{lccc}
\toprule
Prognostic-gene set & Genes & Replication rate & AUC (95\% CI) \\
\midrule
Unadjusted (age only; headline)         & 6{,}616 & 0.63 & 0.82 (0.81--0.83) \\
Subtype-adjusted (grade + IDH)          & 412  & 0.27 & 0.67 (0.61--0.73) \\
Low grade-correlation ($|r|<0.2$)       & 996  & 0.51 & 0.76 (0.73--0.79) \\
High grade-correlation ($|r|\geq0.2$)   & 5{,}620 & 0.65 & 0.85 (0.84--0.86) \\
\bottomrule
\end{tabular}
\end{table*}

### 7.10 Scope: a well-matched second tumor type

To test whether the framework is a *universal* replication predictor or a targeted detector of measurement-artifact discordance, we applied it to a second, unrelated tumor type: breast cancer, using TCGA-BRCA [23] (discovery, RNA-seq) and METABRIC [22] (replication, microarray), with ER status and basal/non-basal subtype as biological anchors and expression from a common processed resource [20, 21]. Here the positive control did **not** meaningfully predict replication: ER-anchor AUC $=0.60$ (95% CI $0.44$–$0.75$; includes chance) for TCGA-BRCA$\to$METABRIC and $0.56$ ($0.52$–$0.59$) for the reverse direction; the basal anchor was $\leq0.54$. Two non-exclusive reasons make breast diagnostic of *when the method helps*: the two cohorts measure genes highly concordantly (genome-wide agreement of expression–ER correlations $r=0.84$), leaving little measurement artifact to flag; and in ER-positive breast almost any gene set is prognostic [25], so cross-cohort non-replication is driven by weak, proliferation-dominated signal and limited power (TCGA-BRCA has 67 events) rather than by measurement fidelity. We cannot fully separate these causes here, but both predict that a measurement-fidelity gate should add little—exactly what we observe. This contrasts sharply with glioma, where CGGA systematically mis-measures a subset of genes (§7.5) and the gate is highly predictive (§7.9). The framework is therefore best understood not as a general replication oracle but as a **targeted detector of measurement-fidelity–driven discordance**: its predictive value scales with how much cross-cohort measurement heterogeneity exists, and it correctly adds little when both cohorts are clean.

---

## 8. Discussion

### 8.1 Interpretation

EGFR expression is a subtype artifact in TCGA: it is elevated in IDH-wildtype tumors, and once *IDH* or subtype is modeled it carries no independent prognostic information and adds no discrimination. This conclusion is robust to VST normalization and to an independent raw-read reprocessing (Figure 6, Tables 5, 9). Crucially, the null extends *within* the IDH-wildtype subtype: neither EGFR expression nor high-level amplification predicts survival there (Table 9). EGFR is therefore best understood as a between-subtype marker of the IDH-wildtype state, with no prognostic signal—by expression or copy number—inside it. This is a more complete characterization than "confounded by IDH": it rules out a residual within-subtype effect that grade-only adjustment could have missed.

### 8.2 The external-validation attempt and its lesson

Our attempt to validate in CGGA is instructive precisely because it initially appeared to contradict TCGA. Positive controls showed that CGGA EGFR expression does not track grade or the canonical IDH relationship, while CGGA clinical signals are strong—and this failure reproduced across three CGGA datasets spanning two platforms with a verified data join, so it is a genuine property of the CGGA expression data, not an analyst error. The apparent IDH-independent effect is therefore an artifact of the expression matrix, not a biological discordance. We guarded against the obvious objection—that we simply discarded the cohort that disagreed—by specifying the positive-control criterion from established biology in advance, testing the between-cohort difference formally, and replicating the failure across independent batches and platforms. The lesson generalizes, and we quantified it: across 6,512 prognostic genes, a gene's positive-control concordance predicts whether its association replicates with AUC $0.82$–$0.93$, across three cohorts, two platforms, and two biological anchors (grade and *IDH*), while a detectability baseline is uninformative (§7.9). This reframes cross-cohort validation as a two-part question—*does the association reproduce, and does the replication cohort measure the gene faithfully?*—with the second part checkable a priori from biology. Practically, gating on the positive control would have prevented the spurious EGFR "discordance" (and 30% of prognostic genes carry the same risk in CGGA-693). The safeguard is cheap, requires only a known biological covariate, and is agnostic to the specific gene, cohort, or platform. Its benefit is nonetheless *scoped* (§7.10): in a second tumor type where both cohorts measure genes concordantly (breast, TCGA-BRCA$\leftrightarrow$METABRIC), the gate had little artifact to flag and did not predict replication. This is the expected behavior of a measurement-fidelity check, and it clarifies the method's remit—flagging artifact-driven discordance—rather than promising to predict all non-replication (much of which is power- or biology-driven). Reporting the naive CGGA result as a genuine discordance—as an earlier version of this analysis did—would have been a false positive that the framework flags.

### 8.3 Strengths

The analysis is transparent and re-runnable, uses standardized comparable effect sizes, confronts the dominant confounder, reports discrimination and PH diagnostics, tests normalization robustness, and—rather than trusting an external cohort—verifies it with positive controls. Negative results are reported without selection.

---

## 9. Limitations

First, the primary finding rests on a single cohort (TCGA); because the CGGA expression matrix could not measure EGFR faithfully (below), no successful external validation of EGFR was obtained, and the conclusion is therefore internally robust and biologically consistent but externally unreplicated for this gene. Second, the CGGA EGFR positive-control failure was shown to be reproducible (three datasets, two platforms) and not a data-handling error (verified join, independent recomputation), but its mechanistic cause—normalization, gene model, or platform—was not identified and would require reprocessing CGGA from raw reads. Third, adjusted models used complete cases without imputation, though the key nested comparison was refit on a common sample. Fourth, the pooled TCGA analysis mixes two separately processed studies, and grade IV is collinear with the glioblastoma study, so batch and grade cannot be separated. Fifth, PH violations were present for some strongly prognostic covariates (the EGFR term itself satisfied PH in the within-subtype model, $p=0.08$). Sixth, the systematic method was demonstrated against a no-gating baseline and a detectability baseline but not against a broad panel of alternative QC metrics, and "grade-discordant" flags a candidate—not proven—measurement artifact. Seventh, and conceptually important, the framework assesses cross-cohort *reproducibility*, not *validity*: it can bless a reproducible association that is nonetheless confounded (as EGFR's is by *IDH*), so it complements, rather than replaces, causal/confounding scrutiny. Eighth, part of the genome-scale AUC reflects shared grade signal between the "prognostic" definition and the grade anchor; the effect persists but is weaker for subtype-independent associations (Table 12). IDH1 *expression* is distinct from the prognostic *IDH* mutation [7].

---

## 10. Future Work

An independent raw-read reprocessing of TCGA (recount3) confirmed the finding is pipeline-robust, but true new-patient external validation still requires an independent cohort whose EGFR expression passes the positive controls. Reprocessing CGGA from raw reads would be the natural test, but its raw reads are controlled-access; an openly available, faithfully quantified glioma RNA-seq cohort with *IDH*, grade, and survival is therefore needed. Integrating DNA-level EGFR amplification with mRNA expression, within molecular subtypes, would further test whether any reproducible EGFR signal exists beyond subtype membership. For the framework, the informative next step is to evaluate additional cohort pairs that *do* differ in measurement fidelity across several tumor types—the regime where the gate is designed to help—and to test better-powered discovery cohorts, since the breast pair (§7.10) was both well-matched and survival-underpowered.

---

## 11. Conclusion

Cross-cohort prognostic replication in transcriptomics is unreliable, and we show that a simple biological positive control detects the subset of discordances that stem from faulty measurement: across 6,512 genes, three glioma cohorts, two platforms, and two anchors, a gene's expression–covariate concordance forecasts whether its prognostic association replicates (AUC $0.82$–$0.93$), while naive detectability does not. Its benefit is appropriately scoped—in a well-matched breast pair (TCGA-BRCA$\leftrightarrow$METABRIC) where both cohorts measure genes concordantly, the gate had little to flag and did not predict replication—so it is best understood as a targeted detector of measurement-fidelity–driven discordance rather than a universal replication oracle. This gives a cheap safeguard: verify that both cohorts measure the gene faithfully before trusting (or distrusting) a cross-cohort comparison. The framework is motivated by EGFR, for which naive analysis manufactures a spurious discordance: in TCGA, EGFR mRNA expression is prognostic only until *IDH* is modeled (no added discrimination; robust to normalization and to an independent raw-read reprocessing), it carries no within-IDH-wildtype signal by expression or amplification, and its "IDH-independent" effect in CGGA is a positive-control failure the framework flags—so EGFR is a marker of the IDH-wildtype state, not an independent prognostic factor. The study thus delivers a reusable methodological safeguard for multi-cohort prognostic research, grounded in a definitive negative characterization of a widely studied marker.

---

## References

1. Brennan CW, Verhaak RGW, McKenna A, et al. The somatic genomic landscape of glioblastoma. *Cell*. 2013;155(2):462–477.
2. Cancer Genome Atlas Research Network. Comprehensive, integrative genomic analysis of diffuse lower-grade gliomas. *N Engl J Med*. 2015;372(26):2481–2498.
3. Grossman RL, Heath AP, Ferretti V, et al. Toward a shared vision for cancer genomic data. *N Engl J Med*. 2016;375(12):1109–1112.
4. Colaprico A, Silva TC, Olsen C, et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic Acids Res*. 2016;44(8):e71.
5. Louis DN, Perry A, Wesseling P, et al. The 2021 WHO Classification of Tumors of the Central Nervous System: a summary. *Neuro Oncol*. 2021;23(8):1231–1251.
6. Verhaak RGW, Hoadley KA, Purdom E, et al. Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and NF1. *Cancer Cell*. 2010;17(1):98–110.
7. Yan H, Parsons DW, Jin G, et al. IDH1 and IDH2 mutations in gliomas. *N Engl J Med*. 2009;360(8):765–773.
8. Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*. 2013;29(1):15–21.
9. Kaplan EL, Meier P. Nonparametric estimation from incomplete observations. *J Am Stat Assoc*. 1958;53(282):457–481.
10. Mantel N. Evaluation of survival data and two new rank order statistics arising in its consideration. *Cancer Chemother Rep*. 1966;50(3):163–170.
11. Cox DR. Regression models and life-tables. *J R Stat Soc Series B*. 1972;34(2):187–220.
12. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B*. 1995;57(1):289–300.
13. Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model*. New York: Springer; 2000. (R package `survival` v3.8-6.)
14. Kassambara A, Kosinski M, Biecek P. survminer: drawing survival curves using 'ggplot2'. R package v0.5.2; 2024.
15. R Core Team. *R: A Language and Environment for Statistical Computing*. Vienna, Austria: R Foundation for Statistical Computing; 2026.
16. Harrell FE Jr, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Stat Med*. 1996;15(4):361–387.
17. Zhao Z, Zhang KN, Wang Q, et al. Chinese Glioma Genome Atlas (CGGA): a comprehensive resource with functional genomic data from Chinese glioma patients. *Genomics Proteomics Bioinformatics*. 2021;19(1):1–12.
18. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*. 2014;15(12):550.
19. Wilks C, Zheng SC, Chen FY, et al. recount3: summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biol*. 2021;22(1):323.
20. Processed TCGA-BRCA and METABRIC datasets (Moanna). Zenodo; doi:10.5281/zenodo.4326602.
21. Cerami E, Gao J, Dogrusoz U, et al. The cBio Cancer Genomics Portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discov*. 2012;2(5):401–404.
22. Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature*. 2012;486(7403):346–352.
23. Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature*. 2012;490(7418):61–70.
24. Michiels S, Koscielny S, Hill C. Prediction of cancer outcome with microarrays: a multiple random validation strategy. *Lancet*. 2005;365(9458):488–492.
25. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol*. 2011;7(10):e1002240.

---

## Supplementary Material

\begin{table*}[t]
\centering
{\textbf{Table S1.} Naive CGGA external-validation models: EGFR HR per 1 SD in primary tumors. Retained for transparency only; these estimates are \emph{not} interpretable, because CGGA EGFR expression fails the positive control (Table 6, Table 7, Figure 6). The apparent IDH-adjusted and treatment-adjusted effects are artifacts of EGFR mis-measurement in CGGA, not evidence of a prognostic effect.}\\[4pt]
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
