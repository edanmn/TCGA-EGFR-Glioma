---
title: "Biological positive controls for cross-cohort validation of prognostic gene expression, with a case study of EGFR in glioma"
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
  - \usepackage{enumitem}
  - \setlist[enumerate]{itemsep=-1.5pt,topsep=2pt,parsep=0pt}
  - \makeatletter\def\fps@figure{H}\makeatother
---

## Abstract

Prognostic gene associations from one cohort frequently fail to replicate in another, and it is hard to know in advance which will. We propose a simple safeguard—a *biological positive control*—and show it predicts replicability at genome scale. The idea: a cross-cohort prognostic comparison is trustworthy only if both cohorts measure the gene faithfully, checkable by whether its expression tracks a known biological covariate (here, tumor grade or *IDH* status) consistently across cohorts. Using The Cancer Genome Atlas (TCGA) as discovery and three Chinese Glioma Genome Atlas (CGGA) datasets (two RNA-seq batches, one microarray) for replication, we scored 6,512 TCGA-prognostic genes. Positive-control concordance predicted replication with **AUC $0.82$ (CGGA-693), $0.93$ (CGGA-325), and $0.92$ (microarray-301)** under the grade anchor and $0.77$–$0.84$ under the *IDH* anchor, whereas detectability ($0.36$–$0.53$) and replication-cohort estimation precision ($0.44$–$0.59$) baselines were at or near chance. In a second tumor type (breast; TCGA-BRCA$\leftrightarrow$METABRIC), where both cohorts measure genes concordantly ($r=0.84$), the gate had little to flag: AUC $0.60$ ($0.44$–$0.75$) forward and $0.56$ ($0.52$–$0.59$) in the better-powered reverse direction—above chance but far below glioma's—**indicating it targets measurement-artifact discordance specifically rather than all non-replication**. The framework is motivated by EGFR, for which naive analysis manufactures a spurious discordance: in TCGA, EGFR expression is prognostic only until *IDH* is modeled (HR $1.31\to1.13$, $p=0.16$; no added discrimination; robust to normalization and to an independent raw-read reprocessing), and it carries no signal *within* IDH-wildtype glioma by expression (HR $1.00$, $0.89$–$1.13$) or amplification. Its apparent IDH-independent effect in CGGA is a positive-control failure—EGFR fails to track grade in all three CGGA datasets ($p<10^{-3}$) despite sound clinical data. EGFR is therefore a marker of the IDH-wildtype state rather than an independent prognostic factor, and biological positive controls provide a cheap gate against measurement-artifact-driven discordance.


## 1. Introduction

Diffuse gliomas vary widely in outcome, and their prognosis is defined primarily by molecular features—*IDH* mutation status and 1p/19q codeletion—rather than histology alone [5]. Against this backdrop, individual genes are frequently proposed as transcriptomic prognostic markers. EGFR is a natural candidate: a defining oncogene of the classical glioblastoma subtype and one of the most frequently altered genes in glioma [1, 6].

Any single-gene prognostic claim in glioma faces a specific threat: because molecular subtype governs both a tumor's transcriptional program and its outcome, a gene whose expression tracks subtype will appear prognostic even if it carries no information beyond subtype. A credible claim must therefore (i) survive adjustment for *IDH* status and molecular subtype, (ii) be robust to expression normalization, and (iii) replicate in an independent cohort—with the replication cohort itself verified to measure the gene faithfully.

Our primary contribution is **methodological**: a biological positive-control framework for cross-cohort prognostic validation, evaluated at genome scale (6,512 genes) across three replication cohorts and two platforms, in which a gene's expression–covariate concordance predicts whether its prognostic association replicates (AUC $0.82$–$0.93$), outperforming detectability and precision baselines and generalizing across two anchors (grade, *IDH*). We delineate its **scope** in a second tumor type (breast, TCGA-BRCA$\leftrightarrow$METABRIC), where concordant measurement leaves the gate little to flag. The framework is motivated and stress-tested by a **case study of EGFR**, for which we additionally provide (i) nested Cox models with discrimination and PH diagnostics showing the association vanishes after *IDH* adjustment, robust to normalization and to an independent raw-read reprocessing; (ii) a within-IDH-wildtype analysis of expression and copy number showing no within-subtype signal; and (iii) a demonstration that EGFR's apparent cross-cohort discordance is a positive-control failure, not biology.


## 2. Background and Related Work

TCGA molecularly characterized thousands of tumors with matched survival data [1, 2, 3]; its glioma studies—glioblastoma [1] and lower-grade glioma [2]—include curated molecular classifications (*IDH* status, 1p/19q codeletion, subtype) [2]. CGGA is an independent resource providing RNA-sequencing and clinical data, including *IDH* status and treatment, for a large Chinese glioma cohort [17]. Under the 2021 WHO classification adult diffuse gliomas are defined molecularly and *IDH* mutation is strongly favorable [5, 7], so "lower-grade glioma" spans prognostically distinct entities. EGFR amplification, mutation (EGFRvIII), and overexpression are recurrent in glioblastoma and mark the classical subtype [1, 6]; EGFR expression is characteristically higher in higher-grade and IDH-wildtype tumors—the fact we exploit as a positive control. Throughout we use Kaplan–Meier estimation with the log-rank test [9, 10], Cox regression [11], Harrell's C-index [16], PH testing via scaled Schoenfeld residuals [13], Benjamini–Hochberg FDR [12], and DESeq2 variance-stabilizing transformation [18].

Cross-study irreproducibility of transcriptomic prognostic signatures is well documented: Michiels et al. [24] found prognostic classifiers unstable under resampling across seven microarray studies, and Venet et al. [25] showed that in ER-positive breast cancer almost any gene set—including biologically irrelevant ones—is associated with outcome, because a dominant proliferation axis correlates with most genes. These motivate two distinctions our framework makes explicit: *reproducibility is not validity* (a reproducible association may still be confounded, as with EGFR and *IDH*), and *non-replication has multiple causes*—limited power, population differences, and measurement artifact. Prior cross-platform quality efforts target dataset-level concordance; we instead ask a per-gene, per-comparison question—does each cohort measure *this* gene faithfully, judged against known biology—and show the answer predicts replication specifically where measurement fidelity varies.


## 3. Problem Statement

For TCGA-LGG, TCGA-GBM, and (for validation) CGGA, let $z_{\text{EGFR}}$ be per-SD standardized EGFR expression. We ask whether $z_{\text{EGFR}}$ is **associated** with overall survival within each cohort; whether any association is **independent** of age, grade, and *IDH*/subtype; whether EGFR improves **discrimination** (C-index) beyond established covariates; whether conclusions are **robust to normalization** (VST rather than TPM); and whether the TCGA conclusions **reproduce externally** in CGGA—given that CGGA must first be shown to measure EGFR faithfully enough to test them.


## 4. Methodology

### 4.1 Data acquisition and expression

TCGA RNA-sequencing and clinical data were retrieved from the Genomic Data Commons [3] with `TCGAbiolinks` [4] (harmonized data, GDC Data Release 45.0, 2025-12-04), using STAR-Counts quantifications [8]; molecular annotations came from the TCGA marker-paper fields [1, 2]. Three CGGA datasets were downloaded from the CGGA portal [17] with their clinical tables (all release 20200506): the mRNAseq_693 and mRNAseq_325 RSEM gene-level matrices and the mRNA-array_301 gene-level microarray matrix. TCGA expression used the TPM assay, $\log_2(\text{TPM}+1)$; the CGGA RSEM matrices used $\log_2(\text{RSEM}+1)$; the microarray matrix is distributed on a log-ratio scale and was used untransformed (§4.6). Expression was standardized to unit SD ($z$-scores) within each analysis sample, so the per-SD unit is defined by the patients entering a given model; the common-sample refit (§6.2) therefore uses a marginally different SD from Table 1. Genes were matched by HGNC symbol; in CGGA, EGFR was verified to be a unique feature. We retained primary solid tumors, one sample per patient. Survival models additionally require positive overall-survival (OS) time (time to death as event, else censored at last follow-up); the expression positive controls, which involve no survival model, use all such tumors with a recorded grade. Grade is handled asymmetrically: TCGA-LGG uses the recorded histologic grade (available for 454 of 511), while all TCGA-GBM patients are treated as WHO grade IV by study definition (281 of 282 carry a recorded grade). For CGGA, WHO grade II/III was the lower-grade analog and grade IV the glioblastoma analog.

### 4.2 Statistical analysis

Cox regression modeled $h(t\mid\mathbf{x}) = h_0(t)\exp(\boldsymbol{\beta}^{\top}\mathbf{x})$. For EGFR we fit nested models adding covariates stepwise: unadjusted; $+$age; $+$grade; $+$*IDH* status; $+$full IDH/1p19q subtype (TCGA). We report the EGFR HR per SD ($\text{HR}>1$ = worse survival), Harrell's C-index [16] with 95% CIs, and the PH test for the EGFR term [13]. So that attenuation across the nested models reflects the added covariate rather than the sample change induced by complete-case fitting, the models were also refit on a single common *IDH*-complete subset; EGFR's incremental value was assessed by a likelihood-ratio test and the change in C-index. Kaplan–Meier curves (median split) are descriptive [9, 10]. We treat *IDH* status as a confounder; because *IDH* mutation is an upstream driver and EGFR expression is partly downstream of the IDH-wildtype state, the adjustment could also remove a mediated path—either way, EGFR carrying no signal after conditioning on *IDH* means it adds no information beyond subtype. Ten EGFR/RTK–PI3K-axis genes (EGFR, IDH1, TP53, CDKN2A, NF1, PTEN, PDGFRA, PIK3R1, PIK3CA, RB1) were each entered per SD into age-adjusted and (in TCGA-LGG) age-plus-*IDH*-adjusted Cox models, with Benjamini–Hochberg FDR correction [12].

### 4.3 Normalization sensitivity and external validation

We recomputed the key TCGA EGFR models on DESeq2 VST-transformed STAR counts [18]. For CGGA we repeated the nested EGFR models (adding radiotherapy and chemotherapy status), and—before interpreting—assessed positive controls. The criterion was **specified a priori from established biology**: because EGFR expression is higher in higher-grade and IDH-wildtype glioma, a cohort that measures EGFR faithfully must reproduce a positive EGFR expression–grade correlation, as TCGA does. Clinical controls (IDH and grade must predict survival) verify the survival and annotation data independently of expression, and the between-cohort difference was tested formally (Fisher $r$-to-$z$). To distinguish a genuine CGGA property from an analyst error, we verified the expression–clinical join, reproduced the result independently, and repeated the controls in CGGA-325 and the array-301 microarray. To test pipeline robustness we re-obtained TCGA RNA-seq from recount3 [19]—an independent Monorail reprocessing of the raw reads—and repeated the nested models.

### 4.4 Within-subtype and copy-number analysis

To test whether EGFR is prognostic within IDH-wildtype glioma, we restricted TCGA (LGG+GBM) to IDH-wildtype tumors and fit Cox models for EGFR expression (per SD) and, separately, for gene-level EGFR copy number obtained from the GDC (ASCAT3 gene-level calls). High-level amplification was pre-specified as copy number $\geq 6$, with sensitivity at CN $\geq 5$ and $\geq 7$ and a continuous analysis on $\log_2$ copy number (per SD). Grade throughout is the original histologic grade recorded by TCGA; because under the WHO 2021 scheme EGFR amplification is itself a grade-4 criterion [5], grade-adjusted amplification analyses carry a partial circularity that the pre-2021 grading mitigates. For the amplification analysis we report a power/precision analysis (detectable HR at 80% power via the Schoenfeld relation) alongside the point estimates.

### 4.5 Positive-control-gated screens

We first generalized the positive control to the 500 most variable genes shared by TCGA and CGGA, ranked by the standard deviation of $\log_2$ expression in TCGA. For each gene we recorded the age-adjusted Cox coefficient and the expression–grade Spearman correlation in each cohort; a gene was "QC-concordant" if its grade correlation had the same sign and (pre-specified) magnitude $>0.1$ in both cohorts, and "QC-discordant" otherwise. Among genes prognostic in TCGA (FDR $<0.05$), we computed cross-cohort replication (same sign, $p<0.05$ in CGGA) and sign-flip rates, stratified by QC status, with binomial 95% CIs and a Fisher exact test on the QC $\times$ replication table. The criteria are deliberately asymmetric—FDR-controlled for discovery, nominal $p<0.05$ for replication—the conventional, more permissive choice, which yields a high replication base rate; the QC contrast is a comparison *within* that base rate, not a claim about its level. We repeated the analysis in a second, independent cohort pair (TCGA$\to$CGGA-325), varied the QC magnitude threshold ($|r|>0.05$, $0.1$, $0.2$), and compared against a no-gating baseline. Throughout we use *QC-concordant* and *QC-discordant* for genes passing and failing this check; when the anchor is grade specifically we also call the latter *grade-discordant*, a deliberately neutral term, since the flag marks a candidate measurement artifact but does not by itself prove mis-measurement—a gene's grade relationship could differ biologically across populations.

To evaluate the framework at genome scale, we extended it to the genes shared by TCGA and all three CGGA cohorts, retaining those with TCGA expression SD $>0.5$ and, for runtime, capping the universe at the 8,000 most variable; of these, 6,512 ($81\%$) were prognostic in TCGA at FDR $<0.05$ and form the evaluation set. Such a large majority qualifies because survival in pooled LGG+GBM is grade-dominated, so "prognostic" is a weak label here; the circularity control (§6.8) tests whether the gate survives removal of that shared grade signal. For each gene we scored positive-control concordance as the cross-cohort product of expression–grade Spearman correlations (grade anchor) and, as an alternative, of expression–*IDH* correlations (IDH anchor), and measured how well each score predicts replication by the area under the ROC curve (Mann–Whitney estimator). We benchmarked these biological anchors against two baselines: a *detectability* baseline (minimum cross-cohort mean expression) and a *precision* baseline (the standard error of the replication-cohort Cox coefficient), the latter to test whether the anchor merely proxies how precisely a gene's effect is estimated in the replication cohort. Replication cohorts were CGGA mRNAseq_693, mRNAseq_325, and the mRNA-array_301 microarray; because these datasets are drawn from one resource, we quantified patient overlap between them and repeated the array-301 evaluation after removing shared patients (§6.8). To probe scope in a second tumor type, we repeated the framework in breast cancer using TCGA-BRCA (discovery) and METABRIC (replication); standardized expression and ER/subtype labels came from a common processed resource [20], overall survival from the respective cohort clinical files [21], and ER status and basal/non-basal subtype served as the biological anchors (in place of grade).

### 4.6 Reproducibility and multiplicity

The pipeline is deterministic (random seeds not applicable). Data provenance is pinned: TCGA via the GDC **Data Release 45.0 (2025-12-04)**, and CGGA datasets mRNAseq_693, mRNAseq_325, and mRNA-array_301 (all release 20200506). CGGA data were used under the CGGA data-use terms, which permit academic reuse with citation [17]. All positive-control and amplification thresholds ($|r|>0.1$; CN $\geq 6$) were pre-specified from biology, not tuned to outcomes. Expression matrices were transformed on their native scales: $\log_2(x+1)$ for the TCGA TPM and CGGA RSEM matrices, and no transformation for the CGGA mRNA-array_301 matrix, which is already log-ratio scaled and contains negative values (applying $\log_2(x+1)$ to it discards every value $\leq-1$). Analysis scripts (steps 00–22, the last two independent verification suites that reconstruct every reported sample size and refit every headline model), their console outputs, and full session information are available at https://github.com/edanmn/TCGA-EGFR-Glioma; the workflow is summarized in Figure S1 and software versions in the Supplementary Methods. The study reports one prespecified primary hypothesis (EGFR $\to$ overall survival, adjusted) evaluated across cohorts and robustness checks; these primary $p$-values are nominal and should be read as a small prespecified family. The pathway and systematic screens are explicitly exploratory and FDR-controlled within each analysis; the single nominally significant secondary observation (the inverse EGFR trend in IDH-wildtype glioblastoma, HR $0.87$, $p=0.034$; Table 3) is not corrected and is not interpreted as confirmatory.


## 5. Experimental Setup

After quality control, TCGA-LGG comprised 511 patients (125 deaths) and TCGA-GBM 282 (227 deaths); molecular composition differs sharply (Table S9). CGGA-693 contains 422 primary tumors with a recorded WHO grade—the set used for the expression positive controls (Table 2)—of which 404 also have evaluable overall survival and enter the survival models: 271 lower-grade (WHO II/III, 99 deaths) and 133 glioblastoma (WHO IV, 110 deaths). Analyses used R 4.6.0 [15] with `TCGAbiolinks` 2.40.0 [4], `SummarizedExperiment` 1.42.0, `DESeq2` [18], `survival` 3.8-6 [13], and `survminer` 0.5.2 [14].


## 6. Results

### 6.1 Unadjusted association in TCGA (descriptive)

Under a median split, high EGFR expression was associated with shorter survival in TCGA-LGG (median 2,988 vs.\ 2,282 days; log-rank $p=0.040$) and pooled ($p=4.85\times10^{-3}$), but not TCGA-GBM ($p=0.99$) (Table S1; Figures S2–S4).

### 6.2 In TCGA, EGFR is not independent of IDH status

The nested Cox models (Table 1) show that in TCGA-LGG the EGFR association—strong when unadjusted (HR $=1.59$ per SD) and after age and grade (HR $=1.31$; $p=8.2\times10^{-3}$)—was **abolished after adjustment for *IDH* status** (HR $=1.13$; $p=0.16$) and molecular subtype (HR $=1.14$; $p=0.10$). EGFR alone discriminated poorly (C-index $0.59$, 95% CI $0.51$–$0.67$) versus $0.83$ ($0.79$–$0.87$) for age, grade, and *IDH*.

To confirm that the attenuation reflects *IDH* rather than the sample change from complete-case fitting, we refit the models on the identical *IDH*-complete subset ($n=452$, 105 events): EGFR remained significant after age and grade (HR $=1.308$; $1.072$–$1.595$) and lost significance only when *IDH* was added (HR $=1.13$; $0.95$–$1.34$; $p=0.16$). A likelihood-ratio test confirmed that EGFR adds nothing beyond age, grade, and *IDH* ($\chi^2=1.98$, df $=1$, $p=0.16$), and it did not improve discrimination: C-index $0.834$ ($0.79$–$0.87$) without EGFR versus $0.830$ ($0.79$–$0.87$) with it.

The EGFR term violates proportional hazards here (PH $p=0.004$), so we checked the null does not depend on a constant-hazard summary. A time-varying effect gave a nominal interaction with $\log$ time ($p=0.045$) and no significant main effect ($p=0.080$); splitting follow-up at three years gave HR $=0.96$ ($0.81$–$1.14$, $p=0.65$) before and $1.54$ ($0.95$–$2.51$, $p=0.079$) after. EGFR is thus not prognostic in either period, but the effect is not constant and a late-emerging association cannot be excluded.

In TCGA-GBM, EGFR was not associated when unadjusted; a weak inverse association in the age- and *IDH*-adjusted model (HR $=0.87$; $p=0.030$; $n=265$/214—a whole-cohort model distinct from the IDH-wildtype-restricted analysis of §6.6) was of borderline significance, in the opposite direction, and not robust across model choices, so we do not interpret it as a real effect. The pooled association (HR $=1.35$; $p=1.98\times10^{-7}$) vanished after age adjustment, demonstrating confounding.

\begin{table*}[t]
\centering
\footnotesize
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

### 6.3 Subtype drives the association; pathway screen and normalization checks

In TCGA-LGG survival separates strongly by molecular subtype (Figure S5) and EGFR expression is elevated in the IDH-wildtype group, so its unadjusted association reflects subtype. In the age-adjusted pathway screen (Table S2) six genes were significant; after adding *IDH*, EGFR, PIK3R1, and PTEN lost significance while TP53, RB1, IDH1, CDKN2A, and PIK3CA remained, consistent with subtype-driven signal (IDH1 *expression* is distinct from the *IDH* mutation). No gene was significant in TCGA-GBM. Recomputing on DESeq2 VST counts reproduced the conclusion (Table S3): the lower-grade EGFR HR fell from $1.50$ to $1.11$ ($p=0.21$) after age, grade, and *IDH*, so the pattern is not a TPM artifact.

### 6.4 External validation is inconclusive: CGGA fails EGFR positive controls

A naive CGGA analysis appeared to show the opposite of TCGA: EGFR was null unadjusted (HR $=1.17$; $p=0.16$) but *significant after* adjusting for age, grade, and *IDH* (HR $=1.44$; $p=4.6\times10^{-4}$), and after treatment (HR $=1.37$; $p=6.8\times10^{-3}$) (Table S8). Before interpreting this as a genuine discordance, we tested positive controls.

The CGGA **clinical** controls are strong and correct: IDH-mutant status is highly protective (HR $=0.23$) and grade predicts survival (grade IV vs.\ II HR $=8.0$), confirming the survival, grade, and *IDH* annotations and the sample–clinical join are sound. But the **EGFR-expression** controls fail. EGFR should rise with grade and be higher in IDH-wildtype tumors, as in TCGA; in CGGA it does neither (Table 2; Figure 1). EGFR expression is essentially flat across WHO grade (Spearman $r=-0.03$, vs.\ $+0.18$ in TCGA) and is *higher* in IDH-mutant tumors (mean $\log_2$ $4.17$ vs.\ $3.63$; the reverse of TCGA's $6.43$ vs.\ $7.44$). Expression–grade correlations are systematically attenuated across the whole gene panel in CGGA relative to TCGA (Table S4), indicating a general weakening of expression–phenotype signal in this matrix rather than an EGFR-specific mislabeling.

Two checks establish this is a genuine property of CGGA's EGFR quantification rather than a processing error. First, the expression–clinical join is complete—all 693 samples match a clinical record, of which 422 are primary tumors with a recorded grade—and the flat EGFR–grade relationship was reproduced by an independent computation. Second, the failure **reproduced in two further CGGA datasets spanning a different platform**: a second RNA-seq batch (mRNAseq_325, $n=229$) and a microarray (mRNA-array_301, $n=264$). In all three datasets the EGFR–grade correlation differed significantly from TCGA (Fisher $r$-to-$z$ $p<10^{-3}$; Table 2), while clinical controls passed in every one (IDH-mutant HR $0.22$–$0.32$; high-grade HR $7$–$9$). The criterion was specified a priori from biology, not after seeing the disagreement, and the failure holds across two RNA-seq batches and a microarray—so it is a property of EGFR measurement in this resource, not one pipeline's quantification.

The apparent "IDH-independent EGFR effect" therefore cannot be read as genuine, and CGGA does not provide a valid external test for this gene; a naive reading would have been a spurious positive that only the positive controls exposed. Identifying the mechanism would require raw-read reprocessing, but consistency across two independent batches points to a systematic processing choice rather than random noise.

\begin{table*}[t]
\centering
\footnotesize
\caption{EGFR positive control across TCGA and three CGGA datasets spanning two platforms, all computed by one method: Spearman correlation of EGFR expression with WHO grade among primary tumors with a recorded grade; Fisher $r$-to-$z$ against the TCGA reference; univariable IDH-mutant HR; and grade IV vs.\ grade II HR (grade II and IV patients only). The correlation is positive in TCGA but null/negative in every CGGA dataset, whereas clinical controls are strong and correct throughout. The expression-level failure is reproducible across batches and platforms; the clinical data are sound.}
\begin{tabular}{lccccc}
\toprule
Cohort & $r$(EGFR, grade) & $n$ & $p$ vs.\ TCGA & IDH-mutant HR & Grade IV vs.\ II HR \\
\midrule
TCGA (reference)        & $+0.179$ & 739 & --- & 0.11 & 14.8 \\
CGGA-693 (RNA-seq)      & $-0.033$ & 422 & $4.5\times10^{-4}$ & 0.23 & 8.0 \\
CGGA-325 (RNA-seq)      & $-0.097$ & 229 & $2.5\times10^{-4}$ & 0.22 & 9.3 \\
CGGA-301 (microarray)   & $-0.086$ & 264 & $2.0\times10^{-4}$ & 0.32 & 7.3 \\
\bottomrule
\end{tabular}
\end{table*}

\begin{center}\includegraphics[width=0.68\linewidth]{figures/positive_control_EGFR_grade.png}\end{center}

**Figure 1.** Positive control: standardized EGFR expression by WHO grade. EGFR rises with grade in TCGA (Spearman $r=+0.18$; canonical biology) but is flat in CGGA ($r=-0.03$), indicating that CGGA EGFR expression does not reflect expected EGFR biology. Points are means, error bars ± standard error. The between-cohort difference is significant (Fisher $r$-to-$z$ $p=4.5\times10^{-4}$). TCGA $n=739$ (grade II: 216, III: 241, IV: 282); CGGA-693 $n=422$ (II: 138, III: 144, IV: 140).

### 6.5 The TCGA result is robust to the processing pipeline

To test whether the finding depends on the GDC pipeline, we re-obtained TCGA from recount3 [19], an independent Monorail reprocessing of the raw reads. EGFR again passed the positive controls—expression rose with grade and was higher in IDH-wildtype tumors ($8.92$ vs.\ $7.99$ in lower-grade; $9.20$ vs.\ $7.63$ in glioblastoma)—and the confounding result replicated: the lower-grade EGFR HR fell from $1.54$ (unadjusted) to $1.27$ (age$+$grade) to $1.11$ ($0.94$–$1.31$; $p=0.22$) after *IDH* (Table S5), matching the GDC-pipeline (Table 1) and VST (Table S3) estimates. The result is thus not a processing artifact, and a faithfully reprocessed cohort recovers the EGFR biology CGGA lacks. CGGA raw reads are controlled-access, so the equivalent reprocessing there—the only route to new-patient external validation—remains the key open step.

### 6.6 EGFR is not prognostic *within* IDH-wildtype glioma

A stronger test than confounder adjustment is whether EGFR carries any prognostic signal *within* the aggressive IDH-wildtype subtype, where subtype confounding no longer applies. Among IDH-wildtype gliomas (TCGA LGG+GBM; $n=338$, 253 deaths), EGFR expression was not prognostic (Table 3): unadjusted HR $=1.00$ ($0.89$–$1.13$), age-adjusted $0.93$ ($p=0.22$), and $+$grade $0.90$ ($p=0.10$), with a weak inverse—not adverse—trend in IDH-wildtype glioblastoma (HR $=0.87$, $p=0.034$, not corrected for multiplicity). This is a well-powered null: with 253 events, 80% power extends to HR $\geq1.19$ per SD, and the observed interval excludes HR $>1.13$. The PH assumption held for the EGFR term ($p=0.084$). Grade here is the *original histologic* grade recorded by TCGA (pre-2021), not the molecular WHO 2021 grade.

Gene-level copy number for the same tumors tells the same story. Of the primary tumors with both copy-number calls and *IDH* status, high-level EGFR amplification was strongly enriched in IDH-wildtype disease ($55\%$ of 167, vs.\ $2.3\%$ of 222 IDH-mutant), confirming the copy-number data are faithful—a positive control the CGGA expression data failed. The 166 IDH-wildtype tumors that also have evaluable survival form the Cox models below. Within IDH-wildtype glioma we found *no evidence* that amplification predicts survival: the continuous effect on $\log_2$ copy number was null (HR $=1.00$ per SD, $0.84$–$1.20$, $p=0.98$, age-adjusted), and binary high-level amplification was null across thresholds (CN $\geq 5$/$6$/$7$: HR $1.09$/$1.00$/$0.92$, all $p>0.6$; Table 3). The binary tests are underpowered—with 127 events and $55\%$ amplified, 80% power (two-sided $\alpha=0.05$) extends only to HR $\geq 1.65$—so we exclude a large but not a modest amplification effect. The continuous analysis is better powered (80% power for HR $\geq1.28$ per SD) and its interval excludes that value; and the expression analysis, which captures amplification's downstream impact, shows a precise null. EGFR therefore carries no detectable within-subtype prognostic signal.

\begin{table*}[t]
\centering
\footnotesize
\caption{EGFR is not a within-subtype prognostic factor in IDH-wildtype glioma (TCGA). Expression (per 1 SD) shows a well-powered null; copy-number amplification shows no evidence of an effect (continuous, reasonably powered) and null across binary thresholds, though binary tests are underpowered (80\% power only for HR $\geq 1.65$). Grade is the pre-2021 histologic grade.}
\begin{tabular}{lllll}
\toprule
Marker & Model & HR (95\% CI) & $p$ & $n$ / events \\
\midrule
Expression   & Unadjusted                 & 1.00 (0.89--1.13) & 0.99 & 338 / 253 \\
Expression   & $+$ age                    & 0.93 (0.82--1.05) & 0.22 & 338 / 253 \\
Expression   & $+$ age, grade             & 0.90 (0.80--1.02) & 0.10 & 330 / 246 \\
Expression   & IDH-wt GBM, $+$ age        & 0.87 (0.76--0.99) & 0.034 & 244 / 203 \\
Amplification& $\log_2$ copy number (per SD), $+$age & 1.00 (0.84--1.20) & 0.98 & 166 / 127 \\
Amplification& CN $\geq 5$ vs.\ not, $+$age & 1.09 (0.76--1.57) & 0.63 & 166 / 127 \\
Amplification& CN $\geq 6$ vs.\ not, $+$age & 1.00 (0.70--1.42) & 0.98 & 166 / 127 \\
Amplification& CN $\geq 7$ vs.\ not, $+$age & 0.92 (0.65--1.32) & 0.66 & 166 / 127 \\
\bottomrule
\end{tabular}
\end{table*}

### 6.7 A positive-control-gated framework for cross-cohort validation

The CGGA episode motivates a general safeguard: before trusting a cross-cohort comparison, verify both cohorts measure the gene faithfully. We formalized this on the 500 most-variable genes shared by TCGA and CGGA, with the expression–grade correlation in each cohort as a *pre-specified* positive control. Of the genes prognostic in TCGA (FDR $<0.05$), replication in CGGA-693 was $92.8\%$ (95% CI $89.9$–$95.0$) among QC-concordant genes versus $24.3\%$ ($11.8$–$41.2$) among QC-discordant genes (Fisher exact $p=4.9\times10^{-21}$, OR $=39.2$); both sign-flips observed in the panel fell among QC-discordant genes, though with only two such events this is suggestive rather than decisive (Table S6; Figure 2). The contrast reproduced in a second cohort pair with no shared patients (TCGA$\to$CGGA-325: $95.4\%$ vs.\ $19.2\%$; Fisher $p=1.1\times10^{-20}$) and exceeded the no-gating baseline ($87.3\%$, $84.0$–$90.2$); $88\%$ of the 500-gene panel passed QC. The contrast persisted across QC thresholds but narrowed as the threshold became more stringent, since a stricter bar moves well-measured genes into the failing set (QC-concordant vs.\ discordant replication: $91.0\%$ vs.\ $5.0\%$ at $|r|>0.05$; $92.8\%$ vs.\ $24.3\%$ at $0.1$; $96.2\%$ vs.\ $53.6\%$ at $0.2$). For EGFR the artifact interpretation is independently supported by its canonical grade/IDH biology (§6.4). Gating thus concentrates non-replicating associations into a small, flaggable minority.

\begin{center}\includegraphics[width=0.64\linewidth]{figures/systematic_method.png}\end{center}

**Figure 2.** Cross-cohort expression QC for TCGA-prognostic genes. Each point is a gene; axes are its expression–grade correlation in TCGA (x) and CGGA (y). QC-concordant genes (blue) lie along the diagonal (measured consistently); QC-discordant genes (red), including EGFR, fall off it and are where cross-cohort prognostic replication breaks down.

### 6.8 The positive control predicts replication genome-wide

To test whether the positive control is a general safeguard rather than an EGFR-specific observation, we scaled it to **6,512 genes prognostic in TCGA** (FDR $<0.05$, age-adjusted) and asked, for each, whether its association replicated in three CGGA cohorts (Table 4; Figure 3). The grade anchor predicted replication with **AUC $0.82$ (CGGA-693), $0.93$ (CGGA-325), and $0.92$ (microarray-301)**; the *IDH* anchor gave AUC $0.77$–$0.84$. A simple grade-concordance flag identified non-replicating associations with 70% precision and 60% recall in CGGA-693, and $32\%$ of the 6,512 prognostic genes were flagged.

Both baselines are uninformative by comparison. Detectability (mean expression) ran at or near chance (AUC $0.36$–$0.53$), so highly expressed genes are not more replicable. More importantly, *estimation precision* in the replication cohort—the standard error of its Cox coefficient—also failed to predict replication (AUC $0.44$–$0.59$), addressing the concern that the anchor merely proxies how noisily a gene is measured: were that so, precision would predict replication at least as well, and it does not.

The three datasets come from one resource and are not fully independent: CGGA-693 and CGGA-325 share no patients, but of the 264 analysed array-301 patients, 61 (23%) also appear in CGGA-325 and 6 in CGGA-693. Removing all 67 leaves 197 patients and barely moves the discrimination—grade-anchor AUC $0.91$ ($0.91$–$0.92$) versus $0.92$ ($0.92$–$0.93$)—though the replication rate falls from $0.68$ to $0.58$, as expected from the smaller sample. The cross-platform finding therefore does not depend on shared patients, and EGFR's positive-control failure persists in the de-duplicated set ($r=-0.084$; Fisher $p=6.6\times10^{-4}$). EGFR falls in the flagged, grade-discordant minority, confirming that its cross-cohort behavior (§6.4) is a measurement artifact the framework catches, not a biological discordance.

**Controlling for circularity.** Because survival in pooled glioma is grade-driven, the prognostic set is enriched for grade-correlated genes, and the grade-anchor QC could predict replication partly *because* both prognostic status and replication share the grade signal it measures. We tested this directly (Table S7). Restricting to associations that remain prognostic *after adjustment for grade and IDH*—i.e., subtype-independent signal—the QC still predicted replication above chance, though this set is small and the interval correspondingly wide (183 genes; AUC $=0.65$, 95% CI $0.55$–$0.74$). It also predicted replication for genes with low grade-correlation ($|r|<0.2$; AUC $=0.75$, $0.72$–$0.78$). The AUC is higher for strongly grade-correlated genes ($0.84$) than for low-grade-correlation genes ($0.75$), so part of the headline value does reflect shared grade signal; but the effect persists for subtype-independent associations, indicating the positive control captures measurement fidelity *beyond* the anchor covariate itself rather than acting tautologically.

\begin{table*}[t]
\centering
\footnotesize
\caption{Genome-scale evaluation: a gene's biological positive-control concordance predicts whether its TCGA prognostic association replicates in a separate cohort (AUC). The three replication datasets come from one resource and are not fully independent (§6.8 quantifies and controls for patient overlap). Grade and IDH anchors are strongly predictive across three cohorts and two platforms; detectability and precision baselines are at or near chance. $n=6,512$ TCGA-prognostic genes (FDR $<0.05$); AUCs are given to three decimals here and rounded to two in the text. Bootstrap 95\% CIs are tight given $n\approx6,500$ (CGGA-693 grade anchor: $0.815$, $0.81$–$0.83$); Table S7 gives the circularity control.}
\begin{tabular}{lccccc}
\toprule
Replication cohort & Replication rate & AUC (grade) & AUC (IDH) & AUC (detectability) & AUC (precision) \\
\midrule
CGGA-693 (RNA-seq)    & 0.622 & 0.815 & 0.771 & 0.358 & 0.591 \\
CGGA-325 (RNA-seq)    & 0.701 & 0.931 & 0.839 & 0.366 & 0.439 \\
CGGA-301 (microarray) & 0.684 & 0.924 & 0.842 & 0.525 & 0.505 \\
\bottomrule
\end{tabular}
\end{table*}

\begin{center}\includegraphics[width=0.64\linewidth]{figures/method_auc.png}\end{center}

**Figure 3.** Biological positive controls predict cross-cohort replication of prognostic associations. AUC for predicting replication of TCGA-prognostic genes in each CGGA cohort, by metric; the grade and IDH anchors far exceed both the detectability and the estimation-precision baselines (dashed line = chance).

### 6.9 Scope: a well-matched second tumor type

To test whether the framework is a *universal* replication predictor or a targeted detector of measurement-artifact discordance, we applied it to a second, unrelated tumor type: breast cancer, using TCGA-BRCA [23] (discovery, RNA-seq) and METABRIC [22] (replication, microarray), with ER status and basal/non-basal subtype as biological anchors, expression from a common processed resource [20] and survival from the cohort clinical files [21]. Here the positive control was at most weakly predictive. In the TCGA-BRCA$\to$METABRIC direction the ER-anchor AUC was $0.60$ (95% CI $0.44$–$0.75$), an interval that includes chance and rests on only 59 prognostic genes; in the better-powered reverse direction (2,568 genes) it was $0.56$ ($0.52$–$0.59$), which *excludes* chance but is far below the $0.82$–$0.93$ seen in glioma. The basal anchor reached at most $0.54$. We therefore do not claim the gate is uninformative in breast, only that its effect is small where glioma's is large. Three non-exclusive reasons make breast diagnostic of *when the method helps*: the cohorts measure genes highly concordantly (genome-wide agreement of expression–ER correlations $r=0.84$), leaving little artifact to flag; in ER-positive breast almost any gene set is prognostic [25], so non-replication reflects weak, proliferation-dominated signal and limited power (TCGA-BRCA has 67 events) rather than measurement fidelity; and the breast anchors (ER, basal subtype) differ from the glioma anchors (grade, *IDH*), confounding anchor informativeness with cohort concordance. One cohort pair cannot separate these, but all three predict the gate should add little—as observed, in contrast to glioma (§6.4, §6.8). The framework is therefore best understood not as a general replication oracle but as a **targeted detector of measurement-fidelity–driven discordance**: its predictive value scales with how much cross-cohort measurement heterogeneity exists, and it correctly adds little when both cohorts are clean.


## 7. Discussion

EGFR expression is a subtype artifact in TCGA: elevated in IDH-wildtype tumors, and once *IDH* or subtype is modeled it carries no independent prognostic information and adds no discrimination—robust to VST normalization and to an independent raw-read reprocessing (Tables S3, S5). Crucially the null extends *within* the IDH-wildtype subtype, where neither expression nor amplification predicts survival (Table 3). EGFR is therefore a between-subtype marker of the IDH-wildtype state—a more complete characterization than "confounded by IDH," since it rules out a residual within-subtype effect that grade-only adjustment could have missed.

Our attempt to validate in CGGA is instructive precisely because it initially appeared to contradict TCGA. Positive controls showed that CGGA EGFR expression tracks neither grade nor the canonical IDH relationship while CGGA clinical signals are strong, and this failure reproduced across three datasets and two platforms, with the CGGA-693 expression–clinical join verified sample by sample. We guarded against the obvious objection—that we discarded the cohort that disagreed—by specifying the criterion from biology in advance, testing the difference formally, and replicating the failure across batches. The lesson generalizes: across 6,512 prognostic genes, positive-control concordance predicts replication with AUC $0.82$–$0.93$ while detectability and precision baselines sit at or near chance ($0.36$–$0.59$) (§6.8). This reframes cross-cohort validation as a two-part question—*does the association reproduce, and does the replication cohort measure the gene faithfully?*—with the second part checkable a priori from biology. Gating would have prevented the spurious EGFR "discordance," and 32% of prognostic genes are QC-flagged in CGGA-693, carrying the same risk. The safeguard is cheap and agnostic to gene, cohort, or platform, but *scoped* (§6.9): where both cohorts measure genes concordantly it has little to flag.

For cross-cohort prognostic studies we recommend: (1) identify a covariate known to associate with the gene (grade for oncogenes, hormone receptor status for breast cancer genes, subtype for pathway members); (2) verify it reproduces with consistent sign and magnitude in both cohorts before interpreting prognostic comparisons; and (3) flag discordant genes (correlations differing by Fisher $r$-to-$z$ $p<0.05$, opposite signs, or $|r|<0.1$ in either cohort) for review or exclusion. The gate needs one known covariate and negligible computation.


## 8. Limitations

First, the primary finding rests on a single cohort: because the CGGA expression matrix could not measure EGFR faithfully, no successful external validation was obtained, so the conclusion is internally robust and biologically consistent but externally unreplicated for this gene. Second, the failure is reproducible (three datasets, two platforms) and not an error in our handling of the CGGA-693 RSEM matrix—join verified, result independently recomputed—but its cause within CGGA (normalization, gene model, or platform) was not identified. Third, adjusted models used complete cases without imputation, though the key nested comparison was refit on a common sample. Fourth, the pooled TCGA analysis mixes two separately processed studies, and grade IV is collinear with the glioblastoma study, so batch and grade cannot be separated. Fifth, the EGFR term violates PH in the *IDH*-adjusted lower-grade model, and the global PH test is violated in that model as well ($p=3.6\times10^{-5}$), driven by grade and *IDH* rather than EGFR; period-specific estimates are null before three years and non-significantly elevated afterwards (§6.2), so a late-emerging effect cannot be excluded, though the EGFR term did satisfy PH in the within-subtype model ($p=0.084$). Sixth, the systematic method was benchmarked against no-gating, detectability, and estimation-precision baselines but not a broad panel of alternative QC metrics, and "grade-discordant" flags a candidate—not proven—measurement artifact. Seventh, the three replication datasets come from a single resource and are not fully independent (61 patients shared between array-301 and CGGA-325 and 6 with CGGA-693, quantified and controlled for in §6.8), and all replication is against one consortium in one tumor type. Eighth, and conceptually important, the framework assesses cross-cohort *reproducibility*, not *validity*: it can bless a reproducible association that is nonetheless confounded (as EGFR's is by *IDH*), so it complements rather than replaces causal scrutiny. Ninth, part of the genome-scale AUC reflects shared grade signal between the "prognostic" definition and the grade anchor; the effect persists but is weaker for subtype-independent associations (Table S7).


## 9. Future Work

True new-patient external validation still requires an independent cohort whose EGFR expression passes the positive controls. Reprocessing CGGA from raw reads would be the natural test, but its raw reads are controlled-access; an openly available, faithfully quantified glioma RNA-seq cohort with *IDH*, grade, and survival is therefore needed. For the framework, the informative next step is to evaluate cohort pairs that *do* differ in measurement fidelity across several tumor types—the regime the gate targets—and better-powered discovery cohorts, since the breast pair (§6.9) was both well-matched and underpowered.


## 10. Conclusion

Cross-cohort prognostic replication in transcriptomics is unreliable, and we show that a simple biological positive control detects the subset of discordances that stem from faulty measurement: across 6,512 genes, three glioma cohorts, two platforms, and two anchors, a gene's expression–covariate concordance forecasts whether its prognostic association replicates (AUC $0.82$–$0.93$), while detectability and precision baselines do not. Its benefit is scoped—in a well-matched breast pair the gate had little to flag—so it is best understood as a targeted detector of measurement-fidelity–driven discordance rather than a universal replication oracle. The safeguard is cheap: verify that both cohorts measure the gene faithfully before trusting a cross-cohort comparison. The framework is motivated by EGFR, for which naive analysis manufactures a spurious discordance—prognostic in TCGA only until *IDH* is modeled, carrying no within-IDH-wildtype signal, and failing the positive control in CGGA. EGFR is therefore a marker of the IDH-wildtype state, not an independent prognostic factor.


## References

\footnotesize

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
13. Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model*. New York: Springer; 2000.
14. Kassambara A, Kosinski M, Biecek P. survminer: survival curves using 'ggplot2'. R package v0.5.2; 2024.
15. R Core Team. *R: A Language and Environment for Statistical Computing*. Vienna, Austria: R Foundation for Statistical Computing; 2026.
16. Harrell FE Jr, Lee KL, Mark DB. Multivariable prognostic models: evaluating assumptions and adequacy, and measuring and reducing errors. *Stat Med*. 1996;15(4):361–387.
17. Zhao Z, Zhang KN, Wang Q, et al. Chinese Glioma Genome Atlas (CGGA): a resource with functional genomic data from Chinese glioma patients. *Genomics Proteomics Bioinformatics*. 2021;19(1):1–12.
18. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*. 2014;15(12):550.
19. Wilks C, Zheng SC, Chen FY, et al. recount3: summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biol*. 2021;22(1):323.
20. Processed TCGA-BRCA and METABRIC datasets (Moanna). Zenodo; doi:10.5281/zenodo.4326602.
21. Cerami E, Gao J, Dogrusoz U, et al. The cBio Cancer Genomics Portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discov*. 2012;2(5):401–404.
22. Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature*. 2012;486(7403):346–352.
23. Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature*. 2012;490(7418):61–70.
24. Michiels S, Koscielny S, Hill C. Prediction of cancer outcome with microarrays: a multiple random validation strategy. *Lancet*. 2005;365(9458):488–492.
25. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol*. 2011;7(10):e1002240.
