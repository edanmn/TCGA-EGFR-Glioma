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

Prognostic gene associations from one cohort frequently fail to replicate in another, and it is hard to know in advance which will. A natural safeguard is a *biological positive control*: trust a cross-cohort comparison only if both cohorts measure the gene faithfully, checked by whether its expression tracks a known covariate consistently. We ask two questions about such a gate---what it actually detects, and how its performance should be scored---and answer both with controlled experiments rather than observational agreement. **Under injected ground truth**, with corruption applied to known genes at known doses in one half of TCGA, detectability is strongly mode-dependent: compression of between-stratum signal is caught at AUC $0.85$--$0.64$, permutation at $0.63$--$0.86$, additive noise at $0.57$--$0.78$, and **flooring essentially not at all ($\leq0.60$)**---a documented blind spot. **Against a matched null**---a held-out, same-population cohort in which no artifact can exist---the score such gates are usually assessed with, the anchor product $r_D r_R$, still reaches AUC $0.79$--$0.90$, so its distance from chance carries no information; the discovery effect size, a baseline the design is apt to omit, reaches $0.74$--$0.78$ alone. Scoring disagreement directly, as $-|r_D-r_R|$, drops the null to $0.53$--$0.61$ and is positive in all three glioma pairs, though between-split variance supplies $77$--$84\%$ of total uncertainty and only one pair is conclusive once that is counted. Applied to EGFR in glioma, the gate flags a gene that is prognostic in TCGA only until *IDH* is modeled (HR $1.31\to1.13$, $p=0.16$), carries no within-IDH-wildtype signal (HR $1.00$, $0.89$--$1.13$), and in three CGGA datasets shows a $2$--$4\times$ collapse in the variance explained by grade and *IDH* that cohort composition does not account for.


**Keywords:** cross-cohort validation; prognostic gene expression; positive controls; calibration; glioma; EGFR; *IDH*; TCGA; CGGA.


## 1. Introduction

Cross-study irreproducibility of transcriptomic prognostic signatures is well documented [24, 25], and validating a discovery cohort's finding in a second cohort is now routine---TCGA$\to$CGGA is among the most used such pairs in glioma [29]. That design carries an assumption rarely tested: that both cohorts measure each gene faithfully enough for disagreement to be informative. A *biological positive control* makes the assumption checkable. Because tumor grade and *IDH* status have known relationships to expression for many genes, a cohort that measures a gene faithfully should reproduce that relationship; one that does not is a cohort whose disagreement means little.

This paper asks whether such a gate works and how one would know. Our contributions are:

1. **A ground-truth evaluation** (§4.1). We inject five mechanistically distinct measurement failures at known doses into one half of TCGA and score the gate against the genes actually corrupted. It detects compression, permutation, noise and contamination dose-dependently, and misses flooring almost entirely. No observational analysis can establish either fact.
2. **A calibration of the metric** (§4.2). The natural way to validate such a gate---its AUC for predicting replication---is confounded. A same-population null reaches $0.79$--$0.90$, so the gap from chance is uninformative; the confound lies in the score's functional form, and a difference form largely removes it.
3. **An honest uncertainty accounting** (§4.2). Between-split variance supplies $77$--$84\%$ of the total in this design, so single-split intervals understate uncertainty $3$--$8\times$. We report what survives.
4. **A worked example** (§4.4): EGFR in glioma, where naive analysis manufactures a cross-cohort discordance the gate correctly flags.


## 2. Background and Related Work

TCGA molecularly characterized thousands of tumors with matched survival [1, 2, 3]; its glioma studies [1, 2] include curated *IDH* status, 1p/19q codeletion and subtype. CGGA provides independent RNA-sequencing and clinical data for a Chinese glioma cohort [17]. Under the 2021 WHO classification adult diffuse gliomas are defined molecularly and *IDH* mutation is strongly favorable [5, 7]. EGFR amplification and overexpression are recurrent in glioblastoma and mark the classical subtype [1, 6]; EGFR expression is characteristically higher in higher-grade and IDH-wildtype tumors---the fact we exploit as a positive control. Throughout we use Kaplan--Meier estimation [9, 10], Cox regression [11], Harrell's C-index [16], PH testing via scaled Schoenfeld residuals [13], Benjamini--Hochberg FDR [12], and DESeq2 VST [18].

Four bodies of work delimit what is new here. **sigQC** [26] standardizes quality control for gene *signatures* across datasets---expression, variability, autocorrelation, scoring stability---but operates on gene sets and does not anchor individual genes to a covariate. **Removing Unwanted Variation** [27] uses *negative* control genes to estimate and subtract unwanted variation; we use the mirror construct, positive controls, and gate rather than correct. Third, the replication-prediction literature outside genomics has established that a discovery study's **effect size** predicts replication well (AUC $\approx0.77$ [28]); this is precisely the baseline a positive-control evaluation is apt to omit, and §4.2 shows it recovers much of the apparent performance. Fourth, Michiels et al. [24] and Venet et al. [25] motivate two distinctions we keep explicit: *reproducibility is not validity*, and *non-replication has multiple causes*.

Against that background our contribution is less the gate than the demonstration of how one must be scored: against a same-population null rather than chance, using a difference rather than a product statistic, and validated against injected ground truth rather than observational agreement.


## 3. Methods

**Cohorts.** TCGA RNA-sequencing and clinical data came from the Genomic Data Commons [3] via `TCGAbiolinks` [4] (GDC Data Release 45.0, 2025-12-04), STAR-Counts quantifications [8], with molecular annotations from the marker papers [1, 2]. Three CGGA datasets [17] (all release 20200506) served as replication cohorts: mRNAseq_693, mRNAseq_325 and the mRNA-array_301 microarray. TCGA used $\log_2(\text{TPM}+1)$ and the CGGA RSEM matrices $\log_2(\text{RSEM}+1)$; the microarray is distributed log-ratio scaled and was used untransformed, since $\log_2(x+1)$ discards every value $\leq-1$. We retained primary solid tumors, one sample per patient, with positive overall-survival time. TCGA-LGG comprised 511 patients (125 deaths) and TCGA-GBM 282 (227 deaths); CGGA-693 contributed 422 primary tumors with a recorded grade, 404 with evaluable survival. Grade is handled asymmetrically: TCGA-LGG uses recorded histologic grade (454 of 511), while all TCGA-GBM are WHO grade IV by study definition. A second tumor type used TCGA-BRCA [23] and METABRIC [22], with expression and ER/subtype labels from a common processed resource [20] and survival from the cohort clinical files [21].

**The gate.** For a gene, the anchor correlation is its Spearman correlation with the anchor covariate (grade, or *IDH*) within a cohort. The *product* score is $r_D r_R$ and the *difference* score is $-|r_D-r_R|$, where $D$ and $R$ index discovery and replication. A gene is *prognostic* if its age-adjusted Cox coefficient is significant at BH-FDR $<0.05$ in discovery, and *replicated* if the replication cohort's coefficient has the same sign at $p<0.05$. Performance is the AUC for predicting replication (Mann--Whitney estimator).

**Ground truth (§4.1).** TCGA was split in half and a random $20\%$ of genes corrupted in one half only, at three doses each of five modes: permutation of a fraction of values; additive Gaussian noise at $k$ SD; compression of the between-stratum component toward the grand mean by $\lambda$ (lowering between-stratum variance while preserving within-stratum variance); flooring below the $q$th quantile; and contamination, replacing a fraction $w$ of the gene by a scale-matched *donor* drawn at random from the same 2,000-gene set. The donor is a real gene and therefore carries its own (grade, *IDH*) structure, which is why contamination leaves the stratum-explained variance share near baseline while destroying the gene's identity (Table 1); it models a mixed or mis-mapped feature rather than dilution by noise. Donors are drawn from the working matrix as corruption proceeds, so a donor may itself already be corrupted; at a 20\% corruption rate this affects about one donor in five and does not bear on the detection scores, which are computed against the known corruption assignment. The gate was scored against the genes actually corrupted (2,000 most-variable genes, 4 replicates per cell, no composition shift). A separate arm biased the split so grade and *IDH* couple more tightly in one half ($\delta=0$--$0.4$).

**Matched null (§4.2).** Within each cohort pair, one discovery set was evaluated against two replication cohorts of identical size: a disjoint, held-out subsample of the discovery cohort, and the real replication cohort subsampled to match. Five settings were run (TCGA$\to$CGGA-693, -325, array-301; TCGA-BRCA$\leftrightarrow$METABRIC) with five quality metrics and $|\beta_D|$ as a negative control, on the gene universe of §4.3. Intervals are paired bootstraps over genes; a cluster bootstrap over 200 co-expression clusters tests the independence assumption; and $K=6$ splits per setting combine within- and between-split variance by $\mathrm{Var}_{\text{total}} = \overline{\mathrm{Var}}_{\text{within}} + (1+1/K)\,\mathrm{Var}_{\text{between}}$.

**Genome-scale evaluation (§4.3).** Genes shared by TCGA and all three CGGA cohorts with TCGA SD $>0.5$, capped at the 8,000 most variable; 6,512 ($81\%$) were prognostic in TCGA at FDR $<0.05$. Anchors were benchmarked against detectability (minimum cross-cohort mean expression), replication-cohort estimation precision, and discovery effect size, with incremental value assessed by adding the anchor to a logistic model already containing effect size.

**EGFR analysis (§4.4).** Nested Cox models added covariates stepwise (unadjusted; $+$age; $+$grade; $+$*IDH*; $+$subtype), reporting HR per SD, C-index [16] and the PH test [13], refit on a common *IDH*-complete subset. Within IDH-wildtype tumors we fit expression and gene-level copy number (ASCAT3; high-level amplification pre-specified at CN $\geq6$). Robustness used DESeq2 VST [18] and an independent Monorail reprocessing from recount3 [19]. Cross-cohort discordance was decomposed into a composition component and a residual by transporting TCGA's within-(grade, *IDH*)-stratum means and variances onto each replication cohort's stratum weights.

**Reproducibility.** Analyses used R 4.6.0 [15] with `TCGAbiolinks` 2.40.0 [4], `survival` 3.8-6 [13], `survminer` 0.5.2 [14] and `DESeq2` [18]. The pipeline is deterministic except for the resampling analyses, which are seeded. All thresholds ($|r|>0.1$; CN $\geq6$) were pre-specified from biology. Analysis scripts (steps 00--40), console outputs and session information are at https://github.com/edanmn/TCGA-EGFR-Glioma, on branch `revision/calibration-and-ground-truth`; steps 23--40, which produce every result in §4.1--4.3, are on that branch and not on the default one. Five verification suites share no code with the analysis pipeline: they reconstruct every reported sample size, refit every headline model, assert every figure transcribed into the manuscript against its source output, and — since a value silently *dropped* from the manuscript would pass all of those — check that each audited value still appears in the manuscript text. `run_all.R` stops non-zero if any fails.


## 4. Results

### 4.1 What the gate detects, under known ground truth

Detectability differs sharply by failure mode (Table 1). *Compression* of the between-stratum signal is the best detected (AUC $0.848$, $0.730$, $0.639$ at $\lambda=0.25$, $0.50$, $0.75$). *Permutation* ($0.630$--$0.858$), *additive noise* ($0.567$--$0.779$) and *contamination* ($0.543$--$0.734$) are also detected dose-dependently. **Flooring---censoring low values, as a detection limit or saturating assay would---is a blind spot**: AUC $0.498$, $0.528$, $0.599$ across doses, barely above chance even where corrupted genes' replication falls to $0.82$ against $0.95$ for clean genes. A gate of this kind should not be relied on where saturation is the suspected failure.

Compression is also the mode that matches EGFR's behavior in CGGA (§4.4), though identifying it needs two statistics rather than one. On variance share alone the match is not unique: at $\lambda=0.50$ compression leaves $0.136$ against the uncorrupted baseline of $0.339$---a ratio of $0.40$, inside EGFR's observed $0.23$--$0.42$ band---but permutation at $50\%$ ($0.28$) and noise at $k=2$ ($0.23$) also fall in that band. What separates them is *within*-stratum variance, which compression preserves by construction while permutation and noise inflate it; EGFR's is unchanged from TCGA's ($1.84$ vs.\ $1.88$, §4.4). On that two-statistic match the gate would have flagged an EGFR-like failure at AUC $0.730$.

Injecting a *confounding-structure shift* instead does not raise the gate's AUC; it lowers it, from $0.943$ at no shift to $0.609$ at the largest, because replication itself collapses ($0.957\to0.450$) as the coupling difference grows from $0.045$ to $0.796$.

\begin{table*}[t]
\centering
\footnotesize
\caption{Controlled experiment with known ground truth: five mechanistically distinct measurement failure modes, injected into a random 20\% of genes in one half of TCGA only (2,000 most-variable genes, 4 replicates per cell, no composition shift). \textbf{All fifteen mode $\times$ dose cells are shown}; none is omitted. AUC is scored against the genes actually corrupted. ``Variance share'' is the fraction of a corrupted gene's variance explained by the six (grade, \emph{IDH}) strata afterwards; the uncorrupted baseline is 0.339. Compression---the signature EGFR shows in CGGA---is the best-detected mode; flooring is a blind spot. Contamination is the one mode whose variance share barely moves (0.28--0.33 against the 0.339 baseline) while its replication collapses: the donor is another gene from the same variable set, which carries its own (grade, \emph{IDH}) structure, so the stratum-explained \emph{quantity} is preserved while the gene's \emph{identity} is destroyed. Its row should therefore not be read as a variance-share signature (\S3). Replicate SDs of the AUC are 0.004--0.026 throughout.}
\begin{tabular}{llcccc}
\toprule
Failure mode & Mechanism & Dose & Corrupt / clean replication & Variance share & AUC (detect corrupted) \\
\midrule
Compression   & shrink between-stratum signal by $\lambda$ & $\lambda=0.25$ & 0.43 / 0.93 & 0.046 & \textbf{0.848} \\
Compression   &                                            & $\lambda=0.50$ & 0.80 / 0.95 & 0.136 & \textbf{0.730} \\
Compression   &                                            & $\lambda=0.75$ & 0.91 / 0.95 & 0.236 & 0.639 \\
Permutation   & shuffle a fraction of values               & 25\%           & 0.92 / 0.97 & 0.197 & 0.630 \\
Permutation   &                                            & 50\%           & 0.67 / 0.95 & 0.098 & 0.728 \\
Permutation   &                                            & 75\%           & 0.24 / 0.96 & 0.035 & 0.858 \\
Additive noise& $x+N(0,k\,\mathrm{sd})$                    & $k=0.5$        & 0.95 / 0.96 & 0.268 & 0.567 \\
Additive noise&                                            & $k=1.0$        & 0.82 / 0.94 & 0.172 & 0.653 \\
Additive noise&                                            & $k=2.0$        & 0.54 / 0.93 & 0.079 & 0.779 \\
Contamination & mix in a donor gene at weight $w$          & $w=0.25$       & 0.89 / 0.93 & 0.308 & 0.543 \\
Contamination &                                            & $w=0.50$       & 0.62 / 0.95 & 0.284 & 0.667 \\
Contamination &                                            & $w=0.75$       & 0.49 / 0.94 & 0.329 & 0.734 \\
Flooring      & censor below the $q$th quantile            & $q=0.25$       & 0.94 / 0.95 & 0.341 & \textbf{0.498} \\
Flooring      &                                            & $q=0.50$       & 0.90 / 0.94 & 0.263 & \textbf{0.528} \\
Flooring      &                                            & $q=0.75$       & 0.82 / 0.94 & 0.184 & \textbf{0.599} \\
\bottomrule
\end{tabular}
\end{table*}

### 4.2 Calibrating the metric

An AUC is interpretable only against the right null, and chance is not it. Giving one discovery set two replication cohorts of identical size---one a disjoint, held-out subsample of the discovery cohort itself, in which no cross-cohort artifact can exist---shows the anchor *product* is **severely confounded**: against that null it still reaches AUC $0.79$--$0.90$ in glioma and $0.72$ in the one evaluable breast direction (Table 2). It does not follow that nothing was detected: for CGGA-693 the null exceeds the real cohort ($-0.075$), but for CGGA-325 and array-301 the real cohort exceeds the null ($+0.044$, $+0.132$). Tightening the null arm's threshold until its replication rate matches the real arm's changes nothing ($0.931\to0.921$), so this is not a base-rate artifact.

**The confound is in the functional form.** The product $r_D r_R$ is large when *both* correlations are large, so it rewards anchor magnitude---which tracks effect size---as much as agreement. The difference $-|r_D-r_R|$ scores only disagreement, and its null AUC is $0.53$--$0.61$ rather than $0.79$--$0.90$. The argument for preferring it is mechanical rather than empirical---a product of two correlations cannot separate *agreement* from *magnitude*, whichever data it is applied to---but the choice was nonetheless made after seeing the product form's null on these three cohort pairs, and its advantage is then reported on the same three. We had no fourth glioma pair to hold out. Readers should treat the size of the difference form's margin over the product form as an in-sample quantity; what does not depend on the selection is the null level of each form, which is a property of the statistic. On a single split its point estimate exceeds the null in all three glioma pairs ($+0.138$, $+0.185$, $+0.148$), recovering the setting the product form gets backwards. The anchor is also not redundant with effect size: adding it to a logistic model already containing $|\beta_D|$ raises the model AUC by $+0.035$ to $+0.199$ (all likelihood-ratio $p<10^{-37}$). The discovery effect size, which by construction cannot carry cross-cohort information, favors the null in every evaluable setting ($-0.035$ to $-0.077$)---the expected behavior of a negative control, and what licenses reading the difference form's separation as genuine.

**The intervals must be much wider than they first appear.** Two sources of variability sit outside a single-split gene bootstrap. Co-expression first: resampling genes independently treats correlated observations as independent, and a cluster bootstrap over 200 co-expression clusters (median 26--29 genes) widens intervals by a median $2.2\times$ (range $1.6$--$2.7$). Split variability second, and it dominates: over $K=6$ splits per setting, **between-split variance is $77$--$84\%$ of the total**, combined intervals are $3.1$--$8.2\times$ wider (median $6.5\times$), and difference-form estimates range $+0.030$ to $+0.283$ across splits. Under that accounting the difference form stays positive in all three glioma pairs ($+0.063$, $+0.198$, $+0.170$, and in all 18 split-level estimates) but only CGGA-325 excludes zero ($+0.091$ to $+0.306$); array-301 is on the boundary ($-0.008$ to $+0.348$) and CGGA-693 spans it. The product form is conclusive only in CGGA-693, where it is inverted ($-0.183$ to $-0.011$). **We therefore claim a consistent direction, not three significant results.** With a few hundred patients per arm, a matched null cannot resolve differences of this size.

An easy error is worth naming. Our first analysis compared arms with a $t$-test across 8 resampling replicates. Replicates resample one fixed cohort, so their spread is Monte-Carlo error shrinking as $1/\sqrt{\text{replicates}}$---re-running with 3, 4, $\ldots$, 8 moves $p$ from $0.32$ to $0.010$ while the effect moves only $+0.030$ to $+0.037$. Such a $p$-value measures how long the simulation ran, not how strong the evidence is.

**Robust to thresholds, sensitive to the anchor.** Varying the discovery and replication thresholds over $\{0.01, 0.05, 0.10\}$ each gives 9 combinations per cell; **all 54 grade-anchor cells reproduce the sign reported here**. Substituting *IDH* for grade, only 39 of 54 agree and CGGA-325's product form inverts in all 9. The calibration is established for the grade anchor specifically; Table 3's *IDH* column is reported uncalibrated.

\begin{table*}[t]
\centering
\footnotesize
\caption{Matched-null calibration, on \textbf{Table 3's own gene universe} (four-cohort intersection, TCGA SD $>0.5$, capped at the 8,000 most variable). Within each row the two arms share a discovery set, a gene universe, and a replication-cohort size by construction; ``null'' is a held-out, disjoint subsample of the discovery cohort, in which no cross-cohort artifact can exist. Intervals are paired bootstraps over \emph{genes} (2,000 resamples), the sampling unit of an AUC; FDR is Benjamini--Hochberg within metric. A positive difference means the metric registers genuine cross-cohort discordance. The anchor \emph{product} is heavily confounded (null 0.79--0.90) and changes sign across settings; the \emph{difference} form is far less so (null 0.53--0.61) and its point estimate exceeds the null in all three glioma pairs, though once split variability is included only CGGA-325 remains conclusive (\S6.9). The effect-size control, which cannot carry cross-cohort information, is negative throughout, as a negative control must be. Splitting the discovery cohort halves its size, so absolute AUCs are not interchangeable with Table 3's; the informative comparison is within a row. Restricting to the 2,000 most variable genes moves point estimates but changes no glioma conclusion except CGGA-693's difference form, which is positive but interval-spanning zero there ($+0.051$; $-0.019$ to $+0.121$). The intervals shown are within-split only: they resample genes from one split. A co-expression cluster bootstrap widens them by a median $2.2\times$, and combining six independent splits widens them a further $3.1$–$8.2\times$, after which only CGGA-325's difference form excludes zero (§4.2). Quote the combined intervals, not these.}
\begin{tabular}{llccrlc}
\toprule
Setting (discovery $n$ / arm $n$ / genes) & Metric & null & real & diff & 95\% CI & FDR \\
\midrule
Glioma: TCGA$\to$CGGA-693  & product    & 0.900 & 0.824 & $-0.075$ & $(-0.093, -0.058)$ & 0.001 \\
(368 / 367 / 5,495)        & difference & 0.530 & 0.669 & $+0.138$ & $(+0.102, +0.176)$ & 0.001 \\
                           & effect size& 0.820 & 0.785 & $-0.035$ & $(-0.056, -0.013)$ & 0.002 \\
\midrule
Glioma: TCGA$\to$CGGA-325  & product    & 0.893 & 0.937 & $+0.044$ & $(+0.032, +0.055)$ & 0.001 \\
(513 / 222 / 6,037)        & difference & 0.611 & 0.796 & $+0.185$ & $(+0.159, +0.210)$ & 0.001 \\
                           & effect size& 0.844 & 0.767 & $-0.077$ & $(-0.093, -0.060)$ & 0.001 \\
\midrule
Glioma: TCGA$\to$array-301 & product    & 0.793 & 0.924 & $+0.132$ & $(+0.120, +0.144)$ & 0.001 \\
(487 / 248 / 6,381)        & difference & 0.613 & 0.760 & $+0.148$ & $(+0.129, +0.167)$ & 0.001 \\
                           & effect size& 0.768 & 0.733 & $-0.036$ & $(-0.051, -0.020)$ & 0.001 \\
\midrule
Breast: METABRIC$\to$BRCA  & product    & 0.716 & 0.566 & $-0.150$ & $(-0.207, -0.094)$ & 0.001 \\
(1,066 / 622 / 1,021)      & difference & 0.544 & 0.487 & $-0.057$ & $(-0.118, +0.004)$ & 0.062 \\
                           & effect size& 0.634 & 0.566 & $-0.068$ & $(-0.125, -0.011)$ & 0.019 \\
\midrule
\multicolumn{7}{l}{\footnotesize Breast: TCGA-BRCA$\to$METABRIC is \emph{not evaluable} under this design at either universe: splitting 622 patients with} \\
\multicolumn{7}{l}{\footnotesize 67 events leaves 26 prognostic genes at 2,000 genes and \textbf{none at all} at 8,000, too few for an AUC in either arm.} \\
\bottomrule
\end{tabular}
\end{table*}

### 4.3 Genome-scale evaluation, correctly read

Scored across 6,512 TCGA-prognostic genes, the grade anchor predicts replication with AUC $0.815$ (CGGA-693), $0.931$ (CGGA-325) and $0.924$ (array-301), and the *IDH* anchor $0.771$--$0.842$ (Table 3, Figure 1). Detectability runs at or near chance ($0.358$--$0.525$) and replication-cohort estimation precision likewise ($0.439$--$0.591$). These two do not, however, license the conclusion they appear to. Both measure *noise*, whereas the objection they were meant to answer is about *signal*: a gene with a strong anchor relationship has a large effect, and large effects replicate. The correct competitor is the discovery effect size, which reaches $0.738$--$0.783$, and the correct null is the same-population one of §4.2, not chance.

A gated screen over the 500 most-variable shared genes shows what the flag buys in practice. Among genes prognostic in TCGA, replication in CGGA-693 was $92.8\%$ ($89.9$--$95.0$) for QC-concordant genes versus $24.3\%$ ($11.8$--$41.2$) for flagged genes (Fisher $p=4.9\times10^{-21}$, OR $=39.2$), against a no-gating baseline of $87.3\%$; $88\%$ of the panel passed QC, and both sign-flips fell among flagged genes. The contrast reproduced in a second pair sharing no patients (TCGA$\to$CGGA-325: $95.4\%$ vs.\ $19.2\%$; Figure S7 plots the panel). Because survival in pooled glioma is grade-dominated, part of this reflects shared grade signal; restricting to associations prognostic *after* grade and *IDH* adjustment leaves the gate above chance but on a small set (183 genes; AUC $0.65$, $0.55$--$0.74$; Table S7).

\begin{table*}[t]
\centering
\footnotesize
\caption{Genome-scale evaluation: a gene's biological positive-control concordance ranks whether its TCGA prognostic association replicates in a separate cohort (AUC). \textbf{Chance is not the appropriate null for these AUCs}: §4.2 shows the same product-form score reaches $0.79$--$0.90$ against a same-population null in which no cross-cohort artifact can exist, and the discovery effect size alone reaches $0.74$--$0.78$. Neither of those figures may be subtracted from this table's, because the matched-null design halves the discovery cohort and so changes both the prognostic gene set and the AUC scale; the licensed comparison is within a row of Table 2, and what this table's numbers require is the \emph{qualitative} correction that a value near $0.9$ here is not evidence of detection. The detectability and precision columns measure noise rather than signal and are not the informative comparison. The three replication datasets come from one resource and are not fully independent (§4.3 quantifies and controls for patient overlap). $n=6,512$ TCGA-prognostic genes (FDR $<0.05$); AUCs are given to three decimals here and rounded to two in the text. Bootstrap 95\% CIs are tight given $n\approx6,500$ (CGGA-693 grade anchor: $0.815$, $0.81$--$0.83$); Table S7 gives the circularity control, Table 2 the matched-null calibration, and Table 1 the controlled experiment.}
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

**Figure 1.** Anchor scores for predicting replication of TCGA-prognostic genes in each CGGA cohort, by metric. The grade and IDH anchors exceed both the detectability and estimation-precision baselines and the dotted chance line at $0.5$—but chance is not the appropriate null, which is the point the figure exists to make. Red segments mark the same-population matched null for the product-form anchor (§4.2): $0.900$, $0.893$, $0.793$. **They are drawn as a reference level, not as a subtractable baseline**: the matched-null design splits the discovery cohort, so its absolute AUCs are on a different scale from the bars (Table 2). What the figure shows is that the null sits in the same range as the reported performance—not that any particular bar clears or fails it. The within-row contrasts of Table 2, and the combined intervals of §4.2, are where that question is settled.

### 4.4 Worked example: EGFR in glioma

In TCGA-LGG the EGFR association is strong unadjusted (HR $=1.59$ per SD) and after age and grade (HR $=1.31$; $p=8.2\times10^{-3}$) but is **abolished by *IDH* adjustment** (HR $=1.13$; $p=0.16$) and by subtype (HR $=1.14$; $p=0.10$). Refitting on a common *IDH*-complete subset ($n=452$) reproduces this, a likelihood-ratio test confirms EGFR adds nothing ($\chi^2=1.98$, $p=0.16$), and discrimination does not improve (C-index $0.834\to0.830$) (Table S10). The EGFR term violates PH in that model ($p=0.004$); period-specific estimates are null before three years (HR $0.96$) and non-significantly elevated after (HR $1.54$), so a late-emerging effect cannot be excluded. Conclusions are unchanged under DESeq2 VST and under an independent recount3 reprocessing [19] (Tables S3, S5).

A stronger test is whether EGFR carries signal *within* the IDH-wildtype subtype, where subtype confounding no longer applies. It does not (Table S12): among IDH-wildtype gliomas ($n=338$, 253 deaths) the unadjusted HR is $1.00$ ($0.89$--$1.13$), age-adjusted $0.93$, and $+$grade $0.90$. This is a well-powered null---with 253 events, 80% power extends to HR $\geq1.19$. Copy number agrees: high-level amplification is strongly enriched in IDH-wildtype disease ($55\%$ of 167 vs.\ $2.3\%$ of 222 IDH-mutant), confirming the copy-number data are faithful, yet within IDH-wildtype glioma amplification does not predict survival (continuous HR $1.00$ per SD, $0.84$--$1.20$; binary null across CN $\geq5/6/7$).

**The gate flags CGGA, correctly but for a specific reason.** CGGA's clinical controls are sound (IDH-mutant HR $0.22$--$0.32$; grade IV vs.\ II HR $7$--$9$), but EGFR's expression controls fail: its grade correlation is $-0.033$, $-0.097$ and $-0.086$ across the three datasets against $+0.179$ in TCGA (Fisher $r$-to-$z$ $p<10^{-3}$; Table S11, Figure S6). Composition explains part of this and should be stated precisely: CGGA's grade strata are less *IDH*-differentiated ($18\%$, $13\%$, $13\%$ IDH-mutant at grade IV vs.\ TCGA's $8\%$), holding *IDH* fixed the EGFR--grade slopes agree ($+0.124$ vs.\ $+0.115$), and within IDH-wildtype tumors the correlation is *positive* in all three CGGA datasets. Genome-wide, composition accounts for only $8$--$25\%$ of mean anchor discordance. What composition does not explain is the collapse in the share of EGFR's variance attributable to grade and *IDH*---from $0.106$ in TCGA to $0.045$, $0.030$ and $0.024$---with within-stratum SD unchanged ($1.88$ vs.\ $1.84$). Reweighting cannot produce this: the transport calculation predicts a *higher* correlation than TCGA's own ($+0.24$ vs.\ $+0.22$) against an observed $+0.035$. Additive noise would inflate within-stratum SD, which it does not. The pattern is signal compression---the mode §4.1 shows is detectable. EGFR is therefore a marker of the IDH-wildtype state rather than an independent prognostic factor, and CGGA cannot adjudicate that for this gene.

### 4.5 Scope: a well-matched second tumor type

In breast cancer (TCGA-BRCA$\leftrightarrow$METABRIC, ER and basal subtype as anchors) the gate has little to flag. The forward ER-anchor AUC is $0.60$ ($0.44$--$0.75$) on only 59 prognostic genes; the better-powered reverse direction gives $0.56$ ($0.52$--$0.59$), where the discovery effect-size baseline ($0.57$) in fact exceeds it, so we cannot claim the anchor adds anything in this pair. The basal anchor reaches at most $0.54$. Three non-exclusive reasons make breast diagnostic of *when* the method helps: the cohorts measure genes highly concordantly (genome-wide agreement of expression--ER correlations $r=0.84$), leaving little artifact to flag; in ER-positive breast almost any gene set is prognostic [25]; and the anchors differ from glioma's. Under a matched null the forward direction is not evaluable at all---splitting 622 patients with 67 events leaves 26 prognostic genes at a 2,000-gene universe and none at 8,000.


## 5. Discussion

The gate does what it was designed to do, and the evidence for that is the controlled experiment rather than any observational agreement: with corruption injected at known doses it identifies compressed genes at AUC $0.73$--$0.85$, permuted and noised genes at $0.63$--$0.86$, and floored genes not at all. That last result is as useful as the first---it tells a practitioner where not to rely on this check.

The methodological lesson is sharper than the one we set out to draw. The gate's apparent genome-scale performance survives comparison with the baselines we first chose but not with the ones we should have chosen: a same-population null reaches $0.79$--$0.90$ with the same score, and discovery effect size alone reaches $0.74$--$0.78$. The diagnosis is the product form's rewarding of anchor magnitude; the remedy is the difference form, whose null is $0.53$--$0.61$. How strong the remaining evidence is depends on how uncertainty is counted, and counting it properly is expensive: split variability supplies $77$--$84\%$ of the total variance, and once included only CGGA-325 excludes zero. What is general here is the *requirement*, not the verdict---split the discovery cohort, build a same-population replication set of matched size, report the AUC against that rather than $0.5$, and judge each setting on what that comparison shows.

The EGFR case is instructive because it initially appeared to contradict TCGA. A second objection---that CGGA simply contains a different patient mix---turns out to be partly right, and stating it precisely matters: composition explains a real share of the marginal grade discordance and, holding *IDH* fixed, the slopes agree. What it does not explain is the collapse in biologically-explained variance with within-stratum variance unchanged. Our recommendation to practitioners is correspondingly specific: run the gate in its difference form, calibrate against a split-half same-population null rather than chance, and do not rely on it where saturation is the suspected failure mode.


## 6. Limitations

The five that bear on the main claim; the remainder, including cohort-specific caveats, are in the Supplementary Discussion.

1. **Uncertainty is dominated by the split, not the genes.** Between-split variance is $77$--$84\%$ of the total, and under the combined accounting only one of three glioma pairs is conclusive (§4.2). This is a limitation of cohort size as much as of method: only a materially larger discovery cohort resolves it.
2. **The calibration is a grade-anchor result.** It is fully robust to the discovery and replication thresholds (54 of 54 grid cells) but not to the anchor: with *IDH*, 39 of 54 agree and CGGA-325's product form inverts. Table 3's *IDH* column is uncalibrated.
3. **The ground-truth experiment covers five modes, not all of them,** and its doses are not calibrated to any real assay. It establishes a blind spot for flooring that we have not characterized further, and the inference that CGGA's EGFR resembles the compression arm rests on a two-statistic signature match, not a demonstrated cause.
4. **Replication rests on four evaluable cohort pairs, three from one consortium,** in two tumor types, one of which shows nothing for any metric. The three CGGA datasets are not fully independent (61 patients shared between array-301 and CGGA-325, 6 with CGGA-693; de-duplication barely moves the discrimination, $0.91$ vs.\ $0.92$).
5. **The framework assesses cross-cohort *reproducibility*, not *validity*.** It can bless a reproducible association that is nonetheless confounded---as EGFR's is by *IDH*---so it complements rather than replaces causal scrutiny.


## 7. Conclusion

Cross-cohort prognostic replication in transcriptomics is unreliable, and a biological positive control does detect the subset of discordances that stem from faulty measurement---but establishing that required discarding the evidence we first assembled for it. Under injected ground truth the gate identifies compressed genes at AUC $0.73$--$0.85$ and floored genes not at all. Against a same-population null the score usually used for such gates reaches $0.79$--$0.90$ where nothing is wrong, so its distance from chance is uninformative; scoring disagreement directly drops that null to $0.53$--$0.61$ and is positive in all three glioma pairs, though with cohorts of this size only one is conclusive once split variability is counted. Run the gate in its difference form, calibrate it against a split-half null rather than against chance, and do not rely on it where saturation is the suspected failure. Applied to EGFR, it correctly flags a comparison that cannot adjudicate the gene: EGFR is a marker of the IDH-wildtype state, not an independent prognostic factor.


## Declarations

**Data availability.** All data analysed here are public and previously
published. TCGA RNA-sequencing, clinical and copy-number data were obtained
from the Genomic Data Commons (GDC Data Release 45.0, 2025-12-04) via
`TCGAbiolinks`; the three CGGA datasets (mRNAseq_693, mRNAseq_325,
mRNA-array_301, release 20200506) from http://www.cgga.org.cn; TCGA-BRCA and
METABRIC expression and subtype labels from the processed deposit of ref. 20,
with survival from cBioPortal (ref. 21); and the recount3 reprocessing from
ref. 19. No new data were generated.

**Code availability.** Analysis scripts, console logs and session information
are at https://github.com/edanmn/TCGA-EGFR-Glioma (branch
`revision/calibration-and-ground-truth`). Steps 23--40, which produce §4.1--4.3,
are on that branch only.

**Ethics.** This is a secondary analysis of de-identified, publicly released
human subjects data obtained under the data-use terms of each resource. No
institutional review board approval was required and no participants were
recruited. The methods are retrospective and are not validated for clinical use.

**Funding.** This research received no specific grant from any funding agency in
the public, commercial, or not-for-profit sectors.

**Competing interests.** The author declares no competing interests.

**Author contributions.** R.K. designed the study, implemented the pipeline and
the verification suites, performed all analyses, and wrote the manuscript.

**Corresponding author.** Rishik Kondadadi (konda052@umn.edu).


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
20. Lupat R, Loi S, Li J. Processed TCGA BRCA and METABRIC datasets used in the Moanna manuscript (v0.1.0). *Zenodo*. 2020. doi:10.5281/zenodo.4326602.
21. Cerami E, Gao J, Dogrusoz U, et al. The cBio Cancer Genomics Portal: an open platform for exploring multidimensional cancer genomics data. *Cancer Discov*. 2012;2(5):401–404.
22. Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature*. 2012;486(7403):346–352.
23. Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature*. 2012;490(7418):61–70.
24. Michiels S, Koscielny S, Hill C. Prediction of cancer outcome with microarrays: a multiple random validation strategy. *Lancet*. 2005;365(9458):488–492.
25. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol*. 2011;7(10):e1002240.
26. Dhawan A, Barberis A, Cheng WC, et al. Guidelines for using sigQC for systematic evaluation of gene signatures. *Nat Protoc*. 2019;14(5):1377–1400.
27. Gagnon-Bartsch JA, Speed TP. Using control genes to correct for unwanted variation in microarray data. *Biostatistics*. 2012;13(3):539–552.
28. Altmejd A, Dreber A, Forsell E, et al. Predicting the replicability of social science lab experiments. *PLoS One*. 2019;14(12):e0225826.
29. Zhao Z, Meng F, Wang W, et al. Comprehensive RNA-seq transcriptomic profiling in the malignant progression of gliomas. *Sci Data*. 2017;4:170024.
