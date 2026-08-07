# TCGA Brain Cancer: EGFR Expression & Survival

An R + [TCGAbiolinks](https://bioconductor.org/packages/TCGAbiolinks) pipeline testing
whether **EGFR expression predicts overall survival** in TCGA brain-cancer cohorts.

**Question:** In lower-grade glioma (LGG), do patients with high EGFR expression
survive shorter than those with low expression — and does that difference disappear
in glioblastoma (GBM)?

The interesting result is the *contrast*: EGFR expression tracks the aggressive,
IDH-wildtype tumors in LGG, so a median split usually separates cleanly there,
while in GBM (where outcomes are uniformly poor) the same split tends to flatten.
The Cox step then asks the honest question — does EGFR still matter **after**
adjusting for age and tumor grade, or is it just a proxy for grade?

## Cohorts
| Study | Cancer | Grade | ~Samples (STAR counts) |
|-------|--------|-------|------------------------|
| `TCGA-LGG` | Lower Grade Glioma | II / III | ~500 |
| `TCGA-GBM` | Glioblastoma | IV | ~170 |
| `LGG+GBM` | Combined | II–IV | ~670 |

## Pipeline
| Script | What it does | Key outputs |
|--------|--------------|-------------|
| `R/00_setup.R` | Install/verify packages (auto-falls back to a Bioconductor mirror) | — |
| `R/01_download.R` | Download RNA-seq + clinical from the GDC, cache as `data/*.rds` | `data/TCGA-*_expr_se.rds` |
| `R/02_survival_km.R` | Kaplan–Meier + log-rank, EGFR median split, per cohort | `results/km_logrank_summary.csv` (Table S1) |
| `R/03_cox_model.R` | Cox models: EGFR alone, then adjusted for age + grade | `results/cox_EGFR_summary.csv` |
| `R/04_pathway_screen.R` | Screen 10 EGFR/RTK–PI3K genes, BH/FDR-corrected | `results/pathway_screen*.csv` (Table S2) |
| `R/05_adjusted_analysis.R` | Nested confounder-adjusted Cox models, C-index, PH tests | `results/cox_EGFR_adjusted.csv` (**Table 1**) |
| `R/06_vst_sensitivity.R` | Re-run key models on DESeq2 VST expression | `results/vst_sensitivity.csv` (Table S3) |
| `R/07_cgga_validation.R` | Naive CGGA external-validation models | `results/cgga_validation.csv` (Table S8) |
| `R/08_forest.R` | Exploratory forest plot — **not used in the manuscript** | `figures/forest_EGFR.png` (orphan) |
| `R/09_cgga_diagnostics.R` | Positive-control QC of CGGA (clinical + expression controls) | `results/cgga_positive_controls.csv` |
| `R/10_review_fixes.R` | Common-sample nested models, C-index CIs, LR test, Fisher r-to-z | `results/common_sample_cox.csv` |
| `R/11_recount3_validation.R` | Independent Monorail raw-read reprocessing check | `results/recount3_validation.csv` (Table S5) |
| `R/12_within_idhwt.R` | EGFR expression within IDH-wildtype glioma | `results/within_idhwt.csv` (**Table 3**) |
| `R/13_systematic_method.R` | 500-gene positive-control-gated screen | `results/systematic_screen.csv`, `figures/systematic_method.png` (**Figure 2**) |
| `R/14_egfr_amplification.R` | EGFR copy number vs survival within IDH-wildtype | `results/egfr_amplification.csv` |
| `R/15_revisions.R` | Amplification power/thresholds, replication CIs + Fisher, QC sensitivity | `results/revisions_summary.txt` |
| `R/16_method_framework.R` | Genome-scale framework: anchors vs detectability/precision baselines | `results/method_auc.csv` (**Table 4**), `figures/method_auc.png` (**Figure 3**) |
| `R/17_breast_generalization.R` | Scope test in TCGA-BRCA ↔ METABRIC | `results/breast_method_auc.csv` |
| `R/18_circularity_control.R` | Does the gate work beyond shared grade signal? | `results/circularity_control.csv` (Table S7) |
| `R/19_overlap_and_baselines.R` | CGGA patient overlap + de-duplicated replication; precision baseline; PH sensitivity | `results/review_fixes.txt`, `results/dedup_array301.csv` |
| `R/20_control_table.R` | Positive-control table by one method for all cohorts; ten-gene panel | `results/control_table.txt` (**Table 2**), `figures/positive_control_EGFR_grade.png` (**Figure 1**) |
| `R/21_assertions.R` | **Verification:** rebuilds every cohort independently, asserts 28 reported counts | `results/assertions.txt` |
| `R/22_stat_recompute.R` | **Verification:** refits every headline model independently, checks 14 statistics | `results/stat_recompute.txt` |
| `R/_helpers.R` | Shared: SummarizedExperiment → tidy survival table | (sourced) |

Scripts 21 and 22 share no code with the analysis pipeline; they exist to catch
silent data-coercion errors (wrong transform, missing values coerced to a valid
level) that produce plausible-looking numbers. Both should pass cleanly:
`28 passed, 0 failed` and `14 passed, 0 failed`.

## Run it
From this directory:

```bash
Rscript run_all.R
```

Or step by step:

```bash
Rscript R/00_setup.R
Rscript R/01_download.R
Rscript R/02_survival_km.R
Rscript R/03_cox_model.R
Rscript R/04_pathway_screen.R
```

The download step is slow (hundreds of MB per cohort) and runs once — later runs
reuse the cached `data/*.rds`. The GDC download also creates a `GDCdata/` folder here.

## How to read the results
- **`km_logrank_summary.csv`** — group sizes, median survival (days) per group, and
  the log-rank p (p < 0.05 = the High/Low curves differ significantly).
- **`cox_EGFR_summary.csv`** — hazard ratio for EGFR. HR is per 1-unit increase in
  `log2(TPM+1)`; **HR > 1 means higher EGFR → worse survival**. Compare the
  `univariate` row to the `adjusted(age+grade)` row: if the HR shrinks toward 1 and
  loses significance after adjustment, EGFR was largely a stand-in for grade.
- **`pathway_screen.csv`** — one age-adjusted Cox per gene; `p_BH` is the FDR-corrected
  p within each cohort. Significant = `p_BH < 0.05`.

## Notes / caveats
- Expression uses the **TPM** assay from STAR counts, `log2(TPM+1)`, primary tumors
  only (`sample_type == "01"`), one sample per patient.
- Overall survival = days_to_death (if dead) or days_to_last_follow_up (if alive).
- Grade parsing is best-effort from the clinical fields present; GBM is treated as
  grade IV. Inspect `df$grade` if a cohort's grade looks off.
- A significant association is **not** causation — it links EGFR level to outcome,
  nothing more.
