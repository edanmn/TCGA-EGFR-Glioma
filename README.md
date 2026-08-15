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
| `R/23_adjusted_replication.R` | 2×2 decomposition of the circularity control (discovery set × replication criterion) | `results/adjusted_replication.csv` |
| `R/24_cgga_composition.R` | Cohort composition vs measurement as explanations for CGGA's EGFR failure | `results/cgga_composition.csv` |
| `R/25_controlled_shift.R` | Controlled experiment: confounding-structure shift and permutation, known doses | `results/controlled_shift.csv` |
| `R/26_transport_anchor.R` | Transport-adjusted anchor; decomposes discordance into composition vs residual | `results/transport_anchor.csv`, `results/pergene_cache.rds` |
| `R/27_effectsize_baseline.R` | The effect-size baseline missing from the glioma evaluation; EGFR variance decomposition | `results/effectsize_baseline.csv` |
| `R/28_matched_null.R` | **Superseded** (replicate t-tests; see `R/33`). Retained for provenance | `results/matched_null.csv` |
| `R/29_breast_verify.R` | Reproduces every §6.10 breast number in both directions | `results/breast_verify.csv` |
| `R/30_multisetting_null.R` | Matched null × 5 settings × 6 quality metrics, with rate-matching control | `results/multisetting_null.csv` |
| `R/31_corruption_modes.R` | Five measurement failure modes injected at known doses | `results/corruption_modes.csv` (**Table 6**) |
| `R/32_paper_audit.R` | **Verification:** asserts every revision figure against its source output | `results/paper_audit.txt` |
| `R/33_gene_level_inference.R` | Gene-level paired bootstrap (2,000-gene robustness run) | `results/gene_level_inference.csv` |
| `R/34_null_full_universe.R` | Matched null on Table 4's gene universe | `results/null_full_universe.csv` (**Table 5**) |
| `R/35_cluster_bootstrap.R` | Co-expression cluster bootstrap: are the gene-level intervals valid? | `results/cluster_bootstrap.csv` |
| `R/_helpers.R` | Shared: SummarizedExperiment → tidy survival table | (sourced) |

Scripts 21, 22 and 32 share no code with the analysis pipeline; they exist to catch
silent data-coercion errors (wrong transform, missing values coerced to a valid
level) that produce plausible-looking numbers, and to catch a manuscript figure
drifting from the code that produced it. All three should pass cleanly:
`28 passed, 0 failed`, `14 passed, 0 failed`, and `177 passed, 0 failed`.
`run_all.R` stops non-zero if any of them reports a failure.

**Scripts 23–35 were added in revision.** They exist because the original
genome-scale claim (Table 4) was benchmarked only against baselines that measure
noise. `R/30`/`R/34` calibrate it against a same-population null, `R/27` supplies
the effect-size baseline it lacked, `R/31` establishes what the gate detects under
injected ground truth, and `R/33`/`R/35` correct the inference (replicate t-tests →
gene bootstrap → co-expression cluster bootstrap). See `QA_REVIEW.md` for the full
audit trail, including the errors found in these scripts themselves.

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
