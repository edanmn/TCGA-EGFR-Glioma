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
| Script | What it does | Outputs |
|--------|--------------|---------|
| `R/00_setup.R` | Install/verify packages (auto-falls back to a Bioconductor mirror) | — |
| `R/01_download.R` | Download RNA-seq + clinical from the GDC, cache as `data/*.rds` | `data/TCGA-*_expr_se.rds` |
| `R/02_survival_km.R` | Kaplan–Meier + log-rank, EGFR median split, per cohort | `figures/KM_EGFR_*.pdf`, `results/km_logrank_summary.csv` |
| `R/03_cox_model.R` | Cox models: EGFR alone, then adjusted for age + grade | `results/cox_EGFR_summary.csv` |
| `R/04_pathway_screen.R` | Screen EGFR pathway genes, BH/FDR-corrected | `results/pathway_screen.csv` |
| `R/_helpers.R` | Shared: SummarizedExperiment → tidy survival table | (sourced) |

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
