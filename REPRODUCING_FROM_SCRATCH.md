# Reproducing this analysis from an empty `data/`

**Written 2026-09-02 after deleting `data/` and rebuilding it from the GDC.**
Everything below was measured, not estimated.

The short version: **the manuscript's cohort definitions reproduce exactly; four
of its GBM regression coefficients do not**, and the reason is upstream of this
repository. Three separate gaps were found, one of which is now fixed.

## What reproduces from a clean checkout

| suite | checks | result |
|---|---:|---|
| `R/32_paper_audit.R` | 221 | **all pass** |
| `R/37_supp_audit.R` | 137 | **all pass** |
| `R/40_manuscript_coverage.R` | 70 | **all pass** |
| `R/21_assertions.R` | 20 of 28 | **all 20 pass**, then halts (see CGGA below) |
| `R/22_stat_recompute.R` | 6 of 14 | LGG exact; **8 GBM/pooled rows differ** |

Every cohort count is exact: LGG 511 patients / 125 deaths, GBM 282 / 227,
pooled 793 / 352, and all category and nested-model sample sizes. The cohorts are
not the problem.

## Gap 1 — copy number could never be rebuilt (FIXED)

`R/14_egfr_amplification.R` and `R/21_assertions.R` both read
`data/<PROJECT>_cnv.rds`. **No script in this repository created it.**
`R/01_download.R` fetches expression only. The original copy-number download was
run ad hoc and only its logs were committed (`results/_cnv_download.log`,
`results/_cnv_prepare.log`), so three assertions depended on a file the pipeline
could not produce. This stayed invisible for as long as the `.rds` sat on disk.

`R/01b_download_cnv.R` now closes it. The query is reconstructed from those two
logs, which record the genome build (hg38), the data type (Gene Level Copy
Number) and the assay names. The one genuinely ambiguous field was
`workflow.type`, resolved by counting files against the logged total:

| workflow.type | files for TCGA-GBM |
|---|---:|
| **ASCAT3** | **511** — matches the log exactly |
| ASCAT2 | 542 |
| AscatNGS | 12 |
| ABSOLUTE LiftOver | 572 |

Rebuilt, the copy-number assertions land at 166/55.4% and 221/2.3% against a
published 167/55.1% and 222/2.3%, and the survival subset reproduces exactly at
166/127. **One patient short in each group** — the ASCAT3 file set has changed
slightly since the original run.

## Gap 2 — GBM sample selection was never deterministic

`R/22_stat_recompute.R` and `R/21_assertions.R` both reduce to one row per
patient with `d[!duplicated(d$patient), ]`, which keeps whichever aliquot appears
**first in column order**. Column order comes from the GDC query response and is
not guaranteed stable between downloads.

| cohort | primary aliquots | patients | patients with >1 aliquot |
|---|---:|---:|---:|
| TCGA-LGG | 516 | 516 | **0** |
| TCGA-GBM | 372 | 284 | **86** |

LGG has no duplicates, so the operation is deterministic there — which is exactly
why all five LGG rows in `R/22` still reproduce to three decimals. GBM has 86,
so up to 86 patients can silently carry a different aliquot's expression than
they did originally: **right patient count, right event count, different values,
exit code 0.**

This has not been changed, because changing it would move the published numbers a
second way on top of Gap 3. It is recorded here as a known defect. If it is
fixed, the rule must be explicit (for example ordering by barcode before
deduplicating) and the manuscript's affected values must be regenerated together.

## Gap 3 — the expression data itself has been re-quantified (NOT FIXABLE HERE)

With cohorts identical and aliquot choice controlled, the GBM coefficients still
differ:

| model | rebuilt | published |
|---|---|---|
| GBM unadjusted | HR 0.95, p 0.411, C 0.530 | HR 0.96, p 0.552, C 0.523 |
| GBM + age, IDH | HR 0.85, p 0.015, C 0.642 | HR 0.87, p 0.030, C 0.639 |
| Pooled unadjusted | HR 1.32, p 1.11e-06, C 0.548 | HR 1.35, p 1.98e-07, C 0.553 |
| Pooled + age, grade, IDH | HR 0.91, p 0.045, C 0.838 | HR 0.92, p 0.078, C 0.837 |

Four deterministic aliquot-selection rules were tried (as-returned, barcode
ascending, barcode descending, patient+barcode). All four give n=282, events=227,
and HR 0.95 with C 0.526–0.530 — **none recovers the published 0.96 / 0.523.**
Same patients, same events, different expression values. GDC has re-processed
TCGA-GBM STAR quantification since the original download.

The effect is small and sits near the null (HR ≈ 0.95, C ≈ 0.53), so the
estimates barely move while the p-values swing. **No conclusion in the manuscript
changes**: the GBM association is null before and after. But the printed
coefficients cannot be reproduced from today's GDC.

`R/01_download.R` queries the GDC live and pins no data release. That is the root
cause and it is the same defect flagged for the docking project's ChEMBL and PDB
queries.

## Gap 4 — CGGA cannot be downloaded by any script here

Eight files under `data/cgga/` are read by 26 scripts:

```
CGGA.mRNA_array_301_gene_level.20200506.txt   cgga_clinical.tsv
CGGA.mRNAseq_325.RSEM-genes.20200506.txt      cgga_genes.tsv
CGGA.mRNAseq_693.RSEM-genes.20200506.txt      cgga325_clinical.tsv
cgga_array_clinical.tsv                       cgga325_genes.tsv
```

No script fetches them and **no URL is recorded anywhere in the repository**.
CGGA requires manual download from cgga.org.cn. `R/21_assertions.R` halts on the
first missing file, which is why it reports 20 checks rather than 28.

## What to do about it

1. **Pin the GDC data release.** Record the file UUIDs from the original
   download in a manifest and query by UUID, so the expression matrix is fixed.
   Without this, Gap 3 recurs on every future rebuild.
2. **Make the aliquot rule explicit** (Gap 2), and regenerate the affected
   manuscript values in the same pass.
3. **Commit a CGGA manifest** — at minimum the download URL, file names and
   checksums, since the data itself is licence-gated.
4. Until 1–3 are done, the honest reproduction claim is **428 of 470 checks from
   a clean checkout**, not 470.

## Note on the cohort table in `README.md`

It lists TCGA-GBM as `~170` samples. A fresh download returns 391 samples and
284 primary-tumour patients. The manuscript's analysis cohort of 282 is correct;
the README's raw-sample figure is stale.
