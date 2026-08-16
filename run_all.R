# run_all.R -- run the whole pipeline end to end from the project root.
#   Rscript run_all.R
# Each step is cached: the download is skipped once data/*.rds exist.
#
# The final two steps are VERIFICATION SUITES. They share no code with the
# analysis pipeline and re-derive every reported sample size and headline
# statistic from the cached objects. If either reports a failure this script
# stops with a non-zero exit: the manuscript and the code have diverged.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(file_arg)) setwd(dirname(normalizePath(file_arg)))

steps <- c(
  "R/00_setup.R",                 # install/verify packages
  "R/01_download.R",              # download LGG + GBM (slow, one-time)
  "R/02_survival_km.R",           # Kaplan-Meier + log-rank            -> Table S1
  "R/03_cox_model.R",             # Cox, adjusted for age + grade
  "R/04_pathway_screen.R",        # pathway screen with FDR correction -> Table S2
  "R/05_adjusted_analysis.R",     # nested confounder-adjusted models  -> Table 1
  "R/06_vst_sensitivity.R",       # VST normalization sensitivity      -> Table S3
  "R/07_cgga_validation.R",       # naive CGGA external validation     -> Table S8
  "R/08_forest.R",                # exploratory forest plot (not used in the paper)
  "R/09_cgga_diagnostics.R",      # CGGA positive-control QC
  "R/10_review_fixes.R",          # common-sample models, C-index CIs, LR test
  "R/11_recount3_validation.R",   # independent raw-read reprocessing  -> Table S5
  "R/12_within_idhwt.R",          # EGFR within IDH-wildtype           -> Table 3
  "R/13_systematic_method.R",     # 500-gene gated screen              -> Figure 2
  "R/14_egfr_amplification.R",    # copy number vs survival
  "R/15_revisions.R",             # power, thresholds, replication CIs
  "R/16_method_framework.R",      # genome-scale framework             -> Table 4, Figure 3
  "R/17_breast_generalization.R", # scope test in breast
  "R/18_circularity_control.R",   # circularity control                -> Table S7
  "R/19_overlap_and_baselines.R", # CGGA overlap, de-duplication, precision baseline
  "R/20_control_table.R",         # positive-control table             -> Table 2, Figure 1
  "R/23_adjusted_replication.R",  # 2x2 circularity decomposition (slow, ~4 min)
  "R/24_cgga_composition.R",      # composition vs measurement for the CGGA failure
  "R/25_controlled_shift.R",      # composition sweep + permutation ground truth
  "R/26_transport_anchor.R",      # transport-adjusted anchor; writes pergene_cache.rds
  "R/27_effectsize_baseline.R",   # the missing effect-size baseline    -> §6.9
  "R/28_matched_null.R",          # matched same-population null        -> §6.9
  "R/29_breast_verify.R",         # reproduces every §6.10 number
  "R/30_multisetting_null.R",     # matched null x 5 settings x 6 metrics -> Table 5
  "R/31_corruption_modes.R",      # five failure modes, ground truth      -> Table 6
  "R/33_gene_level_inference.R",  # gene-level bootstrap, 2,000-gene robustness run
  "R/34_null_full_universe.R",    # matched null on Table 4 universe      -> Table 5
  "R/35_cluster_bootstrap.R",     # co-expression-aware interval check
  "R/36_threshold_sensitivity.R", # FDR x replication-p x anchor grid
  "R/38_split_variance.R",        # split + gene variance combination     -> §6.9
  "R/39_refresh_figures.R"        # redraw Figures 1-3 to match the text
)

for (s in steps) {
  cat(sprintf("\n=== %s ===\n", s))
  source(s)
}

## ---------------- verification gate ----------------
verify <- function(script, results_file, label) {
  cat(sprintf("\n=== %s (verification) ===\n", script))
  source(script)
  txt <- paste(readLines(results_file, warn = FALSE), collapse = "\n")
  hit <- regmatches(txt, regexpr("[0-9]+(?=\\s+FAILED)", txt, perl = TRUE, ignore.case = TRUE))
  if (!length(hit))
    stop(sprintf("%s: could not parse a pass/fail count from %s", label, results_file), call. = FALSE)
  failed <- as.integer(hit[1])
  if (!is.na(failed) && failed > 0)
    stop(sprintf("%s: %d check(s) FAILED -- see %s. The manuscript and the code have diverged.",
                 label, failed, results_file), call. = FALSE)
  cat(sprintf("%s: all checks passed.\n", label))
}
verify("R/21_assertions.R",     "results/assertions.txt",     "Cohort assertions")
verify("R/22_stat_recompute.R", "results/stat_recompute.txt", "Statistical recomputation")
verify("R/32_paper_audit.R",    "results/paper_audit.txt",    "Revision-figure audit")
verify("R/37_supp_audit.R",     "results/supp_audit.txt",     "Supplementary-table audit")

cat("\nPipeline complete, verification passed. See results/ and figures/.\n")
