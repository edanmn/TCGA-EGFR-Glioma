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
  "R/20_control_table.R"          # positive-control table             -> Table 2, Figure 1
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

cat("\nPipeline complete, verification passed. See results/ and figures/.\n")
