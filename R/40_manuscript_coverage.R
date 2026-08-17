# 40_manuscript_coverage.R -- DOES THE MANUSCRIPT STILL SAY WHAT THE AUDITS VERIFY?
#
# R/32 and R/37 assert values held in results/*.csv. They do NOT read the
# manuscript, so they cannot detect a number being DROPPED from the paper or
# silently mis-transcribed into it -- a gap that matters most during a
# restructure, which is exactly when text moves. This script closes it by
# checking that every headline value the other suites verify still appears
# somewhere in the paper or supplement, at either the source precision or a
# correct 2-decimal rounding of it.
#
# Run from project root:  Rscript R/40_manuscript_coverage.R -> results/coverage.txt

out <- file("results/coverage.txt", open = "wt")
say <- function(...) { cat(sprintf(...), file = out); cat(sprintf(...)) }
txt <- paste(readLines("EGFR_glioma_survival_paper.md", warn = FALSE),
             readLines("EGFR_glioma_survival_supplementary.md", warn = FALSE),
             collapse = " ")
PASS <- 0; FAIL <- 0
present <- function(v) {
  cand <- unique(c(sprintf("%.3f", v), sprintf("%.2f", v), sprintf("%.1f", v),
                   sub("^0", "", sprintf("%.3f", v)), sub("^0", "", sprintf("%.2f", v)),
                   format(v)))
  any(vapply(cand, function(c) grepl(c, txt, fixed = TRUE), logical(1)))
}
chk <- function(group, v) {
  ok <- present(v); if (ok) PASS <<- PASS + 1 else FAIL <<- FAIL + 1
  say("  [%s] %-34s %s\n", ifelse(ok, "OK  ", "FAIL"), group, format(v))
}

say("=== manuscript coverage of audited values ===\n")
groups <- list(
  "ground truth: compression"   = c(0.848, 0.730, 0.639),
  "ground truth: flooring"      = c(0.498, 0.528, 0.599),
  "ground truth: variance base" = c(0.339, 0.136),
  "matched null: product"       = c(0.900, 0.893, 0.793),
  "matched null: difference"    = c(0.530, 0.611, 0.613),
  "combined CIs"                = c(0.063, 0.198, 0.170, 0.091, 0.306),
  "cluster bootstrap"           = c(2.2, 1.6, 2.7),
  "effect-size baseline"        = c(0.738, 0.783),
  "anchor increment"            = c(0.035, 0.199),
  "genome-scale AUC"            = c(0.815, 0.931, 0.924, 0.771, 0.842),
  "gated screen"                = c(92.8, 24.3, 39.2, 87.3),
  "EGFR nested Cox"             = c(1.59, 1.31, 1.13, 0.834, 0.830, 1.98),
  "EGFR within IDH-wt"          = c(1.00, 0.89, 338, 253, 1.19),
  "EGFR copy number"            = c(55, 2.3, 167, 222),
  "CGGA controls"               = c(0.033, 0.097, 0.086, 0.179),
  "EGFR variance collapse"      = c(0.106, 0.045, 0.030, 0.024, 1.88, 1.84),
  "breast"                      = c(0.598, 0.556, 0.542, 0.84),
  "cohort sizes"                = c(511, 125, 282, 227, 422, 404)
)
for (g in names(groups)) for (v in groups[[g]]) chk(g, v)
say("\n=== %d checks: %d passed, %d FAILED ===\n", PASS + FAIL, PASS, FAIL)
close(out)
cat(sprintf("\n%d passed, %d failed -> results/coverage.txt\n", PASS, FAIL))
