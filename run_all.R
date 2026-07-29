# run_all.R -- run the whole pipeline end to end from the project root.
#   Rscript run_all.R
# Each step is cached: the download is skipped once data/*.rds exist.

# Make sure we're at the project root (this file's directory).
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(file_arg)) setwd(dirname(normalizePath(file_arg)))

source("R/00_setup.R")        # install/verify packages
source("R/01_download.R")     # download LGG + GBM (slow, one-time)
source("R/02_survival_km.R")  # Kaplan-Meier + log-rank
source("R/03_cox_model.R")    # Cox, adjusted for age + grade
source("R/04_pathway_screen.R") # pathway screen with FDR correction

cat("\nPipeline complete. See results/ and figures/.\n")
