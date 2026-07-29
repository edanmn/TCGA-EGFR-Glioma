# 00_setup.R -- install and load everything the pipeline needs.
# Run this once. It is safe to re-run; installed packages are skipped.

repos <- "https://cloud.r-project.org"
options(repos = c(CRAN = repos))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = repos)

# The main Bioconductor host is occasionally unreachable (TLS connects but the
# HTTP layer stalls). Probe it; if it doesn't answer quickly, use a mirror.
bioc_reachable <- function(url) {
  ok <- tryCatch({
    con <- url(url, open = "rb"); on.exit(close(con))
    length(readBin(con, "raw", 1L)) > 0
  }, error = function(e) FALSE, warning = function(e) FALSE)
  isTRUE(ok)
}
ver <- as.character(BiocManager::version())
main_pkgs <- sprintf("https://bioconductor.org/packages/%s/bioc/src/contrib/PACKAGES", ver)
if (!bioc_reachable(main_pkgs)) {
  mirror <- "https://ftp.gwdg.de/pub/misc/bioconductor"
  message("Main Bioconductor host unreachable; using mirror ", mirror)
  options(BioC_mirror = mirror)
}

# Bioconductor packages
bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment")
need_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need_bioc))
  BiocManager::install(need_bioc, update = FALSE, ask = FALSE)

# CRAN packages
cran_pkgs <- c("survival", "survminer", "dplyr", "tidyr", "ggplot2", "readr")
need_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need_cran))
  install.packages(need_cran, repos = repos)

message("Setup complete. Installed/verified: ",
        paste(c(bioc_pkgs, cran_pkgs), collapse = ", "))
