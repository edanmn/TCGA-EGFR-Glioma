# make_figures_png.R -- render manuscript figures as 300-dpi PNGs.
#   Fig 2-4: EGFR median-split KM (LGG, GBM, pooled)  [unadjusted, descriptive]
#   Fig 5:   LGG overall survival by IDH/1p19q molecular subtype (the confounder)
# Run from project root:  Rscript R/make_figures_png.R

suppressPackageStartupMessages({
  library(survival); library(survminer); library(dplyr); library(SummarizedExperiment)
})
source("R/_helpers.R")

GENE <- "EGFR"
dir.create("figures", showWarnings = FALSE)
YLAB <- "Overall survival probability"

save_km <- function(g, out, h = 7.5) {
  png(out, width = 7, height = h, units = "in", res = 300); print(g); dev.off()
  message("wrote ", out)
}

# ---- Fig 2-4: EGFR median split ----
lgg <- build_analysis("TCGA-LGG", GENE)
gbm <- build_analysis("TCGA-GBM", GENE)
combined <- bind_rows(lgg, gbm); combined$cohort <- "LGG+GBM"
for (nm in c("TCGA-LGG","TCGA-GBM","LGG+GBM")) {
  df <- add_median_split(list("TCGA-LGG"=lgg,"TCGA-GBM"=gbm,"LGG+GBM"=combined)[[nm]], GENE)
  df$grp <- df[[paste0(GENE, "_grp")]]
  fit <- survfit(Surv(time, event) ~ grp, data = df)
  g <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE,
                  legend.labs = c(paste0(GENE," Low"), paste0(GENE," High")),
                  palette = c("#2c7fb8","#d95f0e"),
                  xlab = "Days", ylab = YLAB,
                  title = paste0(nm, ": overall survival by ", GENE, " expression"))
  save_km(g, file.path("figures", paste0("KM_",GENE,"_",gsub("[^A-Za-z0-9]","",nm),".png")))
}

# ---- Fig 5: LGG survival by IDH/1p19q molecular subtype ----
se <- readRDS("data/TCGA-LGG_expr_se.rds")
cd <- as.data.frame(colData(se))
sub <- trimws(as.character(cd$paper_IDH.codel.subtype)); sub[sub %in% c("","NA")] <- NA
sub <- gsub("IDHmut-codel.*","IDHmut-codel", sub)
sub <- gsub("IDHmut-non-codel.*","IDHmut-non-codel", sub)
sub <- gsub("IDHwt.*","IDHwt", sub)
smap <- data.frame(barcode = rownames(cd), subtype = sub, stringsAsFactors = FALSE)
ds <- lgg %>% left_join(smap, by = "barcode") %>% filter(!is.na(subtype))
ds$subtype <- factor(ds$subtype, levels = c("IDHmut-codel","IDHmut-non-codel","IDHwt"))
fit <- survfit(Surv(time, event) ~ subtype, data = ds)
g <- ggsurvplot(fit, data = ds, pval = TRUE, risk.table = TRUE,
                legend.labs = levels(ds$subtype),
                palette = c("#1b9e77","#7570b3","#d95f02"),
                xlab = "Days", ylab = YLAB,
                title = "TCGA-LGG survival by IDH/1p19q subtype")
save_km(g, "figures/KM_LGG_subtype.png", h = 8)
