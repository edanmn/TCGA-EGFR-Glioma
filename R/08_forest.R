# 08_forest.R -- forest plot of EGFR HR (per SD) across TCGA and CGGA
# lower-grade cohorts under increasing adjustment, visualizing the discordance.
# Run from project root:  Rscript R/08_forest.R  -> figures/forest_EGFR.png

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

tc <- read.csv("results/cox_EGFR_adjusted.csv", stringsAsFactors=FALSE) %>%
  filter(cohort=="TCGA-LGG", model %in% c("1. unadjusted","3. + age + grade","4. + age + grade + IDH")) %>%
  mutate(dataset="TCGA-LGG",
         adj=recode(model, "1. unadjusted"="unadjusted",
                    "3. + age + grade"="+ age + grade",
                    "4. + age + grade + IDH"="+ age + grade + IDH"))
cg <- read.csv("results/cgga_validation.csv", stringsAsFactors=FALSE) %>%
  filter(cohort=="CGGA LGr(II/III)", model %in% c("1. unadjusted","2. + age + grade","3. + age + grade + IDH")) %>%
  mutate(dataset="CGGA (II/III)",
         adj=recode(model, "1. unadjusted"="unadjusted",
                    "2. + age + grade"="+ age + grade",
                    "3. + age + grade + IDH"="+ age + grade + IDH"))
d <- bind_rows(tc, cg)
d$adj <- factor(d$adj, levels=c("unadjusted","+ age + grade","+ age + grade + IDH"))
d$dataset <- factor(d$dataset, levels=c("TCGA-LGG","CGGA (II/III)"))

p <- ggplot(d, aes(x=HR, y=adj, colour=dataset)) +
  geom_vline(xintercept=1, linetype=2, colour="grey50") +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0.18,
                 position=position_dodge(width=0.5), linewidth=0.7) +
  geom_point(position=position_dodge(width=0.5), size=2.6) +
  scale_colour_manual(values=c("TCGA-LGG"="#1f78b4","CGGA (II/III)"="#e31a1c")) +
  scale_x_log10() +
  labs(x="EGFR hazard ratio per 1 SD (log scale, 95% CI)", y=NULL, colour="Cohort",
       title="EGFR prognostic effect in lower-grade glioma:\nTCGA vs CGGA under increasing adjustment") +
  theme_bw(base_size=12) + theme(legend.position="top",
       plot.title=element_text(size=12, face="bold"))
ggsave("figures/forest_EGFR.png", p, width=7, height=4.2, dpi=300)
message("wrote figures/forest_EGFR.png")
