# PURPOSE: Create Figure 4 in the manuscript, displaying the relationship between
#          the AdaPT q-values in 2014 with the new 2018 p-values

# Author: Ron Yurko

# Load necessary packages:
library(tidyverse)
library(cowplot)
library(latex2exp)

# ------------------------------------------------------------------------------

# Load the SCZ BrainVar data:
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Load the final model results:
scz_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")

scz_brainvar_adapt_qvals_new_pvals <- 
  tibble(adapt_qvals = scz_with_bd_z_eqtl_slopes_wgcna$qvals,
         updated_scz_pvals = bip_scz_brainvar_data$new_p_scz_18,
         adapt_disc = as.numeric(scz_with_bd_z_eqtl_slopes_wgcna$qvals <= .05)) %>%
  mutate(
    # Impute 1 for non-finite values since they were never in the rejection set:
    adapt_qvals = ifelse(is.infinite(adapt_qvals), 1, adapt_qvals)) %>%
  ggplot(aes(x = -log10(adapt_qvals),
             y = -log10(updated_scz_pvals))) +
  geom_point(alpha = 0.35, color = "gray") + #, aes(color = as.factor(adapt_disc))) +
  geom_segment(x = -log10(.05), xend = -log10(.05),
               y = 0, yend = Inf, linetype = "dotted", color = "darkblue") +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "darkred") +
  annotate("rect", xmin = -log10(.05), xmax = Inf,
           ymin = -log10(.05), ymax = Inf, alpha = 0.15,
           fill = "darkblue") +
  annotate("rect", xmin = -log10(.05), xmax = Inf,
           ymax = -log10(.05), ymin = 0, alpha = 0.15,
           fill = "darkred") +
  annotate("label", x = 1.75, y = 12.5,
           label = TeX('AdaPT discoveries at $\\alpha = 0.05$'),
           color = "darkblue", fill = "white", alpha = 0.75) +
  annotate("label", x = 1.75, y = -log10(.25),
           label = TeX("Replication threshold: -log$_{10}(0.05)$"),
           color = "darkred", alpha = 0.75, fill = "white") +
  annotate("label", x = 1.15, y = -log10(10^(-5)),
           label = "Smooth relationship between new 2018 SCZ p-values and AdaPT q-values",
           color = "black", fill = "white", alpha = 0.75) +
  geom_smooth(color = "black", se = FALSE) +
  theme_bw() +
  labs(x = TeX('-log_{10}(AdaPT q-values from 2014 studies)'),
       y = TeX('-log_{10}(SCZ p-values from new 2018 studies)'))


save_plot(scz_brainvar_adapt_qvals_new_pvals,
          filename = "figures/f4_scz_brainvar_replication.pdf",
          base_aspect_ratio = 2.5)
save_plot(scz_brainvar_adapt_qvals_new_pvals,
          filename = "nonpdf_figures/f4_scz_brainvar_replication.jpg",
          base_aspect_ratio = 2.5)



