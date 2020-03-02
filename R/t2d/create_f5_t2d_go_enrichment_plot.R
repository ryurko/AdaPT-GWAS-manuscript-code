# PURPOSE: Create Figure 5 displaying the T2D GO enrichment results

# AUTHOR: Ron Yurko

library(tidyverse)
library(ggrepel)
library(latex2exp)

# Load in the table of the GO results
t2d_go_results <- read_tsv("data/t2d/t2d_adjpval_low_s0_go_results.txt")
# Modified the file deleting all lines up until the column headers

# Generate the chart with the GO-term enrichment and also the fold enrichment values
t2d_go_enrich_plot <- t2d_go_results %>%
  filter(`upload_1 (over/under)` == "+") %>%
  arrange(desc(`upload_1 (fold Enrichment)`)) %>%
  dplyr::slice(1:10) %>%
  mutate(neglog_fdr = -log10(`upload_1 (FDR)`),
         # Remove the GO part:
         go_label = str_remove(`GO biological process complete`,
                               "\\(GO:[:digit:]{7}\\)"),
         go_label = str_remove(go_label, " process"),
         go_label = fct_reorder(go_label, neglog_fdr, .desc = TRUE)) %>%
  ggplot() +
  geom_bar(aes(x = go_label,
               y = neglog_fdr,
               fill = `upload_1 (fold Enrichment)`),
           width = 0.2, stat = "identity") +
  geom_point(aes(x = go_label,
                 y = neglog_fdr,
                 color = `upload_1 (fold Enrichment)`),
             size = 4) +
  labs(x = "Biological process",
       y = TeX('GO-term enrichment -log_{10}(FDR q-value)'),
       fill = "Fold enrichment",
       color = "Fold enrichment") +
  coord_flip() +
  scale_fill_gradient(low = "darkblue", high = "darkorange") +
  scale_color_gradient(low = "darkblue", high = "darkorange") +
  theme_bw() +
  theme(legend.position = c(.7,.65))

save_plot("figures/f5_t2d_go_enrichment_wide.pdf",
          t2d_go_enrich_plot + theme(legend.direction = "horizontal",
                                     legend.position = c(.65,.75),
                                     legend.text = element_text(size = 8),
                                     axis.text = element_text(size = 8),
                                     axis.title = element_text(size = 10),
                                     legend.title = element_text(size = 8)), ncol = 1, nrow = 1,
          base_width = 6, base_height = 2)
save_plot("nonpdf_figures/f5_t2d_go_enrichment_wide.jpg",
          t2d_go_enrich_plot + theme(legend.direction = "horizontal",
                                     legend.position = c(.65,.75),
                                     legend.text = element_text(size = 8),
                                     axis.text = element_text(size = 8),
                                     axis.title = element_text(size = 10),
                                     legend.title = element_text(size = 8)), ncol = 1, nrow = 1,
          base_width = 6, base_height = 2)
