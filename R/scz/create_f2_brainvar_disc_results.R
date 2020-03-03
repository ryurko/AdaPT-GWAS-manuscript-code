# PURPOSE: Create Figure 2 in the manuscript displaying a comparison of the
#          number of discoveries by the different sets of variables for the 
#          AdaPT BrainVar results

# Author: Ron Yurko

# Load necessary packages:
library(tidyverse)
library(data.table)
library(cowplot)
library(latex2exp)
library(UpSetR)
library(ggupset)

# ------------------------------------------------------------------------------

# Load the full set of original BD and SCZ SNPs that were already matched that
# do not include any pre-processing on the BD 2018 studies:
bip_scz_data_14_18 <- fread("data/bip_schz_data/bip_scz_data_14_18_snps.csv")
bip_scz_data_14_18_filtered <- bip_scz_data_14_18[bip_14_INFO > 0.6 & 
                                                    scz_14_INFO > 0.6,]
rm(bip_scz_data_14_18)

# Next the GTEx eSNPs (check before that the eSNPs are the exact same without 
# doing any pre-processing to account for the BD 2018 data which is good)
bip_scz_gtex_esnps <- fread("data/bip_schz_data/bip_scz_data_14_18_gtex_eqtls_cortical_wgcna.csv")

# Next the BrainVar eSNPs:
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Next proceed to generate the comparison of SCZ enrichment qq-plots using 
# random sampling of the higher p-values after quantiles have been stored for
# the sake of memory when saving the PDF version of the figure:
scz_enrich_plot <- bip_scz_data_14_18_filtered %>%
  as.data.frame() %>%
  dplyr::select(scz_14_P) %>%
  mutate(type = "all") %>%
  bind_rows({(
    bip_scz_gtex_esnps %>%
      as.data.frame() %>%
      dplyr::select(scz_14_P) %>%
      mutate(type = "gtex") 
  )},
  {(
    bip_scz_brainvar_data %>%
      as.data.frame() %>%
      dplyr::select(scz_14_P) %>%
      mutate(type = "brainvar") 
  )}) %>%
  group_by(type) %>%
  arrange(scz_14_P) %>%
  mutate(neglog_p = -log10(scz_14_P),
         exp_neglog_p = -log10(ppoints(n())),
         subset_sample = ifelse(ppoints(n()) > .02,
                                rbernoulli(n(), p = 0.2),
                                1)) %>%
  ungroup() %>%
  filter(subset_sample == 1) %>%
  mutate(type = fct_relevel(type, "brainvar", "gtex", "all")) %>%
  ggplot(aes(x = exp_neglog_p, y = neglog_p,
             color = type)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("darkorange",
                                "darkred",
                                "darkblue"),
                     labels = c("BrainVar eSNPs", 
                                "GTEx eSNPs",
                                "All SNPs")) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5,
              linetype = "dotted", color = "black") +
  theme_bw() +
  theme(legend.direction = "vertical",
        legend.position = c(.75, .25),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = TeX('Expected -log$_{10}$(SCZ p-value)'),
       y = TeX('Observed -log$_{10}$(SCZ p-value)'),
       shape = "SNP type", color = "SNP type") +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Remove these datasets:
rm(bip_scz_data_14_18_filtered)
rm(bip_scz_gtex_esnps)

# ------------------------------------------------------------------------------

# Now will generate the various figures displaying the results for the discoveries
# by the different set of variables - and options

# Load each of the separate model results:
scz_with_bd_z_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_only_s05_2cv.rds")
scz_with_bd_z_eqtl_slopes <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_s05_2cv.rds")
scz_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")
scz_with_wgcna_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_wgcna_only_s05_2cv.rds")
scz_intercept_only <- 
  readRDS("data/bip_schz_data/brainvar_results/scz_intercept_only_s05.rds")
scz_with_no_int <-
  readRDS("data/bip_schz_data/brainvar_results/scz_all_vars_no_int_s05.rds")

# Now create the list of the alpha = 0.05 results:
scz_comp_disc_list <- list("Intercept-only" = which(scz_intercept_only$qvals <= 0.05),
                           "BD z-stats" = which(scz_with_bd_z_only$qvals <= 0.05),
                           "BD z-stats + eQTL slopes" = which(scz_with_bd_z_eqtl_slopes$qvals <= 0.05),
                           "BD z-stats + eQTL slopes + WGCNA (w/ interactions)" = which(scz_with_bd_z_eqtl_slopes_wgcna$qvals <= 0.05),
                           "BD z-stats + eQTL slopes + WGCNA (w/o interactions)" = which(scz_with_no_int$qvals <= 0.05),
                           "WGCNA" = which(scz_with_wgcna_only$qvals <= .05))

# Convert this into the 0-1 membership matrix using the UpSetR package and then
# get it in the data format for ggupset instead:
tidy_scz_comp_disc <- scz_comp_disc_list %>%
  fromList() %>%
  as_tibble() %>%
  mutate(snp = 1:n()) %>%
  gather(disc_set, member, -snp) %>%
  filter(member == 1) %>%
  select(-member) %>%
  group_by(snp) %>%
  summarize(disc_sets = list(disc_set))

# Next create the ggupset plot
scz_disc_upset <- tidy_scz_comp_disc %>%
  ggplot(aes(x = disc_sets)) +
  geom_bar(fill = "darkblue") +
  labs(y = "Intersection size") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_upset(expand = c(0, 0)) +
  theme_combmatrix(combmatrix.panel.point.color.fill = "darkblue",
                   combmatrix.label.make_space = FALSE,
                   combmatrix.panel.line.size = .5,
                   combmatrix.label.text = element_blank(),
                   combmatrix.panel.margin = unit(c(0, 0), "pt"),
                   combmatrix.label.height = unit(95, "pt"))

# Shrink the version of the enrichment plot:
small_scz_enrich_plot <- scz_enrich_plot + 
  theme(legend.direction = "vertical",
        legend.position = c(.75, .25),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  labs(x = TeX('Expected -log$_{10}$(SCZ p-value)'),
       y = TeX('Observed -log$_{10}$(SCZ p-value)'),
       shape = "SNP type", color = "SNP type") +
  guides(color = guide_legend(override.aes = list(size = 4)))

# Combine the enrichment plot with the upset plot:
scz_disc_upset_plus <-
  ggdraw() +
  draw_plot(scz_disc_upset) +
  draw_plot(plot_grid(small_scz_enrich_plot, 
                      labels = "A",
                      label_fontface = "plain", vjust = 1), 
            x = 0.35, y = .375, width = .6, height = .6)

# Create a chart that just counts how many discoveries each type has:
scz_count_disc_chart <- scz_comp_disc_list %>%
  fromList() %>%
  as_tibble() %>%
  mutate(snp = 1:n()) %>%
  gather(disc_set, member, -snp) %>%
  filter(member == 1) %>%
  group_by(disc_set) %>%
  summarize(n_disc = n()) %>%
  mutate(disc_set = fct_reorder(disc_set, n_disc)) %>%
  ggplot(aes(x = disc_set, y = n_disc)) +
  geom_bar(stat = "identity", fill = "darkblue",
           width = .75) +
  geom_text(aes(x = disc_set, y = n_disc / 2,
                label = n_disc),
            color = "white", size = 4) +
  scale_x_discrete(position = "top") +
  scale_y_reverse(position = "right") +
  labs(y = "Number of discoveries") +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 0.5, size = 10, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())



# ------------------------------------------------------------------------------
# Next generate the Manhattan plots of the intercept-only and all variables:
# First create the chromosome axis data for both Manhattan plots:
brainvar_manhattan_axis_df <- bip_scz_brainvar_data %>% 
  group_by(bip_18_CHR) %>% 
  summarize(center = (max(bip_18_BP_cum) + min(bip_18_BP_cum)) / 2)

# Create columns that denote which tests are discoveries for the two different
# versions of results at target alpha = 0.05:
bip_scz_brainvar_data <- bip_scz_brainvar_data %>%
  mutate(int_only_disc = as.numeric(scz_intercept_only$qvals <= 0.05),
         boosted_disc = as.numeric(scz_with_bd_z_eqtl_slopes_wgcna$qvals <= 0.05))

# Now will also make a column for random sampling of the tests with p-values greater
# than 0.1 for saving memory with the PDF images:
bip_scz_brainvar_data <- bip_scz_brainvar_data %>%
  mutate(manhattan_sample = ifelse(scz_14_P <= 0.1,
                                   TRUE,
                                   rbernoulli(n(), p = 0.5)))

# Create the intercept-only Manhattan plot:
int_only_manhattan <- bip_scz_brainvar_data %>%
  filter(manhattan_sample) %>%
  ggplot(aes(x = bip_18_BP_cum, y = -log10(scz_14_P))) +
  # Show all points
  geom_point(aes(color = as.factor(bip_18_CHR)), alpha = 0.5, size = 0.5) +
  geom_point(data = filter(bip_scz_brainvar_data, int_only_disc == 1),
             aes(x = bip_18_BP_cum, y = -log10(scz_14_P)),
             color = "darkorange", alpha = 0.75, size = 0.75) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  # custom x-axis to handle the chromosome labels:
  scale_x_continuous(label = brainvar_manhattan_axis_df$bip_18_CHR, 
                     breaks = brainvar_manhattan_axis_df$center) +
  labs(x = "Chromosome", y = TeX('-log$_{10}$(SCZ p-value)'),
       title = TeX('Manhattan plot of AdaPT intercept-only discoveries for $\\alpha = 0.05$')) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 6))
rm(scz_intercept_only)

# Next using all covariates via gradient boosted trees:
boosted_manhattan <- bip_scz_brainvar_data %>%
  filter(manhattan_sample) %>%
  ggplot(aes(x = bip_18_BP_cum, y = -log10(scz_14_P))) +
  # Show all points
  geom_point(aes(color = as.factor(bip_18_CHR)), alpha = 0.5, size = 0.5) +
  geom_point(data = filter(bip_scz_brainvar_data, boosted_disc == 1),
             aes(x = bip_18_BP_cum, y = -log10(scz_14_P)),
             color = "darkorange", alpha = 0.75, size = 0.75) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  # custom x-axis to handle the chromosome labels:
  scale_x_continuous(label = brainvar_manhattan_axis_df$bip_18_CHR, 
                     breaks = brainvar_manhattan_axis_df$center) +
  labs(x = "Chromosome", y = TeX('-log$_{10}$(SCZ p-value)'),
       title = TeX('Manhattan plot of AdaPT gradient boosted discoveries for $\\alpha = 0.05$')) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 6))

# Now arrange these in a vertical grid without duplicate labels for the
# chromosome axes:
scz_manhattan_plots <-
  plot_grid(int_only_manhattan,
            boosted_manhattan,
            ncol = 1, labels = c("B", "C"), 
            label_fontface = "plain", align = "hv", vjust = 1.2)

# Finally combine everything together:
scz_brainvar_disc_plot_plus <- 
  plot_grid(plot_grid(scz_manhattan_plots,
                      scz_count_disc_chart,
                      ncol = 1, rel_widths = c(1, 1),
                      rel_heights = c(4, 1.75), labels = c("", "D"),
                      label_fontface = "plain"),
            scz_disc_upset_plus, ncol = 2, rel_widths = c(1, 1),
            rel_heights = c(1, 1), labels = c("", "E"),
            label_fontface = "plain")

# Save the figure:
# save_plot("figures/f2_scz_brainvar_results.pdf",
#           scz_brainvar_disc_plot_plus, ncol = 3, nrow = 2,
#           base_width = 5, base_height = 3)
# save_plot("nonpdf_figures/f2_brainvar_results.jpg",
#           scz_brainvar_disc_plot_plus, ncol = 3, nrow = 2,
#           base_width = 5, base_height = 3)
