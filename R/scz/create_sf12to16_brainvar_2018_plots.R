# PURPOSE: Generate the supplementary figures for the SCZ results using all 2018 studies

# AUTHOR: Ron Yurko

# Load necessary packages:
library(tidyverse)
library(data.table)
library(cowplot)
library(latex2exp)
library(xgboost)
library(pdp)

# ------------------------------------------------------------------------------

# Load the SCZ BrainVar data:
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Load each of the separate model results:
scz18_with_bd_z_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz18_with_bd_z_only_s05_2cv.rds")
scz18_with_bd_z_eqtl_slopes <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz18_with_bd_z_eqtl_slopes_s05_2cv.rds")
scz18_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz18_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")
scz18_with_wgcna_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz18_with_wgcna_only_s05_2cv.rds")
scz18_intercept_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz18_intercept_only_s05.rds")
scz18_with_no_int <-
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz18_all_vars_no_int_s05.rds")


# Now create the list of the alpha = 0.05 results:
scz18_comp_disc_list <- list("Intercept-only" = which(scz18_intercept_only$qvals <= 0.05),
                             "BD z-stats" = which(scz18_with_bd_z_only$qvals <= 0.05),
                             "BD z-stats + eQTL slopes" = which(scz18_with_bd_z_eqtl_slopes$qvals <= 0.05),
                             "BD z-stats + eQTL slopes + WGCNA (w/ interactions)" = which(scz18_with_bd_z_eqtl_slopes_wgcna$qvals <= 0.05),
                             "BD z-stats + eQTL slopes + WGCNA (w/o interactions)" = which(scz18_with_no_int$qvals <= 0.05),
                             "WGCNA" = which(scz18_with_wgcna_only$qvals <= .05))

# Convert this into the 0-1 membership matrix using the UpSetR package and then
# get it in the data format for ggupset instead:
tidy_scz18_comp_disc <- scz18_comp_disc_list %>%
  fromList() %>%
  as_tibble() %>%
  mutate(snp = 1:n()) %>%
  gather(disc_set, member, -snp) %>%
  filter(member == 1) %>%
  select(-member) %>%
  group_by(snp) %>%
  summarize(disc_sets = list(disc_set))

# Next create the ggupset plot
scz18_disc_upset <- tidy_scz18_comp_disc %>%
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


# Create a chart that just counts how many discoveries each type has:
scz18_count_disc_chart <- scz18_comp_disc_list %>%
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
  mutate(int_only_disc = as.numeric(scz18_intercept_only$qvals <= 0.05),
         boosted_disc = as.numeric(scz18_with_bd_z_eqtl_slopes_wgcna$qvals <= 0.05))

# Now will also make a column for random sampling of the tests with p-values greater
# than 0.1 for saving memory with the PDF images:
bip_scz_brainvar_data <- bip_scz_brainvar_data %>%
  mutate(manhattan_sample = ifelse(scz_18_P <= 0.1,
                                   TRUE,
                                   rbernoulli(n(), p = 0.5)))

# Create the intercept-only Manhattan plot:
int_only_manhattan <- bip_scz_brainvar_data %>%
  filter(manhattan_sample) %>%
  ggplot(aes(x = bip_18_BP_cum, y = -log10(scz_18_P))) +
  # Show all points
  geom_point(aes(color = as.factor(bip_18_CHR)), alpha = 0.5, size = 0.5) +
  geom_point(data = filter(bip_scz_brainvar_data, int_only_disc == 1),
             aes(x = bip_18_BP_cum, y = -log10(scz_18_P)),
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
rm(scz18_intercept_only)

# Next using all covariates via gradient boosted trees:
boosted_manhattan <- bip_scz_brainvar_data %>%
  filter(manhattan_sample) %>%
  ggplot(aes(x = bip_18_BP_cum, y = -log10(scz_18_P))) +
  # Show all points
  geom_point(aes(color = as.factor(bip_18_CHR)), alpha = 0.5, size = 0.5) +
  geom_point(data = filter(bip_scz_brainvar_data, boosted_disc == 1),
             aes(x = bip_18_BP_cum, y = -log10(scz_18_P)),
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
scz18_manhattan_plots <-
  plot_grid(int_only_manhattan,
            boosted_manhattan,
            ncol = 1, labels = c("A", "B"), 
            label_fontface = "plain", align = "hv", vjust = 1.2)

# Finally combine everything together:
scz18_brainvar_disc_plot_plus <- 
  plot_grid(plot_grid(scz18_manhattan_plots,
                      scz18_count_disc_chart,
                      ncol = 1, rel_widths = c(1, 1),
                      rel_heights = c(4, 1.75), labels = c("", "C"),
                      label_fontface = "plain"),
            scz18_disc_upset, ncol = 2, rel_widths = c(1, 1),
            rel_heights = c(1, 1), labels = c("", "D"),
            label_fontface = "plain")

# Save the figure:
save_plot("figures/sf12_scz18_brainvar_results.pdf",
          scz18_brainvar_disc_plot_plus, ncol = 3, nrow = 2,
          base_width = 5, base_height = 3)
save_plot("nonpdf_figures/sf12_scz18_brainvar_results.jpg",
          scz18_brainvar_disc_plot_plus, ncol = 3, nrow = 2,
          base_width = 5, base_height = 3)

# ------------------------------------------------------------------------------

# Next create the variable importance plots for the model with all covariates
# and interactions:
scz18_pi_model_i <- seq(1, length(scz18_with_bd_z_eqtl_slopes_wgcna$model_fit), by = 2)
scz18_pi_search_importance_data <- map_dfr(1:length(scz18_pi_model_i),
                                           function(step_i) {
                                             pi_model_step_i <- scz18_pi_model_i[step_i]
                                             xgb.importance(model = scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[pi_model_step_i]]) %>%
                                               as.data.frame() %>%
                                               mutate(adapt_step = step_i)
                                           })
# Will highlight the top variables from the final model:
scz18_final_top_pi_importance <- scz18_pi_search_importance_data %>%
  filter(adapt_step == length(scz18_pi_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:5)

# Now create a display of the importance over the search with the top variable
# highlighted:
scz18_pi_importance_search_plot <- scz18_pi_search_importance_data %>%
  filter(Feature %in% scz18_final_top_pi_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
           str_replace_all("ave abs ", "Average |") %>%
           str_replace_all("beta", "beta|") %>%
           str_replace_all("brainvar any gene", "WGCNA module:") %>%
           str_replace_all("grey", "gray") %>%
           str_replace_all("z bip 14", "BD z-statistics")) %>%
  ggplot() +
  geom_line(data = {(
    scz18_pi_search_importance_data %>%
      filter(!(Feature %in% scz18_final_top_pi_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    scz18_pi_search_importance_data %>%
      filter(Feature %in% scz18_final_top_pi_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_replace_all("ave abs ", "Average |") %>%
               str_replace_all("beta", "beta|") %>%
               str_replace_all("brainvar any gene", "WGCNA module:") %>%
               str_replace_all("grey", "gray") %>%
               str_replace_all("z bip 14", "BD z-statistics")) %>%
      filter(adapt_step == length(scz18_pi_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_manual(values = c(
    rep("goldenrod4", 3),
    "darkblue", "gray")) +  
  scale_linetype_manual(guide = FALSE,
                        values = c(c("dotdash", "dotted", "dashed"), rep("solid", 2))) +
  scale_x_continuous(limits = c(1, 25),
                     breaks = seq(1, length(scz18_pi_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\pi_1$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Next repeat but for the \mu models:
scz18_mu_model_i <- seq(2, length(scz18_with_bd_z_eqtl_slopes_wgcna$model_fit), by = 2)
scz18_mu_search_importance_data <- map_dfr(1:length(scz18_mu_model_i),
                                           function(step_i) {
                                             mu_model_step_i <- scz18_mu_model_i[step_i]
                                             xgb.importance(model = scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[mu_model_step_i]]) %>%
                                               as.data.frame() %>%
                                               mutate(adapt_step = step_i)
                                           })
# Will highlight the top variables from the final model:
scz18_final_top_mu_importance <- scz18_mu_search_importance_data %>%
  filter(adapt_step == length(scz18_mu_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:5)
scz18_mu_importance_search_plot <- scz18_mu_search_importance_data %>%
  filter(Feature %in% scz18_final_top_mu_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
           str_replace_all("ave abs ", "Average |") %>%
           str_replace_all("beta", "beta|") %>%
           str_replace_all("brainvar any gene", "WGCNA module:") %>%
           str_replace_all("grey", "gray") %>%
           str_replace_all("z bip 14", "BD z-statistics")) %>%
  ggplot() +
  geom_line(data = {(
    scz18_mu_search_importance_data %>%
      filter(!(Feature %in% scz18_final_top_mu_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    scz18_mu_search_importance_data %>%
      filter(Feature %in% scz18_final_top_mu_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_replace_all("ave abs ", "Average |") %>%
               str_replace_all("beta", "beta|") %>%
               str_replace_all("brainvar any gene", "WGCNA module:") %>%
               str_replace_all("grey", "gray") %>%
               str_replace_all("z bip 14", "BD z-statistics")) %>%
      filter(adapt_step == length(scz18_mu_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_manual(values = c(
    rep("goldenrod4", 3),
    "darkblue", "cyan")) +  
  scale_linetype_manual(guide = FALSE,
                        values = c(c("dotdash", "dotted", "dashed"), rep("solid", 2))) +
  scale_x_continuous(limits = c(1, 25),
                     breaks = seq(1, length(scz18_mu_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\mu$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Arrange these in a grid and save:
scz18_var_imp_grid <-
  plot_grid(scz18_pi_importance_search_plot,
            scz18_mu_importance_search_plot, ncol = 1,
            align = "hv", labels = c("A", "B"),
            label_fontface = "plain")
save_plot("figures/sf13_scz18_var_importance.pdf",
          scz18_var_imp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf13_scz18_var_importance.jpg",
          scz18_var_imp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------

# Next create the change in the partial dependence plot over the course of the 
# search with AdaPT for the \pi_1 model - will do so at quantiles:
scz18_pi_quantiles_search_pdp_data <- map_dfr(1:length(scz18_pi_model_i),
                                              function(step_i) {
                                                pi_model_step_i <- scz18_pi_model_i[step_i]
                                                pdp::partial(scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[pi_model_step_i]], 
                                                             pred.var = "z_bip_14", 
                                                             ice = FALSE, 
                                                             prob = TRUE,
                                                             center = FALSE, 
                                                             plot = FALSE, 
                                                             quantiles = TRUE,
                                                             probs = seq(0, 1, by = .025),
                                                             train = data.matrix(
                                                               bip_scz_brainvar_data[,scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[1]]$feature_names])) %>%
                                                  as.data.frame() %>%
                                                  mutate(adapt_step = step_i)
                                              })
scz18_pi_quantiles_pdp_search_plot <- scz18_pi_quantiles_search_pdp_data %>%
  ggplot() +
  geom_line(aes(x = z_bip_14, y = yhat, 
                color = adapt_step,
                group = as.factor(adapt_step))) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_cowplot() +
  scale_color_gradient(low = "darkblue", high = "darkorange") +
  geom_rug(data = bip_scz_brainvar_data,
           aes(x = z_bip_14), y = rep(1, nrow(bip_scz_brainvar_data)),
           sides = "b",
           color = "black", alpha = 0.15) +
  geom_vline(xintercept = 1.96, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = -1.96, linetype = "dashed", color = "darkred") +
  guides(color = guide_colourbar(barwidth = 10, barheight = .5)) +
  labs(x = "BD z-statistic",
       y = TeX('$\\pi_1$'),
       color = "AdaPT model fitting iteration",
       subtitle = "Dashed red lines indicate z-statistics equal to +/- 1.96",
       title = TeX('Change in partial dependence for $\\pi_1$ and BD z-statistics')) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# Next the change in non-null effect size pdp
scz18_mu_quantiles_search_pdp_data <- map_dfr(1:length(scz18_mu_model_i),
                                              function(step_i) {
                                                mu_model_step_i <- scz18_mu_model_i[step_i]
                                                pdp::partial(scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[mu_model_step_i]], 
                                                             pred.var = "z_bip_14", 
                                                             ice = FALSE, 
                                                             center = FALSE, 
                                                             plot = FALSE, 
                                                             quantiles = TRUE,
                                                             probs = seq(0, 1, by = .025),
                                                             train = data.matrix(
                                                               bip_scz_brainvar_data[,scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[1]]$feature_names])) %>%
                                                  as.data.frame() %>%
                                                  mutate(adapt_step = step_i)
                                              })

scz18_mu_quantiles_pdp_search_plot <- scz18_mu_quantiles_search_pdp_data %>%
  ggplot() +
  geom_line(aes(x = z_bip_14, y = yhat, 
                color = adapt_step,
                group = as.factor(adapt_step))) +
  scale_y_continuous(limits = c(0, 4)) +
  theme_cowplot() +
  scale_color_gradient(low = "darkblue", high = "darkorange") +
  geom_rug(data = bip_scz_brainvar_data,
           aes(x = z_bip_14), y = rep(1, nrow(bip_scz_brainvar_data)),
           sides = "b",
           color = "black", alpha = 0.15) +
  geom_vline(xintercept = 1.96, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = -1.96, linetype = "dashed", color = "darkred") +
  guides(color = guide_colourbar(barwidth = 10, barheight = .5)) +
  labs(x = "BD z-statistic",
       y = TeX('$\\mu$'),
       color = "AdaPT model fitting iteration",
       subtitle = "Dashed red lines indicate z-statistics equal to +/- 1.96",
       title = TeX('Change in partial dependence for $\\mu$ and BD z-statistics')) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# Save the plot together:
scz18_bd_pdp_grid <- plot_grid(
  plot_grid(scz18_pi_quantiles_pdp_search_plot + theme(legend.position = "none"),
            scz18_mu_quantiles_pdp_search_plot + theme(legend.position = "none"), 
            ncol = 1,
            align = "hv", labels = c("A", "B"),
            label_fontface = "plain"),
  get_legend(scz18_pi_quantiles_pdp_search_plot), 
  ncol = 1, rel_heights = c(2, .5))
save_plot("figures/sf14_scz18_bd_pdp.pdf",
          scz18_bd_pdp_grid, ncol = 1, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf14_scz18_bd_pdp.jpg",
          scz18_bd_pdp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------

# Next make pdp plots for the eQTL slopes - for both \pi_1 and \mu models

# First set up a list of the eQTL slope variables:
brainvar_eqtl_slope_list <- list("ave_abs_pre_beta" = "Average |pre beta|",
                                 "ave_abs_post_beta" = "Average |post beta|",
                                 "ave_abs_comp_beta" = "Average |comp beta|")

# Now generate the \pi_1 plots for each of these:
eqtl_slope_pi_plots <- lapply(1:length(brainvar_eqtl_slope_list),
                              function(x) {
                                var_name <- names(brainvar_eqtl_slope_list)[x]
                                var_axis <- brainvar_eqtl_slope_list[[x]]
                                
                                var_pi_quantiles_search_pdp_data <- map_dfr(1:length(scz18_pi_model_i),
                                                                            function(step_i) {
                                                                              pi_model_step_i <- scz18_pi_model_i[step_i]
                                                                              pdp::partial(scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[pi_model_step_i]], 
                                                                                           pred.var = var_name, 
                                                                                           ice = FALSE, 
                                                                                           prob = TRUE,
                                                                                           center = FALSE, 
                                                                                           plot = FALSE, 
                                                                                           quantiles = TRUE,
                                                                                           probs = seq(0, 1, by = .025),
                                                                                           train = data.matrix(
                                                                                             bip_scz_brainvar_data[,scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[1]]$feature_names])) %>%
                                                                                as.data.frame() %>%
                                                                                mutate(adapt_step = step_i)
                                                                            })
                                var_pi_quantiles_search_pdp_data %>%
                                  mutate(adapt_step_factor = as.factor(adapt_step)) %>%
                                  ggplot() +
                                  geom_line(aes_string(x = var_name, 
                                                       y = "yhat", 
                                                       color = "adapt_step",
                                                       group = "adapt_step_factor")) +
                                  scale_y_continuous(limits = c(0, 1)) +
                                  theme_cowplot() +
                                  scale_color_gradient(low = "darkblue", high = "darkorange") +
                                  geom_rug(data = bip_scz_brainvar_data,
                                           aes_string(x = var_name), 
                                           y = rep(1, nrow(bip_scz_brainvar_data)),
                                           sides = "b",
                                           color = "black", alpha = 0.15) +
                                  guides(color = guide_colourbar(barwidth = 10, barheight = .5)) +
                                  labs(x = var_axis,
                                       y = TeX('$\\pi_1$'),
                                       color = "AdaPT model fitting iteration",
                                       title = TeX(paste0('Change in partial dependence for $\\pi_1$ and ',
                                                          var_axis))) +
                                  theme(legend.position = "bottom")
                              })

eqtl_slope_mu_plots <- lapply(1:length(brainvar_eqtl_slope_list),
                              function(x) {
                                var_name <- names(brainvar_eqtl_slope_list)[x]
                                var_axis <- brainvar_eqtl_slope_list[[x]]
                                
                                var_mu_quantiles_search_pdp_data <- map_dfr(1:length(scz18_mu_model_i),
                                                                            function(step_i) {
                                                                              mu_model_step_i <- scz18_mu_model_i[step_i]
                                                                              pdp::partial(scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[mu_model_step_i]], 
                                                                                           pred.var = var_name, 
                                                                                           ice = FALSE, 
                                                                                           center = FALSE, 
                                                                                           plot = FALSE, 
                                                                                           quantiles = TRUE,
                                                                                           probs = seq(0, 1, by = .025),
                                                                                           train = data.matrix(
                                                                                             bip_scz_brainvar_data[,scz18_with_bd_z_eqtl_slopes_wgcna$model_fit[[1]]$feature_names])) %>%
                                                                                as.data.frame() %>%
                                                                                mutate(adapt_step = step_i)
                                                                            })
                                
                                var_mu_quantiles_search_pdp_data %>%
                                  mutate(adapt_step_factor = as.factor(adapt_step)) %>%
                                  ggplot() +
                                  geom_line(aes_string(x = var_name, 
                                                       y = "yhat", 
                                                       color = "adapt_step",
                                                       group = "adapt_step_factor")) +
                                  scale_y_continuous(limits = c(0, 3)) +
                                  theme_cowplot() +
                                  scale_color_gradient(low = "darkblue", high = "darkorange") +
                                  geom_rug(data = bip_scz_brainvar_data,
                                           aes_string(x = var_name), 
                                           y = rep(1, nrow(bip_scz_brainvar_data)),
                                           sides = "b",
                                           color = "black", alpha = 0.15) +
                                  guides(color = guide_colourbar(barwidth = 10, barheight = .5)) +
                                  labs(x = var_axis,
                                       y = TeX('$\\mu$'),
                                       color = "AdaPT model fitting iteration",
                                       title = TeX(paste0('Change in partial dependence for $\\mu$ and ',
                                                          var_axis))) +
                                  theme(legend.position = "bottom")
                              })

# Next arrange in a grid with a shared legend and save:
eqtl_slope_plot_grid <- plot_grid(
  plot_grid(plot_grid(plotlist = 
                        lapply(eqtl_slope_pi_plots,
                               function(eqtl_plot) eqtl_plot + theme(legend.position = "none", plot.title = element_text(size = 12))),
                      ncol = 3, align = "hv", labels = c("A", "B", "C"), label_fontface = "plain"),
            plot_grid(plotlist = 
                        lapply(eqtl_slope_mu_plots,
                               function(eqtl_plot) eqtl_plot + theme(legend.position = "none", plot.title = element_text(size = 12))),
                      ncol = 3, align = "hv",
                      labels = c("D", "E", "F"), label_fontface = "plain"), 
            ncol = 1, align = "hv"),
  get_legend(eqtl_slope_pi_plots[[1]]), ncol = 1,
  rel_heights = c(2, .5))

save_plot("figures/sf15_scz18_brainvar_eqtl_slopes_pdp.pdf",
          eqtl_slope_plot_grid, ncol = 3, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf15_scz18_brainvar_eqtl_slopes_pdp.jpg",
          eqtl_slope_plot_grid, ncol = 3, nrow = 2,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------

# Finally display the WGCNA enrichment plots 
brainvar_any_cols <- colnames(bip_scz_brainvar_data)[
  str_detect(colnames(bip_scz_brainvar_data), "brainvar_any_gene_")]

# Make an example plot to get the legend with appropriate sizing:
darkgreen_pval_distr <- bip_scz_brainvar_data %>%
  ggplot(aes(x = scz_18_P, 
             fill = as.factor(brainvar_any_gene_darkgreen))) +
  geom_density(color = "black", alpha = 0.5) +
  geom_rug(aes(color = as.factor(brainvar_any_gene_darkgreen))) +
  scale_fill_manual(values = c("darkblue", "darkorange"),
                    labels = c("No", "Yes")) +
  scale_color_manual(values = c("darkblue", "darkorange"), guide = FALSE) +
  theme_bw() +
  labs(x = "SCZ p-value", y = "Density",
       fill = "Any cis-eQTL in module?") +
  guides(fill = guide_legend(override.aes = list(size = 12))) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 28))

# Now generate the histograms for each module:
wgcna_member_pval_distr_hist_list <- lapply(brainvar_any_cols,
                                            function(x) {
                                              module_color <- str_remove(x,
                                                                         "brainvar_any_gene_") %>%
                                                str_replace_all("grey", "gray")
                                              # Create the plot - without the legend:
                                              plot_data <- bip_scz_brainvar_data[, c("scz_18_P", x)]
                                              colnames(plot_data) <- c("pvals", "wgcna")
                                              plot_data %>%
                                                ggplot(aes(x = pvals, 
                                                           fill = as.factor(wgcna),
                                                           color = as.factor(wgcna))) +
                                                geom_histogram(breaks = seq(0, 1, by = 0.05),
                                                               alpha = 0.25, position = "identity") +
                                                scale_fill_manual(values = c("darkblue", "darkorange"),
                                                                  labels = c("No", "Yes"), guide = FALSE) +
                                                scale_color_manual(values = c("darkblue", "darkorange"), guide = FALSE) +
                                                theme_bw() + 
                                                facet_wrap(~as.factor(wgcna), ncol = 1,
                                                           scales = "free_y") +
                                                theme(strip.text = element_blank(),
                                                      strip.background = element_blank(),
                                                      plot.title = element_text(size = 20),
                                                      axis.title = element_text(size = 16),
                                                      axis.text = element_text(size = 14)) +
                                                labs(x = "SCZ p-value", y = "Count",
                                                     title = paste0(module_color, " module"))
                                            })
# Now generate all the plots using cowplot:
scz_wgcna_module_enrichment_grid <-
  plot_grid(plot_grid(plotlist =  wgcna_member_pval_distr_hist_list, ncol = 5),
            get_legend(darkgreen_pval_distr), ncol = 1, rel_heights = c(1, .05))

save_plot("figures/sf16_scz18_brainvar_wgcna.pdf",
          scz_wgcna_module_enrichment_grid, ncol = 5, nrow = 4,
          base_width = 5, base_height = 4)

save_plot("nonpdf_figures/sf16_scz18_brainvar_wgcna.jpg",
          scz_wgcna_module_enrichment_grid, ncol = 5, nrow = 4,
          base_width = 5, base_height = 4)










