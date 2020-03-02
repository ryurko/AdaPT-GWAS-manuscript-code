# PURPOSE: Create the supplementary figure displaying the variable importance
#          and partial dependence plots for the SCZ GTEx results

# AUTHOR: Ron Yurko

# Load necessary packages:
library(tidyverse)
library(data.table)
library(cowplot)
library(latex2exp)
library(xgboost)
library(pdp)

# ------------------------------------------------------------------------------
# Now load the GTEx version of results:
scz_gtex_snps <- read_csv("data/bip_schz_data/bip_scz_data_14_18_gtex_eqtls_cortical_wgcna.csv")
scz_gtex_adapt <- readRDS("data/bip_schz_data/gtex_results/cv_tune_results/scz_gtex_cortical_s05_2cv.rds")

# ------------------------------------------------------------------------------

# Start with the pi_1 models (odd numbers) - stacking the importance matrices
# together for each of the models:
gtex_pi_model_i <- seq(1, length(scz_gtex_adapt$model_fit), by = 2)
gtex_pi_search_importance_data <- map_dfr(1:length(gtex_pi_model_i),
                                          function(step_i) {
                                            pi_model_step_i <- gtex_pi_model_i[step_i]
                                            xgb.importance(model = scz_gtex_adapt$model_fit[[pi_model_step_i]]) %>%
                                              as.data.frame() %>%
                                              mutate(adapt_step = step_i)
                                          })
# Will highlight the top variables from the final model:
gtex_final_top_pi_importance <- gtex_pi_search_importance_data %>%
  filter(adapt_step == length(gtex_pi_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:4)

# Now create a display of the importance over the search with the top variable
# highlighted:
gtex_pi_importance_search_plot <- gtex_pi_search_importance_data %>%
  filter(Feature %in% gtex_final_top_pi_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
           str_remove("brain") %>%
           str_replace("ba9", "BA9") %>%
           str_replace("ba24", "BA24") %>%
           str_replace_all("ave abs eqtl slope",
                           "average |eQTL beta|") %>%
           str_replace_all("ave abs cortical eqtl slope",
                           "cortical average |eQTL beta|") %>%
           str_replace_all("z bip 14", "BD z-statistics")) %>%
  ggplot() +
  geom_line(data = {(
    gtex_pi_search_importance_data %>%
      filter(!(Feature %in% gtex_final_top_pi_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    gtex_pi_search_importance_data %>%
      filter(Feature %in% gtex_final_top_pi_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_remove("brain") %>%
               str_replace("ba9", "BA9") %>%
               str_replace("ba24", "BA24") %>%
               str_replace_all("ave abs eqtl slope",
                               "average |eQTL beta|") %>%
               str_replace_all("ave abs cortical eqtl slope",
                               "cortical average |eQTL beta|") %>%
               str_replace_all("z bip 14", "BD z-statistics")) %>%
      filter(adapt_step == length(gtex_pi_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(1, 25),
                     breaks = seq(1, length(gtex_pi_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\pi_1$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Next repeat but for the \mu models:
gtex_mu_model_i <- seq(2, length(scz_gtex_adapt$model_fit), by = 2)
gtex_mu_search_importance_data <- map_dfr(1:length(gtex_mu_model_i),
                                          function(step_i) {
                                            mu_model_step_i <- gtex_mu_model_i[step_i]
                                            xgb.importance(model = scz_gtex_adapt$model_fit[[mu_model_step_i]]) %>%
                                              as.data.frame() %>%
                                              mutate(adapt_step = step_i)
                                          })
# Will highlight the top variables from the final model:
gtex_final_top_mu_importance <- gtex_mu_search_importance_data %>%
  filter(adapt_step == length(gtex_mu_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:4)
gtex_mu_importance_search_plot <- gtex_mu_search_importance_data %>%
  filter(Feature %in% gtex_final_top_mu_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
           str_remove("brain") %>%
           str_replace("ba9", "BA9") %>%
           str_replace("ba24", "BA24") %>%
           str_replace_all("ave abs eqtl slope",
                           "average |eQTL beta|") %>%
           str_replace_all("ave abs cortical eqtl slope",
                           "cortical average |eQTL beta|") %>%
           str_replace_all("z bip 14", "BD z-statistics")) %>%
  ggplot() +
  geom_line(data = {(
    gtex_mu_search_importance_data %>%
      filter(!(Feature %in% gtex_final_top_mu_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    gtex_mu_search_importance_data %>%
      filter(Feature %in% gtex_final_top_mu_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_remove("brain") %>%
               str_replace("ba9", "BA9") %>%
               str_replace("ba24", "BA24") %>%
               str_replace_all("ave abs eqtl slope",
                               "average |eQTL beta|") %>%
               str_replace_all("ave abs cortical eqtl slope",
                               "cortical average |eQTL beta|") %>%
               str_replace_all("z bip 14", "BD z-statistics")) %>%
      filter(adapt_step == length(gtex_mu_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(1, 25),
                     breaks = seq(1, length(gtex_mu_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\mu$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Arrange these in a grid and save:
scz_gtex_var_imp_grid <-
  plot_grid(gtex_pi_importance_search_plot,
            gtex_mu_importance_search_plot, ncol = 1,
            align = "hv", labels = c("A", "B"),
            label_fontface = "plain")
save_plot("figures/sf10_scz_gtex_var_importance.pdf",
          scz_gtex_var_imp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf10_scz_gtex_var_importance.jpg",
          scz_gtex_var_imp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------

#  Generate the change in partial dependence plots
gtex_pi_quantiles_search_pdp_data <- map_dfr(1:length(gtex_pi_model_i),
                                             function(step_i) {
                                               pi_model_step_i <- gtex_pi_model_i[step_i]
                                               pdp::partial(scz_gtex_adapt$model_fit[[pi_model_step_i]], 
                                                            pred.var = "z_bip_14", 
                                                            ice = FALSE, 
                                                            prob = TRUE,
                                                            center = FALSE, 
                                                            plot = FALSE, 
                                                            quantiles = TRUE,
                                                            probs = seq(0, 1, by = .025),
                                                            train = data.matrix(
                                                              scz_gtex_snps[,scz_gtex_adapt$model_fit[[1]]$feature_names])) %>%
                                                 as.data.frame() %>%
                                                 mutate(adapt_step = step_i)
                                             })
gtex_pi_quantiles_pdp_search_plot <- gtex_pi_quantiles_search_pdp_data %>%
  ggplot() +
  geom_line(aes(x = z_bip_14, y = yhat, 
                color = adapt_step,
                group = as.factor(adapt_step))) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_cowplot() +
  scale_color_gradient(low = "darkblue", high = "darkorange") +
  geom_rug(data = scz_gtex_snps,
           aes(x = z_bip_14), y = rep(1, nrow(scz_gtex_snps)),
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

gtex_mu_quantiles_search_pdp_data <- map_dfr(1:length(gtex_mu_model_i),
                                             function(step_i) {
                                               mu_model_step_i <- gtex_mu_model_i[step_i]
                                               pdp::partial(scz_gtex_adapt$model_fit[[mu_model_step_i]], 
                                                            pred.var = "z_bip_14", 
                                                            ice = FALSE, 
                                                            center = FALSE, 
                                                            plot = FALSE, 
                                                            quantiles = TRUE,
                                                            probs = seq(0, 1, by = .025),
                                                            train = data.matrix(
                                                              scz_gtex_snps[,scz_gtex_adapt$model_fit[[1]]$feature_names])) %>%
                                                 as.data.frame() %>%
                                                 mutate(adapt_step = step_i)
                                             })

gtex_mu_quantiles_pdp_search_plot <- gtex_mu_quantiles_search_pdp_data %>%
  ggplot() +
  geom_line(aes(x = z_bip_14, y = yhat, 
                color = adapt_step,
                group = as.factor(adapt_step))) +
  scale_y_continuous(limits = c(0, 5)) +
  theme_cowplot() +
  scale_color_gradient(low = "darkblue", high = "darkorange") +
  geom_rug(data = scz_gtex_snps,
           aes(x = z_bip_14), y = rep(1, nrow(scz_gtex_snps)),
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
scz_gtex_pdp_grid <- plot_grid(
  plot_grid(gtex_pi_quantiles_pdp_search_plot + theme(legend.position = "none"),
            gtex_mu_quantiles_pdp_search_plot + theme(legend.position = "none"), 
            ncol = 1,
            align = "hv", labels = c("A", "B"),
            label_fontface = "plain"),
  get_legend(gtex_pi_quantiles_pdp_search_plot), 
  ncol = 1, rel_heights = c(2, .5))
save_plot("figures/sf11_scz_gtex_bd_pdp.pdf",
          scz_gtex_pdp_grid, ncol = 1, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf11_scz_gtex_bd_pdp.jpg",
          scz_gtex_pdp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)


