# PURPOSE: Create the supplementary figures displaying the partial dependence
#          plots for only using BD z-statistics and no interactions

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

# Load the BD z-stats only and no interaction models
scz_with_bd_z_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_only_s05_2cv.rds")
scz_with_no_int <-
  readRDS("data/bip_schz_data/brainvar_results/scz_all_vars_no_int_s05.rds")

# Set up the model indices
bd_z_pi_model_i <- seq(1, length(scz_with_bd_z_only$model_fit), by = 2)
no_int_pi_model_i <- seq(1, length(scz_with_no_int$model_fit), by = 2)

# Create the PDP for the BD z-stats only for \pi_1
bd_z_pi_quantiles_search_pdp_data <- map_dfr(1:length(bd_z_pi_model_i),
                                             function(step_i) {
                                               pi_model_step_i <- bd_z_pi_model_i[step_i]
                                               pdp::partial(scz_with_bd_z_only$model_fit[[pi_model_step_i]], 
                                                            pred.var = "z_bip_14", 
                                                            ice = FALSE, 
                                                            prob = TRUE,
                                                            center = FALSE, 
                                                            plot = FALSE, 
                                                            quantiles = TRUE,
                                                            probs = seq(0, 1, by = .025),
                                                            train = data.matrix(
                                                              bip_scz_brainvar_data[,scz_with_bd_z_only$model_fit[[1]]$feature_names])) %>%
                                                 as.data.frame() %>%
                                                 mutate(adapt_step = step_i)
                                             })
bd_z_pi_quantiles_pdp_search_plot <- bd_z_pi_quantiles_search_pdp_data %>%
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
       title = TeX('Change in partial dependence BD z-statistics only model')) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# Next for no-interactions:
no_int_pi_quantiles_search_pdp_data <- map_dfr(1:length(no_int_pi_model_i),
                                               function(step_i) {
                                                 pi_model_step_i <- no_int_pi_model_i[step_i]
                                                 pdp::partial(scz_with_no_int$model_fit[[pi_model_step_i]], 
                                                              pred.var = "z_bip_14", 
                                                              ice = FALSE, 
                                                              prob = TRUE,
                                                              center = FALSE, 
                                                              plot = FALSE, 
                                                              quantiles = TRUE,
                                                              probs = seq(0, 1, by = .025),
                                                              train = data.matrix(
                                                                bip_scz_brainvar_data[,scz_with_no_int$model_fit[[1]]$feature_names])) %>%
                                                   as.data.frame() %>%
                                                   mutate(adapt_step = step_i)
                                               })
no_int_pi_quantiles_pdp_search_plot <- no_int_pi_quantiles_search_pdp_data %>%
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
       title = TeX('Change in partial dependence for model without interactions')) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# Arrange the grid and save plots:
bd_z_other_pdp_grid <- plot_grid(
  plot_grid(bd_z_pi_quantiles_pdp_search_plot + theme(legend.position = "none"),
            no_int_pi_quantiles_pdp_search_plot + theme(legend.position = "none"), 
            ncol = 2,
            align = "hv", labels = c("A", "B"),
            label_fontface = "plain"),
  get_legend(bd_z_pi_quantiles_pdp_search_plot), 
  ncol = 1, rel_heights = c(2, .5))
save_plot("figures/sf9_scz_bd_z_other_pdp.pdf",
          bd_z_other_pdp_grid, ncol = 2, nrow = 1,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf9_scz_bd_z_other_pdp.jpg",
          bd_z_other_pdp_grid, ncol = 2, nrow = 1,
          base_width = 6, base_height = 4)











