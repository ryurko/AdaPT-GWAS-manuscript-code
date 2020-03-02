# PURPOSE: Create Figure 3 in the manuscript, displaying the variable importance
#          plots over the entire AdaPT search for the \pi_1 model and the 
#          corresponding change in partial dependence plot. Still include a module
#          specific enrichment plot.

# Author: Ron Yurko

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

# Load the final model results:
scz_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")

# First create a display of the change in variable importance over the AdaPT search

# Start with the pi_1 models (odd numbers) - stacking the importance matrices
# together for each of the models:
pi_model_i <- seq(1, length(scz_with_bd_z_eqtl_slopes_wgcna$model_fit), by = 2)
pi_search_importance_data <- map_dfr(1:length(pi_model_i),
                                     function(step_i) {
                                       pi_model_step_i <- pi_model_i[step_i]
                                       xgb.importance(model = scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[pi_model_step_i]]) %>%
                                         as.data.frame() %>%
                                         mutate(adapt_step = step_i)
                                     })
# Will highlight the top variables from the final model:
final_top_pi_importance <- pi_search_importance_data %>%
  filter(adapt_step == length(pi_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:8)

# Now create a display of the importance over the search with the top variable
# highlighted:
pi_importance_search_plot <- pi_search_importance_data %>%
  filter(Feature %in% final_top_pi_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>%
           str_replace_all("ave abs ", "Average |") %>%
           str_replace_all("beta", "beta|") %>%
           str_replace_all("brainvar any gene", "WGCNA module:") %>%
           str_replace_all("grey", "gray") %>%
           str_replace_all("z bip 14", "BD z-statistics")) %>%
  ggplot() +
  geom_line(data = {(
    pi_search_importance_data %>%
      filter(!(Feature %in% final_top_pi_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    pi_search_importance_data %>%
      filter(Feature %in% final_top_pi_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_replace_all("ave abs ", "Average |") %>%
               str_replace_all("beta", "beta|") %>%
               str_replace_all("brainvar any gene", "WGCNA module:") %>%
               str_replace_all("grey", "gray") %>%
               str_replace_all("z bip 14", "BD z-statistics")) %>%
      filter(adapt_step == length(pi_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = c(#"#E41A1C", "cyan", "#4DAF4A",
    rep("goldenrod4", 3),
    "darkblue",
    "black", "brown",
    "gray50", "salmon1")) +
  scale_linetype_manual(guide = FALSE,
                        values = c("dotdash", "dotted", 
                                   "longdash", #rep("longdash", 3),
                                   rep("solid", 5))) +
  scale_x_continuous(limits = c(1, 25),
                     breaks = seq(1, length(pi_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\pi_1$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Next create the change in the partial dependence plot over the course of the 
# search with AdaPT for the \pi_1 model - will do so at quantiles:
pi_quantiles_search_pdp_data <- map_dfr(1:length(pi_model_i),
                                        function(step_i) {
                                          pi_model_step_i <- pi_model_i[step_i]
                                          pdp::partial(scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[pi_model_step_i]], 
                                                       pred.var = "z_bip_14", 
                                                       ice = FALSE, 
                                                       prob = TRUE,
                                                       center = FALSE, 
                                                       plot = FALSE, 
                                                       quantiles = TRUE,
                                                       probs = seq(0, 1, by = .025),
                                                       train = data.matrix(
                                                         bip_scz_brainvar_data[,scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[1]]$feature_names])) %>%
                                            as.data.frame() %>%
                                            mutate(adapt_step = step_i)
                                        })
pi_quantiles_pdp_search_plot <- pi_quantiles_search_pdp_data %>%
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

# Finally, create a plot revealing the enrichment for the salmon module for SCZ
scz_brainvar_salmon_plot <- bip_scz_brainvar_data %>% 
  ggplot(aes(x = scz_14_P, 
             fill = as.factor(brainvar_any_gene_salmon),
             color = as.factor(brainvar_any_gene_salmon))) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), alpha = 0.5) +
  scale_fill_manual(values = c("darkblue", "darkorange"),
                    labels = c("No", "Yes")) +
  scale_color_manual(values = c("darkblue", "darkorange"), guide = FALSE) +
  theme_bw() + 
  facet_wrap(~as.factor(brainvar_any_gene_salmon), ncol = 1,
             scales = "free_y") +
  theme(strip.text = element_blank(),
        plot.title = element_text(size = 12),
        strip.background = element_blank(), 
        legend.position = c(0.5, 0.3),
        legend.direction = "vertical",
        legend.background = element_blank()) +
  labs(x = "SCZ p-values", y = "Count",
       title = "Distribution of SCZ p-values by salmon module assignment",
       fill = "Any cis-eQTL gene in salmon module?")


# Now create the grid of plots for figure 4 in the paper:
scz_brainvar_variables_f3 <- 
  plot_grid(pi_importance_search_plot, 
            plot_grid(pi_quantiles_pdp_search_plot,
                      scz_brainvar_salmon_plot,
                      ncol = 2, labels = c("B", "C"),
                      label_fontface = "plain", rel_widths = c(1, 1),
                      rel_heights = c(1, 1)),
            labels = c("A", ""), label_fontface = "plain",
            ncol = 1)

# Save
save_plot("figures/f3_scz_variable_plots.pdf",
          scz_brainvar_variables_f3, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/f3_scz_variable_plots.jpg",
          scz_brainvar_variables_f3, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)



