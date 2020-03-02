# PURPOSE: Create supplementary figures for the non-null effect size and partial
#          dependence plots to display, with additional figures conveying the
#          importance of allowing interactions in the model.

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

# Load the final model results:
scz_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")

# ------------------------------------------------------------------------------

# First create the change in variable importance for the \mu model over the AdaPT search
mu_model_i <- seq(2, length(scz_with_bd_z_eqtl_slopes_wgcna$model_fit), by = 2)
mu_search_importance_data <- map_dfr(1:length(mu_model_i),
                                     function(step_i) {
                                       mu_model_step_i <- mu_model_i[step_i]
                                       xgb.importance(model = scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[mu_model_step_i]]) %>%
                                         as.data.frame() %>%
                                         mutate(adapt_step = step_i)
                                     })
# Will highlight the top variables from the final model:
final_top_mu_importance <- mu_search_importance_data %>%
  filter(adapt_step == length(mu_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:8)

# Now create a display of the importance over the search with the top variable
# highlighted:
mu_importance_search_plot <- mu_search_importance_data %>%
  filter(Feature %in% final_top_mu_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>%
           str_replace_all("ave abs ", "Average |") %>%
           str_replace_all("beta", "beta|") %>%
           str_replace_all("brainvar any gene", "WGCNA module:") %>%
           str_replace_all("grey", "gray") %>%
           str_replace_all("z bip 14", "BD z-statistics")) %>%
  ggplot() +
  geom_line(data = {(
    mu_search_importance_data %>%
      filter(!(Feature %in% final_top_mu_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    mu_search_importance_data %>%
      filter(Feature %in% final_top_mu_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_replace_all("ave abs ", "Average |") %>%
               str_replace_all("beta", "beta|") %>%
               str_replace_all("brainvar any gene", "WGCNA module:") %>%
               str_replace_all("grey", "gray") %>%
               str_replace_all("z bip 14", "BD z-statistics")) %>%
      filter(adapt_step == length(mu_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_manual(values = c(
    rep("goldenrod4", 3),
    "darkblue",
    "black", "cyan",
    "gray50", "salmon1")) +
  scale_linetype_manual(guide = FALSE,
                        values = c("dotdash", "dotted", 
                                   "longdash", #rep("longdash", 3),
                                   rep("solid", 5))) +
  scale_x_continuous(limits = c(1, 25),
                     breaks = seq(1, length(mu_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\mu$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Save this plot:
save_plot("figures/sf5_scz_mu_var_importance.pdf",
          mu_importance_search_plot, ncol = 2, nrow = 1,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf5_scz_mu_var_importance.jpg",
          mu_importance_search_plot, ncol = 2, nrow = 1,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------

# Next the change in non-null effect size pdp

mu_quantiles_search_pdp_data <- map_dfr(1:length(mu_model_i),
                                        function(step_i) {
                                          mu_model_step_i <- mu_model_i[step_i]
                                          pdp::partial(scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[mu_model_step_i]], 
                                                       pred.var = "z_bip_14", 
                                                       ice = FALSE, 
                                                       center = FALSE, 
                                                       plot = FALSE, 
                                                       quantiles = TRUE,
                                                       probs = seq(0, 1, by = .025),
                                                       train = data.matrix(
                                                         bip_scz_brainvar_data[,scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[1]]$feature_names])) %>%
                                            as.data.frame() %>%
                                            mutate(adapt_step = step_i)
                                        })

mu_quantiles_pdp_search_plot <- mu_quantiles_search_pdp_data %>%
  ggplot() +
  geom_line(aes(x = z_bip_14, y = yhat, 
                color = adapt_step,
                group = as.factor(adapt_step))) +
  scale_y_continuous(limits = c(0, 3)) +
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

# Save this plot:
save_plot("figures/sf6_scz_mu_pdp_search.pdf",
          mu_quantiles_pdp_search_plot, ncol = 1, nrow = 1,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf6_scz_mu_pdp_search.jpg",
          mu_quantiles_pdp_search_plot, ncol = 1, nrow = 1,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------

# Next make pdp plots for the eQTL slopes - for both \pi_1 and \mu models

# First set up a list of the eQTL slope variables:
brainvar_eqtl_slope_list <- list("ave_abs_pre_beta" = "Average |pre beta|",
                                 "ave_abs_post_beta" = "Average |post beta|",
                                 "ave_abs_comp_beta" = "Average |comp beta|")
pi_model_i <- seq(1, length(scz_with_bd_z_eqtl_slopes_wgcna$model_fit), by = 2)

# Now generate the \pi_1 plots for each of these:
eqtl_slope_pi_plots <- lapply(1:length(brainvar_eqtl_slope_list),
                              function(x) {
                                var_name <- names(brainvar_eqtl_slope_list)[x]
                                var_axis <- brainvar_eqtl_slope_list[[x]]
                                
                                var_pi_quantiles_search_pdp_data <- map_dfr(1:length(pi_model_i),
                                                                            function(step_i) {
                                                                              pi_model_step_i <- pi_model_i[step_i]
                                                                              pdp::partial(scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[pi_model_step_i]], 
                                                                                           pred.var = var_name, 
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
                                
                                var_mu_quantiles_search_pdp_data <- map_dfr(1:length(mu_model_i),
                                                                            function(step_i) {
                                                                              mu_model_step_i <- mu_model_i[step_i]
                                                                              pdp::partial(scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[mu_model_step_i]], 
                                                                                           pred.var = var_name, 
                                                                                           ice = FALSE, 
                                                                                           center = FALSE, 
                                                                                           plot = FALSE, 
                                                                                           quantiles = TRUE,
                                                                                           probs = seq(0, 1, by = .025),
                                                                                           train = data.matrix(
                                                                                             bip_scz_brainvar_data[,scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[1]]$feature_names])) %>%
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

save_plot("figures/sf7_scz_brainvar_eqtl_slopes_pdp.pdf",
          eqtl_slope_plot_grid, ncol = 3, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf7_scz_brainvar_eqtl_slopes_pdp.jpg",
          eqtl_slope_plot_grid, ncol = 3, nrow = 2,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------

# Finally display the WGCNA enrichment plots 

brainvar_any_cols <- colnames(bip_scz_brainvar_data)[
  str_detect(colnames(bip_scz_brainvar_data), "brainvar_any_gene_")]

# Make an example plot to get the legend with appropriate sizing:
darkgreen_pval_distr <- bip_scz_brainvar_data %>%
  ggplot(aes(x = scz_14_P, 
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
                                              plot_data <- bip_scz_brainvar_data[, c("scz_14_P", x)]
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

save_plot("figures/sf8_scz_brainvar_wgcna.pdf",
          scz_wgcna_module_enrichment_grid, ncol = 5, nrow = 4,
          base_width = 5, base_height = 4)

save_plot("nonpdf_figures/sf8_scz_brainvar_wgcna.jpg",
          scz_wgcna_module_enrichment_grid, ncol = 5, nrow = 4,
          base_width = 5, base_height = 4)






