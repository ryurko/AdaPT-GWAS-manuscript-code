# PURPOSE: Create the supplementary figures 17-19 for the T2D results, displaying
#          the enrichment, results with the adjusted p-values, and AdaPT variables

# AUTHOR: Ron Yurko

# Load necessary packages:
library(tidyverse)
library(data.table)
library(cowplot)
library(latex2exp)
library(xgboost)
library(pdp)

# ------------------------------------------------------------------------------

# Load the raw GWAS
t2d_data <- fread("unzip -p data/t2d/Mahajan.NatGenet2018b.T2D.European.zip") 

# Load the T2D eSNPs:
t2d_esnps_data <- readr::read_csv("data/t2d/t2d_esnps_eqtl_slopes_wgcna.csv")

# Create the plot random filtering
t2d_enrichment_points_filter <- t2d_data %>%
  as.data.frame() %>%
  dplyr::select(Pvalue) %>%
  mutate(type = "all",
         Pvalue = as.numeric(Pvalue)) %>%
  bind_rows({(
    t2d_esnps_data %>%
      as.data.frame() %>%
      dplyr::select(Pvalue) %>%
      mutate(type = "gtex") 
  )}) %>%
  group_by(type) %>%
  arrange(Pvalue) %>%
  mutate(neglog_p = -log10(Pvalue),
         exp_neglog_p = -log10(ppoints(n())),
         ppoints_vals = ppoints(n()),
         subset_sample = case_when(
           type == "gtex" & ppoints_vals > .0001 & 
             ppoints_vals <= .001 ~ rbernoulli(n(), p = 0.75),
           type == "gtex" & ppoints_vals > .001 ~ rbernoulli(n(), p = 0.1),
           type == "all" & ppoints_vals > .00001 &
             ppoints_vals <= .0001 ~ rbernoulli(n(), p = 0.05),
           type == "all" & ppoints_vals > .0001 ~ rbernoulli(n(), p = 0.0075),
           TRUE ~ TRUE
         )) %>%
  ungroup() %>%
  filter(subset_sample) %>%
  mutate(type = fct_relevel(type, "gtex", "all")) %>%
  ggplot(aes(x = exp_neglog_p, y = neglog_p,
             color = type)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("darkred",
                                "darkblue"),
                     labels = c("GTEx eSNPs",
                                "All SNPs")) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5,
              linetype = "dotted", color = "black") +
  theme_bw() +
  theme(legend.direction = "vertical",
        legend.position = c(.25, .3),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = TeX('Expected -log$_{10}$(T2D p-value)'),
       y = TeX('Observed -log$_{10}$(T2D p-value)'),
       shape = "SNP type", color = "SNP type") +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Save this plot:
save_plot("figures/sf5_t2d_gtex_enrich.pdf",
          t2d_enrichment_points_filter, ncol = 1, nrow = 1)
save_plot("nonpdf_figures/sf5_t2d_gtex_enrich.jpg",
          t2d_enrichment_points_filter, ncol = 1, nrow = 1)

# Remove the large dataset from the environment:
rm(t2d_data)

# ------------------------------------------------------------------------------
# Next compare the number of discoveries by the adjusted vs not adjusted results

# Load the model results for both the unadjusted and adjusted p-values:
t2d_adj_adapt <- readRDS("data/t2d/adapt_cv_results/t2d_adj_s05_2cv.rds")
t2d_unadj_adapt <- readRDS("data/t2d/adapt_cv_results/t2d_unadj_s05_2cv.rds")

# Create a chart showing the number of discoveries for both the adjusted and 
# unadjusted versions:
t2d_adj_comp_chart <- map_dfr(c(0.01, .05, .1, .15, .2),
                              function(alpha) {
                                data.frame("alpha" = alpha,
                                           "unadj" = length(which(t2d_unadj_adapt$qvals <= alpha)),
                                           "adj" = length(which(t2d_adj_adapt$qvals <= alpha)))
                              }) %>%
  gather(method, n_disc, -alpha) %>%
  mutate(method = fct_relevel(method, "unadj", "adj")) %>%
  ggplot(aes(x = as.factor(alpha), y = n_disc)) +
  geom_bar(aes(fill = method),
           width = 0.1, stat = "identity", position = position_dodge(width = 0.6)) +
  geom_point(aes(color = method),
             size = 4, position = position_dodge(width = 0.6)) +
  theme_bw() +
  ggthemes::scale_fill_colorblind(labels = c("unadjusted", "adjusted")) +
  ggthemes::scale_color_colorblind(labels = c("unadjusted", "adjusted")) +
  labs(x = TeX('$\\alpha$'), y = "Number of discoveries",
       fill = "p-value type", color = "p-value type",
       title = "Comparison of the number of T2D discoveries\nwith or without adjustment to p-values") +
  theme(legend.position = c(.175, .7))

save_plot("figures/sf18_t2d_adj_pval_comp.pdf",
          t2d_adj_comp_chart, ncol = 1, nrow = 1, 
          base_width = 5, base_height = 3)
save_plot("nonpdf_figures/sf18_t2d_adj_pval_comp.jpg",
          t2d_adj_comp_chart, ncol = 1, nrow = 1, 
          base_width = 5, base_height = 3)

# ------------------------------------------------------------------------------

# Finally create the variable importance plots for T2D adjusted models

# First create the \pi_1 variable importance plot
t2d_pi_model_i <- seq(1, length(t2d_adj_adapt$model_fit), by = 2)
t2d_pi_search_importance_data <- map_dfr(1:length(t2d_pi_model_i),
                                         function(step_i) {
                                           pi_model_step_i <- t2d_pi_model_i[step_i]
                                           xgb.importance(model = t2d_adj_adapt$model_fit[[pi_model_step_i]]) %>%
                                             as.data.frame() %>%
                                             mutate(adapt_step = step_i)
                                         })
# Will highlight the top variables from the final model:
t2d_final_top_pi_importance <- t2d_pi_search_importance_data %>%
  filter(adapt_step == length(t2d_pi_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:5)

# Now create a display of the importance over the search with the top variable
# highlighted:
t2d_pi_importance_search_plot <- t2d_pi_search_importance_data %>%
  filter(Feature %in% t2d_final_top_pi_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
           str_replace_all("ave abs ",
                           "average |") %>%
           str_replace_all("slope", "beta|") %>%
           str_replace_all("pancreas any gene", "Pancreas WGCNA module:")) %>%
  ggplot() +
  geom_line(data = {(
    t2d_pi_search_importance_data %>%
      filter(!(Feature %in% t2d_final_top_pi_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    t2d_pi_search_importance_data %>%
      filter(Feature %in% t2d_final_top_pi_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_replace_all("ave abs ",
                               "average |") %>%
               str_replace_all("slope", "beta|") %>%
               str_replace_all("pancreas any gene", "Pancreas WGCNA module:")) %>%
      filter(adapt_step == length(t2d_pi_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_manual(values = c(
    rep("goldenrod4", 2),
    "darkblue", "darkred", "turquoise")) +  
  scale_linetype_manual(guide = FALSE,
                        values = c(c("dotdash", "dotted"), rep("solid", 3))) +
  scale_x_continuous(limits = c(1, 20),
                     breaks = seq(1, length(t2d_pi_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\pi_1$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Next repeat but for the \mu models:
t2d_mu_model_i <- seq(2, length(t2d_adj_adapt$model_fit), by = 2)
t2d_mu_search_importance_data <- map_dfr(1:length(t2d_mu_model_i),
                                         function(step_i) {
                                           mu_model_step_i <- t2d_mu_model_i[step_i]
                                           xgb.importance(model = t2d_adj_adapt$model_fit[[mu_model_step_i]]) %>%
                                             as.data.frame() %>%
                                             mutate(adapt_step = step_i)
                                         })
# Will highlight the top variables from the final model:
t2d_final_top_mu_importance <- t2d_mu_search_importance_data %>%
  filter(adapt_step == length(t2d_mu_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:5)
t2d_mu_importance_search_plot <- t2d_mu_search_importance_data %>%
  filter(Feature %in% t2d_final_top_mu_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
           str_replace_all("ave abs ",
                           "average |") %>%
           str_replace_all("slope", "beta|") %>%
           str_replace_all("pancreas any gene", "Pancreas WGCNA module:")) %>%
  ggplot() +
  geom_line(data = {(
    t2d_mu_search_importance_data %>%
      filter(!(Feature %in% t2d_final_top_mu_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    t2d_mu_search_importance_data %>%
      filter(Feature %in% t2d_final_top_mu_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_replace_all("ave abs ",
                               "average |") %>%
               str_replace_all("slope", "beta|") %>%
               str_replace_all("pancreas any gene", "Pancreas WGCNA module:")) %>%
      filter(adapt_step == length(t2d_mu_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_manual(values = c(
    rep("goldenrod4", 3),
    "darkblue", "darkred")) +  
  scale_linetype_manual(guide = FALSE,
                        values = c(c("dotdash", "dotted", "dashed"), rep("solid", 2))) +
  scale_x_continuous(limits = c(1, 20),
                     breaks = seq(1, length(t2d_mu_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\mu$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Arrange these in a grid and save:
t2d_var_imp_grid <-
  plot_grid(t2d_pi_importance_search_plot,
            t2d_mu_importance_search_plot, ncol = 1,
            align = "hv", labels = c("A", "B"),
            label_fontface = "plain")
save_plot("figures/sf19_t2d_var_importance.pdf",
          t2d_var_imp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf19_t2d_var_importance.jpg",
          t2d_var_imp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)




