# PURPOSE: Create the supplementary figures 20-24 for the BMI results, displaying
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

# First load both the full BMI data and the eSNPs
bmi_data <- read_csv("data/bmi/bmi_18_new_studies.csv",
                     progress = FALSE)
bmi_esnps_data <- read_csv("data/bmi/bmi_esnps_eqtl_slopes_wgcna.csv")

# Create the plot random filtering
bmi_enrichment_points_filter <- bmi_data %>%
  as.data.frame() %>%
  dplyr::select(bmi15_p) %>%
  mutate(type = "all") %>%
  bind_rows({(
    bmi_esnps_data %>%
      as.data.frame() %>%
      dplyr::select(bmi15_p) %>%
      mutate(type = "gtex") 
  )}) %>%
  group_by(type) %>%
  arrange(bmi15_p) %>%
  mutate(neglog_p = -log10(bmi15_p),
         exp_neglog_p = -log10(ppoints(n())),
         ppoints_vals = ppoints(n()),
         subset_sample = case_when(
           type == "gtex" & ppoints_vals > .001 & 
             ppoints_vals <= .01 ~ rbernoulli(n(), p = 0.75),
           type == "gtex" & ppoints_vals > .001 ~ rbernoulli(n(), p = 0.1),
           type == "all" & ppoints_vals > .0001 &
             ppoints_vals <= .001 ~ rbernoulli(n(), p = 0.05),
           type == "all" & ppoints_vals > .001 ~ rbernoulli(n(), p = 0.025),
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
  labs(x = TeX('Expected -log$_{10}$(BMI p-value)'),
       y = TeX('Observed -log$_{10}$(BMI p-value)'),
       shape = "SNP type", color = "SNP type") +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Save this plot:
save_plot("figures/sf20_bmi_gtex_enrich.pdf",
          bmi_enrichment_points_filter, ncol = 1, nrow = 1)
save_plot("nonpdf_figures/sf20_bmi_gtex_enrich.jpg",
          bmi_enrichment_points_filter, ncol = 1, nrow = 1)

# ------------------------------------------------------------------------------
# Next compare the number of discoveries for the adjusted and unadjusted
# p-values BMI AdaPT results

bmi_adj_adapt <- readRDS("data/bmi/adapt_cv_results/bmi_adj15_s05_2cv.rds")
bmi_unadj_adapt <- readRDS("data/bmi/adapt_cv_results/bmi_unadj15_s05_2cv.rds")

# Create a chart showing the number of discoveries for both the adjusted and 
# unadjusted versions:
bmi_adj_comp_chart <- map_dfr(c(0.01, .05, .1, .15, .2),
                              function(alpha) {
                                data.frame("alpha" = alpha,
                                           "unadj" = length(which(bmi_unadj_adapt$qvals <= alpha)),
                                           "adj" = length(which(bmi_adj_adapt$qvals <= alpha)))
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
       title = "Comparison of the number of BMI discoveries\nwith or without adjustment to p-values") +
  theme(legend.position = c(.175, .7))

save_plot("figures/sf21_bmi_adj_pval_comp.pdf",
          bmi_adj_comp_chart, ncol = 1, nrow = 1, 
          base_width = 5, base_height = 3)
save_plot("nonpdf_figures/sf21_bmi_adj_pval_comp.jpg",
          bmi_adj_comp_chart, ncol = 1, nrow = 1, 
          base_width = 5, base_height = 3)

# ------------------------------------------------------------------------------

# Next create the BMI replication chart

# What is the nominal replication rate?
length(which(bmi_esnps_data$bmi18_new_p[which(bmi_adj_adapt$qvals <= 0.05)] <= 0.05))
# [1] 1150

# So the nominal replication rate is just:
length(which(bmi_esnps_data$bmi18_new_p[which(bmi_adj_adapt$qvals <= 0.05)] <= 0.05)) /
  length(which(bmi_adj_adapt$qvals <= 0.05))
# [1] 0.8315257

# First create a chart of the relationship between the q-values from this model
# and the new 2018 p-values:
bmi_gtex_adapt_qvals_new_pvals <- 
  tibble(adapt_qvals = bmi_adj_adapt$qvals,
         updated_scz_pvals = bmi_esnps_data$bmi18_new_p,
         adapt_disc = as.numeric(bmi_adj_adapt$qvals <= .05)) %>%
  mutate(
    # Impute 1 for non-finite values since they were never in the rejection set:
    adapt_qvals = ifelse(is.infinite(adapt_qvals), 1, adapt_qvals)) %>%
  ggplot(aes(x = -log10(adapt_qvals),
             y = -log10(updated_scz_pvals))) +
  geom_point(alpha = 0.35, color = "gray") + 
  geom_segment(x = -log10(.05), xend = -log10(.05),
               y = 0, yend = Inf, linetype = "dotted", color = "darkblue") +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "darkred") +
  annotate("rect", xmin = -log10(.05), xmax = Inf,
           ymin = -log10(.05), ymax = Inf, alpha = 0.15,
           fill = "darkblue") +
  annotate("rect", xmin = -log10(.05), xmax = Inf,
           ymax = -log10(.05), ymin = 0, alpha = 0.15,
           fill = "darkred") +
  annotate("label", x = 2.15, y = 50,
           label = TeX('AdaPT discoveries at $\\alpha = 0.05$'),
           color = "darkblue", fill = "white", alpha = 0.75) +
  annotate("label", x = 2.15, y = -log10(.15),
           label = TeX("Replication threshold: -log$_{10}(0.05)$"),
           color = "darkred", alpha = 0.75, fill = "white") +
  annotate("label", x = 1.15, y = -log10(10^(-9)),
           label = "Smooth relationship between new BMI p-values and AdaPT q-values",
           color = "black", fill = "white", alpha = 0.75) +
  geom_smooth(color = "black", se = FALSE) +
  theme_bw() +
  labs(x = TeX('-log_{10}(AdaPT q-values from 2015 studies)'),
       y = TeX('-log_{10}(BMI p-values from new 2018 studies)'))

save_plot(bmi_gtex_adapt_qvals_new_pvals,
          filename = "figures/sf22_bmi_gtex_replication.pdf",
          base_aspect_ratio = 2.5)
save_plot(bmi_gtex_adapt_qvals_new_pvals,
          filename = "nonpdf_figures/sf22_bmi_gtex_replication.jpg",
          base_aspect_ratio = 2.5)

# ------------------------------------------------------------------------------

# Next create the variable importance and partial dependence plots for the 
# BMI AdaPT results (using the adjusted p-values)


# First create the \pi_1 variable importance plot
bmi_pi_model_i <- seq(1, length(bmi_adj_adapt$model_fit), by = 2)
bmi_pi_search_importance_data <- map_dfr(1:length(bmi_pi_model_i),
                                         function(step_i) {
                                           pi_model_step_i <- bmi_pi_model_i[step_i]
                                           xgb.importance(model = bmi_adj_adapt$model_fit[[pi_model_step_i]]) %>%
                                             as.data.frame() %>%
                                             mutate(adapt_step = step_i)
                                         })
# Will highlight the top variables from the final model:
bmi_final_top_pi_importance <- bmi_pi_search_importance_data %>%
  filter(adapt_step == length(bmi_pi_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:5)

# Now create a display of the importance over the search with the top variable
# highlighted:
bmi_pi_importance_search_plot <- bmi_pi_search_importance_data %>%
  filter(Feature %in% bmi_final_top_pi_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
           str_replace_all("ave abs ",
                           "average |") %>%
           str_replace_all("slope", "beta|") %>%
           str_replace_all("cerebellum any gene", "Cerebellar hemisphere WGCNA module:") %>%
           str_replace_all("brain any gene", "Non-cerebellar hemisphere WGCNA module:") %>%
           str_replace_all("adipose any gene", "Adipose WGCNA module:") %>%
           str_replace_all("whr15 z score", "WHR z-statistics")) %>%
  ggplot() +
  geom_line(data = {(
    bmi_pi_search_importance_data %>%
      filter(!(Feature %in% bmi_final_top_pi_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    bmi_pi_search_importance_data %>%
      filter(Feature %in% bmi_final_top_pi_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_replace_all("ave abs ",
                               "average |") %>%
               str_replace_all("slope", "beta|") %>%
               str_replace_all("cerebellum any gene", "Cerebellar hemisphere WGCNA module:") %>%
               str_replace_all("brain any gene", "Non-cerebellar hemisphere WGCNA module:") %>%
               str_replace_all("adipose any gene", "Adipose WGCNA module:") %>%
               str_replace_all("whr15 z score", "WHR z-statistics")) %>%
      filter(adapt_step == length(bmi_pi_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_manual(values = c(
    "darkred", rep("goldenrod4", 2),
    "turquoise", "darkblue")) +  
  scale_linetype_manual(guide = FALSE,
                        values = c("solid", c("dotdash", "dotted"), rep("solid", 2))) +
  scale_x_continuous(limits = c(1, 25),
                     breaks = seq(1, length(bmi_pi_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\pi_1$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Next repeat but for the \mu models:
bmi_mu_model_i <- seq(2, length(bmi_adj_adapt$model_fit), by = 2)
bmi_mu_search_importance_data <- map_dfr(1:length(bmi_mu_model_i),
                                         function(step_i) {
                                           mu_model_step_i <- bmi_mu_model_i[step_i]
                                           xgb.importance(model = bmi_adj_adapt$model_fit[[mu_model_step_i]]) %>%
                                             as.data.frame() %>%
                                             mutate(adapt_step = step_i)
                                         })
# Will highlight the top variables from the final model:
bmi_final_top_mu_importance <- bmi_mu_search_importance_data %>%
  filter(adapt_step == length(bmi_mu_model_i)) %>%
  arrange(desc(Gain)) %>%
  dplyr::slice(1:5)
bmi_mu_importance_search_plot <- bmi_mu_search_importance_data %>%
  filter(Feature %in% bmi_final_top_mu_importance$Feature) %>%
  mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
           str_replace_all("ave abs ",
                           "average |") %>%
           str_replace_all("slope", "beta|") %>%
           str_replace_all("cerebellum any gene", "Cerebellar hemisphere WGCNA module:") %>%
           str_replace_all("brain any gene", "Non-cerebellar hemisphere WGCNA module:") %>%
           str_replace_all("adipose any gene", "Adipose WGCNA module:") %>%
           str_replace_all("whr15 z score", "WHR z-statistics")) %>%
  ggplot() +
  geom_line(data = {(
    bmi_mu_search_importance_data %>%
      filter(!(Feature %in% bmi_final_top_mu_importance$Feature))
  )}, aes(x = adapt_step, y = Gain, 
          group = Feature), color = "gray90", alpha = 0.5, size = 1) +
  geom_line(aes(x = adapt_step, y = Gain, 
                color = feature_label,
                group = feature_label,
                linetype = feature_label),
            size = 1.5, alpha = 0.8) +
  geom_label_repel(data = {(
    bmi_mu_search_importance_data %>%
      filter(Feature %in% bmi_final_top_mu_importance$Feature) %>%
      mutate(feature_label = str_replace_all(Feature, "_", " ") %>% 
               str_replace_all("ave abs ",
                               "average |") %>%
               str_replace_all("slope", "beta|") %>%
               str_replace_all("cerebellum any gene", "Cerebellar hemisphere WGCNA module:") %>%
               str_replace_all("brain any gene", "Non-cerebellar hemisphere WGCNA module:") %>%
               str_replace_all("adipose any gene", "Adipose WGCNA module:") %>%
               str_replace_all("whr15 z score", "WHR z-statistics")) %>%
      filter(adapt_step == length(bmi_mu_model_i))
  )}, aes(x = adapt_step, y = Gain, label = feature_label,
          color = feature_label),
  direction = "y", nudge_x = .75, segment.size = 0.05,
  hjust = 0) +
  scale_color_manual(values = c(
    rep("goldenrod4", 3),
    "darkred", "darkblue")) +  
  scale_linetype_manual(guide = FALSE,
                        values = c(c("dotdash", "dotted", "dashed"), rep("solid", 2))) +
  scale_x_continuous(limits = c(1, 25),
                     breaks = seq(1, length(bmi_mu_model_i),
                                  by = 1)) +
  theme_cowplot() +
  labs(x = "AdaPT model fitting iteration",
       y = "Importance",
       color = "Variable",
       title = TeX('Change in variable importance across $\\mu$ models in AdaPT search with top variables in final model highlighted')) +
  theme(legend.position = "none")

# Arrange these in a grid and save:
bmi_var_imp_grid <-
  plot_grid(bmi_pi_importance_search_plot,
            bmi_mu_importance_search_plot, ncol = 1,
            align = "hv", labels = c("A", "B"),
            label_fontface = "plain")
save_plot("figures/sf23_bmi_var_importance.pdf",
          bmi_var_imp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf23_bmi_var_importance.jpg",
          bmi_var_imp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------

# Next create the partial dependence plots for the WHR z-statistics:
bmi_pi_quantiles_search_pdp_data <- map_dfr(1:length(bmi_pi_model_i),
                                            function(step_i) {
                                              pi_model_step_i <- bmi_pi_model_i[step_i]
                                              pdp::partial(bmi_adj_adapt$model_fit[[pi_model_step_i]], 
                                                           pred.var = "whr15_z_score", 
                                                           ice = FALSE, 
                                                           prob = TRUE,
                                                           center = FALSE, 
                                                           plot = FALSE, 
                                                           quantiles = TRUE,
                                                           probs = seq(0, 1, by = .025),
                                                           train = data.matrix(
                                                             bmi_esnps_data[,bmi_adj_adapt$model_fit[[1]]$feature_names])) %>%
                                                as.data.frame() %>%
                                                mutate(adapt_step = step_i)
                                            })
bmi_pi_quantiles_pdp_search_plot <- bmi_pi_quantiles_search_pdp_data %>%
  ggplot() +
  geom_line(aes(x = whr15_z_score, y = yhat, 
                color = adapt_step,
                group = as.factor(adapt_step))) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_cowplot() +
  scale_color_gradient(low = "darkblue", high = "darkorange") +
  geom_rug(data = bmi_esnps_data,
           aes(x = whr15_z_score), y = rep(1, nrow(bmi_esnps_data)),
           sides = "b",
           color = "black", alpha = 0.15) +
  geom_vline(xintercept = 1.96, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = -1.96, linetype = "dashed", color = "darkred") +
  guides(color = guide_colourbar(barwidth = 10, barheight = .5)) +
  labs(x = "WHR z-statistic",
       y = TeX('$\\pi_1$'),
       color = "AdaPT model fitting iteration",
       subtitle = "Dashed red lines indicate z-statistics equal to +/- 1.96",
       title = TeX('Change in partial dependence for $\\pi_1$ and WHR z-statistics')) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# Next the change in non-null effect size pdp
bmi_mu_quantiles_search_pdp_data <- map_dfr(1:length(bmi_mu_model_i),
                                            function(step_i) {
                                              mu_model_step_i <- bmi_mu_model_i[step_i]
                                              pdp::partial(bmi_adj_adapt$model_fit[[mu_model_step_i]], 
                                                           pred.var = "whr15_z_score", 
                                                           ice = FALSE, 
                                                           center = FALSE, 
                                                           plot = FALSE, 
                                                           quantiles = TRUE,
                                                           probs = seq(0, 1, by = .025),
                                                           train = data.matrix(
                                                             bmi_esnps_data[,bmi_adj_adapt$model_fit[[1]]$feature_names])) %>%
                                                as.data.frame() %>%
                                                mutate(adapt_step = step_i)
                                            })

bmi_mu_quantiles_pdp_search_plot <- bmi_mu_quantiles_search_pdp_data %>%
  ggplot() +
  geom_line(aes(x = whr15_z_score, y = yhat, 
                color = adapt_step,
                group = as.factor(adapt_step))) +
  scale_y_continuous(limits = c(0, 6)) +
  theme_cowplot() +
  scale_color_gradient(low = "darkblue", high = "darkorange") +
  geom_rug(data = bmi_esnps_data,
           aes(x = whr15_z_score), y = rep(1, nrow(bmi_esnps_data)),
           sides = "b",
           color = "black", alpha = 0.15) +
  geom_vline(xintercept = 1.96, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = -1.96, linetype = "dashed", color = "darkred") +
  guides(color = guide_colourbar(barwidth = 10, barheight = .5)) +
  labs(x = "WHR z-statistic",
       y = TeX('$\\mu$'),
       color = "AdaPT model fitting iteration",
       subtitle = "Dashed red lines indicate z-statistics equal to +/- 1.96",
       title = TeX('Change in partial dependence for $\\mu$ and WHR z-statistics')) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# Save the plot together:
bmi_whr_pdp_grid <- plot_grid(
  plot_grid(bmi_pi_quantiles_pdp_search_plot + theme(legend.position = "none"),
            bmi_mu_quantiles_pdp_search_plot + theme(legend.position = "none"), 
            ncol = 1,
            align = "hv", labels = c("A", "B"),
            label_fontface = "plain"),
  get_legend(bmi_pi_quantiles_pdp_search_plot), 
  ncol = 1, rel_heights = c(2, .5))
save_plot("figures/sf24_bmi_whr_pdp.pdf",
          bmi_whr_pdp_grid, ncol = 1, nrow = 2,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf24_bmi_whr_pdp.jpg",
          bmi_whr_pdp_grid, ncol = 2, nrow = 2,
          base_width = 6, base_height = 4)




