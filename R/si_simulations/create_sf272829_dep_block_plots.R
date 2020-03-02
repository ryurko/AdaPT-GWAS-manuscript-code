# PURPOSE: Generate the figures summarizing the FDP and power of the dependent
#          block simulations

# AUTHOR: Ron Yurko

# Access packages
library(tidyverse)
library(cowplot)
library(latex2exp)

# Load the simulation results
sim_results <- read_csv("data/si_simulations/dep_sim_alt_blocks.csv") %>%
  mutate(fdp = ifelse(is.na(fdp), 0, fdp),
         fdp_blocks = ifelse(is.na(fdp_blocks), 0 ,fdp_blocks),
         fdp_test_blocks = ifelse(is.na(fdp_test_blocks), 0, fdp_test_blocks))

# First define the function that will be used for displaying the power and FDP
# facet labels correctly:
# taken from https://stackoverflow.com/questions/14181234/facet-labels-involving-a-greek-symbol
my_label_bquote <- function(expr1 = (beta[1] == .(x)),
                            expr2 = (gamma[1] == .(x))) {
  quoted1<- substitute(expr1)
  quoted2 <- substitute(expr2)
  function(variable, value) {
    value <- as.character(value)
    if(variable == 'beta_1')
      lapply(value, function(x)
        eval(substitute(bquote(expr1, list(x = x)),list(expr1 = quoted1))))
    else
      lapply(value, function(x) 
        eval(substitute(bquote(expr2, list(x = x)),list(expr2 = quoted2))))
  }
}


# Just focus on settings with non-zero coefficients
# ---------------------------------
# mu_floor == 0.5

# power plot
power_mu5_plot <- sim_results %>%
  filter(alpha %in% c(0.05),
         mu_floor == .5,
         pi_beta1 == pi_beta2,
         pi_beta1 > 0,
         mu_beta1 == mu_beta2,
         mu_beta1 > 0) %>%
  mutate(method = fct_relevel(method, "bh", "adapt_int_only",
                              "adapt_regular"),
         beta_1 = pi_beta1,
         gamma_1 = mu_beta1) %>% 
  ggplot(aes(x = as.factor(block_rho), y = power, 
             fill = method, color = method)) +
  stat_summary(fun.y = mean, geom = "point",
               size = 1, position = position_dodge(width = 0.7)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               size = .5, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .75) + 
  labs(x = TeX('Block correlation $\\rho$'),
       y = "Observed power",
       color = "Method", fill = "Method") +
  ggthemes::scale_fill_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  ggthemes::scale_color_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  facet_grid(beta_1 ~ gamma_1, scales = "free_y", labeller = my_label_bquote()) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank()) 

# FDP plot
fdp_mu5_plot <- sim_results %>%
  filter(alpha %in% c(0.05),
         mu_floor == .5,
         pi_beta1 == pi_beta2,
         pi_beta1 > 0,
         mu_beta1 == mu_beta2,
         mu_beta1 > 0) %>%
  mutate(method = fct_relevel(method, "bh", "adapt_int_only",
                              "adapt_regular"),
         beta_1 = pi_beta1,
         gamma_1 = mu_beta1) %>% 
  ggplot(aes(x = as.factor(block_rho), y = fdp, 
             fill = method, color = method)) +
  geom_hline(aes(yintercept = alpha), color = "darkred", linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point",
               size = 1, position = position_dodge(width = 0.7)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               size = .5, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .75) + 
  labs(x = TeX('Block correlation $\\rho$'),
       y = "Observed FDP",
       color = "Method", fill = "Method") +
  ggthemes::scale_fill_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  ggthemes::scale_color_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  facet_grid(beta_1 ~ gamma_1, scales = "free_y", labeller = my_label_bquote()) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank()) 


# ---------------------------------
# mu_floor == 1

# power plot
power_mu1_plot <- sim_results %>%
  filter(alpha %in% c(0.05),
         mu_floor == 1,
         pi_beta1 == pi_beta2,
         pi_beta1 > 0,
         mu_beta1 == mu_beta2,
         mu_beta1 > 0) %>%
  mutate(method = fct_relevel(method, "bh", "adapt_int_only",
                              "adapt_regular"),
         beta_1 = pi_beta1,
         gamma_1 = mu_beta1) %>% 
  ggplot(aes(x = as.factor(block_rho), y = power, 
             fill = method, color = method)) +
  stat_summary(fun.y = mean, geom = "point",
               size = 1, position = position_dodge(width = 0.7)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               size = .5, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .75) + 
  labs(x = TeX('Block correlation $\\rho$'),
       y = "Observed power",
       color = "Method", fill = "Method") +
  ggthemes::scale_fill_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  ggthemes::scale_color_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  facet_grid(beta_1 ~ gamma_1, scales = "free_y", labeller = my_label_bquote()) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank()) 

# FDP plot
fdp_mu1_plot <- sim_results %>%
  filter(alpha %in% c(0.05),
         mu_floor == 1,
         pi_beta1 == pi_beta2,
         pi_beta1 > 0,
         mu_beta1 == mu_beta2,
         mu_beta1 > 0) %>%
  mutate(method = fct_relevel(method, "bh", "adapt_int_only",
                              "adapt_regular"),
         beta_1 = pi_beta1,
         gamma_1 = mu_beta1) %>% 
  ggplot(aes(x = as.factor(block_rho), y = fdp, 
             fill = method, color = method)) +
  geom_hline(aes(yintercept = alpha), color = "darkred", linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point",
               size = 1, position = position_dodge(width = 0.7)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               size = .5, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .75) + 
  labs(x = TeX('Block correlation $\\rho$'),
       y = "Observed FDP",
       color = "Method", fill = "Method") +
  ggthemes::scale_fill_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  ggthemes::scale_color_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  facet_grid(beta_1 ~ gamma_1, scales = "free_y", labeller = my_label_bquote()) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank()) 


# ---------------------------------
# mu_floor == 1.5

# power plot
power_mu15_plot <- sim_results %>%
  filter(alpha %in% c(0.05),
         mu_floor == 1.5,
         pi_beta1 == pi_beta2,
         pi_beta1 > 0,
         mu_beta1 == mu_beta2,
         mu_beta1 > 0) %>%
  mutate(method = fct_relevel(method, "bh", "adapt_int_only",
                              "adapt_regular"),
         beta_1 = pi_beta1,
         gamma_1 = mu_beta1) %>% 
  ggplot(aes(x = as.factor(block_rho), y = power, 
             fill = method, color = method)) +
  stat_summary(fun.y = mean, geom = "point",
               size = 1, position = position_dodge(width = 0.7)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               size = .5, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .75) + 
  labs(x = TeX('Block correlation $\\rho$'),
       y = "Observed power",
       color = "Method", fill = "Method") +
  ggthemes::scale_fill_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  ggthemes::scale_color_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  facet_grid(beta_1 ~ gamma_1, scales = "free_y", labeller = my_label_bquote()) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank()) 

# FDP plot
fdp_mu15_plot <- sim_results %>%
  filter(alpha %in% c(0.05),
         mu_floor == 1.5,
         pi_beta1 == pi_beta2,
         pi_beta1 > 0,
         mu_beta1 == mu_beta2,
         mu_beta1 > 0) %>%
  mutate(method = fct_relevel(method, "bh", "adapt_int_only",
                              "adapt_regular"),
         beta_1 = pi_beta1,
         gamma_1 = mu_beta1) %>% 
  ggplot(aes(x = as.factor(block_rho), y = fdp, 
             fill = method, color = method)) +
  geom_hline(aes(yintercept = alpha), color = "darkred", linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point",
               size = 1, position = position_dodge(width = 0.7)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               size = .5, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .75) + 
  labs(x = TeX('Block correlation $\\rho$'),
       y = "Observed FDP",
       color = "Method", fill = "Method") +
  ggthemes::scale_fill_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  ggthemes::scale_color_colorblind(labels = c("BH", "AdaPT: intercept-only", "AdaPT: x1 + x2")) +
  facet_grid(beta_1 ~ gamma_1, scales = "free_y", labeller = my_label_bquote()) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank()) 

# ----------------------------------------------------------------------

# Next create the grid of these plots - just stack them two at a at time, 1 column
# this is the supplement after so who cares how much space they take up

# First mu_floor 0.5
mu5_grid <- plot_grid(plot_grid(fdp_mu5_plot + theme(legend.position = "none"),
                                power_mu5_plot + theme(legend.position = "none"), 
                                ncol = 1, align = "hv", 
                                label_fontface = "plain", labels = c("A", "B")),
                      get_legend(power_mu5_plot), ncol = 1, rel_heights = c(2, .25))

save_plot("figures/sf27_dep_pval_sims_05.pdf", 
          mu5_grid, ncol = 1, nrow = 2)
save_plot("nonpdf_figures/sf27_dep_pval_sims_05.jpg", 
          mu5_grid, ncol = 1, nrow = 2)


# Next for mu_floor 1
mu1_grid <- plot_grid(plot_grid(fdp_mu1_plot + theme(legend.position = "none"),
                                power_mu1_plot + theme(legend.position = "none"), 
                                ncol = 1, align = "hv", 
                                label_fontface = "plain", labels = c("A", "B")),
                      get_legend(power_mu1_plot), ncol = 1, rel_heights = c(2, .25))

save_plot("figures/sf28_dep_pval_sims_1.pdf", 
          mu1_grid, ncol = 1, nrow = 2)
save_plot("nonpdf_figures/sf28_dep_pval_sims_1.jpg", 
          mu1_grid, ncol = 1, nrow = 2)

# Next for mu_floor 1.5
mu15_grid <- plot_grid(plot_grid(fdp_mu15_plot + theme(legend.position = "none"),
                                 power_mu15_plot + theme(legend.position = "none"), 
                                 ncol = 1, align = "hv", 
                                 label_fontface = "plain", labels = c("A", "B")),
                       get_legend(power_mu15_plot), ncol = 1, rel_heights = c(2, .25))

save_plot("figures/sf29_dep_pval_sims_15.pdf", 
          mu15_grid, ncol = 1, nrow = 2)
save_plot("nonpdf_figures/sf29_dep_pval_sims_15.jpg", 
          mu15_grid, ncol = 1, nrow = 2)
