# PURPOSE: Create the figure displaying the difference in power for selecting
#          the number of CV steps and initial rejection threshold

# AUTHOR: Ron Yurko

#Load necessary packages:
library(tidyverse)
library(cowplot)
library(latex2exp)

# Load the simulation results
sim_s0_n_cv_results <- read_csv("data/si_simulations/s0_n_cv_sims.csv")

# Create the figure showing the difference in power for the different thresholds
# and number of cross-validation steps (only show for alpha <= 0.1)

power_diff_s0_n_cv <- sim_s0_n_cv_results %>%
  filter(alpha <= 0.1 & alpha > .01) %>%
  mutate(alpha = paste0("alpha == ", alpha)) %>%
  group_by(alpha, sim_index, n_cv_steps) %>%
  summarize(disc_power_05_25_diff = disc_power[which(s0_value == .05)] - 
              disc_power[which(s0_value == .25)],
            disc_power_25_45_diff = disc_power[which(s0_value == .25)] - 
              disc_power[which(s0_value == .45)],
            disc_power_05_45_diff = disc_power[which(s0_value == .05)] - 
              disc_power[which(s0_value == .45)]) %>%
  gather(type, disc_power_diff, -alpha, -sim_index, -n_cv_steps) %>%
  ggplot(aes(y = disc_power_diff,
             x = as.factor(n_cv_steps),
             color = type)) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point",
               size = 5, position = position_dodge(width = 0.7)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", alpha = 0.5,
               size = 1, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .3) +
  labs(x = "Number of CV steps",
       color = "Threshold difference",
       y = "Observed difference in power",
       title = "Difference in power between starting thresholds by number of CV steps",
       subtitle = "Points denote averages with intervals for +/- two standard errors") +
  ggthemes::scale_color_colorblind(labels = c("0.05 - 0.25",
                                              "0.05 - 0.45",
                                              "0.25 - 0.45")) +
  facet_wrap(~alpha, labeller = label_parsed,
             ncol = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom") 

# Save
#save_plot("figures/sf25_power_diff_sims_s0_n_cv.pdf",
#          power_diff_s0_n_cv, ncol = 2, nrow = 1, base_width = 4, base_height = 3.5)
#save_plot("nonpdf_figures/sf25_power_diff_sims_s0_n_cv.jpg",
#          power_diff_s0_n_cv, ncol = 2, nrow = 1, base_width = 4, base_height = 3.5)

# Next plot the difference in power for s_0 = 0.05 by the number of CV steps  
power_diff_n_cv_at_s05 <- sim_s0_n_cv_results %>%
  filter(alpha <= 0.1 & alpha > .01,
         s0_value == 0.05) %>%
  group_by(alpha, sim_index) %>%
  summarize(disc_power_2_5_diff = disc_power[which(n_cv_steps == 2)] - 
              disc_power[which(n_cv_steps == 5)],
            disc_power_1_2_diff = disc_power[which(n_cv_steps == 1)] - 
              disc_power[which(n_cv_steps == 2)],
            disc_power_1_5_diff = disc_power[which(n_cv_steps == 1)] - 
              disc_power[which(n_cv_steps == 5)]) %>%
  gather(type, disc_power_diff, -alpha, -sim_index) %>%
  ungroup() %>%
  mutate(alpha = paste0("alpha == ", alpha)) %>%
  ggplot(aes(y = disc_power_diff,
             x = type)) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point",
               size = 5, position = position_dodge(width = 0.7)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", alpha = 0.5,
               size = 1, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .3) +
  labs(x = "CV step difference",
       y = "Observed power difference",
       title = TeX('Difference in power between number of CV steps for $s_0 = 0.05$'),
       subtitle = "Points denote averages with intervals for +/- two standard errors") +
  scale_x_discrete(labels = c("1 - 2",
                              "1 - 5",
                              "2 - 5")) +
  facet_wrap(~alpha, labeller = label_parsed,
             ncol = 3) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom") 

# Save
#save_plot("figures/sf26_power_diff_sims_n_cv_at_s05.pdf",
#          power_diff_n_cv_at_s05, ncol = 2, nrow = 1, base_width = 4, base_height = 3)
#save_plot("nonpdf_figures/sf26_power_diff_sims_n_cv_at_s05.jpg",
#          power_diff_n_cv_at_s05, ncol = 2, nrow = 1, base_width = 4, base_height = 3)

