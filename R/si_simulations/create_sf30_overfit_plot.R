# PURPOSE: Create overfitting simulation figures

# Access packages
library(tidyverse)
library(cowplot)
library(latex2exp)

# Load the simulation results
overfit_sims_data <- read_csv("data/si_simulations/overfit_sims.csv")

overfit_fdp <- overfit_sims_data %>%
  filter(alpha %in% c(0.05, 0.1)) %>%
  mutate(alpha_label = paste0("alpha == ", alpha)) %>%
  ggplot(aes(x = as.factor(n_trees), y = true_fdp)) +
  geom_violin(fill = "gray", alpha = 0.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", alpha = 0.5,
               size = 1, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .3) +
  geom_hline(aes(yintercept = alpha), color = "darkred", linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point", color = "black",
               size = 3) + 
  labs(x = "Number of trees",
       y = "Observed FDP",
       title = TeX('Distribution of simulation FDP by number of trees and target $\\alpha$'),
       subtitle = "Points denote averages with +/- two standard error intervals") +
  facet_wrap(~alpha_label, labeller = label_parsed,
             ncol = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom") 


overfit_power <- overfit_sims_data %>%
  filter(alpha %in% c(0.05, 0.1)) %>%
  mutate(alpha_label = paste0("alpha == ", alpha)) %>%
  ggplot(aes(x = as.factor(n_trees), y = disc_power)) +
  geom_violin(fill = "gray", alpha = 0.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", alpha = 0.5,
               size = 1, position = position_dodge(width = 0.7),
               fun.args = list(mult = 2), width = .3) +
  stat_summary(fun.y = mean, geom = "point", color = "black",
               size = 3) + 
  labs(x = "Number of trees",
       y = "Observed power",
       title = TeX('Distribution of simulation power by number of trees and target $\\alpha$'),
       subtitle = "Points denote averages with +/- two standard error intervals") +
  facet_wrap(~alpha_label, labeller = label_parsed,
             ncol = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom") 
# beautiful!

# Create the grid of plots:
overfit_grid <-
  plot_grid(overfit_fdp, overfit_power,
            labels = c("A", "B"), align = "hv",
            label_fontface = "plain",
            ncol = 1)

# Save
#save_plot("figures/sf30_overfit_sims.pdf",
#          overfit_grid, ncol = 2, nrow = 2)
#save_plot("nonpdf_figures/sf30_overfit_sims.jpg",
#          overfit_grid, ncol = 2, nrow = 2)

