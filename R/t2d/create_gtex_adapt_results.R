# PURPOSE: Generate the AdaPT CV results for T2D using both the observed
#          and adjusted p-values as a comparison. Plus, generate the intercept
#          only results.

# AUTHOR: Ron Yurko

# Access necessary packages:
library(tidyverse)
library(future)
# Use the development version of furrr:
# devtools::install_github("DavisVaughan/furrr")
library(furrr)
# Use the modified version of AdaPT:
# devtools::install_github("ryurko/adaptMT")
library(adaptMT)

# Load the T2D eSNPs:
t2d_esnps_data <- readr::read_csv("data/t2d/t2d_esnps_eqtl_slopes_wgcna.csv")

# What are the Pancreas WGCNA columns:
wgcna_cols <- colnames(t2d_esnps_data)[which(str_detect(colnames(t2d_esnps_data), 
                                                                  "pancreas_any_gene_"))]

# Create the list of arguments to search over
args_search <- list("nrounds100md2" = list("nrounds" = 100,
                                           "max_depth" = 2,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 10),
                    "nrounds150md2" = list("nrounds" = 150,
                                           "max_depth" = 2,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 10),
                    "nrounds100md3" = list("nrounds" = 100,
                                           "max_depth" = 3,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 10),
                    "nrounds150md3" = list("nrounds" = 150,
                                           "max_depth" = 3,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 10))

# Next with the others as well:
t2d_all_variables <- c("pancreas_ave_abs_eqtl_slope", wgcna_cols,
                       "liver_ave_abs_eqtl_slope", 
                       "adipose_ave_abs_eqtl_slope",
                       "adipose_subcutaneous_ave_abs_eqtl_slope",
                       "adipose_visceral_omentum_ave_abs_eqtl_slope")

# Generate the results with the observed p-values
t2d_unadj_adapt <- adapt_xgboost_cv(as.matrix(
  t2d_esnps_data[,t2d_all_variables]),
  t2d_esnps_data$Pvalue,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = args_search,
  muargs = args_search,
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2),
  s0 = rep(0.05, nrow(t2d_esnps_data)))

# Save these results:
# saveRDS(t2d_unadj_adapt,
#         "data/t2d/adapt_cv_results/t2d_unadj_s05_2cv.rds")

# Now the adjusted p-values
t2d_adj_adapt <- adapt_xgboost_cv(as.matrix(
  t2d_esnps_data[,t2d_all_variables]),
  t2d_esnps_data$adj_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = args_search,
  muargs = args_search,
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2),
  s0 = rep(0.05, nrow(t2d_esnps_data)))

# Save these results:
# saveRDS(t2d_adj_adapt,
#        "data/t2d/adapt_cv_results/t2d_adj_s05_2cv.rds")

# Now the intercept-only results:
adapt_intercept_results <- adapt_glm(x = t2d_esnps_data,
                                     pvals = t2d_esnps_data$adj_pval,
                                     pi_formulas = "1",
                                     mu_formulas = "1",
                                     verbose = list(print = FALSE, 
                                                    fit = FALSE,
                                                    ms = FALSE),
                                     s0 = rep(0.05, nrow(t2d_esnps_data)))
# saveRDS(adapt_intercept_results, 
#         "data/t2d/adapt_cv_results/t2d_intercept_only_s05.rds")


