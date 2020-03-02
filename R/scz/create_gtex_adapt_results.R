# PURPOSE: Generate the AdaPT CV results using the GTEx eSNPs with all variables
#          and then the intercept-only results

# Author: Ron Yurko

# Access necessary packages:
library(tidyverse)
library(future)
# Use the development version of furrr:
# devtools::install_github("DavisVaughan/furrr")
library(furrr)
# Use the modified version of AdaPT:
# devtools::install_github("ryurko/adaptMT")
library(adaptMT)

# Load the data:
bip_scz_gtex_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_gtex_eqtls_cortical_wgcna.csv")

# Also a vector of just BD with the cortical variables:
cortical_variables <- c("z_bip_14", "ave_abs_cortical_eqtl_slope", 
                        "brain_anterior_cingulate_cortex_ba24_ave_abs_eqtl_slope",
                        "brain_frontal_cortex_ba9_ave_abs_eqtl_slope",
                        colnames(bip_scz_gtex_data)[which(str_detect(colnames(bip_scz_gtex_data), "cortical_any_"))])

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

# Generate the AdaPT CV results with only cortical variables
scz_cortical_adapt_cv_results <- adapt_xgboost_cv(as.matrix(
  bip_scz_gtex_data[,cortical_variables]),
  bip_scz_gtex_data$scz_14_P,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = args_search,
  muargs = args_search,
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2),
  s0 = rep(.05, nrow(bip_scz_gtex_data)))

# Save these results:
# saveRDS(scz_cortical_adapt_cv_results,
#        "data/bip_schz_data/gtex_results/cv_tune_results/scz_gtex_cortical_s05_2cv.rds")

# Next generate intercept-only results:
adapt_intercept_results <- adapt_glm(x = bip_scz_gtex_data,
                                     pvals = bip_scz_gtex_data$scz_14_P,
                                     pi_formulas = "1",
                                     mu_formulas = "1",
                                     verbose = list(print = FALSE, 
                                                    fit = FALSE,
                                                    ms = FALSE),
                                     s0 = rep(0.05, nrow(bip_scz_gtex_data)))
# saveRDS(adapt_intercept_results, 
#         "data/bip_schz_data/gtex_results/scz_gtex_intercept_only_s05.rds")


