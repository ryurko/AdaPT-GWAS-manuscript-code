# PURPOSE: Generate the AdaPT CV results for BMI using different sets of variables
#          along with intercept-only results, and for all 2018 studies.
#          Additionally, need to generate the results using both the observed
#          and adjusted p-values as a comparison

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

# Now load the complete model data with the any WGCNA results:
bmi_esnps_data <- readr::read_csv("data/bmi/bmi_esnps_eqtl_slopes_wgcna.csv")

# Create a vector of average slopes:
ave_slope_vars <-  c("adipose_subcutaneous_ave_abs_eqtl_slope",
                     "adipose_visceral_omentum_ave_abs_eqtl_slope",
                     "brain_amygdala_ave_abs_eqtl_slope",
                     "brain_anterior_cingulate_cortex_ba24_ave_abs_eqtl_slope",
                     "brain_caudate_basal_ganglia_ave_abs_eqtl_slope",
                     "brain_cerebellar_hemisphere_ave_abs_eqtl_slope",
                     "brain_frontal_cortex_ba9_ave_abs_eqtl_slope",
                     "brain_hippocampus_ave_abs_eqtl_slope",
                     "brain_hypothalamus_ave_abs_eqtl_slope", 
                     "brain_nucleus_accumbens_basal_ganglia_ave_abs_eqtl_slope",
                     "brain_putamen_basal_ganglia_ave_abs_eqtl_slope",
                     "brain_spinal_cord_cervical_c-1_ave_abs_eqtl_slope",
                     "brain_substantia_nigra_ave_abs_eqtl_slope",
                     "adipose_ave_abs_eqtl_slope",
                     "brain_ave_abs_eqtl_slope")

# Get the WGCNA module variables:
wgcna_vars <- c(colnames(bmi_esnps_data)[which(str_detect(colnames(bmi_esnps_data), "adipose_any_"))],
                colnames(bmi_esnps_data)[which(str_detect(colnames(bmi_esnps_data), "brain_any_"))],
                colnames(bmi_esnps_data)[which(str_detect(colnames(bmi_esnps_data), "cerebellum_any_"))])

# Now a vector of all variables considered:
bmi_variables <- c("whr15_z_score", ave_slope_vars, wgcna_vars)

args_search <- list("nrounds100md2" = list("nrounds" = 100,
                                           "max_depth" = 2,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 5),
                    "nrounds150md2" = list("nrounds" = 150,
                                           "max_depth" = 2,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 5),
                    "nrounds100md3" = list("nrounds" = 100,
                                           "max_depth" = 3,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 5),
                    "nrounds150md3" = list("nrounds" = 150,
                                           "max_depth" = 3,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 5))

# Generate the AdaPT CV results for BMI unajdusted p-values
bmi_adapt_cv_unadj_results <- adapt_xgboost_cv(as.matrix(
  bmi_esnps_data[,bmi_variables]),
  bmi_esnps_data$bmi15_p,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = args_search,
  muargs = args_search,
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2),
  s0 = rep(0.05, nrow(bmi_esnps_data)))

# Save these results:
# saveRDS(bmi_adapt_cv_unadj_results,
#         "data/bmi/adapt_cv_results/bmi_unadj15_s05_2cv.rds")
# 
# Now the adjusted results
bmi_adapt_cv_adj_results <- adapt_xgboost_cv(as.matrix(
  bmi_esnps_data[,bmi_variables]),
  bmi_esnps_data$adj_bmi15_p,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = args_search,
  muargs = args_search,
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2),
  s0 = rep(0.05, nrow(bmi_esnps_data)))

# Save these results:
# saveRDS(bmi_adapt_cv_adj_results,
#         "data/bmi/adapt_cv_results/bmi_adj15_s05_2cv.rds")
# 
whr_only_args_search <- list("nrounds100md1" = list("nrounds" = 100,
                                                    "max_depth" = 1,
                                                    "min_child_weight" = 1,
                                                    "verbose" = 0,
                                                    "nthread" = 1),
                             "nrounds150md1" = list("nrounds" = 150,
                                                    "max_depth" = 1,
                                                    "min_child_weight" = 1,
                                                    "verbose" = 0,
                                                    "nthread" = 1),
                             "nrounds100md6" = list("nrounds" = 100,
                                                    "max_depth" = 6,
                                                    "min_child_weight" = 1,
                                                    "verbose" = 0,
                                                    "nthread" = 1),
                             "nrounds150md6" = list("nrounds" = 150,
                                                    "max_depth" = 6,
                                                    "min_child_weight" = 1,
                                                    "verbose" = 0,
                                                    "nthread" = 1))

# Generate BMI results only using WHR z-scores:
bmi_adapt_whr_only_cv_adj_results <- adapt_xgboost_cv(as.matrix(
  bmi_esnps_data[,"whr15_z_score"]),
  bmi_esnps_data$adj_bmi15_p,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = whr_only_args_search,
  muargs = whr_only_args_search,
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2),
  s0 = rep(0.05, nrow(bmi_esnps_data)))

# Save these results:
# saveRDS(bmi_adapt_whr_only_cv_adj_results,
#         "data/bmi/adapt_cv_results/bmi_adj15_whr_only_s05_2cv.rds")

# Now with the eQTL slopes
whr_slopes_args_search <- list("nrounds100md3" = list("nrounds" = 100,
                                                      "max_depth" = 3,
                                                      "min_child_weight" = 1,
                                                      "verbose" = 0,
                                                      "nthread" = 5),
                               "nrounds150md3" = list("nrounds" = 150,
                                                      "max_depth" = 3,
                                                      "min_child_weight" = 1,
                                                      "verbose" = 0,
                                                      "nthread" = 5),
                               "nrounds100md6" = list("nrounds" = 100,
                                                      "max_depth" = 6,
                                                      "min_child_weight" = 1,
                                                      "verbose" = 0,
                                                      "nthread" = 5),
                               "nrounds150md6" = list("nrounds" = 150,
                                                      "max_depth" = 6,
                                                      "min_child_weight" = 1,
                                                      "verbose" = 0,
                                                      "nthread" = 5))

bmi_adapt_whr_slopes_only_cv_adj_results <- adapt_xgboost_cv(as.matrix(
  bmi_esnps_data[,c("whr15_z_score", ave_slope_vars)]),
  bmi_esnps_data$adj_bmi15_p,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = whr_slopes_args_search,
  muargs = whr_slopes_args_search,
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2),
  s0 = rep(0.05, nrow(bmi_esnps_data)))
 
# Save these results:
# saveRDS(bmi_adapt_whr_slopes_only_cv_adj_results,
#         "data/bmi/adapt_cv_results/bmi_adj15_whr_slopes_only_s05_2cv.rds")


# Results for the new 2018 studies:
bmi18_adapt_cv_results <- adapt_xgboost_cv(as.matrix(
  bmi_esnps_data[,bmi_variables]),
  bmi_esnps_data$bmi18_p,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = args_search,
  muargs = args_search,
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2),
  s0 = rep(0.05, nrow(bmi_esnps_data)))
 
# Save these results:
# saveRDS(bmi18_adapt_cv_results,
#         "data/bmi/adapt_cv_results/bmi18_s05_2cv.rds")

# Next generate intercept-only results:
adapt_intercept_results <- adapt_glm(x = bmi_esnps_data,
                                     pvals = bmi_esnps_data$adj_bmi15_p,
                                     pi_formulas = "1",
                                     mu_formulas = "1",
                                     verbose = list(print = FALSE, 
                                                    fit = FALSE,
                                                    ms = FALSE),
                                     s0 = rep(0.05, nrow(bmi_esnps_data)))
# saveRDS(adapt_intercept_results, 
#         "data/bmi/adapt_cv_results/bmi_intercept_only_s05.rds")


