# PURPOSE: Generate the AdaPT CV results using the BrainVar eSNPs with all
#          models returned, using the GWAS results with all 2018 studies. 
#          Additionally generate the AdaPT intercept-only results

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
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Get the variables:
ave_columns <- stringr::str_subset(colnames(bip_scz_brainvar_data),
                                   "ave_")
# Now each of the wgcna_label columns:
wgcna_cols <- stringr::str_subset(colnames(bip_scz_brainvar_data),
                                  "brainvar_any_gene_")

# Create the list of variables to use for SCZ:
scz_variable_list <- list("bd_z_only" = c("z_bip_14"),
                          "bd_z_eqtl_slopes" = c("z_bip_14", ave_columns),
                          "bd_z_eqtl_slopes_wgcna" = c("z_bip_14", ave_columns,
                                                       wgcna_cols),
                          "wgcna_only" = c(wgcna_cols))


# Initialize the list of XGBoost model choices for different choices of variables:
bd_only_args <- list("nrounds100md1" = list("nrounds" = 100,
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

bd_eqtl_args <-  list("nrounds100md3" = list("nrounds" = 100,
                                             "max_depth" = 3,
                                             "min_child_weight" = 1,
                                             "verbose" = 0,
                                             "nthread" = 1),
                      "nrounds150md3" = list("nrounds" = 150,
                                             "max_depth" = 3,
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

bd_eqtl_wgcna_args <- list("nrounds100md2" = list("nrounds" = 100,
                                                  "max_depth" = 2,
                                                  "min_child_weight" = 1,
                                                  "verbose" = 0,
                                                  "nthread" = 3),
                           "nrounds150md2" = list("nrounds" = 150,
                                                  "max_depth" = 2,
                                                  "min_child_weight" = 1,
                                                  "verbose" = 0,
                                                  "nthread" = 3),
                           "nrounds100md3" = list("nrounds" = 100,
                                                  "max_depth" = 3,
                                                  "min_child_weight" = 1,
                                                  "verbose" = 0,
                                                  "nthread" = 3),
                           "nrounds150md3" = list("nrounds" = 150,
                                                  "max_depth" = 3,
                                                  "min_child_weight" = 1,
                                                  "verbose" = 0,
                                                  "nthread" = 3))

wgcna_only_args <- list("nrounds100md1" = list("nrounds" = 100,
                                               "max_depth" = 1,
                                               "min_child_weight" = 1,
                                               "verbose" = 0,
                                               "nthread" = 3),
                        "nrounds150md1" = list("nrounds" = 150,
                                               "max_depth" = 1,
                                               "min_child_weight" = 1,
                                               "verbose" = 0,
                                               "nthread" = 3),
                        "nrounds100md2" = list("nrounds" = 100,
                                               "max_depth" = 2,
                                               "min_child_weight" = 1,
                                               "verbose" = 0,
                                               "nthread" = 3),
                        "nrounds150md2" = list("nrounds" = 150,
                                               "max_depth" = 2,
                                               "min_child_weight" = 1,
                                               "verbose" = 0,
                                               "nthread" = 3),
                        "nrounds100md3" = list("nrounds" = 100,
                                               "max_depth" = 3,
                                               "min_child_weight" = 1,
                                               "verbose" = 0,
                                               "nthread" = 3),
                        "nrounds150md3" = list("nrounds" = 150,
                                               "max_depth" = 3,
                                               "min_child_weight" = 1,
                                               "verbose" = 0,
                                               "nthread" = 3))

scz_variable_arg_search <- list("bd_z_only" = bd_only_args,
                                "bd_z_eqtl_slopes" = bd_eqtl_args,
                                "bd_z_eqtl_slopes_wgcna" = bd_eqtl_wgcna_args,
                                "wgcna_only" = wgcna_only_args)

# Now create the results in parallel for the three different sets of information
plan(multiprocess(workers = 4))
results <- furrr::future_map_dfr(1:length(scz_variable_list),
                                 function(var_i) {
                                   
                                   # Generate the AdaPT CV results with just 1 step
                                   cv_model_names <- names(scz_variable_arg_search[[var_i]])
                                   adapt_results <- adapt_xgboost_cv(as.matrix(
                                     bip_scz_brainvar_data[,scz_variable_list[[var_i]]]),
                                     bip_scz_brainvar_data$scz_18_P,
                                     verbose = list(print = FALSE, 
                                                    fit = FALSE,
                                                    ms = FALSE),
                                     piargs = scz_variable_arg_search[[var_i]],
                                     muargs = scz_variable_arg_search[[var_i]],
                                     s0 = rep(0.05, nrow(bip_scz_brainvar_data)),
                                     n_folds = 5,
                                     niter_ms = 10,
                                     nms = as.integer(2))
                                   
                                   # Save these results:
                                   # saveRDS(adapt_results,
                                   #         paste0("data/bip_schz_data/brainvar_results/cv_tune_results/scz18_with_",
                                   #                names(scz_variable_list)[var_i], "_s05_2cv.rds"))
                                   
                                   # Now return a dataframe of results over
                                   # target alpha values:
                                   target_alphas <- c(.01, .05, .1, .15, .2)
                                   # Which was the selected model:
                                   best_model_i <- which.max(adapt_results$model_holdout_ll_sums[[2]])
                                   best_model_name <- cv_model_names[best_model_i]
                                   
                                   # Generate a data frame of results for the vector of alpha values:
                                   do.call(rbind,
                                           lapply(target_alphas, 
                                                  function(alpha) {
                                                    # Access the discoveries for alpha:
                                                    adapt_disc <- which(adapt_results$qvals <= alpha)
                                                    
                                                    # Return the fdp and power:
                                                    return(data.frame("method" = "adapt_cv",
                                                                      "variables" = names(scz_variable_list)[var_i],
                                                                      "alpha" = alpha,
                                                                      "n_disc" = length(adapt_disc),
                                                                      "selected_model" = best_model_name))   
                                                  }))
                                   
                                 })

# Save the dataframe summarizing results:
# readr::write_csv(results,
#                  "data/bip_schz_data/brainvar_results/cv_tune_results/scz18_adapt_s05_2cv_summary.csv")

# Next generate intercept-only results:
adapt_intercept_results <- adapt_glm(x = bip_scz_brainvar_data,
                                     pvals = bip_scz_brainvar_data$scz_18_P,
                                     pi_formulas = "1",
                                     mu_formulas = "1",
                                     verbose = list(print = FALSE, 
                                                    fit = FALSE,
                                                    ms = FALSE),
                                     s0 = rep(0.05, nrow(bip_scz_brainvar_data)))
# saveRDS(adapt_intercept_results, 
#         "data/bip_schz_data/brainvar_results/scz18_intercept_only_s05.rds")
