# PURPOSE: Generate simulation results from a fake truth model demonstrating
#          that overfitting the boosting models does not affect the FDR control
#          but only affects the power. Will NOT use the CV here but rather
#          just fix the boosting settings to be constant throughout the search.

# AUTHOR: Ron Yurko

# ------------------------------------------------------------------------------

# Access necessary packages:
library(tidyverse)
library(future)
# Use the development version of furrr:
# devtools::install_github("DavisVaughan/furrr")
library(furrr)
# Use the modified version of AdaPT:
# devtools::install_github("ryurko/adaptMT")
library(adaptMT)

# Load the data
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Load the model generated with all variables:
scz14_model <- readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")

# Grab the WGCNA cols: -------
wgcna_cols <- colnames(bip_scz_brainvar_data)[which(stringr::str_detect(colnames(bip_scz_brainvar_data), 
                                                                        "brainvar_any_gene_"))]
# Get the variables: -------
ave_columns <- colnames(bip_scz_brainvar_data)[which(stringr::str_detect(colnames(bip_scz_brainvar_data), "ave_"))]

# Create a vector of variables to use for SCZ: ------
scz_variables <- c("z_bip_14", ave_columns, wgcna_cols)

# Next proceed to loop through the simulations in parallel (10 at a time)
plan(multiprocess(workers = 10))
sim_results <- future_map_dfr(c(1:50),
                                         function(sim_i) {
                                           # Generate simulated test types and p-values using this model:
                                           # Generate fake p-values from this:
                                           
                                           sim_nonnull_prob <- predict(scz14_model$model_fit[[1]],
                                                                       newdata = as.matrix(bip_scz_brainvar_data[,scz14_model$model_fit[[1]]$feature_names]),
                                                                       type = "response")
                                           
                                           fake_true_test_types <- ifelse(sapply(sim_nonnull_prob,
                                                                                 function(p) rbinom(1, 1, p)) == 1,
                                                                          "h1", "h0")
                                           n_h1 <- length(which(fake_true_test_types == "h1"))
                                           n_h0 <- length(which(fake_true_test_types == "h0"))
                                           
                                           # Now the observed effect size:
                                           sim_alt_means <- predict(scz14_model$model_fit[[2]],
                                                                    newdata = as.matrix(bip_scz_brainvar_data[,scz14_model$model_fit[[2]]$feature_names]),
                                                                    type = "response")
                                           
                                           test_y <- ifelse(fake_true_test_types == "h1", 
                                                            rexp(nrow(bip_scz_brainvar_data), 
                                                                 1/sim_alt_means), 
                                                            rexp(nrow(bip_scz_brainvar_data), 1))
                                           sim_pvals <- exp(-test_y)
                                           
                                           # Generate the updated replicated test p-vals with smaller standard errors:
                                           test_y_2 <- ifelse(fake_true_test_types == "h1", 
                                                              rexp(nrow(bip_scz_brainvar_data), 
                                                                   1/sim_alt_means), 
                                                              rexp(nrow(bip_scz_brainvar_data), 1))
                                           sim_pvals_2 <- exp(-test_y_2)
                                           # First convert sim_pvals_2 to zstats
                                           sim_alt_mean_abs_zstats_2 <- abs(qnorm(sim_pvals_2 / 2))
                                           
                                           # Now multiply the 2014 se / new 2018 se:
                                           sim_alt_mean_new_abs_zstats_2 <- sim_alt_mean_abs_zstats_2 * 
                                             (bip_scz_brainvar_data$scz_14_SE / bip_scz_brainvar_data$new_se_scz_18)
                                           
                                           # Get updated p-values from this - but only for the alternatives:
                                           new_sim_pvals_2 <- ifelse(fake_true_test_types == "h1",
                                                                     2 * pnorm(-sim_alt_mean_new_abs_zstats_2),
                                                                     sim_pvals_2)
                                           # Now loop over the potential number of trees:
                                           vary_ntree_results <- future_map_dfr(seq(100, 900, by = 200),
                                                                                function(n_trees) {
                                                                                  # Create the boosting arguments:
                                                                                  adapt_args = list("nrounds" = n_trees,
                                                                                                    "max_depth" = 3,
                                                                                                    "min_child_weight" = 1,
                                                                                                    "verbose" = 0,
                                                                                                    "nthread" = 2)
                                                                                  # Now get the AdaPT XGBoost model results as before using sim_pvals_1:
                                                                                  sim_adapt_results <- adapt_xgboost(as.matrix(bip_scz_brainvar_data[,scz_variables]),
                                                                                                                     sim_pvals,
                                                                                                                     verbose = list(print = FALSE, 
                                                                                                                                    fit = FALSE, 
                                                                                                                                    ms = FALSE),
                                                                                                                     piargs = adapt_args,
                                                                                                                     muargs = adapt_args,
                                                                                                                     # use low starting threshold
                                                                                                                     s0 = rep(0.05, nrow(bip_scz_brainvar_data)))
                                                                                  
                                                                                  target_alphas <- c(.01, .05, .1, .15, .2)
                                                                                  # Generate a data frame of results for the vector of alpha values:
                                                                                  alpha_adapt_ntree_results_df <- do.call(rbind,
                                                                                                                          lapply(target_alphas, 
                                                                                                                                 function(alpha) {
                                                                                                                                   # Access the discoveries for alpha:
                                                                                                                                   adapt_disc <- which(sim_adapt_results$qvals <= alpha)
                                                                                                                                   n_rej <- length(adapt_disc)
                                                                                                                                   # If there are any rejection then compute the actual
                                                                                                                                   # FDP as well as the replication FDP based 
                                                                                                                                   # on which tests are replicated using 
                                                                                                                                   # the simple nominal significance threshold:
                                                                                                                                   if (n_rej > 0) {
                                                                                                                                     disc_power <- length(which(fake_true_test_types[adapt_disc] == "h1")) / n_h1
                                                                                                                                     actual_fdp <- length(which(fake_true_test_types[adapt_disc] == "h0")) / n_rej
                                                                                                                                     n_rep <- length(which(new_sim_pvals_2[adapt_disc] <= alpha))
                                                                                                                                     rep_prop <- n_rep / n_rej
                                                                                                                                     
                                                                                                                                   } else {
                                                                                                                                     disc_power <- 0
                                                                                                                                     actual_fdp <- 0
                                                                                                                                     n_rep <- 0
                                                                                                                                     rep_prop <- 0
                                                                                                                                   }
                                                                                                                                   
                                                                                                                                   # Return the fdp and power:
                                                                                                                                   data.frame("n_trees" = n_trees,
                                                                                                                                              "alpha" = alpha,
                                                                                                                                              "n_alt" = n_h1,
                                                                                                                                              "n_null" = n_h0,
                                                                                                                                              "n_disc" = n_rej,
                                                                                                                                              "disc_power" = disc_power,
                                                                                                                                              "true_fdp" = actual_fdp,
                                                                                                                                              "n_rep_disc" = n_rep,
                                                                                                                                              "rep_prop" = rep_prop,
                                                                                                                                              "sim_index" = sim_i)
                                                                                                                                 }))
                                                                                  return(alpha_adapt_ntree_results_df)
                                                                                })
                                           
                                           # Return the results:
                                           return(vary_ntree_results)
                                           
                                         })


# Save the results
#readr::write_csv(sim_results,
#                 "data/si_simulations/overfit_sims.csv")





