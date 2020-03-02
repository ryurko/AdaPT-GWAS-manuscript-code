# PURPOSE: Generate simulations varying the starting threshold s0 for AdaPT to see
#          the impact on performance. Does lower s0 lead to better results on 
#          average? How does this vary with the number of CV steps?

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

# Initialize the potential parameters
args_search <- list("nrounds100md2" = list("nrounds" = 100,
                                           "max_depth" = 2,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 3),
                    "nrounds150md3" = list("nrounds" = 150,
                                           "max_depth" = 3,
                                           "min_child_weight" = 1,
                                           "verbose" = 0,
                                           "nthread" = 3))


# Next proceed to loop through the simulations in parallel 
plan(multiprocess(workers = 5)) # Distribute across 5 cores
sim_s0_n_cv_results <- future_map_dfr(1:100, 
                                      function(i) {
                                        
                                        # First generate the probability of being non-null for each test using 
                                        # the first non-null probability model
                                        sim_nonnull_prob <- predict(scz14_model$model_fit[[1]],
                                                                    newdata = as.matrix(bip_scz_brainvar_data[,scz14_model$model_fit[[1]]$feature_names]),
                                                                    type = "response")
                                        
                                        # Next the "truth":
                                        sim_test_type <- ifelse(sapply(sim_nonnull_prob,
                                                                       function(p) rbinom(1, 1, p)) == 1,
                                                                "h1", "h0")
                                        n_h1 <- length(which(sim_test_type == "h1"))
                                        n_h0 <- length(which(sim_test_type == "h0"))
                                        
                                        # Now the observed effect size with the first non-null effect size model:
                                        sim_alt_means <- predict(scz14_model$model_fit[[2]],
                                                                 newdata = as.matrix(bip_scz_brainvar_data[,scz14_model$model_fit[[2]]$feature_names]),
                                                                 type = "response")
                                        
                                        # Now store the test_effects for the both studies using the type of test to dicate whether or not it is alt or null:
                                        sim_test_effect_1 <- ifelse(sim_test_type == "h1", 
                                                                    rexp(nrow(bip_scz_brainvar_data), 
                                                                         1/sim_alt_means), 
                                                                    rexp(nrow(bip_scz_brainvar_data), 1))
                                        sim_test_effect_2 <- ifelse(sim_test_type == "h1", 
                                                                    rexp(nrow(bip_scz_brainvar_data), 
                                                                         1/sim_alt_means), 
                                                                    rexp(nrow(bip_scz_brainvar_data), 1))
                                        
                                        # Now the p-values for the two tests:
                                        sim_pvals_1 <- exp(-sim_test_effect_1)
                                        sim_pvals_2 <- exp(-sim_test_effect_2)
                                        
                                        # Next generate updated alternative p-values for the second study
                                        # reflecting the change in standard errors for the 2018 studies
                                        # as compared to 2014:
                                        
                                        # First convert sim_pvals_2 to zstats
                                        sim_alt_mean_abs_zstats_2 <- abs(qnorm(sim_pvals_2 / 2))
                                        
                                        # Now multiply the 2014 se / new 2018 se:
                                        sim_alt_mean_new_abs_zstats_2 <- sim_alt_mean_abs_zstats_2 * 
                                          (bip_scz_brainvar_data$scz_14_SE / bip_scz_brainvar_data$new_se_scz_18)
                                        
                                        # Get updated p-values from this - but only for the alternatives:
                                        new_sim_pvals_2 <- ifelse(sim_test_type == "h1",
                                                                  2 * pnorm(-sim_alt_mean_new_abs_zstats_2),
                                                                  sim_pvals_2)
                                        
                                        # Generate results over s0 and n_cv steps:
                                        s0_results <- future_map_dfr(c(.05, .25, .45),
                                                                     function(s0_value) {
                                                                       # now the results for different n_cv_steps:
                                                                       future_map_dfr(c(1, 2, 5),
                                                                                      function(n_cv_steps) {
                                                                                        sim_adapt_results <- adapt_xgboost_cv(
                                                                                          as.matrix(bip_scz_brainvar_data[,scz_variables]),
                                                                                          sim_pvals_1,
                                                                                          verbose = list(print = FALSE, 
                                                                                                         fit = FALSE,
                                                                                                         ms = FALSE),
                                                                                          piargs = args_search,
                                                                                          muargs = args_search,
                                                                                          n_folds = 5,
                                                                                          niter_ms = 10,
                                                                                          nms = as.integer(n_cv_steps),
                                                                                          s0 = rep(s0_value, nrow(bip_scz_brainvar_data)))
                                                                                        
                                                                                        best_model_i <- which.max(sim_adapt_results$model_holdout_ll_sums[[n_cv_steps]])
                                                                                        best_model_name <- names(args_search)[best_model_i]
                                                                                        
                                                                                        target_alphas <- c(.01, .05, .1, .15, .2)
                                                                                        # Generate a data frame of results for the vector of alpha values:
                                                                                        alpha_adapt_n_cv_steps_results_df <- do.call(rbind,
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
                                                                                                                                                disc_power <- length(which(sim_test_type[adapt_disc] == "h1")) / n_h1
                                                                                                                                                actual_fdp <- length(which(sim_test_type[adapt_disc] == "h0")) / n_rej
                                                                                                                                                n_rep <- length(which(new_sim_pvals_2[adapt_disc] <= alpha))
                                                                                                                                                rep_prop <- n_rep / n_rej
                                                                                                                                                
                                                                                                                                              } else {
                                                                                                                                                disc_power <- 0
                                                                                                                                                actual_fdp <- 0
                                                                                                                                                n_rep <- 0
                                                                                                                                                rep_prop <- 0
                                                                                                                                              }
                                                                                                                                              
                                                                                                                                              # Return the fdp and power:
                                                                                                                                              data.frame("n_cv_steps" = n_cv_steps,
                                                                                                                                                         "s0_value" = s0_value,
                                                                                                                                                         "alpha" = alpha,
                                                                                                                                                         "n_alt" = n_h1,
                                                                                                                                                         "n_null" = n_h0,
                                                                                                                                                         "n_disc" = n_rej,
                                                                                                                                                         "disc_power" = disc_power,
                                                                                                                                                         "true_fdp" = actual_fdp,
                                                                                                                                                         "n_rep_disc" = n_rep,
                                                                                                                                                         "rep_prop" = rep_prop,
                                                                                                                                                         "final_model_settings" = best_model_name,
                                                                                                                                                         "sim_index" = i)
                                                                                                                                            }))
                                                                                        
                                                                                        return(alpha_adapt_n_cv_steps_results_df)
                                                                                        
                                                                                      }) 
                                                                     })
                                        
                                        return(s0_results)
                                        
                                      })

# Save the results
# readr::write_csv(sim_s0_n_cv_results,
#                 "data/si_simulations/s0_n_cv_sims.csv")
