# PURPOSE: Generate the replication simulation results using the final non-null
#          effect size model from the AdaPT search. 

# AUTHOR: Ron Yurko

# Load necessary packages:
library(tidyverse)

# ------------------------------------------------------------------------------

# Load the SCZ BrainVar data:
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Load the final model results:
scz_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")

# Which of the SNPs are discoveries at alpha = 0.05 in 2014:
scz_14_disc_type <- ifelse(scz_with_bd_z_eqtl_slopes_wgcna$qvals <= 0.05,
                           "h1", "h0")

# Using the final effect size model generate the predicted mu values:
sim_alt_means_14 <- 
  predict(scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[38]],
          newdata = as.matrix(bip_scz_brainvar_data[,scz_with_bd_z_eqtl_slopes_wgcna$model_fit[[38]]$feature_names]),
          type = "response")

# Now simulate with this final effect size model:
start_time <- Sys.time()
sim_change_se_results <- map_dbl(1:10000,
                                 function(sim_i) {
                                   # Generate the simulated observed effect sizes 
                                   # for the data - only modifying the discoveries from 2014
                                   sim_test_effect <- ifelse(scz_14_disc_type == "h1", 
                                                             rexp(nrow(bip_scz_brainvar_data), 
                                                                  1 / sim_alt_means_14), 
                                                             rexp(nrow(bip_scz_brainvar_data), 1))
                                   
                                   # Get the p-vals
                                   sim_pvals <- exp(-sim_test_effect)
                                   # What the z-stats
                                   sim_zstats <- abs(qnorm(sim_pvals / 2))
                                   # Update the standard errors
                                   sim_new_zstats <- sim_zstats * (bip_scz_brainvar_data$scz_14_SE / 
                                                                     bip_scz_brainvar_data$new_se_scz_18)
                                   new_sim_pvals <- 2 * pnorm(-sim_new_zstats)
                                   # Now return how many of the discoveries are less than 0.05?
                                   length(which(new_sim_pvals[which(scz_14_disc_type == "h1")] <= 0.05))
                                 })
end_time <- Sys.time()
end_time - start_time
# Time difference of 1.330613 mins
summary(sim_change_se_results / length(which(scz_14_disc_type == "h1")))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5125  0.5552  0.5670  0.5667  0.5777  0.6275

