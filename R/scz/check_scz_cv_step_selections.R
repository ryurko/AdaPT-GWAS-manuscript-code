# PURPOSE: View the AdaPT CV selection at each step for the SCZ results

# AUTHOR: Ron Yurko

# Load the AdaPT results
scz_with_bd_z_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_only_s05_2cv.rds")
scz_with_bd_z_eqtl_slopes <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_s05_2cv.rds")
scz_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")
scz_with_wgcna_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_wgcna_only_s05_2cv.rds")

# Initialize the vector of arguments:
scz_with_bd_z_eqtl_slopes_wgcna_args <- c("nrounds100md2", "nrounds150md2",
                                          "nrounds100md3", "nrounds150md3")
scz_with_bd_z_only_args <- c("nrounds100md1", "nrounds150md1",
                             "nrounds100md6", "nrounds150md6")
scz_with_bd_z_eqtl_slopes_args <- c("nrounds100md3", "nrounds150md3",
                                    "nrounds100md6", "nrounds150md6")
wgcna_args <- c("nrounds100md1", "nrounds150md1",
                "nrounds100md2", "nrounds150md2",
                "nrounds100md3", "nrounds150md3")

# Results with all variables:
# First step:
scz_with_bd_z_eqtl_slopes_wgcna_args[which.max(scz_with_bd_z_eqtl_slopes_wgcna$model_holdout_ll_sums[[1]])]
# [1] "nrounds150md3"
# Second step:
scz_with_bd_z_eqtl_slopes_wgcna_args[which.max(scz_with_bd_z_eqtl_slopes_wgcna$model_holdout_ll_sums[[2]])]
# [1] "nrounds150md3"

# Results using only BD z-statistics
# First step:
scz_with_bd_z_only_args[which.max(scz_with_bd_z_only$model_holdout_ll_sums[[1]])]
# [1] "nrounds150md1"
# Second step:
scz_with_bd_z_only_args[which.max(scz_with_bd_z_only$model_holdout_ll_sums[[2]])]
# [1] "nrounds150md6"

# Results using BD z-statistics + eQTL slopes
# First step:
scz_with_bd_z_eqtl_slopes_args[which.max(scz_with_bd_z_eqtl_slopes$model_holdout_ll_sums[[1]])]
# [1] "nrounds150md3"
# Second step:
scz_with_bd_z_eqtl_slopes_args[which.max(scz_with_bd_z_eqtl_slopes$model_holdout_ll_sums[[2]])]
# [1] "nrounds150md6"

# Results using WGCNA only
# First step:
wgcna_args[which.max(scz_with_wgcna_only$model_holdout_ll_sums[[1]])]
# [1] "nrounds150md3"
# Second step:
wgcna_args[which.max(scz_with_wgcna_only$model_holdout_ll_sums[[2]])]

# Next load the GTEx SCZ results:
scz_gtex_adapt <- readRDS("data/bip_schz_data/gtex_results/cv_tune_results/scz_gtex_cortical_s05_2cv.rds")
# First step:
scz_with_bd_z_eqtl_slopes_wgcna_args[which.max(scz_gtex_adapt$model_holdout_ll_sums[[1]])]
# [1] "nrounds100md2"
# Second step:
scz_with_bd_z_eqtl_slopes_wgcna_args[which.max(scz_gtex_adapt$model_holdout_ll_sums[[2]])]
# [1] "nrounds150md3"



