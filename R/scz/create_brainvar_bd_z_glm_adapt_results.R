# PURPOSE: Create the AdaPT GLM results using splines

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
library(splines)

# Load the data:
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Generate the results using the original approach by Lei & Fithian (2018)
bd_z_formulas <- paste0("ns(x = z_bip_14, df = ", 2:5, ")")
adapt_bd_z_glm <- adapt_glm(x = bip_scz_brainvar_data,
                            pvals = bip_scz_brainvar_data$scz_14_P,
                            pi_formulas = bd_z_formulas,
                            mu_formulas = bd_z_formulas,
                            s0 = rep(0.05, nrow(bip_scz_brainvar_data)))
# Save these results:
# saveRDS(adapt_bd_z_glm,
#         "data/bip_schz_data/brainvar_results/scz_bd_z_glm_s05.rds")



