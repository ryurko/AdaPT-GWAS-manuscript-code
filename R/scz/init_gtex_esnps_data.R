# PURPOSE: Create the GTEx eSNPs dataset used for AdaPT along with their
#          respective expression/co-expression covariates

# Author: Ron Yurko

# Access necessary packages
library(data.table)
library(tidyverse)

# Load the SNPs BD and SCZ data before any pre-processing takes place, that was
# generated in the R/scz/init_gwas_data.R script:
bip_scz_data_14_18 <- fread("data/bip_schz_data/bip_scz_data_14_18_snps.csv")

# Load the eQTL data that was generated in R/gtex_wgcna/init_eqtl_data.R, which
# has eQTLs for any tissue type 
egenes_eqtl_data <- fread("data/gtex_eqtl/eqtl_data.csv") 

# Filter the joined dataset to only be the SNPs that are eqtls for any tissue
bip_scz_eqtl_14_18 <- bip_scz_data_14_18[SNP %in% unique(egenes_eqtl_data$eqtl_name),]

# Now filter these on the INFO score for both SCZ and BD in 2014:
bip_scz_eqtl_14_18_filtered <- bip_scz_eqtl_14_18[bip_14_INFO > 0.6 & scz_14_INFO > 0.6, ]

# Filter the eQTL data to only be the eSNPs for cortical tissues
cortical_eqtl_data <- egenes_eqtl_data %>%
  filter(eqtl_name %in% bip_scz_eqtl_14_18_filtered$SNP,
         str_detect(region, "(BA24)|(BA9)"))

# Create three different summaries: separately for each tissue, as well as across
cortical_esnp_summary <- cortical_eqtl_data %>%
  mutate(region = tolower(region)) %>%
  # First calculate the average across both types:
  group_by(eqtl_name) %>%
  mutate(ave_abs_cortical_eqtl_slope = mean(abs(eqtl_slope), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(eqtl_name, region) %>%
  summarize(ave_abs_cortical_eqtl_slope = first(ave_abs_cortical_eqtl_slope),
            ave_abs_eqtl_slope = mean(abs(eqtl_slope), na.rm = TRUE)) %>%
  gather(variable, value, -eqtl_name, -region, -ave_abs_cortical_eqtl_slope) %>%
  unite(region_variable, "region", "variable") %>%
  spread(region_variable, value)
cortical_esnp_summary[is.na(cortical_esnp_summary)] <- 0

# Now left-join both of these datasets to the eSNP dataset:
bip_scz_eqtl_14_18_filtered <- bip_scz_eqtl_14_18_filtered %>%
  left_join(cortical_esnp_summary, by = c("SNP" = "eqtl_name"))
bip_scz_eqtl_14_18_filtered[is.na(bip_scz_eqtl_14_18_filtered)] <- 0

# Next join the WGCNA information:
wgcna_cortical_data <- read_csv("data/gtex_eqtl/wcgna_results/cortical_wgcna_data.csv")

# First go through each SNP in the model dataset and construct a list corresponding
# to the vector of genes for each:
esnp_gene_list <- map(bip_scz_eqtl_14_18_filtered$SNP,
                      function(x) {
                        cortical_eqtl_data %>%
                          filter(eqtl_name == x) %>%
                          pull(ensg_gene)
                      })
names(esnp_gene_list) <- bip_scz_eqtl_14_18_filtered$SNP

cortical_any_gene_member <- map_dfc(unique(wgcna_cortical_data$wgcna_label),
                                    function(x) {
                                      
                                      # Get the vector of candidate genes
                                      # based on the module assignment:
                                      candidate_genes <- 
                                        wgcna_cortical_data$gene[
                                          wgcna_cortical_data$wgcna_label == x
                                          ]
                                      
                                      # Now go through each SNP in the data,
                                      # and see if any of the genes are 
                                      # members:
                                      membership <- sapply(esnp_gene_list,
                                                           function(y) {
                                                             as.numeric(any(y %in% candidate_genes))
                                                           })
                                      result <- data.frame(members = membership)
                                      colnames(result) <- paste0("cortical_any_gene_", x)
                                      return(result)
                                    })

# Join these columns
bip_scz_eqtl_14_18_filtered <- bip_scz_eqtl_14_18_filtered %>%
  bind_cols(cortical_any_gene_member)

# Save the data:
# write_csv(bip_scz_eqtl_14_18_filtered,
#           "data/bip_schz_data/bip_scz_data_14_18_gtex_eqtls_cortical_wgcna.csv")

