# PURPOSE: Create the text file of associated cis-eQTL genes with the AdaPT discoveries
#          on the 2014 SCZ dataset

# AUTHOR: Ron Yurko

# Access necessary packages:
library(tidyverse)
library(readxl)

# Load the data:
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Load the final model results:
scz_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")

# Define the function used to read in excel sheets for the BrainVar eQTL data
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# Load the BrainVar eQTL sheets:
brainvar_eqtl_sheets <- read_excel_allsheets("data/brainvar_eqtl/eqtls.xlsx")

# Create a vector of SNPs that are discoveries:
scz_brainvar_14_disc_snps <- bip_scz_brainvar_data$SNP[which(scz_with_bd_z_eqtl_slopes_wgcna$qvals <= 0.05)]

# How many associated genes?
brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05 %>%
  filter(rsID %in% scz_brainvar_14_disc_snps) %>%
  mutate(ensembl_id = str_sub(GeneId, 1, 15)) %>%
  pull(ensembl_id) %>%
  unique() %>% 
  length()
# 136

# Write the list of associated genes:
scz_brainvar_gene_file <- file("data/gene_lists/scz_brainvar_disc_genes.txt")
brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05 %>%
  filter(rsID %in% scz_brainvar_14_disc_snps) %>%
  mutate(ensembl_id = str_sub(GeneId, 1, 15)) %>%
  pull(ensembl_id) %>%
  unique() %>%
  paste(collapse = "\n") %>% 
  writeLines(scz_brainvar_gene_file) 
close(scz_brainvar_gene_file) 


