# PURPOSE: Create the BrainVar eSNPs dataset used for AdaPT along with their
#          respective expression/co-expression covariates

# Author: Ron Yurko

# Access necessary packages
library(readxl)
library(data.table)
library(tidyverse)

# The following function is from code here: 
# https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(sheet_i) readxl::read_excel(filename,
                                                             sheet = sheet_i))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}

# Load the different BrainVar sheets, first the eQTL sheets:
brainvar_eqtl_sheets <- read_excel_allsheets("data/brainvar_eqtl/eqtls.xlsx")

# Next the WGCNA sheets (this is more straightforward to join)
brainvar_wgcna_sheets <- read_excel_allsheets("data/brainvar_eqtl/wgcna.xlsx")

# Load the SNPs BD and SCZ data before any pre-processing takes place, that was
# generated in the R/scz/init_gwas_data.R script:
bip_scz_data_14_18 <- fread("data/bip_schz_data/bip_scz_data_14_18_snps.csv")

# Filter to only include those with INFO scores > 0.6 in both BD and SCZ 2014:
bip_scz_data_14_18 <- bip_scz_data_14_18[bip_14_INFO > 0.6 & 
                                           scz_14_INFO > 0.6,]

# Based on the documentation in the Werling et al. (2019) supplementary materials, there were two
# different methods used to generate two _different_ sets of eQTLs: HCP and SVA.
# The HCP method was considered for the primary analysis. The first sheet of 
# interest in the brainvar_eqtl_sheets is named  HCP_allEQTLs_FDR0.05. 
# The the ColumnKey sheet contains the glossary describing the different columns 
# in the dataset. Will only use the SNPs that are also in the HCP eQTL sheet:
bip_scz_data_14_18 <- bip_scz_data_14_18[SNP %in% 
                                           unique(brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05$rsID),]

# Now join the eQTL columns of interest to use by summarizing the eSNPs in the
# BD and SCZ data:
bd_scz_brainvar_eqtl_summary <- brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05 %>%
  filter(rsID %in% bip_scz_data_14_18$SNP) %>%
  group_by(rsID) %>%
  # Now for each eSNP calculate its average and max betas for each of the 
  # three different betas:
  summarize(ave_abs_pre_beta = mean(abs(Beta_Pre), na.rm = TRUE),
            max_abs_pre_beta = max(abs(Beta_Pre), na.rm = TRUE),
            ave_abs_post_beta = mean(abs(Beta_Post), na.rm = TRUE),
            max_abs_post_beta = max(abs(Beta_Post), na.rm = TRUE),
            ave_abs_comp_beta = mean(abs(Beta_Comp), na.rm = TRUE),
            max_abs_comp_beta = max(abs(Beta_Comp), na.rm = TRUE),
            # Indicators if any of the SNP's base pairs are the best 
            any_gene_best_best = as.numeric(any(Gene_Best_Best == "Yes")),
            any_gene_best_pre = as.numeric(any(Gene_Best_Pre == "Yes")),
            any_gene_best_post = as.numeric(any(Gene_Best_Post == "Yes")),
            any_gene_best_comp = as.numeric(any(Gene_Best_Comp == "Yes")),
            # Next its most popular category for each categorization:
            freq_cat_by_eqtl = names(sort(table(Cat_byEqtl),
                                          decreasing = TRUE)[1]),
            # Note that the column is incorrectly called Cat_byConstanteBest
            freq_cat_by_gene_best = names(sort(table(Cat_byConstanteBest),
                                               decreasing = TRUE)[1]))

# Join this data to the BD and SCZ summary stats data:
bip_scz_data_14_18 <- bip_scz_data_14_18 %>%
  inner_join(bd_scz_brainvar_eqtl_summary,
             by = c("SNP" = "rsID"))

# Next need to join the WGCNA assignments from the WGCNA sheets using the GeneId
# column from the eQTL sheet. First go through each SNP in the model dataset and 
# construct a list corresponding to the vector of genes for each:
snp_gene_list <- map(bip_scz_data_14_18$SNP,
                     function(x) {
                       as.data.table(brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05) %>%
                         .[rsID == x, GeneId]
                     })
names(snp_gene_list) <- bip_scz_data_14_18$SNP

# Any gene-pair member of BrainVar WGCNA modules?
brainvar_any_gene_member <- map_dfc(unique(brainvar_wgcna_sheets$`Module Genes`$Module),
                                    function(x) {
                                      
                                      # Get the vector of candidate genes
                                      # based on the module assignment:
                                      candidate_genes <- 
                                        brainvar_wgcna_sheets$`Module Genes`$gene_id[
                                          brainvar_wgcna_sheets$`Module Genes`$Module == x
                                          ]
                                      
                                      # Now go through each SNP in the data,
                                      # and see if any of the genes are 
                                      # members:
                                      membership <- sapply(snp_gene_list,
                                                           function(y) {
                                                             as.numeric(any(str_sub(y, 1, 15) %in% 
                                                                              candidate_genes))
                                                           })
                                      result <- data.frame(members = membership)
                                      colnames(result) <- paste0("brainvar_any_gene_", x)
                                      return(result)
                                    })

# Join this data and save:
bip_scz_data_14_18 <- bip_scz_data_14_18 %>%
  bind_cols(brainvar_any_gene_member)

# Now finally create the cumulative position of the SNP for each chromosome
# following the same code here https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html

# Compute chromosome size for the 2018 BD values
bip_scz_data_14_18 <- bip_scz_data_14_18 %>%
  # Create an index to return the ordering - because the position columns
  # were added after the fact...
  mutate(snp_index = 1:n()) 
bip_scz_data_14_18 <- bip_scz_data_14_18 %>%
  group_by(bip_18_CHR) %>% 
  summarise(chr_len = max(bip_18_BP, na.rm = TRUE)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(bip_scz_data_14_18, ., 
            by = "bip_18_CHR") %>%
  # Add a cumulative position of each SNP
  arrange(bip_18_CHR, bip_18_BP) %>%
  mutate(bip_18_BP_cum = bip_18_BP + tot) %>%
  dplyr::select(-tot) %>% 
  arrange(snp_index) %>%
  dplyr::select(-snp_index)

# Save the data:
# readr::write_csv(bip_scz_data_14_18, 
#                  path = "data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")


