# PURPOSE: Create the list of associated cis-eQTL genes with the AdaPT discoveries
#          for the T2D GWAS

# Load necessary packages:
library(tidyverse)

# Load the data:
t2d_esnps_data <- readr::read_csv("data/t2d/t2d_esnps_eqtl_slopes_wgcna.csv")

# Load the final adjusted model results:
t2d_adj_adapt <- readRDS("data/t2d/adapt_cv_results/t2d_adj_s05_2cv.rds")

# How many associated genes?
read_csv("data/gtex_eqtl/eqtl_data.csv") %>%
  filter(eqtl_name %in% t2d_esnps_data$rsid[which(t2d_adj_adapt$qvals <= .05)]) %>%
  pull(ensg_gene) %>%
  unique() %>% 
  length()
# 5970

# Now save this list of genes to get the GO enrichment for
t2d_gtex_gene_file <- file("data/gene_lists/t2d_gtex_disc_genes.txt")
read_csv("data/gtex_eqtl/eqtl_data.csv") %>%
  filter(eqtl_name %in% t2d_esnps_data$rsid[which(t2d_adj_adapt$qvals <= .05)]) %>%
  pull(ensg_gene) %>%
  unique() %>% 
  paste(collapse = "\n") %>% 
  writeLines(t2d_gtex_gene_file)  
close(t2d_gtex_gene_file)
