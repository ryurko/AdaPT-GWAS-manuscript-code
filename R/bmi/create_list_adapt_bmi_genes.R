# PURPOSE: Create the list of associated cis-eQTL genes with the AdaPT discoveries
#          for the BMI GWAS

# Load necessary packages:
library(tidyverse)

# Load the data:
bmi_esnps_data <- read_csv("data/bmi/bmi_esnps_eqtl_slopes_wgcna.csv")

# Load the final adjusted model results:
bmi_adj_adapt <- readRDS("data/bmi/adapt_cv_results/bmi_adj15_s05_2cv.rds")

# Now save this list of genes to get the GO enrichment for
bmi_gtex_gene_file <- file("data/gene_lists/bmi_gtex_disc_genes.txt")
read_csv("data/gtex_eqtl/eqtl_data.csv") %>%
  filter(eqtl_name %in% bmi_esnps_data$rsid[which(bmi_adj_adapt$qvals <= .05)]) %>%
  pull(ensg_gene) %>%
  unique() %>% 
  paste(collapse = "\n") %>% 
  writeLines(bmi_gtex_gene_file)  
close(bmi_gtex_gene_file)
