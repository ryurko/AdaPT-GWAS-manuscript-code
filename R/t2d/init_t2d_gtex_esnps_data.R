# PURPOSE: Create the GTEx eSNPs datasets with covariates for AdaPT

# AUTHOR: Ron Yurko

# Access necessary packages
library(data.table)
library(tidyverse)

# Load the eQTL data that was generated before - has eQTLs for any tissue type
# and then see which of the eQTLs are in the datasets
egenes_eqtl_data <- fread("data/gtex_eqtl/eqtl_data.csv") 
# Create a column that is the CHR:POS:
egenes_eqtl_data <- egenes_eqtl_data %>%
  .[, ':=' (chr_pos = paste0(eqtl_chr, ":", eqtl_pos),
            allele_seq = paste0(eqtl_ref, ":", eqtl_alt))] 

# Load the 
t2d_data <- fread("unzip -p data/t2d/Mahajan.NatGenet2018b.T2D.European.zip") %>%
  .[, ':=' (allele_seq1 = paste0(EA, ":", NEA),
            allele_seq2 = paste0(NEA, ":", EA))]

# Now filter downn to the SNPs that are in the GTEx eQTLs:
t2d_eqtl_data <- t2d_data %>%
  .[SNP %in% egenes_eqtl_data$chr_pos,]
rm(t2d_data)

# Get the duplicates:
t2d_duplicated_esnps_i <- which(duplicated(t2d_eqtl_data$SNP))
t2d_duplicated_esnps <- unique(t2d_eqtl_data$SNP[t2d_duplicated_esnps_i])
# OK so there are CHR:POS where there are two different effect alleles - 
# simplest step is determine which match the GTEx data exactly in terms of 
# the effect and reference alleles:

# Loop through each row of the t2d_eqtl_data and find the matching position
# in the eQTL data, and see if it matches the alleles (or is flipped):
allele_match <- sapply(1:nrow(t2d_eqtl_data),
                       function(snp_i) {
                         eqtl_snp_data <- egenes_eqtl_data %>%
                           .[chr_pos == t2d_eqtl_data$SNP[snp_i]]
                         any(eqtl_snp_data$allele_seq == 
                               t2d_eqtl_data$allele_seq1[snp_i]) |
                           any(eqtl_snp_data$allele_seq == 
                                 t2d_eqtl_data$allele_seq2[snp_i])
                       })
length(which(allele_match))
# [1] 176278

# Only grab those that match:
t2d_eqtl_data_match <- t2d_eqtl_data[which(allele_match),]
t2d_duplicated_match_esnps_i <- which(duplicated(t2d_eqtl_data_match$SNP))

# Now for each of these - loop through the eQTL dataset to get the unique rs id:
rs_ids <- lapply(1:nrow(t2d_eqtl_data_match),
                 function(snp_i) {
                   eqtl_snp_data <- egenes_eqtl_data %>%
                     .[chr_pos == t2d_eqtl_data_match$SNP[snp_i]]
                   eqtl_snp_data$eqtl_name[which((eqtl_snp_data$allele_seq == 
                                                    t2d_eqtl_data_match$allele_seq1[snp_i]) |
                                                   (eqtl_snp_data$allele_seq == 
                                                      t2d_eqtl_data_match$allele_seq2[snp_i]))]
                 })

unique_rs_ids_list <- lapply(rs_ids, unique)
t2d_eqtl_data_match$rsid <- unlist(unique_rs_ids_list)

# Only use the eSNPs with non "." RS IDs:
t2d_esnps_data <- t2d_eqtl_data_match %>%
  .[rsid != ".",]

# Load the egenes for any tissue type but SNPs from the T2D eSNPs and the region
# is either Pancreas, Liver, or one of the two Adipose tissues
t2d_egenes_eqtl_data <- fread("data/gtex_eqtl/eqtl_data.csv") %>%
  .[eqtl_name %in% t2d_esnps_data$rsid &
      str_detect(region, "(Adipose)|(Liver)|(Pancreas)"),]

# Next summarize these SNPs across their eQTL gene pairs and region
t2d_esnp_summary <- t2d_egenes_eqtl_data %>%
  mutate(region = tolower(region)) %>%
  group_by(eqtl_name, region) %>%
  summarize(n_genes = n(),
            n_gene_ids = length(unique(ensg_gene)),
            ave_abs_eqtl_slope = mean(abs(eqtl_slope), na.rm = TRUE),
            max_abs_eqtl_slope = max(abs(eqtl_slope), na.rm = TRUE)) %>%
  gather(variable, value, -eqtl_name, -region) %>%
  unite(region_variable, "region", "variable") %>%
  spread(region_variable, value)
t2d_esnp_summary[is.na(t2d_esnp_summary)] <- 0

# Make one that aggregates over the Adipose tissues:
t2d_esnp_adipose_summary <- t2d_egenes_eqtl_data %>%
  filter(str_detect(region, "Adipose")) %>%
  mutate(region = "adipose") %>%
  group_by(eqtl_name, region) %>%
  summarize(n_genes = n(),
            n_gene_ids = length(unique(ensg_gene)),
            ave_abs_eqtl_slope = mean(abs(eqtl_slope), na.rm = TRUE),
            max_abs_eqtl_slope = max(abs(eqtl_slope), na.rm = TRUE)) %>%
  gather(variable, value, -eqtl_name, -region) %>%
  unite(region_variable, "region", "variable") %>%
  spread(region_variable, value)
t2d_esnp_adipose_summary[is.na(t2d_esnp_adipose_summary)] <- 0

# Now left-join both of these datasets to the original dataset:
t2d_esnps_data <- t2d_esnps_data %>%
  as.data.frame() %>%
  left_join(t2d_esnp_summary, by = c("rsid" = "eqtl_name")) %>%
  left_join(t2d_esnp_adipose_summary, by = c("rsid" = "eqtl_name"))
t2d_esnps_data[is.na(t2d_esnps_data)] <- 0

# Load the pancreas WGCNA network data (just the labels for now):
wgcna_pancreas_data <- read_csv("data/gtex_eqtl/wcgna_results/pancreas_wgcna_data.csv")

# First go through each SNP in the model dataset and construct a list corresponding
# to the vector of genes for each:
all_esnp_gene_list <- map(t2d_esnps_data$rsid,
                          function(x) {
                            t2d_egenes_eqtl_data[eqtl_name == x, ensg_gene]
                          })
names(all_esnp_gene_list) <- t2d_esnps_data$rsid

# Any gene-pair member of Pancreas tissue WGCNA modules?
pancreas_any_gene_member <- map_dfc(unique(wgcna_pancreas_data$wgcna_label),
                                    function(x) {
                                      
                                      # Get the vector of candidate genes
                                      # based on the module assignment:
                                      candidate_genes <- 
                                        wgcna_pancreas_data$gene[
                                          wgcna_pancreas_data$wgcna_label == x
                                          ]
                                      
                                      # Now go through each SNP in the data,
                                      # and see if any of the genes are 
                                      # members:
                                      membership <- sapply(all_esnp_gene_list,
                                                           function(y) {
                                                             as.numeric(any(y %in% candidate_genes))
                                                           })
                                      result <- data.frame(members = membership)
                                      colnames(result) <- paste0("pancreas_any_gene_", x)
                                      return(result)
                                    })

t2d_esnps_data <- t2d_esnps_data %>%
  bind_cols(pancreas_any_gene_member)

# Save this dataset:
# write_csv(t2d_esnps_data,
#           "data/t2d/t2d_esnps_eqtl_slopes_wgcna.csv")



