# PURPOSE: Create the GTEx eSNPs dataset used for AdaPT with the appropriate covariates

# AUTHOR: Ron Yurko

# Access necessary packages
library(data.table)
library(tidyverse)

# First load the BMI data:
bmi_data <- read_csv("data/bmi/bmi_18_new_studies.csv",
                     progress = FALSE)

# Load the eqtl_data:
eqtl_data <- read_csv("data/gtex_eqtl/eqtl_data.csv", progress = FALSE)

# Filter the BMI data to be the eQTLs
bmi_esnps_data <- bmi_data %>%
  filter(SNP %in% unique(eqtl_data$eqtl_name))
rm(bmi_data)

# Filter the eQTL data to only be these eSNPs and then either brain or adipose
# tissue eQTLs
adipose_brain_eqtl_data <- eqtl_data %>%
  filter(eqtl_name %in% bmi_esnps_data$SNP,
         str_detect(region, "(Brain)|(Adipose)")) %>%
  # Remove the duplicate Brain tissues - Brain_Cortex and Brain_Cerebellum
  # (since these are not measure the same time as the other brain tissues):
  filter(!(region %in% c("Brain_Cortex", "Brain_Cerebellum")))

# Next summarize these SNPs across their eQTL gene pairs and region
region_esnp_summary <- adipose_brain_eqtl_data %>%
  mutate(region = tolower(region)) %>%
  group_by(eqtl_name, region) %>%
  summarize(n_genes = n(),
            n_gene_ids = length(unique(ensg_gene)),
            ave_abs_eqtl_slope = mean(abs(eqtl_slope), na.rm = TRUE),
            max_abs_eqtl_slope = max(abs(eqtl_slope), na.rm = TRUE)) %>%
  gather(variable, value, -eqtl_name, -region) %>%
  unite(region_variable, "region", "variable") %>%
  spread(region_variable, value)
region_esnp_summary[is.na(region_esnp_summary)] <- 0

# Now versions for only Adipose and Brain (excluding cerebellum) aggregates:
aggregate_esnp_summary <- adipose_brain_eqtl_data %>%
  filter(region != "Brain_Cerebellar_Hemisphere") %>%
  mutate(region = tolower(region),
         region = str_extract(region, "(adipose)|(brain)")) %>%
  group_by(eqtl_name, region) %>%
  summarize(n_genes = n(),
            n_gene_ids = length(unique(ensg_gene)),
            ave_abs_eqtl_slope = mean(abs(eqtl_slope), na.rm = TRUE),
            max_abs_eqtl_slope = max(abs(eqtl_slope), na.rm = TRUE)) %>%
  gather(variable, value, -eqtl_name, -region) %>%
  unite(region_variable, "region", "variable") %>%
  spread(region_variable, value)
aggregate_esnp_summary[is.na(aggregate_esnp_summary)] <- 0

# Now left-join both of these datasets to the eSNP dataset:
bmi_esnps_data <- bmi_esnps_data %>%
  left_join(region_esnp_summary, by = c("SNP" = "eqtl_name")) %>%
  left_join(aggregate_esnp_summary, by = c("SNP" = "eqtl_name"))
bmi_esnps_data[is.na(bmi_esnps_data)] <- 0

# Next need to join the WGCNA information for the cerebellum, adipose, and
# brain (non-cerebellum tissues):
wgcna_cerebellum_data <- read_csv("data/gtex_eqtl/wcgna_results/cerebellum_wgcna_data.csv")
wgcna_adipose_data <- read_csv("data/gtex_eqtl/wcgna_results/adipose_wgcna_data.csv")
wgcna_brain_data <- read_csv("data/gtex_eqtl/wcgna_results/brain_wgcna_data.csv")

# First go through each SNP in the model dataset and construct a list corresponding
# to the vector of genes for each:
esnp_gene_list <- map(bmi_esnps_data$SNP,
                      function(x) {
                        adipose_brain_eqtl_data %>%
                          filter(eqtl_name == x) %>%
                          pull(ensg_gene)
                      })
names(esnp_gene_list) <- bmi_esnps_data$SNP

# make the any gene member datasets for each:

# First for cerebellum:
cerebellum_any_gene_member <- map_dfc(unique(wgcna_cerebellum_data$wgcna_label),
                                      function(x) {
                                        
                                        # Get the vector of candidate genes
                                        # based on the module assignment:
                                        candidate_genes <- 
                                          wgcna_cerebellum_data$gene[
                                            wgcna_cerebellum_data$wgcna_label == x
                                            ]
                                        
                                        # Now go through each SNP in the data,
                                        # and see if any of the genes are 
                                        # members:
                                        membership <- sapply(esnp_gene_list,
                                                             function(y) {
                                                               as.numeric(any(y %in% candidate_genes))
                                                             })
                                        result <- data.frame(members = membership)
                                        colnames(result) <- paste0("cerebellum_any_gene_", x)
                                        return(result)
                                      })

# Next for adipose:
adipose_any_gene_member <- map_dfc(unique(wgcna_adipose_data$wgcna_label),
                                   function(x) {
                                     
                                     # Get the vector of candidate genes
                                     # based on the module assignment:
                                     candidate_genes <- 
                                       wgcna_adipose_data$gene[
                                         wgcna_adipose_data$wgcna_label == x
                                         ]
                                     
                                     # Now go through each SNP in the data,
                                     # and see if any of the genes are 
                                     # members:
                                     membership <- sapply(esnp_gene_list,
                                                          function(y) {
                                                            as.numeric(any(y %in% candidate_genes))
                                                          })
                                     result <- data.frame(members = membership)
                                     colnames(result) <- paste0("adipose_any_gene_", x)
                                     return(result)
                                   })

# The brain non-cerebellum tissues:
brain_any_gene_member <- map_dfc(unique(wgcna_brain_data$wgcna_label),
                                 function(x) {
                                   
                                   # Get the vector of candidate genes
                                   # based on the module assignment:
                                   candidate_genes <- 
                                     wgcna_brain_data$gene[
                                       wgcna_brain_data$wgcna_label == x
                                       ]
                                   
                                   # Now go through each SNP in the data,
                                   # and see if any of the genes are 
                                   # members:
                                   membership <- sapply(esnp_gene_list,
                                                        function(y) {
                                                          as.numeric(any(y %in% candidate_genes))
                                                        })
                                   result <- data.frame(members = membership)
                                   colnames(result) <- paste0("brain_any_gene_", x)
                                   return(result)
                                 })

# Join all of these columns
bmi_esnps_data <- bmi_esnps_data %>%
  bind_cols(cerebellum_any_gene_member) %>%
  bind_cols(adipose_any_gene_member) %>%
  bind_cols(brain_any_gene_member)

# Load the information with the WHR tests to join
bmi_whr_studies <- read_csv("data/bmi/screened_bmi_18_whr_15.csv")

# Join the WHR p-values and z-scores:
bmi_esnps_data <- bmi_esnps_data %>%
  left_join(dplyr::select(bmi_whr_studies, SNP, whr15_p, whr15_z_score),
            by = "SNP")

# Now for this dataset create an adjusted p-value columns modifying the
# values that are exactly 1

# First see what are the maximums for each of these (less than 1):
max(bmi_esnps_data$bmi15_p[which(bmi_esnps_data$bmi15_p < 1)])
# [1] 0.9984
table(bmi_esnps_data$bmi15_p[which(bmi_esnps_data$bmi15_p > .97)])
# 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99    1 
# 1333 1306 1343 1365 1397 1361 1133 1337 1476  640 
# For the data with all SNPs let's use 0.97 as the lower bound

set.seed(1389)
bmi_esnps_data <- bmi_esnps_data %>%
  mutate(adj_bmi15_p = ifelse(bmi15_p == 1, 
                              runif(nrow(bmi_esnps_data),
                                    min = 0.97, 1 - 1e-15), bmi15_p))

# Save this final dataset:
# write_csv(bmi_esnps_data,
#           "data/bmi/bmi_esnps_eqtl_slopes_wgcna.csv")





