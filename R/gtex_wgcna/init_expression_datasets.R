# PURPOSE: Initializes the expression datasets with just the specific
# tissue types for the WGCNA analysis using the reference GTEx file containing
# the tissue type information for each of the samples measured from

# AUTHOR: Ron Yurko

# Load necessary packages
library(data.table)
library(tidyverse)

# See data/gtex_rna_seq/GTEx_v7_Annotations_SampleAttributesDS.txt for all files and see
# data/gtex_rna_seq/GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx for an
# appropriate glossary
ref_file <- fread("data/gtex_rna_seq/GTEx_v7_Annotations_SampleAttributesDS.txt",
                  verbose = FALSE)

ref_file_glossary <- readxl::read_excel(
  "data/gtex_rna_seq/GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx")

# This file initializes the expression datasets with just the specific
# tissue types for the WGCNA analysis using the reference file containing
# the tissue type information for each of the samples measured from

# Read in the reference table and then filtering on the SMTSD field to the desired tissue types:

# 1) Cerebellar Hemisphere
cerebellum_samples <- fread("data/gtex_rna_seq/GTEx_v7_Annotations_SampleAttributesDS.txt",
                            verbose = FALSE) %>%
  filter(str_detect(SMTSD, "Cerebellar Hemisphere")) %>%
  pull(SAMPID)

# 2) Adipose
adipose_samples <- fread("data/gtex_rna_seq/GTEx_v7_Annotations_SampleAttributesDS.txt",
                         verbose = FALSE) %>%
  filter(str_detect(SMTSD, "Adipose")) %>%
  pull(SAMPID)

# 3) Pancrease
pancreas_samples <- fread("data/gtex_rna_seq/GTEx_v7_Annotations_SampleAttributesDS.txt",
                          verbose = FALSE) %>%
  filter(str_detect(SMTSD, "Pancreas")) %>%
  pull(SAMPID)

# 4) Cortical tissue samples: either Brain - Anterior cingulate cortex (BA24) or Brain - Frontal Cortex (BA9) 
cortical_samples <- fread("data/gtex_rna_seq/GTEx_v7_Annotations_SampleAttributesDS.txt",
                          verbose = FALSE) %>%
  filter(str_detect(SMTSD, "(BA24)|(BA9)")) %>% 
  pull(SAMPID)


# Next needd to load the TPM file, only selecting columns matching these
# samples and saving the files:
tpm_counts <- fread(input = paste0("zcat < ", 
                                   "data/gtex_rna_seq/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"),
                    verbose = FALSE)
tpm_colnames <- colnames(tpm_counts)

# Next save each of the respective TPM files by first matching on columns then saving:
# 1) Cerebellar Hemisphere
match_cerebellum_samples <- cerebellum_samples[which(cerebellum_samples %in% tpm_colnames)]
tpm_counts[, c("Name", "Description", match_cerebellum_samples), with = FALSE] %>%
  write_csv("data/gtex_eqtl/expression_datasets/tpm_cerebellum_samples.csv")

# 2) Adipose
match_adipose_samples <- adipose_samples[which(adipose_samples %in% tpm_colnames)]
tpm_counts[, c("Name", "Description", match_adipose_samples), with = FALSE] %>%
  write_csv("data/gtex_eqtl/expression_datasets/tpm_adipose_samples.csv")

# 3) Pancreas
match_pancreas_samples <- pancreas_samples[which(pancreas_samples %in% tpm_colnames)]
tpm_counts[, c("Name", "Description", match_pancreas_samples), with = FALSE] %>%
  write_csv("data/gtex_eqtl/expression_datasets/tpm_pancreas_samples.csv")

# 4) Cortical 
match_cortical_samples <- cortical_samples[which(cortical_samples %in% tpm_colnames)]
tpm_counts[, c("Name", "Description", match_cortical_samples), with = FALSE] %>%
  write_csv("data/gtex_eqtl/expression_datasets/tpm_cortical_samples.csv")

# Create results using all Brain tissues outside of cerebellum - but also don't
# use the duplicate for frontal cortex - only use the BA9 version

brain_samples <- fread("data/gtex_rna_seq/GTEx_v7_Annotations_SampleAttributesDS.txt",
                       verbose = FALSE) %>%
  filter(str_detect(SMTSD, "Brain")) %>%
  # Remove the Cerebellum and duplicate cortex tissues:
  filter(!(SMTSD %in% c("Brain - Cerebellar Hemisphere",
                        "Brain - Cerebellum",
                        "Brain - Cortex"))) %>% 
  pull(SAMPID)

# Now proceed to generate the TPM dataset for the brain samples found:
match_brain_samples <- brain_samples[which(brain_samples %in% tpm_colnames)]
tpm_counts[, c("Name", "Description", match_brain_samples), with = FALSE] %>%
  write_csv("data/gtex_eqtl/expression_datasets/tpm_brain_samples.csv")
