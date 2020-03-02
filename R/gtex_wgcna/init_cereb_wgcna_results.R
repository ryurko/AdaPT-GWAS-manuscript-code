# PURPOSE: Generate WGCNA results for Cerebellum tissues:

# AUTHOR: Ron Yurko

# Access WGCNA
library(WGCNA)

# Enable parallel processing for WGCNA with five threads 
# NOTE enableWGCNAThreads does not work - but this is a contradiction with the
# documentation unfortunately...
allowWGCNAThreads(nThreads = 5)

# Load frontal cortex samples:
cereb_tpm_data <- readr::read_csv("data/gtex_eqtl/expression_datasets/tpm_cerebellum_samples.csv")

# Now access the grex package in order to filter the genes to only be protein-coding
library(grex)

# Also access tidyverse
library(tidyverse)

# Since the gene IDs provided in the raw GTEx data are in the GENCODE form, 
# need to remove the .version part of them which can be done with the 
# cleanid() function in grex:
cereb_tpm_data <- cereb_tpm_data %>%
  mutate(clean_name = cleanid(Name))

# Now use grex to generate the mapping table of the genes (they are same set
# of genes so can just do once):
gtex_gene_table <- grex(cereb_tpm_data$clean_name)

# Get the protein_coding genes:
protein_coding_genes <- gtex_gene_table %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(ensembl_id)

# Filter to only include these genes for the data:
cereb_tpm_data <- cereb_tpm_data %>%
  filter(clean_name %in% protein_coding_genes)

# Let's remove all genes with 0 expression for over 50% of the provided samples - 
# this follows what Bert said (also captures genes with 0s for everything)
zero_half_cereb_expression <- which(apply(dplyr::select(cereb_tpm_data,
                                                           -c(Name, clean_name,
                                                              Description)), 1,
                                             function(x) length(which(x == 0)) / length(x) >= .5))

# Now remove these genes
cereb_tpm_data <- cereb_tpm_data[-zero_half_cereb_expression,]

# Grab only the samples columns and tranpose so each row is a sample, but
# take the log transformation of the expressions (TPMs) + 1:
cereb_log_tpm <- apply(dplyr::select(zero_half_cereb_expression,
                                        -c(Name, clean_name,
                                           Description)),
                          1, function(x) log(x + 1, base = 2))


# Set the column names to be the gene names:
colnames(cereb_log_tpm) <- cereb_tpm_data$Name

# Now need to choose which power to use based on
# the scale-free topology criterion
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))

# We're going to use unsigned now
cereb_power_search <- pickSoftThreshold(cereb_log_tpm, 
                                           dataIsExpr = TRUE, 
                                           powerVector = powers,
                                           corFnc = cor,
                                           corOptions = list(use = 'p'),
                                           networkType = "unsigned")

# Create two plots showing the R-squared and connectivity:
cereb_power_search$fitIndices %>%
  ggplot(aes(x = Power, y = SFT.R.sq)) +
  geom_hline(yintercept = 0.8, color = "darkorange") +
  geom_text(aes(label = Power), color = "darkblue",
            size = 5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Soft threshold power", 
       y = "Scale free topology model fit R^2",
       title = "Scale independence chart for cerebellum tissue samples") +
  theme_bw()

# Will use a power of 7

# Now generate the WGCNA results using the blockwise implementation to run
# locally with the selected powers from above:

# Set the proper directory to the wcgna_results folder
#setwd("data/gtex_eqtl/wcgna_results")
cereb_bwnet <- blockwiseModules(cereb_log_tpm, maxBlockSize = 20000,
                                   power = 7, TOMType = "unsigned", minModuleSize = 30,
                                   numericLabels = TRUE,
                                   saveTOMs = TRUE,
                                   saveTOMFileBase = "cereb_TOM_blockwise",
                                   verbose = 3)
#setwd("data")

# Convert labels to colors for plotting
cereb_mergedColors <- labels2colors(cereb_bwnet$colors)

# Create a dataframes with the gene names and the labels:
cereb_wgcna_data <- data.frame("gene" = colnames(cereb_log_tpm),
                               "wgcna_label" = cereb_mergedColors)

# Save them to the wcgna_results folder
# write_csv(cereb_wgcna_data, 
#          "data/gtex_eqtl/wcgna_results/cerebellum_wgcna_data.csv")
