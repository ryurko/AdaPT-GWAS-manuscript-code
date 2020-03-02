# PURPOSE: Create a single Excel file containing the sheets of discoveries for
#          each of the GWAS results considered (SCZ, T2D, and BMI)
# ------------------------------------------------------------------------------

# Access packages
library(tidyverse)
library(readxl)
library(xlsx)

# Define the function used to read in excel sheets
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# ------------------------------------------------------------------------------
# Initialize a single results workbooks that will store all of the sheets of results:
results_wb <- createWorkbook()

# All discoveries will be at target FDR level $\alpha = 0.05$

# ------------------------------------------------------------------------------
# Create the discovery sheets for the 2014 SCZ BrainVar results:

# Load the SCZ BrainVar data:
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Load the 2014 model results:
scz_brainvar_2014_results <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")

# Create a dataframe of just the 2014 discoveries using the gradient boosted trees
# model with necessary contextual information:
scz_brainvar_14_esnp_discoveries <- bip_scz_brainvar_data[which(scz_brainvar_2014_results$qvals <= 0.05),] %>%
  dplyr::select(SNP, scz_18_CHR, scz_18_BP, 
                scz_14_P, scz_14_OR, scz_14_SE, 
                scz_14_A1, scz_14_A2,
                new_p_scz_18, new_beta_scz_18, new_se_scz_18,
                scz_18_P, scz_18_OR, scz_18_SE) %>%
  mutate(scz_new_18_or = exp(new_beta_scz_18),
         scz_14_adapt_qvalue = scz_brainvar_2014_results$qvals[which(scz_brainvar_2014_results$qvals <= 0.05)]) %>%
  dplyr::select(SNP, scz_18_CHR, scz_18_BP, 
                scz_14_adapt_qvalue,
                scz_14_P, scz_14_OR, scz_14_SE, 
                scz_14_A1, scz_14_A2,
                new_p_scz_18, scz_new_18_or, new_se_scz_18,
                scz_18_P, scz_18_OR, scz_18_SE) %>%
  rename(snp = SNP,
         chr = scz_18_CHR,
         bp = scz_18_BP,
         scz_14_pvalue = scz_14_P,
         scz_14_or = scz_14_OR,
         scz_14_se = scz_14_SE,
         scz_14_effect_allele = scz_14_A1,
         scz_14_noneffect_allele = scz_14_A2,
         scz_new_18_pvalue = new_p_scz_18,
         scz_new_18_se = new_se_scz_18,
         scz_all_18_pvalue = scz_18_P,
         scz_all_18_or = scz_18_OR,
         scz_all_18_se = scz_18_SE)

# Next load the BrainVar eQTL and WGCNA sheets:
brainvar_eqtl_sheets <- read_excel_allsheets("data/brainvar_eqtl/eqtls.xlsx")
brainvar_wgcna_sheets <- read_excel_allsheets("data/brainvar_eqtl/wgcna.xlsx")

# Now create a vector of the eSNPs
scz_brainvar_14_disc_esnps <- scz_brainvar_14_esnp_discoveries %>% pull(snp)

# Next create dataframe of the cis-eQTL genes assoicated with with all of the 
# discovery eSNPs:
scz_brainvar_14_disc_esnp_genes <- brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05 %>%
  filter(rsID %in% scz_brainvar_14_disc_esnps) %>%
  dplyr::select(GeneId, rsID) %>%
  mutate(ensembl_id = str_sub(GeneId, 1, 15),
         gene_symbol = str_sub(GeneId, 17, nchar(GeneId))) %>%
  dplyr::select(-GeneId) %>%
  rename(snp = rsID) %>%
  # Finally left-join the eSNP information:
  left_join(scz_brainvar_14_esnp_discoveries, by = "snp")

# ------------------------------------------------------------------------------
# Next repeat but for the all-2018 discoveries
scz_brainvar_all_18_results <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz18_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")

# Create a dataframe using the gradient boosted trees model with necessary contextual information:
scz_brainvar_all_18_esnp_discoveries <- bip_scz_brainvar_data[which(scz_brainvar_all_18_results$qvals <= 0.05),] %>%
  dplyr::select(SNP, scz_18_CHR, scz_18_BP, 
                scz_14_P, scz_14_OR, scz_14_SE, 
                scz_14_A1, scz_14_A2,
                new_p_scz_18, new_beta_scz_18, new_se_scz_18,
                scz_18_P, scz_18_OR, scz_18_SE) %>%
  mutate(scz_new_18_or = exp(new_beta_scz_18),
         scz_all_18_adapt_qvalue = scz_brainvar_all_18_results$qvals[which(scz_brainvar_all_18_results$qvals <= 0.05)]) %>%
  dplyr::select(SNP, scz_18_CHR, scz_18_BP, 
                scz_all_18_adapt_qvalue,
                scz_14_P, scz_14_OR, scz_14_SE, 
                scz_14_A1, scz_14_A2,
                new_p_scz_18, scz_new_18_or, new_se_scz_18,
                scz_18_P, scz_18_OR, scz_18_SE) %>%
  rename(snp = SNP,
         chr = scz_18_CHR,
         bp = scz_18_BP,
         scz_14_pvalue = scz_14_P,
         scz_14_or = scz_14_OR,
         scz_14_se = scz_14_SE,
         scz_14_effect_allele = scz_14_A1,
         scz_14_noneffect_allele = scz_14_A2,
         scz_new_18_pvalue = new_p_scz_18,
         scz_new_18_se = new_se_scz_18,
         scz_all_18_pvalue = scz_18_P,
         scz_all_18_or = scz_18_OR,
         scz_all_18_se = scz_18_SE)

# Now create a vector of the eSNPs
scz_brainvar_all_18_disc_esnps <- scz_brainvar_all_18_esnp_discoveries %>% pull(snp)

# Next create dataframe of the cis-eQTL genes associated with with all of the 
# discovery eSNPs:
scz_brainvar_all_18_disc_esnp_genes <- brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05 %>%
  filter(rsID %in% scz_brainvar_all_18_disc_esnps) %>%
  dplyr::select(GeneId, rsID) %>%
  mutate(ensembl_id = str_sub(GeneId, 1, 15),
         gene_symbol = str_sub(GeneId, 17, nchar(GeneId))) %>%
  dplyr::select(-GeneId) %>%
  rename(snp = rsID) %>%
  # Finally left-join the eSNP information:
  left_join(scz_brainvar_all_18_esnp_discoveries, by = "snp")

# ------------------------------------------------------------------------------
# Next the SCZ GTEx 2014-only discoveries

# Load the data:
scz_gtex_esnps_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_gtex_eqtls_cortical_wgcna.csv")
scz_gtex_adapt <- readRDS("data/bip_schz_data/gtex_results/cv_tune_results/scz_gtex_cortical_s05_2cv.rds")

# Create a dataframe of just the 2014 discoveries using the gradient boosted trees
# model with necessary contextual information:
scz_gtex_esnp_discoveries <- scz_gtex_esnps_data[which(scz_gtex_adapt$qvals <= 0.05),] %>%
  dplyr::select(SNP, scz_18_CHR, scz_18_BP, 
                scz_14_P, scz_14_OR, scz_14_SE, 
                scz_14_A1, scz_14_A2,
                new_p_scz_18, new_beta_scz_18, new_se_scz_18,
                scz_18_P, scz_18_OR, scz_18_SE) %>%
  mutate(scz_new_18_or = exp(new_beta_scz_18),
         scz_14_adapt_qvalue = scz_gtex_adapt$qvals[which(scz_gtex_adapt$qvals <= 0.05)]) %>%
  dplyr::select(SNP, scz_18_CHR, scz_18_BP, 
                scz_14_adapt_qvalue,
                scz_14_P, scz_14_OR, scz_14_SE, 
                scz_14_A1, scz_14_A2,
                new_p_scz_18, scz_new_18_or, new_se_scz_18,
                scz_18_P, scz_18_OR, scz_18_SE) %>%
  rename(snp = SNP,
         chr = scz_18_CHR,
         bp = scz_18_BP,
         scz_14_pvalue = scz_14_P,
         scz_14_or = scz_14_OR,
         scz_14_se = scz_14_SE,
         scz_14_effect_allele = scz_14_A1,
         scz_14_noneffect_allele = scz_14_A2,
         scz_new_18_pvalue = new_p_scz_18,
         scz_new_18_se = new_se_scz_18,
         scz_all_18_pvalue = scz_18_P,
         scz_all_18_or = scz_18_OR,
         scz_all_18_se = scz_18_SE)

# Next create dataframe of the cis-eQTL genes associated with with all of the 
# discovery eSNPs:
scz_gtex_disc_esnp_egene <- read_csv("data/gtex_eqtl/eqtl_data.csv") %>%
  filter(eqtl_name %in% scz_gtex_esnp_discoveries$snp) %>%
  dplyr::select(eqtl_name, ensg_gene, gene_name) %>%
  rename(snp = eqtl_name,
         ensembl_id = ensg_gene,
         gene_symbol = gene_name) %>%
  left_join(scz_gtex_esnp_discoveries, by = "snp")


# ------------------------------------------------------------------------------
# Next the T2D discoveries

# Load the data:
t2d_esnps_data <- readr::read_csv("data/t2d/t2d_esnps_eqtl_slopes_wgcna.csv")
t2d_low_s0_adj_results <- readRDS("data/t2d/adapt_cv_results/t2d_adj_s05_2cv.rds")

# Create a dataframe of just the 2014 discoveries using the gradient boosted trees
# model with necessary contextual information:
t2d_gtex_esnp_discoveries <- t2d_esnps_data[which(t2d_low_s0_adj_results$qvals <= 0.05),] %>%
  dplyr::select(rsid, Chr, Pos, EA, NEA, Beta, SE, Pvalue, adj_pval) %>%
  mutate(t2d_or = exp(Beta),
         t2d_adj_adapt_qvalue = t2d_low_s0_adj_results$qvals[which(t2d_low_s0_adj_results$qvals <= 0.05)]) %>%
  dplyr::select(rsid, Chr, Pos, EA, NEA, t2d_or, Beta, SE, Pvalue, adj_pval, t2d_adj_adapt_qvalue) %>%
  rename(snp = rsid,
         chr = Chr,
         bp = Pos,
         t2d_effect_allele = EA,
         t2d_noneffect_allele = NEA,
         t2d_beta = Beta,
         t2d_se = SE,
         t2d_unadjusted_pvalue = Pvalue,
         t2d_adjusted_pvalue = adj_pval)

# Create the dataset with the associated cis-eQTL genes:
t2d_gtex_disc_esnp_egene <- read_csv("data/gtex_eqtl/eqtl_data.csv") %>%
  filter(eqtl_name %in% t2d_gtex_esnp_discoveries$snp) %>%
  dplyr::select(eqtl_name, ensg_gene, gene_name) %>%
  rename(snp = eqtl_name,
         ensembl_id = ensg_gene,
         gene_symbol = gene_name) %>%
  left_join(t2d_gtex_esnp_discoveries, by = "snp")


# ------------------------------------------------------------------------------
# Finally for the BMI discoveries
bmi_esnps_data <- read_csv("data/bmi/bmi_esnps_eqtl_slopes_wgcna.csv")
bmi_adj_15_adapt <- readRDS("data/bmi/adapt_cv_results/bmi_adj15_s05_2cv.rds")

# Create a dataframe of the 2015 discoveries
bmi_gtex_esnp_discoveries <- bmi_esnps_data[which(bmi_adj_15_adapt$qvals <= 0.05),] %>%
  dplyr::select(SNP, bmi18_chr, bmi18_pos, bmi15_A1, bmi15_A2, bmi15_b, 
                bmi15_se, bmi15_p, adj_bmi15_p) %>%
  mutate(bmi_15_adj_adapt_qvalue = bmi_adj_15_adapt$qvals[which(bmi_adj_15_adapt$qvals <= 0.05)]) %>%
  dplyr::select(SNP, bmi18_chr, bmi18_pos, bmi15_A1, bmi15_A2, bmi15_b, 
                bmi15_se, bmi15_p, adj_bmi15_p, bmi_15_adj_adapt_qvalue) %>%
  rename(snp = SNP,
         chr = bmi18_chr,
         bp = bmi18_pos,
         bmi_15_effect_allele = bmi15_A1,
         bmi_15_noneffect_allele = bmi15_A2,
         bmi_15_beta = bmi15_b,
         bmi_15_se = bmi15_se,
         bmi_15_unadjusted_pvalue = bmi15_p,
         bmi_15_adjusted_pvalue = adj_bmi15_p)

# Create the dataset with the associated cis-eQTL genes:
bmi_gtex_disc_esnp_egene <- read_csv("data/gtex_eqtl/eqtl_data.csv") %>%
  filter(eqtl_name %in% bmi_gtex_esnp_discoveries$snp) %>%
  dplyr::select(eqtl_name, ensg_gene, gene_name) %>%
  rename(snp = eqtl_name,
         ensembl_id = ensg_gene,
         gene_symbol = gene_name) %>%
  left_join(bmi_gtex_esnp_discoveries, by = "snp")

# ------------------------------------------------------------------------------
# Now finally make a column key for each table and concatenate them on top of
# another:

scz_brainvar_column_key <- data.frame("column_key" = union(colnames(scz_brainvar_14_disc_esnp_genes),
                                                           colnames(scz_brainvar_all_18_disc_esnp_genes))) %>%
  mutate(description = c("dbSNP rsID", "Ensembl ID", "Gene symbol", 
                         "Chromosome", "Genomic position (hg19)",
                         "SCZ 2014-only AdaPT q-value",
                         "SCZ 2014-only p-value", "SCZ 2014-only odds ratio based on effect allele",
                         "SCZ 2014-only standard error", "SCZ effect allele",
                         "SCZ non-effect allele", 
                         "SCZ 2018-only p-value", "SCZ 2018-only odds ratio based on effect allele",
                         "SCZ 2018-only standard error", 
                         "SCZ all-2018 p-value", "SCZ all-2018 odds ratio based on effect allele",
                         "SCZ all-2018 standard error",
                         "SCZ all-2018 AdaPT q-value"),
         sheet = c(rep("SCZ BrainVar (2014-only or all 2018)", 
                       length(union(colnames(scz_brainvar_14_disc_esnp_genes),
                                    colnames(scz_brainvar_all_18_disc_esnp_genes))))))

scz_gtex_column_key <- data.frame("column_key" = colnames(scz_gtex_disc_esnp_egene)) %>%
  mutate(description = c("dbSNP rsID", "Ensembl ID", "Gene symbol",
                         "Chromosome", "Genomic position (hg19)",
                         "SCZ 2014-only AdaPT q-value",
                         "SCZ 2014-only p-value", "SCZ 2014-only odds ratio based on effect allele",
                         "SCZ 2014-only standard error", "SCZ effect allele",
                         "SCZ non-effect allele", 
                         "SCZ 2018-only p-value", "SCZ 2018-only odds ratio based on effect allele",
                         "SCZ 2018-only standard error", 
                         "SCZ all-2018 p-value", "SCZ all-2018 odds ratio based on effect allele",
                         "SCZ all-2018 standard error"),
         sheet = c(rep("SCZ GTEx", ncol(scz_gtex_disc_esnp_egene))))

t2d_gtex_column_key <- data.frame("column_key" = colnames(t2d_gtex_disc_esnp_egene)) %>%
  mutate(description = c("dbSNP rsID", "Ensembl ID", "Gene symbol",
                         "Chromosome", "Genomic position (hg19)",
                         "T2D effect allele",
                         "T2D non-effect allele", 
                         "T2D odds ratio based on effect allele",
                         "T2D log odds ratio based on effect allele",
                         "T2D standard error",
                         "T2D unadjusted p-value",
                         "T2D adjusted p-value",
                         "T2D AdaPT q-value (using adjusted p-values)"),
         sheet = c(rep("T2D GTEx", ncol(t2d_gtex_disc_esnp_egene))))

bmi_gtex_column_key <- data.frame("column_key" = c(colnames(bmi_gtex_disc_esnp_egene))) %>%
  mutate(description = c("dbSNP rsID", "Ensembl ID", "Gene symbol",
                         "Chromosome", "Genomic position (hg19)",
                         "BMI effect allele",
                         "BMI noneffect allele", 
                         "BMI 2015-only beta based on effect allele",
                         "BMI 2015-only standard error",
                         "BMI 2015-only unadjusted p-value",
                         "BMI 2015-only adjusted p-value",
                         "BMI 2015-only AdaPT q-value (using adjusted p-values)"),
         sheet = c(rep("BMI GTEx", ncol(bmi_gtex_disc_esnp_egene))))

# Now stack them and combine duplicates:
all_column_key <- bind_rows(scz_brainvar_column_key, scz_gtex_column_key,
                            t2d_gtex_column_key, bmi_gtex_column_key) %>%
  group_by(column_key, description) %>%
  summarize(sheets = paste0(unique(sheet), collapse = ", "),
            nsheets = length(unique(sheet))) %>%
  arrange(desc(nsheets)) %>%
  dplyr::select(-nsheets) %>%
  ungroup() %>%
  # Add a blank row to the bottom of this and then a message saying that each
  # row in the various results are eSNP - cis-eQTL gene pair, mean each eSNP
  # and each gene can appear more than once, but a pairing can only appear once.
  bind_rows(data.frame("column_key" = c(" ", 
                                        "Each row in each tab of results corresponds to a unique eSNP - associated cis-eQTL gene pair, meaning each eSNP and its test statistics can appear more than once. To view the unique eSNP results returned by AdaPT just examine distinct values for the snp column. The associated cis-eQTL genes for each eSNP are provided in this format for convenience."),
                       "description" = c(" ", " "),
                       "sheets" = c(" ", " ")))

# Now make the five sheets with this information:

# (1) Column Key
column_key_sheet <- createSheet(results_wb, "Column Key")
addDataFrame(as.data.frame(all_column_key), 
             sheet = column_key_sheet, 
             startColumn = 1, row.names = FALSE)

# (2) BrainVar eSNPs for 2014-only results:
scz_14_sheet <- createSheet(results_wb, "SCZ BrainVar 2014-only")
addDataFrame(as.data.frame(scz_brainvar_14_disc_esnp_genes), 
             sheet = scz_14_sheet, 
             startColumn = 1, row.names = FALSE)

# (3) BrainVar eSNPs for all 2018 results:
scz_18_sheet <- createSheet(results_wb, "SCZ BrainVar all 2018")
addDataFrame(as.data.frame(scz_brainvar_all_18_disc_esnp_genes), 
             sheet = scz_18_sheet, 
             startColumn = 1, row.names = FALSE)

# (4) GTEx eSNPs:
scz_gtex_sheet <- createSheet(results_wb, "SCZ GTEx")
addDataFrame(as.data.frame(scz_gtex_disc_esnp_egene), 
             sheet = scz_gtex_sheet, 
             startColumn = 1, row.names = FALSE)

# (5) T2D GTEx eSNPs
t2d_gtex_sheet <- createSheet(results_wb, "T2D GTEx")
addDataFrame(as.data.frame(t2d_gtex_disc_esnp_egene), 
             sheet = t2d_gtex_sheet, 
             startColumn = 1, row.names = FALSE)

# (6) GTEx eSNPs
bmi_gtex_sheet <- createSheet(results_wb, "BMI GTEx")
addDataFrame(as.data.frame(bmi_gtex_disc_esnp_egene), 
             sheet = bmi_gtex_sheet, 
             startColumn = 1, row.names = FALSE)

# Save it:
saveWorkbook(results_wb, 
             "data/adapt_gwas_results.xlsx")

