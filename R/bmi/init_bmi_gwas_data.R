# PURPOSE: Create the initial dataset from the BMI studies in 2015 with 
#          replication results from the more recent studies in 2018

# Author: Ron Yurko

# Access necessary packages
library(tidyverse)
library(data.table)

# Load the three data files:
bmi_18 <- read.table("data/bmi/bmi_giantuk_18.txt", header = TRUE)
bmi_15 <- read.table("data/bmi/bmi_euro_15.txt", header = TRUE)
whr_15 <- read.table("data/bmi/whr_euro_15.txt", header = TRUE)

# Rename the main columns to match:
bmi_18 <- bmi_18 %>% rename(A1 = Tested_Allele, A2 = Other_Allele,
                            b = BETA, se = SE, p = P,
                            chr = CHR, pos = POS, ref_af = Freq_Tested_Allele_in_HRS)
bmi_15 <- bmi_15 %>% rename(ref_af = Freq1.Hapmap)
whr_15 <- whr_15 %>% rename(SNP = MarkerName, A1 = Allele1, A2 = Allele2,
                            chr = Chr, pos = Pos, ref_af = FreqAllele1HapMapCEU)

# Vector SNPs common in the BMI 18 and BMI 15:
snps_bmi <- intersect(unique(str_trim(bmi_18$SNP)),
                      unique(str_trim(bmi_15$SNP)))

# Vector of SNPs common in the BMI 18 and WHR 15:
snps_bmi18_whr <- intersect(unique(str_trim(bmi_18$SNP)),
                            unique(str_trim(whr_15$SNP)))

# Vector of SNPs common across all three:
snps_all <- intersect(snps_bmi18_whr,
                      unique(str_trim(bmi_15$SNP)))

# Compare the betas of the SNPs that are in all three datasets
all_bmi_18 <- bmi_18 %>%
  dplyr::select(SNP, b, p, A1, A2, se, N, chr, pos, ref_af) %>%
  filter(SNP %in% snps_all) %>%
  mutate(dataset = "bmi18",
         z_score = b / se)
all_whr_15 <- whr_15 %>%
  dplyr::select(SNP, b, p, A1, A2, se, N, chr, pos, ref_af) %>%
  filter(SNP %in%  snps_all) %>%
  mutate(dataset = "whr15",
         z_score = b / se)
all_bmi_15 <- bmi_15 %>%
  dplyr::select(SNP, b, p, A1, A2, se, N, ref_af) %>%
  filter(SNP %in%  snps_all) %>%
  mutate(dataset = "bmi15",
         z_score = b / se)

# First need to remove duplicates in each:
bmi_dups <- all_bmi_18 %>%
  group_by(SNP) %>%
  count() %>%
  filter(n > 1) %>%
  pull(SNP) %>%
  as.character()
# there are 9 duplicates

whr_dups <- all_whr_15 %>%
  group_by(SNP) %>% 
  count() %>%
  filter(n > 1) %>%
  pull(SNP) %>%
  as.character()
#  0 duplicates

bmi_dups_15 <- all_bmi_15 %>%
  group_by(SNP) %>%
  count() %>%
  filter(n > 1) %>%
  pull(SNP) %>%
  as.character()
# 9 duplicates

# Two bmi datasets have the same duplicates

# The 9 duplicates for the 2018 data are just repeats
# Only 2 of the 11 duplicates for the 2015 data are flipped tests

flipped_snps <- all_bmi_15 %>%
  filter(SNP %in% bmi_dups_15) %>%
  arrange(SNP) %>%
  group_by(SNP) %>%
  summarise(n_unique_b = length(unique(b))) %>%
  filter(n_unique_b > 1) %>%
  pull(SNP) %>%
  as.character()

# Might be easier just to work with each separately given the size of the data
# and then just rename and join the columns:
rm(bmi_18, whr_15, bmi_15)

clean_all_bmi_18 <- all_bmi_18 %>%
  # Grab distinct rows
  distinct() %>%
  # Mark the duplicates
  mutate(dups = if_else(SNP %in% flipped_snps, 1, 0)) %>%
  # Now group by the SNP and arrange by the SNP and A1:
  group_by(SNP, dataset) %>%
  arrange(SNP, dataset, A1) %>%
  # Which allele number is it?
  mutate(n_allele = 1:n()) %>%
  # Remove any observations that are duplicates and n_allele greater than 1
  filter(dups == 0 | (dups == 1 & n_allele == 1)) %>%
  ungroup() %>%
  dplyr::select(-dups, -n_allele)

clean_all_whr_15 <- all_whr_15 %>%
  # Grab distinct rows
  distinct() %>%
  # Mark the duplicates
  mutate(dups = if_else(SNP %in% flipped_snps, 1, 0)) %>%
  # Now group by the SNP and arrange by the SNP and A1:
  group_by(SNP, dataset) %>%
  arrange(SNP, dataset, A1) %>%
  # Which allele number is it?
  mutate(n_allele = 1:n()) %>%
  # Remove any observations that are duplicates and n_allele greater than 1
  filter(dups == 0 | (dups == 1 & n_allele == 1)) %>%
  ungroup() %>%
  dplyr::select(-dups, -n_allele)

clean_all_bmi_15 <- all_bmi_15 %>%
  # Grab distinct rows
  distinct() %>%
  # Mark the duplicates
  mutate(dups = if_else(SNP %in% flipped_snps, 1, 0)) %>%
  # Now group by the SNP and arrange by the SNP and A1:
  group_by(SNP, dataset) %>%
  arrange(SNP, dataset, A1) %>%
  # Which allele number is it?
  mutate(n_allele = 1:n()) %>%
  # Remove any observations that are duplicates and n_allele greater than 1
  filter(dups == 0 | (dups == 1 & n_allele == 1)) %>%
  ungroup() %>%
  dplyr::select(-dups, -n_allele)

rm(all_bmi_18, all_whr_15, all_bmi_15)

# Next create three datasets with the columns joined together:

# First make a pipeline to tidy the datasets:
tidy_bmi_data <- . %>%
  gather(variable, value, -SNP, -dataset) %>%
  unite(temp, dataset, variable) %>%
  spread(temp, value)

# Apply this to each one:
bmi_18_data <- clean_all_bmi_18 %>% tidy_bmi_data
whr_15_data <- clean_all_whr_15 %>% tidy_bmi_data  
bmi_15_data <- clean_all_bmi_15 %>% tidy_bmi_data
rm(clean_all_bmi_15, clean_all_bmi_18, clean_all_whr_15, testing)

# Now join together the datasets and create the proper columns for storing the
# minor allele frequencies:
all_bmi_whr_18_data <- bmi_18_data %>% 
  inner_join(whr_15_data, by = "SNP") %>%
  mutate_at(vars(matches("(_b)|(_z_score)|(_p)|(_se)|(_N)|(_ref_af)|(_chr)|(_pos)")), as.numeric) %>%
  filter((bmi18_A1 == whr15_A1 & bmi18_A2 == whr15_A2) |
           (bmi18_A1 == whr15_A2 & bmi18_A2 == whr15_A1)) %>%
  mutate(allele_flip = if_else(bmi18_A1 == whr15_A2, TRUE, FALSE),
         bmi18_b = if_else(allele_flip, -1 * bmi18_b, bmi18_b),
         bmi18_z_score = if_else(allele_flip, -1 * bmi18_z_score, bmi18_z_score)) %>%
  # Rename bmi18_maf to be the bmi18_ref_af and then create the actual maf column
  # using whether or not the ref_af is over .5:
  mutate(bmi18_maf = ifelse(bmi18_ref_af >= .5, 1 - bmi18_ref_af, bmi18_ref_af),
         whr15_maf = ifelse(whr15_ref_af >= .5, 1 - whr15_ref_af, whr15_ref_af))

common_bmi_data <- bmi_15_data %>% 
  inner_join(bmi_18_data, by = "SNP") %>%
  mutate_at(vars(matches("(_b)|(_z_score)|(_p)|(_se)|(_N)|(_ref_af)|(_chr)|(_pos)")), as.numeric) %>%
  filter((bmi18_A1 == bmi15_A1 & bmi18_A2 == bmi15_A2) |
           (bmi18_A1 == bmi15_A2 & bmi18_A2 == bmi15_A1)) %>%
  mutate(allele_flip = if_else(bmi18_A1 == bmi15_A2, TRUE, FALSE),
         bmi18_b = if_else(allele_flip, -1 * bmi18_b, bmi18_b),
         bmi18_z_score = if_else(allele_flip, -1 * bmi18_z_score, bmi18_z_score)) %>%
  mutate(bmi18_maf = ifelse(bmi18_ref_af >= .5, 1 - bmi18_ref_af, bmi18_ref_af))


# Now only consider the SNPs where the maf are .05 in both the BMI18 and WHR 15
# ( since these are the same as all three )
screened_bmi_whr_18_data <- all_bmi_whr_18_data %>%
  filter(bmi18_maf >= .05, whr15_maf >= .05)

# Save this dataset to use for the covariates
# write_csv(screened_bmi_whr_18_data, "data/bmi/screened_bmi_18_whr_15.csv")

# make the main dataset
bmi_data <- common_bmi_data %>%
  filter(SNP %in% unique(screened_bmi_whr_18_data$SNP))

# Create the columns separating the new studies in 2018 from the studies in 2015 
# that are included in the total 2018 meta-analysis results

# Now perform the calculations based on the METAL algebra:
bmi_data <- bmi_data[, var_bmi18 := bmi18_se^2]
bmi_data <- bmi_data[, var_bmi15 := bmi15_se^2]
bmi_data <- bmi_data[, bmi18_new_b := (bmi18_b * var_bmi15 - bmi15_b * var_bmi18) / 
                             (var_bmi15 - var_bmi18)]
bmi_data <- bmi_data[, bmi18_new_se := sqrt((var_bmi15 * var_bmi18) / 
                                                    (var_bmi15 - var_bmi18))]
bmi_data <- bmi_data[, bmi18_new_z_score := bmi18_new_b / bmi18_new_se]
bmi_data <- bmi_data[, bmi18_new_p := 2*pnorm(-abs(bmi18_new_z_score))]

# Save this dataset:
#write_csv(bmi_data, 
#          path = "data/bmi/bmi_18_new_studies.csv")

