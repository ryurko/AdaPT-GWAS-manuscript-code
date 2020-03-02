# PURPOSE: Create the initial dataset from the SCZ1 vs BIP1 study in 2014 with 
#          replication results from the more recent studies in 2018.

# Author: Ron Yurko

# Access necessary packages
library(data.table)
library(magrittr)

# Will load the two separate datasets for BIP against controls, and SCZ against controls:
bip_data_14 <- fread("data/bip_schz_data/bip_scz_14/BPSCZ.bp-only.results.txt")
scz_data_14 <- fread("data/bip_schz_data/bip_scz_14/BPSCZ.scz-only.results.txt")

# Modify the columnn names for both datasets to have prefixes specifying
# which ones they are from except for SNP:
colnames(bip_data_14) <- sapply(colnames(bip_data_14),
                                function(x) ifelse(x == "SNP", "SNP",
                                                   paste0("bip_14_", x)))
colnames(scz_data_14) <- sapply(colnames(scz_data_14),
                                function(x) ifelse(x == "SNP", "SNP",
                                                   paste0("scz_14_", x)))

# Set the SNP column to be the keys for joining them together:
setkey(bip_data_14, SNP)
setkey(scz_data_14, SNP)

# Now perform an inner join on SNPs:
bip_scz_data_14 <- bip_data_14[scz_data_14, nomatch = 0]
# Remove the initial datasets:
rm(bip_data_14, scz_data_14)

# Next load the BIP and SCZ 2018 datasets. Since these are massive, will load
# them each individually and in the process only use SNPs from the 2014 data,
# then modify the column names and join. 
bip_data_18 <- fread("gunzip -c data/bip_schz_data/bip_scz_18/BDvsCONT.sumstats.gz") %>%
  .[SNP %in% bip_scz_data_14$SNP,]
colnames(bip_data_18) <- sapply(colnames(bip_data_18),
                                function(x) ifelse(x == "SNP", "SNP",
                                                   paste0("bip_18_", x)))
setkey(bip_data_18, SNP)
bip_scz_data_14_18 <- bip_scz_data_14[bip_data_18, nomatch = 0]
rm(bip_scz_data_14, bip_data_18)
scz_data_18 <- fread("gunzip -c data/bip_schz_data/bip_scz_18/SCZvsCONT.sumstats.gz") %>%
  .[SNP %in% bip_scz_data_14_18$SNP,]
colnames(scz_data_18) <- sapply(colnames(scz_data_18),
                                function(x) ifelse(x == "SNP", "SNP",
                                                   paste0("scz_18_", x)))
setkey(scz_data_18, SNP)
bip_scz_data_14_18 <- bip_scz_data_14_18[scz_data_18, nomatch = 0]
rm(scz_data_18)

# Next thing to do is create the columns denoting the direction of the effects
# ensuring that the same alleles are considered for each. This means checking
# for each of the four datasets what A1 and A2 were, then flipping the sign
# of the effect based on the flipping of these alleles. We treat the bip_14 
# data / alleles as the reference for which all flipping will be examined relative to.

# Now create six columns of two types: 
# (1) indicator if test flipped is relative to bip_14_A1 and bip_14_A2, 
# (2) adjusting the flipped odds ratio:
bip_scz_data_14_18[,':='(scz_14_allele_status = ifelse(scz_14_A1 == bip_14_A2 & 
                                                         scz_14_A2 == bip_14_A1,
                                                       "flip",
                                                       ifelse(scz_14_A1 == bip_14_A1 & 
                                                                scz_14_A2 == bip_14_A2,
                                                              "match", "other")),
                         scz_18_allele_status = ifelse(scz_18_A1 == bip_14_A2 & 
                                                         scz_18_A2 == bip_14_A1,
                                                       "flip",
                                                       ifelse(scz_18_A1 == bip_14_A1 & 
                                                                scz_18_A2 == bip_14_A2,
                                                              "match", "other")),
                         bip_18_allele_status = ifelse(bip_18_A1 == bip_14_A2 & 
                                                         bip_18_A2 == bip_14_A1,
                                                       "flip",
                                                       ifelse(bip_18_A1 == bip_14_A1 & 
                                                                bip_18_A2 == bip_14_A2,
                                                              "match", "other")))]
# View the frequencies:
table(bip_scz_data_14_18$scz_14_allele_status)
#   flip  match 
# 596341 546561 
table(bip_scz_data_14_18$scz_18_allele_status)
#   match   other 
# 1142137     765
table(bip_scz_data_14_18$bip_18_allele_status)
#   match   other 
# 1142137     765
table(bip_scz_data_14_18$bip_18_allele_status[which(bip_scz_data_14_18$scz_18_allele_status == "other")])
# other 
#   765 

# Okay so the 'other' SNPs for the 2018 results match - meaning the same 765 
# have different alleles tested for both BIP and SCZ in 2018. We can remove these
# mismatched SNPs. Then only the 2014 SCZ results need to have their odds ratio
# "flipped" since the remainder of the 2018 SCZ and BIP results match the 2014
# BIP alleles exactly.

# Remove the 765 'other' SNPs:
bip_scz_data_14_18 <- bip_scz_data_14_18[scz_18_allele_status == "match",]

# Next create the new column match_scz_14_OR to match the 2014 allele directions
bip_scz_data_14_18[, match_scz_14_OR := ifelse(scz_14_allele_status == "flip",
                                               exp(-log(scz_14_OR)),
                                               scz_14_OR)]

# Now calculate the effect size as the log of the odds ratio:
bip_scz_data_14_18[, ':='(beta_bip_14 = log(bip_14_OR),
                          match_beta_scz_14 = log(match_scz_14_OR),
                          beta_scz_18 = log(scz_18_OR),
                          beta_bip_18 = log(bip_18_OR))]

# Next create columns with the z-stats using the log(OR) / SE
bip_scz_data_14_18[, ':='(z_bip_14 = beta_bip_14 / bip_14_SE,
                          z_scz_14 = match_beta_scz_14 / scz_14_SE,
                          z_scz_18 = beta_scz_18 / scz_18_SE,
                          z_bip_18 = beta_bip_18 / bip_18_SE)]

# Create absolute value columns just to have:
bip_scz_data_14_18[, ':='(abs_z_bip_14 = abs(z_bip_14),
                          abs_z_scz_14 = abs(z_scz_14),
                          abs_z_scz_18 = abs(z_scz_18),
                          abs_z_bip_18 = abs(z_bip_18))]

# Now create columns separating the 2014 and 2018 studies from one another
# using the fact that both the 2014 and 2018 paper use the fixed-effect 
# inverse-variance based meta-analysis. Will denote these columns with the
# new_ and follow the same code used for separating the new studies for the BMI
# GIANT data, since the same meta-analysis procedure was used.

# First calculate the variances for all four:
bip_scz_data_14_18[, ':='(var_bip_14 = bip_14_SE^2,
                          var_scz_14 = scz_14_SE^2,
                          var_bip_18 = bip_18_SE^2,
                          var_scz_18 = scz_18_SE^2)]

# Next separate the effects and standard errors from the studies independent in 2018:
bip_scz_data_14_18[, ':='(new_beta_bip_18 = (beta_bip_18 * var_bip_14 - 
                                               beta_bip_14 * var_bip_18) /
                            (var_bip_14 - var_bip_18),
                          new_se_bip_18 = sqrt((var_bip_14 * var_bip_18) /
                                                 (var_bip_14 - var_bip_18)),
                          new_beta_scz_18 = (beta_scz_18 * var_scz_14 -
                                               match_beta_scz_14 * var_scz_18) /
                            (var_scz_14 - var_scz_18),
                          new_se_scz_18 = sqrt((var_scz_14 * var_scz_18) /
                                                 (var_scz_14 - var_scz_18)))]
# Using these new effect sizes and betas, get the new z-scores and p-values:
bip_scz_data_14_18[, ':='(new_z_bip_18 = new_beta_bip_18 / new_se_bip_18,
                          new_z_scz_18 = new_beta_scz_18 / new_se_scz_18)]
bip_scz_data_14_18[, ':='(new_p_bip_18 = 2 * pnorm(-abs(new_z_bip_18)),
                          new_p_scz_18 = 2 * pnorm(-abs(new_z_scz_18)))]

# Will only keep the SNPs for which the studies can be separated from standard
# errors decreasing in 2018 compared to 2014, just for SCZ
bip_scz_data_14_18 <- bip_scz_data_14_18[!is.na(new_p_scz_18),]

# Save this dataset prior to filtering on any conditions:
# readr::write_csv(bip_scz_data_14_18, "data/bip_schz_data/bip_scz_data_14_18_snps.csv")
