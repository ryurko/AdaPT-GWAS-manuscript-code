# PURPOSE: Initialize the GTEx eQTL dataset to use

# AUTHOR: Ron Yurko

# Access necessary packages
library(data.table)
library(R.utils)

# Get the list of files
eqtl_files <- list.files("data/gtex_eqtl/GTEx_Analysis_v7_eQTL",
                         pattern = "egenes", full.names = T)


# Now loop through the files to load and stack the datasets together:
eqtl_data <- lapply(eqtl_files,
                    function(x) {
                      region <- gsub("\\.v7\\.egenes\\.txt\\.gz","",
                                     gsub("data/gtex_eqtl/GTEx_Analysis_v7_eQTL/",
                                          "", x))
                      eqtl <- read.table(gzfile(x), header = TRUE,
                                         sep = "\t")
                      eqtl <- eqtl[
                        c("variant_id", "gene_id","gene_name","gene_chr","gene_start","gene_end",
                          "rs_id_dbSNP147_GRCh37p13","chr","pos","ref","alt","maf",
                          "qval", "slope", "slope_se", "pval_beta", "log2_aFC", "log2_aFC_lower",
                          "log2_aFC_upper")]
                      data.frame(region = region, ensg_gene = eqtl$gene_id,
                                 variant_id = eqtl$variant_id,
                                 gene_name = eqtl$gene_name,
                                 gene_chr = eqtl$gene_chr,
                                 gene_start = eqtl$gene_start,
                                 gene_end = eqtl$gene_end,
                                 eqtl_name = eqtl$rs_id_dbSNP147_GRCh37p13,
                                 eqtl_chr = eqtl$chr,
                                 eqtl_pos = eqtl$pos,
                                 eqtl_ref = eqtl$ref,
                                 eqtl_alt = eqtl$alt,
                                 eqtl_maf = eqtl$maf,
                                 eqtl_qval = eqtl$qval,
                                 eqtl_pval_beta = eqtl$pval_beta,
                                 eqtl_slope = eqtl$slope,
                                 eqtl_slope_se = eqtl$slope_se,
                                 eqtl_log2_afc = eqtl$log2_aFC,
                                 eqtl_log2_afc_lower = eqtl$log2_aFC_lower,
                                 eqtl_log2_afc_upper = eqtl$log2_aFC_upper)
                    })

eqtl_data <- as.data.frame(rbindlist(eqtl_data))
egenes <- eqtl_data[eqtl_data$eqtl_qval < .05, ]
# Save this dataset:
readr::write_csv(egenes, "data/gtex_eqtl/eqtl_data.csv")
