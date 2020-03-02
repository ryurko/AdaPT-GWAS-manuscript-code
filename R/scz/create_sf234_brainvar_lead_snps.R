# PURPOSE: Perform post-hoc analysis of discoveries, identifying the lead SNPs
#          among the set of discoveries using LD pruning. Generate supplementary
#          Figures 2 - 4.

# AUTHOR: Ron Yurko

# ------------------------------------------------------------------------------

# Load the SCZ BrainVar data:
bip_scz_brainvar_data <- read_csv("data/bip_schz_data/bip_scz_data_14_18_brainvar_eqtls_hcp_wgcna.csv")

# Load each of the separate model results:
scz_with_bd_z_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_only_s05_2cv.rds")
scz_with_bd_z_eqtl_slopes <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_s05_2cv.rds")
scz_with_bd_z_eqtl_slopes_wgcna <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_bd_z_eqtl_slopes_wgcna_s05_2cv.rds")
scz_with_wgcna_only <- 
  readRDS("data/bip_schz_data/brainvar_results/cv_tune_results/scz_with_wgcna_only_s05_2cv.rds")
scz_intercept_only <- 
  readRDS("data/bip_schz_data/brainvar_results/scz_intercept_only_s05.rds")
scz_with_no_int <-
  readRDS("data/bip_schz_data/brainvar_results/scz_all_vars_no_int_s05.rds")

# Set-up the reference data:

# Get the file paths for the reference data: -----
g1000_eur_fam_path <- "data/reference_data/g1000_eur.fam"
g1000_eur_bim_path <- "data/reference_data/g1000_eur.bim"
g1000_eur_bed_path <- "data/reference_data/g1000_eur.bed"

# Read in the PLINK data -----
ref_snps_plink <- read.plink(g1000_eur_bed_path, 
                             g1000_eur_bim_path,
                             g1000_eur_fam_path)
ref_genotypes <- ref_snps_plink$genotypes
rm(ref_snps_plink)

# Only include the eSNPs
ref_genotypes <- ref_genotypes[, which(colnames(ref_genotypes) %in%
                                         bip_scz_brainvar_data$SNP)]
print(ref_genotypes)
# A SnpMatrix with  503 rows and  24883 columns
# Row names:  HG00096 ... NA20832 
# Col names:  rs1048488 ... rs1129880 

# Only 24883 of 25076 in 1000 Genomes reference data 

# Create the numeric version of this matrix with alleles encoded as 0, 1, or 2:
numeric_ref_genotypes <- as(ref_genotypes, "numeric")

# Define the flip matrix function:
flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

numeric_ref_genotypes <- flip_matrix(numeric_ref_genotypes)
# Save this dataset to the problem_genes folder:
# write.table(numeric_ref_genotypes, 
#             file = "data/reference_data/brainvar_esnps_genotype.csv", 
#             row.names = FALSE, col.names = TRUE)
rm(ref_genotypes)
# numeric_ref_genotypes <- read.table("data/reference_data/brainvar_esnps_genotype.csv",
#                                     header = TRUE)

# Next for each set of discoveries and  - determine the lead SNPs and see which
# of the lead SNPs are in the set of discoveries:
scz_comp_qval_list <- list("Intercept-only" = scz_intercept_only$qvals,
                           "BD z-stats" = scz_with_bd_z_only$qvals,
                           "BD z-stats + eQTL slopes" = scz_with_bd_z_eqtl_slopes$qvals,
                           "BD z-stats + eQTL slopes + WGCNA (w/ interactions)" = scz_with_bd_z_eqtl_slopes_wgcna$qvals,
                           "BD z-stats + eQTL slopes + WGCNA (w/o interactions)" = scz_with_no_int$qvals,
                           "WGCNA" = scz_with_wgcna_only$qvals)
scz_comp_disc_list <- list("Intercept-only" = which(scz_intercept_only$qvals <= 0.05),
                           "BD z-stats" = which(scz_with_bd_z_only$qvals <= 0.05),
                           "BD z-stats + eQTL slopes" = which(scz_with_bd_z_eqtl_slopes$qvals <= 0.05),
                           "BD z-stats + eQTL slopes + WGCNA (w/ interactions)" = which(scz_with_bd_z_eqtl_slopes_wgcna$qvals <= 0.05),
                           "BD z-stats + eQTL slopes + WGCNA (w/o interactions)" = which(scz_with_no_int$qvals <= 0.05),
                           "WGCNA" = which(scz_with_wgcna_only$qvals <= .05))

# Return as a list:
ordered_brainvar_esnps <- bip_scz_brainvar_data %>%
  filter(SNP %in% colnames(numeric_ref_genotypes)) %>%
  arrange(bip_18_CHR, bip_18_BP)
# Get the correct column order for the genotype matrix using this sorted
# chromosome and positional information
geno_snp_order_i <- sapply(colnames(numeric_ref_genotypes),
                           function(snp) which(ordered_brainvar_esnps$SNP == snp))
geno_snp_order_i <- unlist(geno_snp_order_i)
ordered_fake_ref <- snp_fake(nrow(numeric_ref_genotypes), 
                             ncol(numeric_ref_genotypes))
ordered_fake_ref$genotypes[] <- numeric_ref_genotypes[, geno_snp_order_i]
ordered_brainvar_esnps$snp_sort_index <- 1:nrow(ordered_brainvar_esnps)

scz_ind_loci_disc_list <- map(1:length(scz_comp_disc_list),
                              function(adapt_i) {
                                result_name <- names(scz_comp_disc_list)[adapt_i]
                                temp_data <- bip_scz_brainvar_data %>%
                                  mutate(adapt_qvals = scz_comp_qval_list[[adapt_i]],
                                         adapt_qvals = ifelse(is.infinite(adapt_qvals), 1, adapt_qvals),
                                         adapt_disc = as.numeric(adapt_qvals <= 0.05)) %>%
                                  filter(SNP %in% colnames(numeric_ref_genotypes)) %>%
                                  arrange(bip_18_CHR, bip_18_BP)
                                adapt_ld_clumps <- snp_clumping(ordered_fake_ref$genotypes,
                                                                infos.chr = temp_data$bip_18_CHR,  
                                                                infos.pos = temp_data$bip_18_BP,
                                                                S = -log10(temp_data$adapt_qvals),
                                                                thr.r2 = .1, size = 500, ncores = 4)
                                # Print the number of clumps found:
                                print(paste0(result_name, " : ", length(adapt_ld_clumps)))
                                # Return the overlap between the ld clumps and discoveries:
                                temp_data <- temp_data %>%
                                  mutate(adapt_lead_snp = ifelse(1:n() %in% adapt_ld_clumps, 1, 0))
                                which(temp_data$adapt_lead_snp == 1 &
                                        temp_data$adapt_disc == 1) %>%
                                  return
                              })

names(scz_ind_loci_disc_list) <- names(scz_comp_disc_list)

# Now using these list of lead SNPs in each set of discoveries, generate a replica
# of Figure 2 (excluding the enrichment plot):

# Convert this into the 0-1 membership matrix using the UpSetR package and then
# get it in the data format for ggupset instead:
tidy_scz_comp_ind_loci_disc <- scz_ind_loci_disc_list %>%
  fromList() %>%
  as_tibble() %>%
  mutate(snp = 1:n()) %>%
  gather(disc_set, member, -snp) %>%
  filter(member == 1) %>%
  select(-member) %>%
  group_by(snp) %>%
  summarize(disc_sets = list(disc_set))

# Next create the ggupset plot
scz_ind_loci_disc_upset <- tidy_scz_comp_ind_loci_disc %>%
  ggplot(aes(x = disc_sets)) +
  geom_bar(fill = "darkblue") +
  labs(y = "Intersection size") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_upset(expand = c(0, 0)) +
  theme_combmatrix(combmatrix.panel.point.color.fill = "darkblue",
                   combmatrix.label.make_space = FALSE,
                   combmatrix.panel.line.size = .5,
                   combmatrix.label.text = element_blank(),
                   combmatrix.panel.margin = unit(c(0, 0), "pt"),
                   combmatrix.label.height = unit(95, "pt"))


# Create a chart that just counts how many discoveries each type has:
scz_count_ind_loci_disc_chart <- scz_ind_loci_disc_list %>%
  fromList() %>%
  as_tibble() %>%
  mutate(snp = 1:n()) %>%
  gather(disc_set, member, -snp) %>%
  filter(member == 1) %>%
  group_by(disc_set) %>%
  summarize(n_disc = n()) %>%
  mutate(disc_set = fct_reorder(disc_set, n_disc)) %>%
  ggplot(aes(x = disc_set, y = n_disc)) +
  geom_bar(stat = "identity", fill = "darkblue",
           width = .75) +
  geom_text(aes(x = disc_set, y = n_disc / 2,
                label = n_disc),
            color = "white", size = 4) +
  scale_x_discrete(position = "top") +
  scale_y_reverse(position = "right") +
  labs(y = "Number of independent loci in discovery set") +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 0.5, size = 10, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())


# Next create the updated versions of the Manhattan plots with the AdaPT q-values:

# First create the chromosome axis data for both Manhattan plots:
ordered_brainvar_manhattan_axis_df <- ordered_brainvar_esnps %>% 
  group_by(bip_18_CHR) %>% 
  summarize(center = (max(bip_18_BP_cum) + min(bip_18_BP_cum)) / 2)

# Add columns with the intercept-only and all covariates q-values, along with 
# indicators for discoveries and lead SNPs. Add in a random sampling indicator
# to use for helping with saving PDFs of the images:
update_ordered_brainvar_esnps <- bip_scz_brainvar_data %>%
  mutate(int_only_qvals = scz_comp_qval_list$`Intercept-only`,
         boosted_qvals = scz_comp_qval_list$`BD z-stats + eQTL slopes + WGCNA (w/ interactions)`,
         int_only_qvals = ifelse(is.infinite(int_only_qvals), 1, int_only_qvals),
         boosted_qvals = ifelse(is.infinite(boosted_qvals), 1, boosted_qvals),
         int_only_disc = as.numeric(int_only_qvals <= 0.05),
         boosted_disc = as.numeric(boosted_qvals <= 0.05))  %>%
  filter(SNP %in% colnames(numeric_ref_genotypes)) %>%
  arrange(bip_18_CHR, bip_18_BP) %>%
  mutate(snp_sort_index = 1:n(),
         lead_snp_int_only_disc = ifelse(snp_sort_index %in% scz_ind_loci_disc_list$`Intercept-only`,
                                         1, 0),
         lead_snp_boosted_disc = ifelse(snp_sort_index %in% scz_ind_loci_disc_list$`BD z-stats + eQTL slopes + WGCNA (w/ interactions)`,
                                        1, 0),
         int_only_sample = ifelse(int_only_qvals <= 0.4,
                                  TRUE,
                                  rbernoulli(n(), p = 0.5)),
         boosted_sample = ifelse(boosted_qvals <= 0.4,
                                 TRUE,
                                 rbernoulli(n(), p = 0.5)))

# What about a version of the q-values instead:
int_only_qval_manhattan <- update_ordered_brainvar_esnps %>%
  filter(int_only_sample) %>%
  ggplot(aes(x = bip_18_BP_cum, y = -log10(int_only_qvals))) +
  # Show all points
  geom_point(aes(color = as.factor(bip_18_CHR)), alpha = 0.5, size = 0.5) +
  geom_point(data = filter(update_ordered_brainvar_esnps, int_only_disc == 1),
             aes(x = bip_18_BP_cum, y = -log10(int_only_qvals),
                 shape = as.factor(lead_snp_int_only_disc), 
                 size = as.factor(lead_snp_int_only_disc)),
             color = "darkorange", alpha = 0.75) +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "darkred") +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_shape_manual(values = c(1, 4)) +
  scale_size_manual(values = c(1, 3)) +
  # custom x-axis to handle the chromosome labels:
  scale_x_continuous(label = ordered_brainvar_manhattan_axis_df$bip_18_CHR, 
                     breaks = ordered_brainvar_manhattan_axis_df$center) +
  labs(x = "Chromosome", y = TeX('-log$_{10}$(SCZ q-value)'),
       title = TeX('Manhattan q-value plot of AdaPT intercept-only discoveries for $\\alpha = 0.05$')) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 6))

# What about a version of the q-values instead:
boosted_qval_manhattan <- update_ordered_brainvar_esnps %>%
  filter(boosted_sample) %>%
  ggplot(aes(x = bip_18_BP_cum, y = -log10(boosted_qvals))) +
  # Show all points
  geom_point(aes(color = as.factor(bip_18_CHR)), alpha = 0.5, size = 0.5) +
  geom_point(data = filter(update_ordered_brainvar_esnps, boosted_disc == 1),
             aes(x = bip_18_BP_cum, y = -log10(boosted_qvals),
                 shape = as.factor(lead_snp_boosted_disc), 
                 size = as.factor(lead_snp_boosted_disc)),
             color = "darkorange", alpha = 0.75) +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "darkred") +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  scale_shape_manual(values = c(1, 4)) +
  scale_size_manual(values = c(1, 3)) +
  # custom x-axis to handle the chromosome labels:
  scale_x_continuous(label = ordered_brainvar_manhattan_axis_df$bip_18_CHR, 
                     breaks = ordered_brainvar_manhattan_axis_df$center) +
  labs(x = "Chromosome", y = TeX('-log$_{10}$(SCZ q-value)'),
       title = TeX('Manhattan q-value plot of AdaPT gradient boosted discoveries for $\\alpha = 0.05$')) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 6))

# Arrange the Manhattan plots together
scz_qval_manhattan_plots <-
  plot_grid(int_only_qval_manhattan,
            boosted_qval_manhattan,
            ncol = 1, labels = c("A", "B"), 
            label_fontface = "plain", align = "hv", vjust = 1.2)

# Finally combine everything together:
scz_brainvar_ind_loci_disc_plot <- 
  plot_grid(plot_grid(scz_qval_manhattan_plots,
                      scz_count_ind_loci_disc_chart,
                      ncol = 1, rel_widths = c(1, 1),
                      rel_heights = c(4, 1.75), labels = c("", "C"),
                      label_fontface = "plain"),
            scz_ind_loci_disc_upset, ncol = 2, rel_widths = c(1, 1),
            rel_heights = c(1, 1), labels = c("", "D"),
            label_fontface = "plain")

# Save the figure:
save_plot("figures/sf2_scz_ind_loci_results.pdf",
          scz_brainvar_ind_loci_disc_plot, ncol = 3, nrow = 2,
          base_width = 5, base_height = 3)
save_plot("/nonpdf_figures/sf2_ind_loci_results.jpg",
          scz_brainvar_ind_loci_disc_plot, ncol = 3, nrow = 2,
          base_width = 5, base_height = 3)

# ------------------------------------------------------------------------------

# Next plot the number of independent loci discoveries by chromosome:
brainvar_ld_loci_chr_disc_plot <- update_ordered_brainvar_esnps %>%
  filter(lead_snp_int_only_disc == 1) %>%
  group_by(bip_18_CHR) %>%
  summarize(n_disc = sum(int_only_disc)) %>%
  mutate(type = "Intercept-only") %>%
  bind_rows({
    data.frame("bip_18_CHR" = 3,
               "n_disc" = 0,
               "type" = "Intercept-only")
  }) %>%
  bind_rows({
    update_ordered_brainvar_esnps %>%
      filter(lead_snp_boosted_disc == 1) %>%
      group_by(bip_18_CHR) %>%
      summarize(n_disc = sum(boosted_disc)) %>%
      mutate(type = "BD z-stats + eQTL slopes + WGCNA (w/ interactions)")
  }) %>%
  ggplot(aes(x = as.factor(bip_18_CHR),
             y = n_disc, fill = type)) +
  geom_bar(color = "black", stat = "identity", position = "dodge",
           alpha = 0.5) +
  theme_bw() +
  #ggthemes::scale_fill_colorblind(drop = FALSE) +
  scale_fill_manual(values = c("darkorange", "darkblue")) +
  scale_x_discrete(breaks = c(1:22)) +
  labs(x = "Chromosome", y = "Number of independent loci in discovery set",
       fill = "Type",
       title = "Comparison of the number of independent loci in discovery set by chromosome") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

save_plot("figures/sf3_scz_ind_loci_chr_results.pdf",
          brainvar_ld_loci_chr_disc_plot, ncol = 1, nrow = 1,
          base_width = 6, base_height = 4)
save_plot("nonpdf_figures/sf3_ind_loci_chr_results.jpg",
          brainvar_ld_loci_chr_disc_plot, ncol = 1, nrow = 1,
          base_width = 6, base_height = 4)

# ------------------------------------------------------------------------------
# Finally generate a figure of the results with the pruning based on the actual
# p-values instead of the AdaPT q-values:

scz_pval_ind_loci_disc_list <- map(1:length(scz_comp_disc_list),
                                   function(adapt_i) {
                                     #result_name <- names(scz_comp_disc_list)[adapt_i]
                                     temp_data <- bip_scz_brainvar_data %>%
                                       mutate(adapt_qvals = scz_comp_qval_list[[adapt_i]],
                                              adapt_qvals = ifelse(is.infinite(adapt_qvals), 1, adapt_qvals),
                                              adapt_disc = as.numeric(adapt_qvals <= 0.05)) %>%
                                       filter(SNP %in% colnames(numeric_ref_genotypes)) %>%
                                       arrange(bip_18_CHR, bip_18_BP)
                                     pval_ld_clumps <- snp_clumping(ordered_fake_ref$genotypes,
                                                                    infos.chr = temp_data$bip_18_CHR,  
                                                                    infos.pos = temp_data$bip_18_BP,
                                                                    S = -log10(temp_data$scz_14_P),
                                                                    thr.r2 = .1, size = 500, ncores = 4)
                                     
                                     # Return the overlap between the ld clumps and discoveries:
                                     temp_data <- temp_data %>%
                                       mutate(adapt_lead_snp = ifelse(1:n() %in% pval_ld_clumps, 1, 0))
                                     which(temp_data$adapt_lead_snp == 1 &
                                             temp_data$adapt_disc == 1) %>%
                                       return
                                   })
names(scz_pval_ind_loci_disc_list) <- names(scz_comp_disc_list)

# Create a chart that just counts how many discoveries each type has:
scz_count_pval_ind_loci_disc_chart <- scz_pval_ind_loci_disc_list %>%
  fromList() %>%
  as_tibble() %>%
  mutate(snp = 1:n()) %>%
  gather(disc_set, member, -snp) %>%
  filter(member == 1) %>%
  group_by(disc_set) %>%
  summarize(n_disc = n()) %>%
  mutate(disc_set = fct_reorder(disc_set, n_disc)) %>%
  ggplot(aes(x = disc_set, y = n_disc)) +
  geom_bar(stat = "identity", fill = "darkblue",
           width = .75) +
  geom_text(aes(x = disc_set, y = n_disc / 2,
                label = n_disc),
            color = "white", size = 4) +
  scale_x_discrete(position = "top") +
  scale_y_reverse(position = "right") +
  labs(y = "Number of independent loci in discovery set") +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 0.5, size = 10, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())

# Now save:
save_plot("figures/sf4_scz_ld_pval_disc_count.pdf",
          scz_count_pval_ind_loci_disc_chart, ncol = 1, nrow = 1,
          base_width = 8, base_height = 4)
save_plot("nonpdf_figures/sf4_scz_ld_pval_disc_count.jpg",
          scz_count_pval_ind_loci_disc_chart, ncol = 1, nrow = 1,
          base_width = 8, base_height = 4)


