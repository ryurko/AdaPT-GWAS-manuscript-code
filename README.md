# AdaPT-GWAS-manuscript-code

This repository contains the code for reproducing the results and figures in 
[_A selective inference approach for FDR control using multi-omics covariates yields insights into disease risk_](https://www.biorxiv.org/content/early/2019/10/16/806471.full.pdf).

The folders are organized in the following manner:

- [R](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/R) - all scripts for initializing datasets, generating results and figures in manuscript,
- [figures](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/figures) - files for final figures displayed in the manuscript (including supplementary materials),
- [nonpdf_figures](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/nonpdf_figures) - files for non-pdf final figures displayed in the manuscript (including supplementary materials),
- [data](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/data) - folder
containing files necessary for the analysis.

Due to their respective terms and conditions, we are unable to publicly post the initial genome-wide association studies (GWAS) datasets for which the analysis was conducted on. However, all datasets considered in the manuscript are publicly available in the following locations:

- __Body mass index (BMI)__ and __Waist-to-hip ratio (WHR)__ GWAS results are available to download from the [GIANT consortium](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_Consortium_data_files). The _2015-only_ BMI studies correspond to the [_GWAS Anthropometric 2015 BMI Summary Statistics_ for EUR ancestry](https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz) (Locke et al., 2015). The _all 2018_ BMI studies for replication
analysis are from the [updated _BMI and Height GIANT and UK BioBank Meta-analysis Summary Statistics_](https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz) (Yengo et al., 2018). The WHR studies used as covariates in the
BMI analysis are from the [combined EUR ancestry results](https://portals.broadinstitute.org/collaboration/giant/images/5/54/GIANT_2015_WHR_COMBINED_EUR.txt.gz) (Shungin et al., 2015).
- __Type 2 diabetes (T2D)__ GWAS results are available to download from the [Diabetes Genetics Replication And Meta-analysis (DIAGRAM) consortium](https://www.diagram-consortium.org/downloads.html), under the _T2D GWAS meta-analysis - Unadjusted for BMI_ (Mahajan et al., 2018).
- __Schizophrenia (SCZ)__ and __bipolar disorder (BD)__ GWAS results are available to download
from the [Psychiatric Genomics Consortium](https://www.med.unc.edu/pgc/download-results/). The _2014-only_ studies correspond to the results from Ruderfer et al. (2014), while the _all 2018_ studies for replication analysis are from Ruderfer et al. (2018).
- [__Genotype-Tissue Expression (GTEx)__ V7 project datasets](https://gtexportal.org/home/datasets) (GTEx Consortium, 2015) were accessed for identifying eQTLs and generating gene co-expression modules. [The Single-Tissue cis-eQTL Data](https://gtexportal.org/home/datasets#filesetFilesDiv55) file corresponding to "eGene and significant variant-gene associations based on permutations" was accessed for identifying GTEx eQTLs. The [RNA-Seq Data](https://gtexportal.org/home/datasets#filesetFilesDiv54) gene TPM counts were used for creating gene co-expression covariates.
- BrainVar [eQTLs](https://www.biorxiv.org/content/biorxiv/early/2019/03/22/585430/DC5/embed/media-5.xlsx?download=true) and [WGCNA](https://www.biorxiv.org/content/biorxiv/early/2019/03/22/585430/DC4/embed/media-4.xlsx?download=true) results are available to download from the [available
pre-print for Werling et al. (2019) on bioRxiv](https://www.biorxiv.org/content/10.1101/585430v1).

## References

- AE Locke, et al., Genetic studies of body mass index yield new insights for obesity biology. _Nature_ __518__, 197 EP – (2015).
- L Yengo, et al., Meta-analysis of genome-wide association studies for height and body mass index in 700000 individuals of european ancestry. _Hum. Mol. Genet._ __27__, 3641 – 3649 (2018).
- D Shungin, et al., New genetic loci link adipose and insulin biology to body fat distribution. _Nature_ __518__ (2015).
- A Mahajan, et al., Fine-mapping type 2 diabetes loci to single-variant resolution using high-density imputation and islet-specific epigenome maps. _Nat. Genet._ __50__, 1505 – 1513 (2018).
- DM Ruderfer, et al., Polygenic dissection of diagnosis and clinical dimensions of bipolar disorder and schizophrenia. _Mol. psychiatry_ __19__, 1017 – 1024 (2014).
- DM Ruderfer, et al., Genomic dissection of bipolar disorder and schizophrenia, including 28 subphenotypes. _Cell_ __173__, 1705 – 1715.e16 (2018).
- GTEx Consortium, The genotype-tissue expression (gtex) pilot analysis: Multitissue gene regulation in humans. Science __348__, 648 – 660 (2015).
- DM Werling, et al., Whole-genome and rna sequencing reveal variation and transcriptomic coordination in the developing human prefrontal cortex. _bioRxiv_ (2019).





