# Organization of R code

This folder contains the scripts for initializing datasets, generating results
and figures in the manuscript. The code is separated into folders by the respective
type of analysis conducted:

- [scz](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/R/scz) - contains
all code relevant to the main analysis presented in the paper on application of AdaPT
to SCZ GWAS results,
- [bmi](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/R/bmi) - contains
all code relevant to the supplementary analysis on application of AdaPT to BMI GWAS results,
- [t2d](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/R/t2d) - contains
all code relevant to the supplementary analysis on application of AdaPT to T2D GWAS results,
- [gtex_wgcna](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/R/gtex_wgcna) - contains all code necessary for initializing the GTEx eQTL datasets and the WGCNA results using GTEx RNA-Seq data,
- [si_simulations](https://github.com/ryurko/AdaPT-GWAS-manuscript-code/blob/master/R/si_simulations) - contains all code necessary for supplementary simulations to determine the number of CV steps, effect of dependent tests, and overfitting.

All code is presented under the assumption of working in a R Project version of this repository and that the necessary initial datasets are available in the `data` folder. See the README in the `data` folder for more information regarding the required initial datasets.
