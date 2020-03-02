# PURPOSE: View the AdaPT CV selection at each step for BMI results

# AUTHOR: Ron Yurko

# Load the three AdaPT results
bmi_adj15 <- readRDS("data/bmi/adapt_cv_results/bmi_adj15_s05_2cv.rds")
bmi_adj15_whr_only <- readRDS("data/bmi/adapt_cv_results/bmi_adj15_whr_only_s05_2cv.rds")
bmi_adj15_whr_slopes_only <- readRDS("data/bmi/adapt_cv_results/bmi_adj15_whr_slopes_only_s05_2cv.rds")


bmi_adj15_args <- c("nrounds100md2", "nrounds150md2",
                    "nrounds100md3", "nrounds150md3")
bmi_adj15_whr_only_args <- c("nrounds100md1", "nrounds150md1",
                             "nrounds100md6", "nrounds150md6")
bmi_adj15_whr_slopes_only_args <- c("nrounds100md3", "nrounds150md3",
                                    "nrounds100md6", "nrounds150md6")

# BMI results with all variables:
# First step:
bmi_adj15_args[which.max(bmi_adj15$model_holdout_ll_sums[[1]])]
# [1] "nrounds100md2"
# Second step:
bmi_adj15_args[which.max(bmi_adj15$model_holdout_ll_sums[[2]])]
# [1] "nrounds150md3"

# BMI results using only WHR z-statistics
# First step:
bmi_adj15_whr_only_args[which.max(bmi_adj15_whr_only$model_holdout_ll_sums[[1]])]
# [1] "nrounds150md1"
# Second step:
bmi_adj15_whr_only_args[which.max(bmi_adj15_whr_only$model_holdout_ll_sums[[2]])]
# [1] "nrounds150md1"

# BMI results using WHR z-statistics + eQTL slopes
# First step:
bmi_adj15_whr_slopes_only_args[which.max(bmi_adj15_whr_slopes_only$model_holdout_ll_sums[[1]])]
# [1] "nrounds100md3"
# Second step:
bmi_adj15_whr_slopes_only_args[which.max(bmi_adj15_whr_slopes_only$model_holdout_ll_sums[[2]])]
# [1] "nrounds150md3"


