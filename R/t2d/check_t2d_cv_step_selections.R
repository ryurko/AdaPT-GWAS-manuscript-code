# PURPOSE: View the AdaPT CV selection at each step for T2D results

# AUTHOR: Ron Yurko

# Load the final adjusted model results:
t2d_adj_adapt <- readRDS("data/t2d/adapt_cv_results/t2d_adj_s05_2cv.rds")

# Vector of the AdaPT arguments:
args_search <- c("nrounds100md2", "nrounds150md2",
                 "nrounds100md3", "nrounds150md3")

# First step:
args_search[which.max(t2d_adj_adapt$model_holdout_ll_sums[[1]])]
# [1] "nrounds100md2"
# Second step:
args_search[which.max(t2d_adj_adapt$model_holdout_ll_sums[[2]])]
# [1] "nrounds150md3"
