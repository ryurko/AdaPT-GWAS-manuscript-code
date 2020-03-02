# PURPOSE: Generate block correlation simulations with randomly dispersed
#          alternatives to a subset of the number of blocks so there are far 
#          fewer blocks containing alternatives

# Access necessary packages:
library(tidyverse)
library(future)
# Use the development version of furrr:
# devtools::install_github("DavisVaughan/furrr")
library(furrr)
# Use the modified version of AdaPT:
# devtools::install_github("ryurko/adaptMT")
library(adaptMT)
library(Matrix)
library(MASS)

# First a function that returns a correlation matrix for a given correlation value
# and block size
create_single_rho_block <- function(rho, block_dim) {
  block <- matrix(rho, nrow = block_dim, ncol = block_dim)
  diag(block) <- 1
  return(block)
}

# Next a function that takes a vector of true effect sizes, along with a 
# correlation value and returns the simulated vector of dependent block sizes
create_dep_block <- function(effect_sizes, block_cor) {
  return(mvrnorm(1, 
                 mu = matrix(effect_sizes, 
                             ncol = 1), 
                 Sigma = create_single_rho_block(block_cor,
                                                 length(effect_sizes))))
}

# Function that takes in the number of tests, total number of blocks, and then
# randomly assigns which of the blocks to assign tests to:
sim_random_alt_block_structure <- function(n_tests, n_blocks,
                                           block_rho, 
                                           pi_int, pi_beta1, pi_beta2,
                                           mu_floor, mu_beta1, mu_beta2) {
  # Generate two independent uniform covariates between 0 and 1:
  data <- data.frame(x1 = runif(n_tests, 0, 1),
                     x2 = runif(n_tests, 0, 1))
  
  # Now use the logit relationship with the given information to generate
  # the probability of being non-null
  inv_logit <- function(x) {exp(x) / (1 + exp(x))}
  data$pi <- inv_logit(data$x1 * pi_beta1 + data$x2 * pi_beta2 + pi_int)
  
  # Now the test types (later on will explore correlated binary variables...)
  data$test_type <- ifelse(sapply(data$pi,
                                  function(p) rbinom(1, 1, p)) == 1,
                           "h1", "h0")
  # Add a test index column for keeping tracking of the tests to then assign
  # which blocks each belongs to:
  data$test_index <- 1:nrow(data)
  
  # Now depending on the test type generate the true means:
  data$effect_size_center <- pmax(mu_floor, mu_beta1 * data$x1 + mu_beta2 * data$x2)
  data$effect_size_center <- ifelse(data$test_type == "h0", 0,
                                    data$effect_size_center)
  
  # How many true alternatives are there?
  n_alt <- length(which(data$test_type == "h1"))
  # Thus nulls:
  n_null <- n_tests - n_alt
  
  # Create a vector of the alt test indices:
  alt_test_indices <- which(data$test_type == "h1")
  # For the nulls:
  null_test_indices <- which(data$test_type == "h0")
  
  # What's the block size? Assuming n_tests is a multiple of n_blocks (can always
  # add a check for that):
  block_size <- n_tests / n_blocks
  
  # What is the minimum number of blocks required to fill all of the alternative
  # results - round up:
  min_alt_blocks <- ceiling(n_alt / block_size)
  # Now pick the number of alt blocks to be the average between this and
  # the max number of blocks:
  n_alt_blocks <- ceiling(mean(c(min_alt_blocks, n_blocks)))
  
  # Now make a vector of alt_block_labels just going up to this number of blocks,
  # but with the full block sizes:
  possible_alt_block_labels <- unlist(lapply(1:n_alt_blocks, function(x) rep(x, block_size)))
  
  # Randomly assign WITHOUT replacement, these possible block labels to the 
  # alternative tests. Do this by randomly accessing the indices:
  sample_indices <- sample.int(length(possible_alt_block_labels), size = n_alt,
                               replace = FALSE) 
  alt_block_labels <- possible_alt_block_labels[sample_indices]
  remaining_alt_block_labels <- possible_alt_block_labels[-sample_indices]
  
  # Make a vector of remaining null labels joining the ones that weren't sampled
  # for the alt blocks and join with the remaining null blocks:
  remaining_null_block_labels <- c(remaining_alt_block_labels,
                                   unlist(lapply((n_alt_blocks + 1):n_blocks, 
                                                 function(x) rep(x, block_size))))
  
  # Now sample which null test is assigned to which block:
  null_block_labels <- sample(remaining_null_block_labels,
                              n_null, replace = FALSE)
  
  # Now make a column of the block labels:
  data$block_labels <- rep(1, nrow(data))
  data$block_labels[alt_test_indices] <- alt_block_labels
  data$block_labels[null_test_indices] <- null_block_labels
  
  # Now for each block generate a vector of observed effect sizes:
  data$obs_effect_size <- rep(NA, nrow(data))
  for (block_i in 1:n_blocks) {
    block_centers <- data$effect_size_center[which(data$block_labels == block_i)]
    data$obs_effect_size[which(data$block_labels == block_i)] <- create_dep_block(block_centers, block_rho)
  }
  
  # Finally the two-sided p-values:
  data$pvals <- 2*pnorm(-abs(data$obs_effect_size))
  # Return a list with data frame of these results, along with information about
  # the settings:
  return(list("data" = data,
              "n_alt" = n_alt,
              "n_null" = n_null,
              "block_size" = block_size))
}


# Create the AdaPT functions:
# This function takes in a simulated dataset using the sim_dep_pvals_data function
# above, and returns the fdp and power for the intercept-only AdaPT model results
# at the target level alpha vector:
get_adapt_intercept_results <- function(sim_data, alphas) {
  # Generate the intercept only model results:
  adapt_int_only <- adapt_glm(sim_data,
                              sim_data$pvals,
                              pi_formulas = "1",
                              mu_formulas = "1",
                              verbose = list(print = FALSE, 
                                             fit = FALSE,
                                             ms = FALSE),
                              nfits = 10)
  
  # Vector of true nulls:
  true_nulls <- which(sim_data$test_type == "h0")
  # Vector of true alternatives:
  true_alts <- which(sim_data$test_type == "h1")
  
  # Which blocks have alternatives in them:
  alt_block_labels <- unique(sim_data$block_labels[true_alts])
  
  # Which tests are in these blocks?
  tests_in_alt_blocks <- which(sim_data$block_labels %in% alt_block_labels)
  
  # Generate a data frame of results for the vector of alpha values:
  do.call(rbind,
          lapply(alphas, 
                 function(alpha) {
                   # Access the discoveries for alpha:
                   adapt_disc <- which(adapt_int_only$qvals <= alpha)
                   
                   # Discovery blocks:
                   adapt_disc_blocks <- unique(sim_data$block_labels[adapt_disc])
                   
                   # Return the fdp and power:
                   return(data.frame("method" = "adapt_int_only",
                                     "alpha" = alpha,
                                     "n_true_alt" = length(true_alts),
                                     "n_disc" = length(adapt_disc),
                                     "fdp" = length(which(adapt_disc %in% true_nulls)) / length(adapt_disc),
                                     "power" = length(which(adapt_disc %in% true_alts)) / length(true_alts),
                                     "n_true_alt_blocks" = length(alt_block_labels),
                                     "n_disc_blocks" = length(adapt_disc_blocks),
                                     "fdp_blocks" = length(which(!(adapt_disc_blocks %in% alt_block_labels))) / length(adapt_disc_blocks),
                                     "power_blocks" = length(which((adapt_disc_blocks %in% alt_block_labels))) / length(alt_block_labels),
                                     "fdp_test_blocks" = length(which(!(adapt_disc %in% tests_in_alt_blocks))) / length(adapt_disc),
                                     "power_test_blocks" = length(which(adapt_disc %in% tests_in_alt_blocks)) / length(tests_in_alt_blocks)))   
                 }))
  
}

# This function returns BH results for a dataset constructed
# using the function sim_dep_pvals_data function above for target level alphas:
get_bh_results <- function(sim_data, alphas) {
  # Vector of true nulls:
  true_nulls <- which(sim_data$test_type == "h0")
  # Vector of true alternatives:
  true_alts <- which(sim_data$test_type == "h1")
  
  # Which blocks have alternatives in them:
  alt_block_labels <- unique(sim_data$block_labels[true_alts])
  
  # Which tests are in these blocks?
  tests_in_alt_blocks <- which(sim_data$block_labels %in% alt_block_labels)
  
  # Get the BH and BY adjusted p-values:
  bh_pvals <- p.adjust(sim_data$pvals, method = "BH")
  
  # Generate a data frame of results for the vector of alpha values:
  do.call(rbind,
          lapply(alphas, 
                 function(alpha) {
                   # Get the BH discoveries for alpha:
                   bh_disc <- which(bh_pvals <= alpha)
                   
                   # Discovery blocks:
                   bh_disc_blocks <- unique(sim_data$block_labels[bh_disc])
                   
                   # Return the fdp and power:
                   return(data.frame("method" = "bh",
                                     "alpha" = alpha,
                                     "n_true_alt" = length(true_alts), 
                                     "n_disc" = length(bh_disc),
                                     "fdp" = length(which(bh_disc %in% true_nulls)) / length(bh_disc),
                                     "power" = length(which(bh_disc %in% true_alts)) / length(true_alts),
                                     "n_true_alt_blocks" = length(alt_block_labels),
                                     "n_disc_blocks" = length(bh_disc_blocks),
                                     "fdp_blocks" = length(which(!(bh_disc_blocks %in% alt_block_labels))) / length(bh_disc_blocks),
                                     "power_blocks" = length(which((bh_disc_blocks %in% alt_block_labels))) / length(alt_block_labels),
                                     "fdp_test_blocks" = length(which(!(bh_disc %in% tests_in_alt_blocks))) / length(bh_disc),
                                     "power_test_blocks" = length(which(bh_disc %in% tests_in_alt_blocks)) / length(tests_in_alt_blocks)))    
                 }))
}

# This function takes in a simulated dataset using the sim_bivariate_dep_pvals_data
# function above, and returns the fdp and power for the AdaPT model using XGBoost
# with the x1 and x2 covariates:
get_adapt_boosting_results <- function(sim_data, alphas) {
  # Generate the regular covariate results:
  adapt_regular <- adapt_xgboost(as.matrix(sim_data[,c("x1", "x2")]),
                                 sim_data$pvals,
                                 verbose = list(print = FALSE, 
                                                fit = FALSE,
                                                ms = FALSE),
                                 piargs = list("nrounds" = 50,
                                               "max_depth" = 1,
                                               "nthread" = 1,
                                               "verbose" = 0),
                                 muargs = list("nrounds" = 50,
                                               "max_depth" = 1,
                                               "nthread" = 1,
                                               "verbose" = 0),
                                 nfits = 10)
  
  # Vector of true nulls:
  true_nulls <- which(sim_data$test_type == "h0")
  # Vector of true alternatives:
  true_alts <- which(sim_data$test_type == "h1")
  
  # Which blocks have alternatives in them:
  alt_block_labels <- unique(sim_data$block_labels[true_alts])
  
  # Which tests are in these blocks?
  tests_in_alt_blocks <- which(sim_data$block_labels %in% alt_block_labels)
  
  # Generate a data frame of results for the vector of alpha values:
  do.call(rbind,
          lapply(alphas, 
                 function(alpha) {
                   # Access the discoveries for alpha:
                   adapt_reg_disc <- which(adapt_regular$qvals <= alpha)
                   
                   # Discovery blocks:
                   adapt_reg_disc_blocks <- unique(sim_data$block_labels[adapt_reg_disc])
                   
                   # Return the fdp and power:
                   return(data.frame("method" = "adapt_regular",
                                     "alpha" = alpha,
                                     "n_true_alt" = length(true_alts),
                                     "n_disc" = length(adapt_reg_disc), 
                                     "fdp" = length(which(adapt_reg_disc %in% true_nulls)) / length(adapt_reg_disc),
                                     "power" = length(which(adapt_reg_disc %in% true_alts)) / length(true_alts),
                                     "n_true_alt_blocks" = length(alt_block_labels),
                                     "n_disc_blocks" = length(adapt_reg_disc_blocks), 
                                     "fdp_blocks" = length(which(!(adapt_reg_disc_blocks %in% alt_block_labels))) / length(adapt_reg_disc_blocks),
                                     "power_blocks" = length(which((adapt_reg_disc_blocks %in% alt_block_labels))) / length(alt_block_labels),
                                     "fdp_test_blocks" = length(which(!(adapt_reg_disc %in% tests_in_alt_blocks))) / length(adapt_reg_disc),
                                     "power_test_blocks" = length(which(adapt_reg_disc %in% tests_in_alt_blocks)) / length(tests_in_alt_blocks)))   
                 }))
  
}

# Generate simulation settings, and for ease will set the effect sizes for
# both covariates to be the exact same:

sim_n_tests <- c(10000) 
sim_block_rho <- seq(0, 1, by = 0.25)
sim_n_blocks <- 500 
sim_pi_int <- c(-3)
sim_pi_beta <- c(0, 1, 2, 3)
sim_mu_floor <- c(.5, 1, 1.5)
sim_mu_beta <- c(0, 0.5, .75, 1)
sim_index <- c(1:100)

cov_sim_settings_df <- as.data.frame(expand.grid(sim_block_rho, sim_n_tests,
                                                 sim_n_blocks,
                                                 sim_pi_int, 
                                                 sim_pi_beta, sim_pi_beta,
                                                 sim_mu_floor,
                                                 sim_mu_beta, sim_mu_beta,
                                                 sim_index))
colnames(cov_sim_settings_df) <- c("block_rho", "n_tests",
                                   "n_blocks",
                                   "pi_int", "pi_beta1", "pi_beta2",
                                   "mu_floor", "mu_beta1", "mu_beta2",
                                   "sim_index")

# Generate the results ------
plan(tweak(multiprocess, workers = 20)) # This code was distributed on a cluster
results <- furrr::future_map_dfr(1:nrow(cov_sim_settings_df),
                                 function(x) {
                                   
                                   # Grab the simulation settings:
                                   sim_settings_row <- cov_sim_settings_df[x,]
                                   
                                   # Generate the simulated data:
                                   sim_data_list <- 
                                     sim_random_alt_block_structure(sim_settings_row$n_tests,
                                                                    sim_settings_row$n_blocks,
                                                                    sim_settings_row$block_rho, 
                                                                    sim_settings_row$pi_int,
                                                                    sim_settings_row$pi_beta1,
                                                                    sim_settings_row$pi_beta2,
                                                                    sim_settings_row$mu_floor,
                                                                    sim_settings_row$mu_beta1,
                                                                    sim_settings_row$mu_beta2)
                                   
                                   alphas <- c(.01, .05, .1, .15, .2)
                                   # Now return a dataframe of results for the
                                   # AdaPT intercept, boosting, and BH results
                                   # over target alpha values:
                                   method_results <- rbind(get_adapt_intercept_results(sim_data_list$data,
                                                                                       alphas),
                                                           get_adapt_boosting_results(sim_data_list$data,
                                                                                      alphas),
                                                           get_bh_results(sim_data_list$data, alphas))
                                   
                                   
                                   # Join the columns about the simulation settings and return:
                                   final_results <- cbind(method_results,
                                                          do.call(rbind,
                                                                  lapply(1:nrow(method_results),
                                                                         function(x) {
                                                                           data.frame(cbind(sim_settings_row,
                                                                                            data.frame("n_alt" = sim_data_list$n_alt,
                                                                                                       "n_null" = sim_data_list$n_null,
                                                                                                       "block_size" = sim_data_list$block_size)))
                                                                         }
                                                                  )))
                                   
                                   # Return
                                   return(final_results)
                                   
                                 })



# Save the results
# write_csv(results, 
#           "data/si_simulations/dep_sim_alt_blocks.csv")

