################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Test new Cholesky decomposition function / MVNorm
################################################################################

rm(list = ls()); gc()

# Load packages -----------------------------------------------------------

pacman::p_load(
  Matrix, MASS, reshape2,
  tidyverse, magrittr, janitor, here
)

# Load my function --------------------------------------------------------

source(here("source", "fn_mv_normal_matrix.R"))

# Set scalars -------------------------------------------------------------

# Number of simulations
n_sim <- 1000

# Number of draws per simulation
n_draws <- 1e4

# Set variances
var_diag <- c(7, 50, 300)

# Take square roots
sd_diag <- as.matrix(sqrt(var_diag), ncol = 1)

# Set correlations
var_corr_mx <- matrix(
  c(1, 0.2, 0.5,
    0.2, 1, 0.8,
    0.5, 0.8, 1),
  nrow = 3, ncol = 3, byrow = T
  )

# Multiply together to get var-covar matrix
vc_mx <- tcrossprod(sd_diag) * var_corr_mx

# Convert to Matrix form
vc_MX <- Matrix(vc_mx)

# Set mean vector
mean_vec <- rep(0, 3)

# Generate multivariate normal --------------------------------------------

# Create arrays to save results of each method
array_mass <- array_custom <- array(0, dim = c(n_draws, 3, n_sim))

# Generate random draws using 2 methods
set.seed(978)

for(i in 1:n_sim){
  array_mass[,,i] <- mvrnorm(n_draws, mu = mean_vec, Sigma = vc_mx)
  array_custom[,,i] <- as.matrix(
    rmvnorm_Matrix(n_draws, mean = mean_vec, sigma = vc_MX)
  )
}

# Summarize results -------------------------------------------------------

# Take empirical var-covar for each simulation
vc_emp_mass <- simplify2array(apply(array_mass, 3, cov, simplify = F))

vc_emp_custom <- simplify2array(apply(array_custom, 3, cov, simplify = F))

# Take mean of each entry of var-covar matrix
vc_emp_avg_mass <- apply(vc_emp_mass, c(1,2), mean)
vc_emp_avg_custom <- apply(vc_emp_custom, c(1,2), mean)

# Compare to truth
bias_avg_mass <- vc_emp_avg_mass - vc_mx
bias_avg_custom <- vc_emp_avg_custom - vc_mx
