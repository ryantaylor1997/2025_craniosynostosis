################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Compare Joint Model Fits
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load soap film object ("cranio_soap")
load(file = here("data", "cleaned", "soap_object.rda"))

# Load observation-level smooth object ("obs_smooth_list")
load(file = here("data", "cleaned", "obs_level_smooth.rda"))

# Test 1. Check if chol crossprod returns original ------------------------

# Sum up penalty matrices to make matrix to work with
pen_total <- Reduce("+", obs_smooth_list$S)

## Take cross-product of soap film penalty matrix using base R function
chol_base <- chol(pen_total)

# Take cross product
chol_base_sq <- crossprod(chol_base)

# Compare to original matrix
chol_base_diff <- mean(pen_total - chol_base_sq)

## Repeat: base R, allowing columns to be permuted
chol_pivot <- chol(pen_total, pivot = T)

chol_pivot_sq <- crossprod(chol_pivot)

chol_pivot_diff <- mean(pen_total - chol_pivot_sq)

# Also compare to base decomposition
chol_pivot_order <- attr(chol_pivot, "pivot")

chol_pivot_reord <- chol_pivot[, order(chol_pivot_order)]

chol_reord_diff <- mean(chol_base - chol_pivot_reord)

cat("Avg. difference between pivot and no pivot:", chol_reord_diff)

## Try Matrix package

# Convert matrix to Matrix format
Pen_Total <- Matrix(pen_total)

chol_mat <- chol(Pen_Total, pivot = F)

chol_mat_sq <- crossprod(chol_mat)

chol_mat_diff <- mean(Pen_Total - chol_mat_sq)

# Compare to base R
chol_diff_base_mat <- mean(as.matrix(chol_mat) - chol_base)

## Try Matrix package with pivot

chol_matpiv <- chol(Pen_Total, pivot = T)

chol_matpiv_sq <- crossprod(chol_matpiv)

chol_matpiv_diff <- mean(Pen_Total - chol_matpiv_sq)

# Compare to base R
chol_diff_reord_mat <- mean(as.matrix(chol_matpiv) - chol_pivot)


# Try using base R on large matrix ----------------------------------------

## Generate large matrix similar to the one we need to invert

# Take sum of observation-level smooth penalties
obs_pen_total <- Reduce("+", obs_smooth_list$S)

# Take cross-product of observation-level design matrix
obs_design <- crossprod(obs_smooth_list$X)

# Take sum of soap-film smooth penalties
soap_pen_total <- Reduce("+", cranio_soap$S)

# Combine to generate test matrix
large_test_mx <- kronecker(soap_pen_total, obs_design) +
  kronecker(
    diag(x=1, nrow = ncol(soap_pen_total)),
    obs_pen_total
    )

# Benchmark cholesky decompositions
bench_chol <- microbenchmark::microbenchmark(
  "base" = chol(large_test_mx),
  "Matrix" = chol(Matrix(large_test_mx)),
  times = 10
)
