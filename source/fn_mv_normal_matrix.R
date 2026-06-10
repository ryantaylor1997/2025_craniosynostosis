################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Function to draw multivariate random normal using Matrix methods
################################################################################


# New MVNorm random draw generator ----------------------------------------

# Take number of samples, mean vector, variance matrix
# - Based on mvtnorm::rmvnorm cholesky method
# - Innovation: accept Matrix package sparse formats as mean and variance
rmvnorm_Matrix <- function(n,
                           mean = rep(0, ncol(sigma)),
                           sigma = Diagonal(length(mean))){

  ## Take Cholesky decomposition

  if(is.matrix(sigma)){ # First, where sigma is a traditional matrix

    # Cholesky decomp with pivoted rows / columns
    # - Documentation says this is "R s.t. A = R'R"
    sigma_ch <- chol(sigma, pivot = T)

    # Identify order of pivot
    ch_order <- attr(sigma_ch, "pivot")

    # Reorder columns to get real square root matrix
    sigma_sqrt <- sigma_ch[, order(ch_order)]

  } else { # Second, where sigma is a format from Matrix package

    # Take Cholesky decomposition using Matrix package
    # - Documentation says this is "L' s.t. A = LL'"
    sigma_sqrt <- Matrix::chol(sigma, pivot = F)
  }

  ## Convert mean vector to a matrix the size of the random draw


  if(is.null(dim(mean)) | min(dim(mean)) == 1){
    # First, the case when mean is a vector or 1-dim matrix
    mean_mx <- Matrix(rep(mean, n), nrow = n, byrow = TRUE)
  } else if(nrow(mean) == n) {
    # Next, where mean is a matrix and each row is the mean for a separate draw
    mean_mx <- mean
  }

  ## Take random draw of (univariate) standard normal
  rands <- rnorm(n * ncol(sigma))

  # Convert to matrix
  rands <- Matrix(rands, nrow = n, byrow = TRUE)

  ## Transform to multivariate normal draw

  # Multiply by square root of variance matrix
  mvn_draw <- rands %*% sigma_sqrt + mean_mx

  return(mvn_draw)
}
