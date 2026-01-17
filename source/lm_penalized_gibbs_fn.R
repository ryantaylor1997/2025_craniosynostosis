### Write Gibbs sampler for Bayesian linear regression

# Notes:
# - Assumes that each penalty matrix is the full K x K but is block diagonal
# - Assumes at least one set of penalized terms
# - Assumes one distinct smoothing parameter lambda for each S penalty matrix
# - If there are unpenalized terms, assumes one matrix is unpenalized
# - Assumes that the one unpenalized matrix is the first in the list

lm_penalized_gibbs <- function(
    y, X, # Data
    S, has_unpenalized = FALSE, # Penalty matrices
    sigmasq0, lambda0, beta0, # Default parameters
    a_pri, b_pri, c_pri, d_pri, # Prior hyperparameters
    iters = 1e4, burn_pct = 0.1 # MCMC parameters
){

  # Set up time tracking
  time_check <- Sys.time()
  time_pts <- c("start" = time_check)
  timer <- c("setup" = 0,
             "lambda_s" = 0,
             "beta_update_c" = 0,
             "x_inverse" = 0,
             "beta_moments" = 0,
             "beta_draw" = 0,
             "sigma2_update" = 0,
             "lambdas_update" = 0)

  ### Reformat penalty matrix
  if(is.matrix(S)){ S <- list(S) }

  ### Calculate dimensions

  ## Design matrix
  K <- ncol(X)
  N <- nrow(X)

  stopifnot(length(y) == N)

  ## Iterations of burn-in
  iters_burn <- ceiling(iters * burn_pct)

  ## Penalty dimensions

  # number of penalties
  S_blocks <- length(S)
  num_penals <- S_blocks - as.numeric(has_unpenalized)

  # non-zero columns in each penalty matrix
  S_block_cols <- map(S, ~which(apply(., 2, sum) != 0))

  # number of parameters for each penalty
  S_block_Ks <- map_dbl(S_block_cols, length)

  # First number of parameters might not impact any lambda; if so, remove it
  if(has_unpenalized){ S_block_Ks <- S_block_Ks[-1] }

  ### Initialize hyperparameters

  # Assign initial values
  sigmasq <- sigmasq0
  lambda <- rep_len(lambda0, num_penals)
  beta <- rep_len(beta0, K)

  # Initialize beta moments before update
  beta_expect <- beta
  beta_var <- diag(K)
  xtx_plus_pen <- diag(K)

  # Compute hyperparameters that don't update (gamma / Igamma shapes)
  a_post <- N/2 + K/2 + a_pri
  c_post <- S_block_Ks/2 + c_pri

  ### Create objects to save output

  # Posterior samples
  beta_out <- array(0, dim = c(iters - iters_burn, K))
  sigmasq_out <- array(0, dim = c(iters - iters_burn, 1))
  lambda_out <- array(0, dim = c(iters - iters_burn, num_penals))

  # Burn-in period
  beta_burn <- array(0, dim = c(iters_burn, K))
  sigmasq_burn <- array(0, dim = c(iters_burn, 1))
  lambda_burn <- array(0, dim = c(iters_burn, num_penals))

  timer["setup"] <- timer["setup"] +
    difftime(Sys.time(), time_check, units = "secs")
  time_check <- Sys.time()

  ### Run MCMC loop

  for(iter in 1:iters){

    # Output on checkpoints
    if(iter %% iters_burn == 0) cat("at turn", iter, "\n")

    ## (1) Update beta (coefficients)

    # Full penalty matrix including lambdas
    lambda_S_list <- if(has_unpenalized){
      # If there are unpenalized terms, do not multiply 1st matrix by any lambda
      c(
        S[1],
        map2(lambda, S[-1],
             function(l_const, s_mx){
               l_const * s_mx
             })
      )
    } else {
      # If no unpenalized terms, scale by one lambda for each penalty matrix
      map2(lambda, S,
           function(l_const, s_mx){
             l_const * s_mx
           })
    }

    # Sum together these individual (lambda * S) block diagonal matrices
    lambda_S <- Reduce("+", lambda_S_list)

    timer["lambda_s"] <- timer["lambda_s"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    # Call to Cpp function for beta update
    penalized_beta_inplace(
      y = y, X = X, lambda_S = lambda_S, sigmasq = sigmasq,
      beta_expect,
      beta_var,
      xtx_plus_pen
    )

    timer["beta_update_c"] <- timer["beta_update_c"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    # Generate new beta sample
    beta <- t(rmvnorm(n = 1, beta_expect, beta_var))

    timer["beta_draw"] <- timer["beta_draw"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ## (2) Update sigma^2 (residual variance)

    # 2nd parameter of inverse gamma
    b_post <-
      b_pri - crossprod(y, X) %*% beta +
      0.5 * (sum(y^2) + crossprod(beta, xtx_plus_pen) %*% beta)

    # Generate new sigma squared sample
    precision <- rgamma(1, a_post, b_post)
    sigmasq <- 1 / precision

    timer["sigma2_update"] <- timer["sigma2_update"] +
      difftime(Sys.time(), time_check, units = "secs")
    time_check <- Sys.time()

    ## (3) Update lambda (beta function smoothness factor)

    # Calculate posterior gamma shape
    d_post <- map(
      # Exclude any unpenalized terms from posterior
      if(has_unpenalized){ S[-1] } else { S },
      function(S_j){
        (1 / (2 * sigmasq)) * crossprod(beta, S_j) %*% beta + d_pri
      })

    # Update lambdas
    lambda <- map2_dbl(c_post, d_post,
                       function(c_new, d_new){
                         rgamma(1, c_new, d_new)
                       })

    timer["lambdas_update"] <- timer["lambdas_update"] +
      difftime(Sys.time(), time_check, units = "mins")
    time_check <- Sys.time()

    ## Save results
    iter_post <- iter - iters_burn

    if(iter_post > 0){
      beta_out[iter_post,] <- beta
      sigmasq_out[iter_post,] <- sigmasq
      lambda_out[iter_post,] <- lambda
    } else{
      beta_burn[iter,] <- beta
      sigmasq_burn[iter,] <- sigmasq
      lambda_burn[iter,] <- lambda
    }
  }

  time_pts <- c(time_pts, "end" = Sys.time())
  timer["total"] <- difftime(time_pts["end"], time_pts["start"], units = "secs")
  timer <- enframe(timer, name = "step", value = "secs")

  return(list("beta" = beta_out,
              "sigma_sq" = sigmasq_out,
              "lambda" = lambda_out,
              "burn_beta" = beta_burn,
              "burn_sigma_sq" = sigmasq_burn,
              "burn_lambda" = lambda_burn,
              "start-stop" = time_pts,
              "timing" = timer))
}
