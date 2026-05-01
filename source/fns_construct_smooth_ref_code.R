### Write functions to set up design matrix
##  - Given a constructed smooth (with constraint applied), factor variable, and new data:
## .. and generate full design matrix for new data,
## .. assuming "reference coding" with reference spline and difference functions

## 1) Given smooth, subtract column means and drop last column
apply_sweep_and_drop <- function(sm, drop_col = NULL, newdata = NULL){

  # Extract design matrix
  sm_X0 <- sm$X

  # Calculate column means
  sm_means <- colMeans(sm_X0)

  # Identify number of columns
  sm_k <- length(sm_means)

  # Determine whether column to drop is provided
  if(is.null(drop_col)) drop_col <- sm_k

  # Determine whether to sweep and drop from original or new matrix
  if(is.null(newdata)){
    nd_X0 <- sm_X0
  } else {
    nd_X0 <- PredictMat(sm, data = newdata)
  }

  # Subtract these from design matrix AND remove last column
  X_ref <- sweep(nd_X0[, -drop_col], 2, STATS = as.array(sm_means[-drop_col]))

  # Locate full penalty matrix
  sm_S0 <- sm$S[[1]]

  # Filter to subset corresponding to remaining design matrix
  S_ref <- sm_S0[-drop_col, -drop_col]

  # Export
  out <- list("X" = X_ref,
              "S" = S_ref,
              "means" = sm_means,
              "col_dropped" = drop_col)

  out
}

## 2) Given constructed smooth and a factor, make full design and penalty
##    - assume "reference" coding for splines
make_reference_coded_splines <- function(sm, drop_col = NULL, newdata = NULL,
                                         by_var){

  # Apply constraint to smooth object
  constrained_list <- apply_sweep_and_drop(sm = sm,
                                           drop_col = drop_col,
                                           newdata = newdata)

  # Get data with factor variable: either from the smooth or a new data frame
  if(is.null(newdata)){
    nd <- get(attr(sm, "data"))
  } else {
    nd <- newdata
  }

  # Define vector of factor assignments in data
  by_vec <- nd[[by_var]]

  # Define list of effect values for which we need separate splines
  by_effects <- levels(by_vec)[-1]

  # Generate one design matrix for each effect
  X_by <- map(by_effects,
              ~ as.numeric(by_vec == .x) * constrained_list$X)

  # Merge with original spline basis
  X_splines <- Reduce(cbind, x = X_by, init = constrained_list$X)

  out <- list("X" = X_splines)

  # Construct penalty matrix for all splines together
  if(is.null(newdata)){
    # Define indicator matrix list
    # Each matrix is 0 except for a 1 on the diagonal
    # This 1 represents the block where the penalty matrix S is populated
    spline_indicator_list <- map(1:length(levels(by_vec)),
                                 .f = function(j, z = rep(0, length(levels(by_vec)))){
                                   z[j] <- 1
                                   diag(z, nrow = length(z))
                                 })

    # For each matrix in the list, take kronecker product with S from the smooth
    S_splines <- map(spline_indicator_list,
                     ~kronecker(.x, constrained_list$S))

    out[["S"]] <- S_splines
  }
  out
}

## 3) Add parametric covariates
setup_reference_model <- function(sm, drop_col = NULL, by_var, newdata = NULL,
                                  param_formula = NULL){

  # Make spline matrix
  ref_splines <- make_reference_coded_splines(sm = sm, drop_col = drop_col,
                                              by_var = by_var, newdata = newdata)

  # If there is no parametric formula, add unpenalized intercept
  if(is.null(param_formula)){
    design_param <- matrix(1, nrow = nrow(ref_splines$X))

    # Include 1 x 1 matrix for unpenalized intercept
    S_unpenal <- matrix(1, 1, 1)
  } else{
    # Otherwise, create matrix columns from formula

    # Get data with parametric variables: if no new data, get it from smooth
    if(is.null(newdata)){
      nd <- get(attr(sm, "data"))
    } else {
      nd <- newdata
    }

    # Create parametric matrix from this new data
    design_param <- model.matrix(param_formula, nd)
    colnames(design_param) <- NULL

    # Create identity matrix for these unpenalized parameters
    S_unpenal <- diag(nrow = ncol(design_param))
  }

  # Merge model matrix with spline matrix for full design matrix
  X_design <- cbind(design_param, ref_splines$X)

  out <- list("X" = X_design)

  if(is.null(newdata)){
    ## Create full penalty / precision matrix

    # Add 0's to unpenalized matrix to match size of penalties
    S_unpenal_full <- as.matrix(bdiag(S_unpenal,
                                      diag(0, nrow = nrow(ref_splines$S[[1]]))))

    # Add 0's for penalized matrices at beginning of other matrices
    S_full_size <- map(ref_splines$S,
                       ~as.matrix(
                         bdiag(
                           diag(0, nrow = nrow(S_unpenal)),
                           .x)))

    # Append unpenalized chunk to beginning of penalty matrix list
    S_full <- c(list(S_unpenal_full), S_full_size)

    out[["S"]] <- S_full
  }

  out
}

### OR: Use smoothCon to construct design and penalty, then combine
construct_reference_smooth <- function(sm, dat, by_var,
                                       param_formula = NULL){

  # Add "by" variable to pre-specified smooth
  sm$by <- by_var

  # Construct smooth with smoothCon
  smooth_list <- smoothCon(
    sm,
    dat,
    knots = NULL,
    absorb.cons = T, # Apply identifiability constraints
    scale.penalty = FALSE, # Don't multiply penalty scale
    sparse.cons = -1, # If no existing constraint, apply sweep and drop
    apply.by = FALSE # Also return a design matrix not separated by smooths
  )

  out <- list("smooth" = smooth_list)

  # Extract list of design matrices
  X_list <- map(smooth_list, ~.$X)

  # Combine to form total design matrix
  ref_X <- Reduce("cbind", X_list[-1], init = smooth_list[[1]]$X0)

  # Construct penalty matrix for all splines together
  S_list <- map(smooth_list, ~.$S[[1]])

  # Construct penalty matrix for all splines together
  # Define indicator matrix list
  # Each matrix is 0 except for a 1 on the diagonal
  # This 1 represents the block where the penalty matrix S is populated
  spline_indicator_list <- map(1:length(smooth_list),
                               .f = function(j,
                                             z = rep(0, length(smooth_list))){
                                 z[j] <- 1
                                 diag(z, nrow = length(z))
                               })

  # For each matrix in the list, take kronecker product with S from the smooth
  ref_S <- map2(spline_indicator_list, S_list,
                ~kronecker(.x, .y))

  # If there is no parametric formula, add unpenalized intercept
  if(is.null(param_formula)){
    design_param <- matrix(1, nrow = nrow(ref_X))

    # Include 1 x 1 matrix for unpenalized intercept
    S_unpenal <- matrix(1, 1, 1)
  } else{
    # Otherwise, create matrix columns from formula

    # Create parametric matrix from the data
    design_param <- model.matrix(param_formula, dat)
    colnames(design_param) <- NULL

    # Create identity matrix for these unpenalized parameters
    S_unpenal <- diag(nrow = ncol(design_param))
  }

  # Merge model matrix with spline matrix for full design matrix
  X_design <- cbind(design_param, ref_X)

  out[["X"]] <- X_design

  ## Create full penalty / precision matrix

  # Add 0's to unpenalized matrix to match size of penalties
  S_unpenal_full <- as.matrix(bdiag(S_unpenal,
                                    diag(0, nrow = nrow(ref_S[[1]]))))

  # Add 0's for penalized matrices at beginning of other matrices
  S_full_size <- map(ref_S,
                     ~as.matrix(
                       bdiag(
                         diag(0, nrow = nrow(S_unpenal)),
                         .x)))

  # Append unpenalized chunk to beginning of penalty matrix list
  S_full <- c(list(S_unpenal_full), S_full_size)

  out[["S"]] <- S_full

  out
}
