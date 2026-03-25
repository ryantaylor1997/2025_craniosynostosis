# Define function for predicting GAM on new data (in homage to get_predictions)
obtain_predictions <- function(model, newdata,
                               rm_effects = c("~"), # By default match nothing; assume no variable includes tilde
                               ci_crit = 1.96){

  x <- newdata

  # Get design matrix of basis functions
  nd_mat <- predict(model, newdata = x, type = "lpmatrix")

  # Find the columns for effects we want to remove
  rm_cols <- map(rm_effects,
                 ~grep(., colnames(nd_mat), fixed = T, value = T)) %>%
    unlist() %>% unique()

  # For these columns, set (linear predictor) columns to 0
  nd_mat[, rm_cols] <- 0

  # Get coefficients of the model
  m_coef <- coef(model)

  # Calculate predictions
  x$fit <- (nd_mat %*% m_coef) %>% as.vector()

  # Calculate variances (without calculating all n x n covariances)
  x$se <- (nd_mat %*% vcov(model) * nd_mat) %>% rowSums() %>% sqrt()

  # Add confidence interval
  x %<>% mutate(
    ci_ll = fit - ci_crit * se,
    ci_ul = fit + ci_crit * se
  )

  return(x)
}
