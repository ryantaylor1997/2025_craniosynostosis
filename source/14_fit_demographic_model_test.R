
# Load GCV-based coefficients ---------------------------------------------

## Load all coefficients previously fit, along with demographics
gcv_coeffs <- read_csv(here::here("analysis", "intermediate", "individual_soap_coeffs.csv"))

cranio_coeff_df <- cranio_sub %>%
  left_join(gcv_coeffs %>% select(fname, matches("beta")))

# Visualize some soap film basis coefficients by age ----------------------

# Clean coefficient data to be plotted
coeff_toplot <- cranio_coeff_df %>%
  pivot_longer(cols = matches("beta"),
               names_to = "beta_id",
               values_to = "beta_val") %>%
  mutate(beta_num = as.numeric(str_remove_all(beta_id, "beta_")))

# Plot selection of coefficients by age
coeff_plots <- ggplot(coeff_toplot %>%
                        filter(beta_num %in% so_plot_knots),
                      aes(x = age, y = beta_val, color = fusion_type)) +
  geom_point(alpha = 0.2) +
  geom_smooth(alpha = 0.5) +
  geom_rug(inherit.aes = F, aes(x = age), color = "gray50") +
  scale_color_discrete(palette = fusion_color_dict$color) +
  facet_grid(beta_id ~ fusion_type, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Age (Days)", y = "Soap Basis Coefficient")


# Set up coefficient smoother ---------------------------------------------

### Create cubic regression basis splines for age

# Construct smooth without constraint or factor levels
cranio_coeff_smooth <- s(age, bs = "ts", k = cranio_knots_age)

## Create design and penalty matrices using functions we wrote
cranio_coeff_ingredients <- construct_reference_smooth(
  sm = cranio_coeff_smooth, dat = cranio_coeff_df,
  by_var = "fusion_type", param_formula = ~sex + fusion_type
)

# Combine suture fusion penalties so we estimate the same lambda
cranio_coeff_S_combo <- list(cranio_coeff_ingredients$S[[1]],
                             Reduce("+", cranio_coeff_ingredients$S[-1]))

# Adjust penalty matrices for underflow
scale_coeff_penalty <- norm(cranio_soap$S[[1]], type = "O") /
  norm(cranio_coeff_ingredients$smooth[[1]]$S[[1]], type = "O")

cranio_coeff_S_combo[-1] <- map(cranio_coeff_S_combo[-1], ~ .x * scale_coeff_penalty)

# Save intermediate steps to share
save(cranio_coeff_ingredients,
     file = here::here("analysis", "intermediate", "demographics_model_matrices.rda"))

# Run Gibbs sampler for coefficient models --------------------------------


# Set 1 covariate as outcome
# - (plotted above, should be lower for sagittal and higher for metopic)
coeff_y <- cranio_coeff_df$beta_38

# Set default parameters
def_pars_coeff <- c(
  "beta" = 0.0,
  "sigma_sq" = var(coeff_y) / length(coeff_y),
  "lambda" = 1e-3,
  "a" = 1.0,
  "b" = 1.0,
  "c" = 1.0,
  "d" = 1.0
)

# Set iterations
bayes_coeff_iters <- 5e3

# Run function
coeff_bayes1 <- lm_penalized_gibbs(
  y = coeff_y, X = cranio_coeff_ingredients$X,
  S = cranio_coeff_S_combo, has_unpenalized = TRUE,
  beta0 = def_pars_coeff["beta"],
  sigmasq0 = def_pars_coeff["sigma_sq"],
  lambda0 = def_pars_coeff["lambda"],
  a_pri = def_pars_coeff["a"],
  b_pri = def_pars_coeff["b"],
  c_pri = def_pars_coeff["c"],
  d_pri = def_pars_coeff["d"],
  iters = bayes_coeff_iters,
  burn_pct = 0.1
)

# Visualize MCMC diagnostics ----------------------------------------------

# Check timing
mcmc_coeff_time <- coeff_bayes1$timing %>%
  mutate(pct = as.numeric(secs) / last(secs))

# Make accessible all parameters together

coeff_bayes1_tr <- imap(
  coeff_bayes1 %>%
    list_modify("timing" = zap(),
                "start-stop" = zap()),
  function(x, xn){
    x <- as_tibble(x)
    if(ncol(x) == 1){
      names(x) <- xn
    } else {
      names(x) <- paste0(xn, ".", 1:ncol(x))
    }
    x %<>%
      mutate(iter = 1:n(), .before = everything()) %>%
      pivot_longer(-iter,
                   names_pattern = "(burn_)?(.*)",
                   names_to = c("draw", "name"))

    x
  }
) %>%
  unname() %>%
  list_rbind() %>%
  mutate(draw_cat = if_else(draw == "burn_",
                            "Burn-In", "Post-Burn-In"),
         param_num = as.numeric(str_extract(name, "[0-9]+"))) %>%
  mutate(param_num = coalesce(param_num, 1))

# Add default parameters to this table
trace_coeff_default <- coeff_bayes1_tr %>%
  filter(draw_cat == "Burn-In") %>%
  distinct(draw_cat, name, param_num) %>%
  mutate(
    iter = 0,
    value = case_when(
      str_detect(name, "sigma_sq") ~ def_pars_coeff["sigma_sq"],
      str_detect(name, "lambda") ~ def_pars_coeff["lambda"],
      str_detect(name, "beta") ~ def_pars_coeff["beta"],
      T ~ NA
    ))

coeff_bayes1_tr %<>% bind_rows(trace_coeff_default)

# Trace plot of sigma
sigma_trace_coeff <- ggplot(coeff_bayes1_tr %>%
                              filter(str_detect(name, "sigma_sq"))) +
  geom_line(aes(x = iter, y = value)) +
  facet_grid(rows = name ~ draw_cat, scales = "free") +
  labs(title = "Sigma Squared Trace")

# Trace plot of lambdas
lambda_trace_coeff <- ggplot(coeff_bayes1_tr %>%
                               filter(str_detect(name, "lambda"))) +
  geom_line(aes(x = iter, y = value), alpha = 0.5) +
  facet_grid(param_num ~ draw_cat, scales = "free") +
  labs(title = "Lambda Trace")

# Trace plot of betas
beta_trace_coeff <- ggplot(coeff_bayes1_tr %>%
                             filter(str_detect(name, "beta"))) +
  geom_line(aes(x = iter, y = value), alpha = 0.5) +
  facet_grid(param_num ~ draw_cat, scales = "free") +
  labs(title = "Beta Trace")

# Plot demographic test model predictions ---------------------------------

### Illustrate coefficient functions from these splines

## Create new data

# Create data shell
coeff_newdata <- cranio_coeff_df %>%
  mutate(min_age= min(age),
         max_age = max(age)) %>%
  distinct(fusion_type, min_age, max_age) %>%
  rowwise() %>%
  mutate(age = list(seq(min_age, max_age, length.out = 50))) %>%
  ungroup() %>%
  mutate(sex = list(c(0, 1))) %>%
  unnest(age) %>% unnest(sex) %>%
  mutate(sex = factor(sex, levels = levels(cranio_coeff_df$sex)))

# Convert this to design matrix
coeff_new_X_list <- map(cranio_coeff_ingredients$smooth,
                       ~PredictMat(.x, data = coeff_newdata))

# Combine to create Normative and add unpenalized terms
coeff_new_design <- Reduce(
  "cbind", c(
    list(
      model.matrix(~sex + fusion_type, data = coeff_newdata),
      Reduce("+", coeff_new_X_list)),
    coeff_new_X_list[-1]
  ))

## Calculate predictions
coeff_bayes1_preds <- tcrossprod(coeff_new_design, coeff_bayes1$beta)

coeff_bayes1_pred_int <- t(apply(coeff_bayes1_preds, 1,
                                 FUN = function(x) matrix(
                                   c(mean(x),
                                     quantile(x, 0.05),
                                     quantile(x, 0.95)),
                                   nrow = 1
                                 )))

colnames(coeff_bayes1_pred_int) <- c("est", "ll", "ul")

# Add predictions to data
coeff_newdata %<>% bind_cols(as.data.frame(coeff_bayes1_pred_int))

## Plot truth and predictions from Bayesian example
coeff_bayes1_plot <- ggplot(coeff_newdata %>% filter(sex == 0),
                            aes(x = age, y = est, ymin = ll, ymax = ul)) +
  geom_line(aes(color = fusion_type)) +
  geom_point(inherit.aes = F,
             data = enframe(cranio_coeff_smooth$xp),
             aes(x = value, y = 0),
             fill = "yellow", shape = 23) +
  geom_ribbon(aes(fill = fusion_type), alpha = 0.1) +
  facet_wrap(~fusion_type, scales = "free_y") +
  scale_color_discrete(palette = fusion_color_dict$color) +
  scale_fill_discrete(palette = fusion_color_dict$color) +
  labs(x = "Age (Days)", y = "Est. Basis Coefficient",
       color = "Fusion Type",
       fill = "Fusion Type",
       title = "Bayesian Hierarchical Coefficient Functions") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Add y limit so that patterns are more visible
coeff_bayes1_plot_lim <- coeff_bayes1_plot +
  coord_cartesian(ylim = max(coeff_y)*c(-0.5, 1.25)) +
labs(title = "Bayesian Hierarchical Coefficient Functions (Axes Fixed)")
