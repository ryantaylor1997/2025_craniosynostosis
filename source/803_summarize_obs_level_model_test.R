################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Summarize results of test Bayesian GAM on one basis coefficient
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load model fit on 1 coefficient ("coeff_bayes1", "coeff_y", "def_pars_coeff")
load(file = here("results", "obs_level_test_coefficient_model.rda"))

# Load all coefficients fit with GCV, with demographics ("cranio_models_beta")
load(file = here("results", "individual_soap_film_fits.rda"))

# Load observation-level design matrix ("obs_smooth_list")
load(file = here::here("data", "cleaned", "obs_level_smooth.rda"))

# Load color mapping for fusion types ("fusion_color_dict")
load(file = here("data", "intermediate", "fusion_color_mapping.rda"))

# Visualize some soap film basis coefficients by age ----------------------

# Clean coefficient data to be plotted
coeff_toplot <- cranio_models_beta %>%
  pivot_longer(cols = matches("beta"),
               names_to = "beta_id",
               values_to = "beta_val") %>%
  mutate(beta_num = as.numeric(str_remove_all(beta_id, "beta_")))

# Define selection of basis functions to plot
plot_knots <- c(1:5, 15:19, 38:42, 65:69)

# Plot selection of coefficients by age
coeff_plots <- ggplot(coeff_toplot %>%
                        filter(beta_num %in% plot_knots),
                      aes(x = age, y = beta_val, color = fusion_type)) +
  geom_point(alpha = 0.2) +
  geom_smooth(alpha = 0.5) +
  geom_rug(inherit.aes = F, aes(x = age), color = "gray50") +
  scale_color_discrete(palette = fusion_color_dict$color) +
  facet_grid(beta_id ~ fusion_type, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Age (Days)", y = "Soap Basis Coefficient")

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
coeff_newdata <- cranio_models_beta %>%
  mutate(min_age= min(age),
         max_age = max(age)) %>%
  distinct(fusion_type, min_age, max_age) %>%
  rowwise() %>%
  mutate(age = list(seq(min_age, max_age, length.out = 50))) %>%
  ungroup() %>%
  mutate(sex = list(c(0, 1))) %>%
  unnest(age) %>% unnest(sex) %>%
  mutate(sex = factor(sex, levels = levels(cranio_models_beta$sex)))

# Convert this to design matrix
coeff_new_X_list <- map(obs_smooth_list$smooth,
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
             data = enframe(obs_smooth_list[[1]]$xp),
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

