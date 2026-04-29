################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Compare Joint Model Fits
################################################################################

rm(list = ls()); gc()

# Load previous runs

load(file = here::here("analysis", "intermediate",
                       "joint_model_fit_2026-04-21.rda"))

cranio_joint_recent <- cranio_joint

load(file = here::here("analysis", "intermediate",
                       "joint_model_fit 03.14.rda"))

cranio_joint_prev <- cranio_joint

rm(cranio_joint); gc()

# Load soap film object for penalty scales
load(file = here::here("analysis", "intermediate", "soap_object.rda"))

# Load demographic smooth object for scale terms
load(file = here::here("analysis", "intermediate", "demographics_model_matrices.rda"))

# Collect posterior mean values

model_summ_prev <- map(cranio_joint_prev[c("sigma_sq",
                                            "tau_sq",
                                            "lambda_basis",
                                            "lambda_demo")],
                       rowMeans) %>%
  enframe() %>%
  unnest_longer(value) %>%
  group_by(name) %>%
  mutate(num = row_number(), .after = name) %>%
  ungroup() %>%
  mutate(model = "previous")

model_summ_rec <- map(cranio_joint_recent[c("sigma_sq",
                                            "tau_sq",
                                            "lambda_basis",
                                            "lambda_demo")],
                      rowMeans) %>%
  enframe() %>%
  unnest_longer(value) %>%
  group_by(name) %>%
  mutate(num = row_number(), .after = name) %>%
  ungroup() %>%
  mutate(model = "recent")

scales <- c(scale_soap_pen, scale_coeff_penalty)

model_summ_scale <- model_summ_prev %>%
  filter(str_detect(name, "lambda")) %>%
  mutate(model = "previous_scaled") %>%
  mutate(value = value / scales)

comp_tbl <- bind_rows(model_summ_prev, model_summ_rec, model_summ_scale) %>%
  arrange(name, num, model)
