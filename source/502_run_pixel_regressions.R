################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run pixel-level regression models
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Load file of pixel-level data for all images ("cranio_point_clean")
load(file = here("data", "cleaned", "point_data_clean.rda"))

# Run models --------------------------------------------------------------

## Run models for pointwise differences in each pixel
pixel_models <- cranio_point_clean %>%
  mutate(
    model = map(data,
                function(d){
                  d_reg <- d %>% left_join(fusion_dict, by = "fusion_type")

                  mdl <- bam(diff ~
                               fused_Sagittal + fused_Metopic +
                               fused_RCoronal * fused_LCoronal +
                               sex + s(age, bs = "cr"),
                             data = d_reg)

                  return(mdl)
                })) %>%
  mutate(coeffs = map(
    model,
    function(m){

      # Make constrast vector identifying all coronal variables
      is_coronal <- as.numeric(str_detect(names(coef(m)), "Coronal"))

      # Compute variance of bicoronal combination
      bicoronal_v <- crossprod(is_coronal, vcov(m)) %*% is_coronal
      bicoronal_e <- crossprod(is_coronal, coef(m))

      # Add bicoronal as row in tidy summary
      bicoronal_row <- tibble(
        term = "fused_Bicoronal",
        estimate = bicoronal_e[1,1],
        std.error = sqrt(bicoronal_v[1,1])) %>%
        mutate(statistic = estimate / std.error,
               p.value = 2 * pnorm(abs(estimate / std.error),
                                   lower.tail = FALSE))

      # Output coefficient table (smaller than model)
      return(tidy(m, parametric = TRUE) %>%
               bind_rows(bicoronal_row))
    }))

# Save all models
save(pixel_models,
     file = here("results", "models_pixel_regressions.rda"))
