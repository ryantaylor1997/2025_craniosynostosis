################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Clean image/observation-level characteristics file
################################################################################

source(here::here("source", "000_definitions.R"))

# Load files --------------------------------------------------------------

# Read in covariates data
patient_images <- read_csv(file.path(DATA_FOLDER, "merged_covariates.csv"))

# List all files in folder with any extension (so no folders)
matrix_files <- list.files(MATRIX_FOLDER, pattern = "\\.")

# Classify sutures --------------------------------------------------------

# Currently stored in binary columns with names ending in "fused"
# Reshape data to identify all combinations of sutures
fusion_df <- patient_images %>%
  clean_names() %>%
  select(id_fmt, matches("fused")) %>%
  pivot_longer(cols = -id_fmt) %>%
  # Reformat fusion names
  mutate(
    name = name %>%
      str_remove("_fused") %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  ) %>%
  mutate(fusion = if_else(value == 1, name, NA_character_)) %>%
  group_by(id_fmt) %>%
  summarize(num_fusions = sum(value, na.rm = T),
            fusions = paste(na.omit(fusion), collapse = "/")) %>%
  ungroup()

# Clean data --------------------------------------------------------------

cranio <- patient_images %>%
  clean_names() %>%
  ## Check if files in data were found in data folder
  mutate(
    fname = paste0(id_fmt, "_desc-unsigneddist.RData"),
    fname_exists = if_else(fname %in% matrix_files, 1, 0)
  ) %>%
  ## Classify sutures
  # (1) count number of sutures
  # (2) identify bicoronal (the only useful instance of multiple sutures)
  # (3) identify other instances of multiple sutures
  # (4) assign remaining suture fusion types
  left_join(fusion_df, by = "id_fmt") %>%
  mutate(fusion_type = case_when(
    num_fusions == 0 ~ "Normative",
    fusions == "Left Coronal/Right Coronal" ~ "Bicoronal",
    num_fusions > 1 ~ "Multiple",
    T ~ fusions
  )) %>%
  mutate(fusion_type = fct_infreq(fusion_type)) %>%
  mutate(fusion_type = fct_relevel(fusion_type, "Multiple", after = Inf)) %>%
  ## Identify which visit each image is from and how many visits they had
  arrange(id_base, age) %>%
  group_by(id_base) %>%
  mutate(visit_num = 1:n()) %>%
  # Count visits total and visits under age limit (currently 1 year)
  mutate(visits = max(visit_num)) %>%
  ungroup()

save(cranio,
     file = here("data", "intermediate", "obs_data_complete.rda"))

# Create subset data ------------------------------------------------------

# Exclude rare fusion categories: multiple, lambdoid
# Exclude sphenoidal involvement (in Craniosynostosis group without suture)
# Exclude patients with previous surgery (10)
# Exclude patients above age cap of 3 yrs
# Only keep last picture taken
cranio_sub <- cranio %>%
  filter(!fusion_type %in% c("Multiple",
                             "Left Lambdoid", "Right Lambdoid") &
           sphenoidal_involvement %in% c(NA, 0) &
           previous_surgery %in% c(NA, 0) &
           age <= cranio_max_age &
           visit_num == visits) %>%
  # Reduce number of columns
  select(id_base, fname, fusion_type,
         age, sex) %>%
  arrange(id_base, age) %>%
  # Make sex variable categorical; drop unused fusion types
  mutate(sex = factor(sex),
         fusion_type = fct_drop(fusion_type))

save(cranio_sub,
     file = here("data", "intermediate", "obs_data_subset.rda"))

# Check data against folder -----------------------------------------------

# Print a table to check we have all files
fname_summ <- tibble(
  in_data_not_folder = sum(1 - cranio$fname_exists),
  in_folder_not_data = sum(!matrix_files %in% cranio$fname)
)

# Summarize Fusion Types --------------------------------------------------

# Summarize counts of images and patients by fusion type
fusion_summ_sub <- cranio_sub %>%
  group_by(fusion_type) %>%
  summarize(pics = n(),
            patients = n_distinct(id_base)) %>%
  ungroup() %>%
  arrange(desc(pics)) %>%
  mutate(pics_pct = pics / sum(pics),
         patients_pct = patients / sum(patients)) %>%
  adorn_totals("row")

### Summarize frequencies of each suture fusion type
fusion_summ <- cranio %>%
  group_by(fusion_type) %>%
  summarize(pics = n(),
            patients = n_distinct(id_base)) %>%
  ungroup() %>%
  arrange(desc(pics)) %>%
  mutate(pics_pct = pics / sum(pics),
         patients_pct = patients / sum(patients)) %>%
  adorn_totals("row")

# Investigate fusion type -------------------------------------------------

# Sanity check: fusion category against collaborators' "group" column
fusion_group_summ <- tabyl(cranio, fusion_type, group)

# 6 in "Craniosynostosis" group with no sutures all have sphenoidal involvement
sphen_summ <- cranio %>%
  group_by(fusion_type, group, sphenoidal_involvement) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  arrange(desc(n))

## Investigate multiple sutures
fusion_multiple_summ <- cranio %>%
  filter(fusion_type == "Multiple") %>%
  group_by(num_fusions, fusions) %>%
  summarize(pics = n(),
            patients = n_distinct(id_base)) %>%
  ungroup() %>%
  arrange(desc(pics))

max_multiple_count <- max(fusion_multiple_summ$pics)

## Investigate previous surgery
surg_summ <- cranio %>%
  filter(previous_surgery == 1) %>%
  group_by(fusion_type, group, sphenoidal_involvement) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  adorn_totals("row")

# Summarize visits --------------------------------------------------------

# Summarize visits on average and max by fusion group
visits_summ <- cranio %>%
  group_by(fusion_type) %>%
  summarize(patients = n_distinct(id_base),
            pics = n(),
            visits_avg = mean(visits),
            visits_max = max(visits)) %>%
  ungroup() %>%
  arrange(desc(visits_avg))

# Visualize Age Ranges ----------------------------------------------------

# Set colors for fusions to remain consistent in future plots
fusion_color_dict <- distinct(cranio, fusion_type) %>%
  arrange(fusion_type) %>%
  mutate(color = c(palette.colors(palette = "Set1")[-c(6, 9)],
                   palette.colors(palette = "Set1")[c(6, 9)]))

save(fusion_color_dict,
     file = here::here("data", "intermediate",
                       "fusion_color_mapping.rda"))

# Plot age densities
age_plot <- ggplot(cranio) +
  geom_density_ridges(aes(x = age,
                          y =  fusion_type,
                          fill = fusion_type)) +
  scale_fill_discrete(guide = "none",
                      palette = fusion_color_dict$color) +
  labs(x = "Age (Days)", y = "Fusion Type")

# Summarize densities
have_fusions <- cranio %>%
  filter(fusion_type != "Normative")

have_fusions_quant <- have_fusions %$%
  quantile(age, probs = c(0.9, 0.95, 0.975, 0.99, 1)) %>%
  enframe(name = "Percentile (with Fusions)", value = "Age") %>%
  mutate(Years = Age / 365.25)

# Determine what percentile is our preset max age among patients with fusion
max_age_summ <- tibble("Age Cap" = cranio_max_age,
                       "Percentile" = mean(have_fusions$age <= cranio_max_age))

