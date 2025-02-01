# PREPARE SURVEY DATA  ----
source(here::here("R/functions/fxn_xlsx_to_csv.R"))
# Collate and reshape surveys  
processed_data <- fxn_prepare_input_data(
  list_years = 2019:2022,
  project_paths = project_paths,
  lookup_tables = lookup_tables
)
#
# ========================================================== -----
# CALCULATE RICHNESS AND ABUNDANCE  -----
source(here("R/functions/fxn_calculate_metrics.R"))
# Define surveys to use ----
# Exclude years 2016 to 2018 and any inactive plots
input_surveys <- processed_data %>%
  filter(index %in% list_active_surveys)

# Calculate richness ----
rich_frb <- fxn_rich(
  index_data = processed_data,
  index_subset = "frb",
  index_native = "n1"
)
rich_nat <- fxn_rich(
  index_data = processed_data,
  index_subset = "nat",
  index_native = "n1"
)
rich_non <- fxn_rich(
  index_data = processed_data,
  index_subset = "non",
  index_native = "n0"
)

rich_subsets3 <-
  bind_rows(
    rich_non,
    rich_nat,
    rich_frb
  ) %>%
  write_csv(here(
    project_paths$path_in_data_derived, 
    "richness.csv"
  ))

# Calculate abundance ----
abun_frb <- fxn_abun(
  index_data = processed_data,
  index_subset = "frb",
  index_native = "n1"
)
abun_nat <- fxn_abun(
  index_data = processed_data,
  index_subset = "nat",
  index_native = "n1"
)
abun_non <- fxn_abun(
  index_data = processed_data,
  index_subset = "non",
  index_native = "n0"
)
abun_subsets3 <-
  bind_rows(
    abun_non,
    abun_nat,
    abun_frb
  ) %>%
  write_csv(here(
    project_paths$path_in_data_derived,
    "abundance.csv"
  ))

# ========================================================== -----