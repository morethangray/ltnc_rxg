# Load libraries, functions, workflows -----
rm(list = ls())
#
library(tidyverse) ## To manipulate data frames
library(here) ## To manage directories
#
source(here("R/functions/fxn_utilities.R"))
#
# ========================================================== -----
# PREPARE SURVEY DATA  ----
source(here("R/functions/fxn_xlsx_to_csv.R"))
# Collate and reshape surveys  
list_years <- 2019:2022
all_surveys <- fxn_xlsx_to_csv(list_years)
#
# ========================================================== -----
# CALCULATE RICHNESS AND ABUNDANCE  -----
source(here("R/functions/fxn_calculate_metrics.R"))
# Define surveys to use ----
# Exclude years 2016 to 2018 and any inactive plots
input_surveys <- all_surveys %>%
  filter(index %in% list_active_surveys)

# Calculate richness ----
rich_frb <- fxn_rich(
  index_subset = "frb",
  index_native = "n1"
)
rich_nat <- fxn_rich(
  index_subset = "nat",
  index_native = "n1"
)
rich_non <- fxn_rich(
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
    path_in_data_derived, 
    "rich_nat-frb-non_2019-2022.csv"
  ))

# Calculate abundance ----
abun_frb <- fxn_abun(
  index_subset = "frb",
  index_native = "n1"
)
abun_nat <- fxn_abun(
  index_subset = "nat",
  index_native = "n1"
)
abun_non <- fxn_abun(
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
    path_in_data_derived,
    "abun_nat-frb-non_2019-2022.csv"
  ))

# ========================================================== -----