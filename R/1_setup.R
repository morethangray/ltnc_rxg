rm(list = ls())

# Set seed for reproducible simulations during model selection
set.seed(912)

# Manage dependencies
source(here::here("R/functions/fxn_setup_dependencies.R"))

# Define file paths
source(here::here("R/functions/fxn_setup_file_paths.R"))

# Create lookup tables
source(here::here("R/functions/fxn_setup_lookup_tables.R"))

# Prepare data from xlsx 
source(here::here("R/functions/fxn_prepare_input_data.R"))

# Create the processed data 
processed_data <- fxn_prepare_input_data(list_years = 2019:2022, 
                                         project_paths = project_paths, 
                                         lookup_tables = lookup_tables)

# Calculate metrics for richness and abundance
source(here::here("R/functions/fxn_calculate_metrics.R"))

# Load richness and abundance data
source(here::here("R/functions/fxn_load_rich_abun.R"))

# Run function to load richness and abundance data
rich_abun <- fxn_load_rich_abun(project_paths = project_paths)