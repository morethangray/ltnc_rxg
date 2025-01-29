# ========================================================== -----
# Load libraries ----
# File management ---
library(glue)
library(readxl)   ## To read xlsx files
# library(lubridate)   ## To work with dates and times
# library(stringr)   ## To wrangle character variables
# library(testthat)  ## For data checks
# library(openxlsx)   ## To write xlsx files
# library(fs)   ## To manage directories
# library(hms)  ## For working with times
# library(janitor)   ## To clean data tables

#
# Define basic functions ----
#   %nin%: Negate %in% 
"%nin%" <- Negate("%in%")
#   substrRight: Subset character string from the right 
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}
# Define file paths ----
fxn_define_file_paths <- function(){
  #
  # For R scripts and work
  path_r <- here()
  path_in <- here("input")
  #
  # Input files
  path_in_data_raw <<- here(path_in, "data_raw")
  path_in_data_derived <<- here(path_in, "data_derived")
  path_in_lookup <<- here(path_in, "lookup-tables")
  #
  # Output files
  path_out <<- here("output")
  path_out_explore <<- here(path_out, "1_explore")
  path_out_prepare <<- here(path_out, "2_prepare")
  path_out_summary <<- here(path_out, "3_summary")
  
  path_out_prepare <<- here(path_out, "0_prepare")
  path_out_clean <<- here(path_out, "1_clean")
  path_out_explore <<- here(path_out, "2_explore")
  path_out_abundance <<- here(path_out, "4_abundance")
  path_out_richness <<- here(path_out, "5_richness")
  #
}
#
# Execute function to define file paths
fxn_define_file_paths()
# ========================================================== -----
# Create helpers ----
#   lookup_plants ----
lookup_plants <-
  read_excel(here(path_in_lookup, "attributes_plant.xlsx"),
             sheet = "crosswalk"
  ) %>%
  # specify package because distinct() may be blocked
  dplyr::distinct(
    orig_name,
    taxonomic_name,
    genus_species,
    f_native,
    f_forb,
    include_plant,
    has_mult_taxa,
    is_native,
    is_forb
  )
#   lookup_plot_names ----
lookup_plot_names <-
  read_xlsx(here(path_in_lookup, "attributes_plots.xlsx"),
            sheet = "lookup-plot-name"
  )
#   lookup_grazing_interval ----
lookup_grazing_interval <-
  read_xlsx(here(path_in_lookup, "attributes_plots.xlsx"),
            sheet = "grazing-interval"
  ) %>%
  # specify package because distinct() may be blocked
  dplyr::select(
    index,
    starts_with("f_"),
    grazer,
    interval,
    in_v1,
    -f_grazed
  ) %>%
  dplyr::distinct()

#   list_active_surveys ----
# Define active surveys (2019 to 2022)
# Exclude years 2016 to 2018 and any inactive plots
list_active_surveys <-
  read_xlsx(here(path_in_lookup, "attributes_plots.xlsx"),
            sheet = "surveyed-plots"
  ) %>%
  filter(
    include_plot == TRUE,
    include_data == TRUE,
    year > 2018,
  ) %>%
  pull(index)

#   lookup_annotation ----
lookup_annotation <-
  read_xlsx(here(path_in_lookup, "attributes_plots.xlsx"),
            sheet = "plot-attributes"
  ) %>%
  mutate(f_new = ifelse(interval == "New plot", "n1", "n0")) %>%
  filter(
    index %in% list_active_surveys,
  ) %>%
  distinct(
    treatment,
    plot_name,
    plot_type,
    year,
    interval,
    grazer,
    f_year,
    f_break,
    f_new,
    f_one_yr,
    f_two_yr,
    index,
    index_g
  )



# ========================================================== -----
