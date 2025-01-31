# ========================================================== -----
# ---------------------------------------------------------- -----
# !!! UPDATE THESE IF YOU CHANGE XLSX FILE OR SHEET NAMES !!! ----
# ---------------------------------------------------------- -----
# ---------------------------------------------------------- -----
# DEFINE FILE AND SHEET NAMES ----
# Plant attributes xlsx file
xlsx_name_plants <- "attributes_plant.xlsx"
sheet_name_plants <- "crosswalk"
#
# Plot attributes xlsx file
xlsx_name_plots <- "attributes_plots.xlsx"
# Plot names 
sheet_name_plots <- "lookup-plot-name"
# Grazing interval
sheet_name_grazing <- "grazing-interval"
# Active surveys 
sheet_name_surveys <- "surveyed-plots"
# Plot attributes
sheet_name_attributes <- "plot-attributes"
# 
# ========================================================== -----

#' Create Lookup Tables from Excel Sources
#'
#' Reads and processes lookup tables from multiple Excel sheets
#'
#' @param xlsx_name_plants Filename for plants Excel file
#' @param sheet_name_plants Sheet name for plants data
#' @param xlsx_name_plots Filename for plots Excel file
#' @param sheet_name_plots Sheet name for plot names
#' @param sheet_name_grazing Sheet name for grazing intervals
#' @param sheet_name_surveys Sheet name for surveys
#' @param sheet_name_attributes Sheet name for plot attributes
#'
#' @return List of processed lookup tables
#' @import here
#' @import readxl
#' @import dplyr
#' 
#' @export
fxn_create_lookup_tables <- function(
    xlsx_name_plants = "attributes_plant.xlsx", 
    sheet_name_plants = "crosswalk",
    xlsx_name_plots = "attributes_plots.xlsx",
    sheet_name_plots = "lookup-plot-name",
    sheet_name_grazing = "grazing-interval",
    sheet_name_surveys = "surveyed-plots",
    sheet_name_attributes = "plot-attributes"
) {
  # Validate packages
  required_packages <- c("here", "readxl", "dplyr")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Source file path definition 
  source(here::here("R/functions/fxn_define-file-paths.R"))
  
  # Validate file existence
  excel_files <- c(
    file.path(project_paths$path_in_lookup, xlsx_name_plants),
    file.path(project_paths$path_in_lookup, xlsx_name_plots)
  )
  
  for(file_path in excel_files) {
    if(!file.exists(file_path)) {
      stop(paste("File not found:", file_path))
    }
  }
  
  # Create lookup tables
  lookup_plants <- readxl::read_excel(
    here::here(project_paths$path_in_lookup, xlsx_name_plants),
    sheet = sheet_name_plants
  ) %>%
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
  
  lookup_plot_names <- readxl::read_excel(
    here::here(project_paths$path_in_lookup, xlsx_name_plots),
    sheet = sheet_name_plots
  )
  
  lookup_grazing_interval <- readxl::read_excel(
    here::here(project_paths$path_in_lookup, xlsx_name_plots),
    sheet = sheet_name_grazing
  ) %>%
    dplyr::select(
      index,
      starts_with("f_"),
      grazer,
      interval,
      in_v1,
      -f_grazed
    ) %>%
    dplyr::distinct()
  
  list_active_surveys <- readxl::read_excel(
    here::here(project_paths$path_in_lookup, xlsx_name_plots),
    sheet = sheet_name_surveys
  ) %>%
    dplyr::filter(
      include_plot == TRUE,
      include_data == TRUE,
      year > 2018
    ) %>%
    dplyr::pull(index)
  
  lookup_annotation <- readxl::read_excel(
    here::here(project_paths$path_in_lookup, xlsx_name_plots),
    sheet = sheet_name_attributes
  ) %>%
    dplyr::mutate(f_new = dplyr::if_else(interval == "New plot", "n1", "n0")) %>%
    dplyr::filter(
      index %in% list_active_surveys
    ) %>%
    dplyr::distinct(
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
  
  # Return list of lookup tables
  return(list(
    lookup_plants = lookup_plants,
    lookup_plot_names = lookup_plot_names,
    lookup_grazing_interval = lookup_grazing_interval,
    list_active_surveys = list_active_surveys,
    lookup_annotation = lookup_annotation
  ))
}

# Run function with default parameters
lookup_tables <- fxn_create_lookup_tables()
