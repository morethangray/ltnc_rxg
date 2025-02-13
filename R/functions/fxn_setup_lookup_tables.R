#' Create Lookup Tables from Excel Sources
#'
#' Reads and processes lookup tables from multiple Excel sheets
#'
#' @param csv_plants File name for plant attributes csv file
#' @param csv_plots File name for plot attributes csv file
#' @param csv_surveys File name for surveys by year csv file
#' @param csv_plot_names File name for plot names csv file
#' @param csv_subsets File name for species subsets csv file
#' @param csv_terms File name for model term labels csv file
#' @param csv_means File name for marginal mean labels csv file
#' @param csv_contrasts File name for model contrast labels csv file
#'
#' @return List of processed lookup tables
#' @import here
#' @import readxl
#' @import dplyr
#' 
#' @export
fxn_setup_lookup_tables <- function(
    csv_plants = "plant_attributes.csv", 
    csv_plots = "plot_attributes.csv", 
    csv_surveys = "surveys_by_year.csv", 
    csv_plot_names = "plot_names.csv", 
    csv_grazing_interval = "grazing_interval.csv",
    csv_subsets = "species_subsets.csv", 
    csv_terms = "model_terms.csv", 
    csv_means = "model_means.csv", 
    csv_contrasts = "model_contrasts.csv"
) {
  # Validate packages
  required_packages <- c("here", "readxl", "dplyr")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Source file path definition 
  source(here::here("R/functions/fxn_setup_file_paths.R"))
  
  # Validate file existence
  check_files <- c(
    file.path(project_paths$path_in_lookup, csv_plants),
    file.path(project_paths$path_in_lookup, csv_plots), 
    file.path(project_paths$path_in_lookup, csv_surveys), 
    file.path(project_paths$path_in_lookup, csv_plot_names), 
    file.path(project_paths$path_in_lookup, csv_grazing_interval), 
    file.path(project_paths$path_in_lookup, csv_subsets), 
    file.path(project_paths$path_in_lookup, csv_terms), 
    file.path(project_paths$path_in_lookup, csv_means), 
    file.path(project_paths$path_in_lookup, csv_contrasts)
  )
  
  for(file_path in check_files) {
    if(!file.exists(file_path)) {
      stop(paste("File not found:", file_path))
    }
  }
  
  # Create lookup tables
  lookup_plants <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_plants)
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
  
  lookup_plot_names <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_plot_names)
  )
  
  lookup_grazing_interval <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_grazing_interval)
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
  
  list_active_surveys <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_surveys)
  ) %>%
    dplyr::filter(
      include_plot == TRUE,
      include_data == TRUE,
      year > 2018
    ) %>%
    dplyr::pull(index)
  
  lookup_plots <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_plots)
  ) %>%
    dplyr::mutate(f_new = dplyr::if_else(interval == "New plot", "n1", "n0")) %>%
    dplyr::filter(index %in% list_active_surveys) %>%
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
  
  lookup_species_subset <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_subsets)
  ) %>%
    dplyr::select(subset, abbr_subset)
  
  lookup_model_term <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_terms)
  ) %>%
    dplyr::select(term, abbr_term)
  
  lookup_model_mean <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_means)
  ) %>%
    dplyr::select(value, abbr_value)
  
  lookup_model_contrast <- readr::read_csv(
    here::here(project_paths$path_in_lookup, csv_contrasts)
  ) %>%
    dplyr::select(contrast, abbr_contrast)
  

  # Return list of lookup tables
  return(list(
    lookup_plants = lookup_plants,
    lookup_plot_names = lookup_plot_names,
    lookup_grazing_interval = lookup_grazing_interval,
    list_active_surveys = list_active_surveys,
    lookup_plots = lookup_plots, 
    lookup_species_subset = lookup_species_subset, 
    lookup_model_term = lookup_model_term, 
    lookup_model_mean = lookup_model_mean, 
    lookup_model_contrast = lookup_model_contrast
    
  ))
}


# Run function with default parameters
lookup_tables <- fxn_setup_lookup_tables()
