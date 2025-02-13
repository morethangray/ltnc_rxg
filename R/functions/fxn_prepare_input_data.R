#' Prepare Input Data from Excel Files
#'
#' Reads, processes, and combines data from multiple Excel files across specified years.
#' Reshapes data to long format and joins with lookup tables for plant and plot attributes.
#'
#' @param list_years Vector of years to process
#' @param lookup_tables List of lookup tables from fxn_setup_lookup_tables()
#' @param project_paths List of file paths from fxn_setup_file_paths()
#'
#' @details
#' This function requires project paths and lookup tables that must be generated first:
#' ```r
#' project_paths <- fxn_setup_file_paths()
#' lookup_tables <- fxn_setup_lookup_tables()
#' ```
#'
#' @return A data frame of combined and processed yearly data
#' @import here
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import readxl
#' @import glue
#' @export
#' #'
#' @examples
#' \dontrun{
#' # First generate required objects
#' project_paths <- fxn_setup_file_paths()
#' lookup_tables <- fxn_setup_lookup_tables()
#' 
#' # Then process data
#' processed_data <- fxn_prepare_input_data(
#'   list_years = 2019:2022,
#'   project_paths = project_paths,
#'   lookup_tables = lookup_tables
#' )}
fxn_prepare_input_data <- function(list_years, project_paths, lookup_tables) {
  # Validate required packages
  required_packages <- c("here", "dplyr", "tidyr", "purrr", "readxl", "glue")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Create custom not-in operator if it doesn't exist
  `%nin%` <- function(x, y) !(`%in%`(x, y))
  
  # Initialize list for yearly data
  datalist_year <- list()
  
  # index_year = list_years[4]
  for (index_year in list_years) {
    # Read xlsx ----
    path <- here::here(
      project_paths$path_in_data,
      glue::glue("wantrup_{index_year}.xlsx")
    )
    
    # Validate file existence
    if (!file.exists(path)) {
      stop(glue::glue("File not found: {path}"))
    }
    
    input_data <- path %>%
      readxl::excel_sheets() %>%
      purrr::set_names() %>%
      purrr::map_dfr(readxl::read_excel,
                     path = path,
                     .id = "sheet"
      )
    
    # Reshape and annotate survey table ----
    column_count <- ncol(input_data)
    
    # Gather into long format
    long_data <- input_data %>%
      tidyr::gather(orig_plot, value, 3:all_of(column_count)) %>%
      tidyr::drop_na(value) %>%
      dplyr::mutate(
        year = index_year,
        value = as.numeric(value)
      ) %>%
      # Join plant attributes
      dplyr::left_join(lookup_tables$lookup_plants, "orig_name") %>%
      dplyr::filter(include_plant == TRUE) %>%
      dplyr::mutate(
        native = dplyr::if_else(is_native == TRUE, "Native", "Non-native"),
        forb = dplyr::if_else(is_forb == TRUE, "Forb", "Non-forb"),
        f_native = dplyr::if_else(is_native == TRUE, "n1", "n0"),
        f_forb = dplyr::if_else(is_forb == TRUE, "f1", "f0")
      ) %>%
      # Join plot attributes
      dplyr::left_join(lookup_tables$lookup_plot_names, "orig_plot") %>%
      dplyr::mutate(
        index = paste(plot_type, plot_name, year, sep = "_"),
        index_g = paste(index, f_grazed, sep = "_"),
        index_g_n_f = paste(index_g, f_native, f_forb, sep = "_")
      ) %>%
      # Join temporal attributes
      dplyr::left_join(lookup_tables$lookup_grazing_interval, "index") %>%
      # Reorder columns
      dplyr::select(
        plot_type,
        plot_name,
        year,
        month,
        treatment,
        grazer,
        native,
        forb,
        interval,
        genus_species,
        value,
        dplyr::starts_with("f_"),
        dplyr::starts_with("index"),
        taxonomic_name
      ) %>%
      dplyr::distinct()
    
    # Combine counts for taxa with ssp. or var. ----
    list_multi_gs <- long_data %>%
      dplyr::distinct(genus_species, taxonomic_name) %>%
      dplyr::group_by(genus_species) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      dplyr::filter(n > 1) %>%
      dplyr::pull(genus_species)
    
    one_gs <- long_data %>%
      dplyr::filter(genus_species %nin% list_multi_gs) %>%
      dplyr::select(-taxonomic_name)
    
    # Calculate the maximum count for Apr & May surveys ----
    percent_cover <- long_data %>%
      dplyr::filter(genus_species %in% list_multi_gs) %>%
      dplyr::group_by(dplyr::across(c(-value, -taxonomic_name))) %>%
      dplyr::summarize(value = sum(value), .groups = "drop") %>%
      dplyr::bind_rows(one_gs) %>%
      dplyr::distinct() %>%
      dplyr::group_by(dplyr::across(c(-month, -value))) %>%
      dplyr::summarize(value = max(value), .groups = "drop") %>%
      dplyr::relocate(value, .after = genus_species)
    
    datalist_year[[index_year]] <- percent_cover
  }
  
  bind_years <- do.call(dplyr::bind_rows, datalist_year) %>%
    dplyr::filter(index %in% lookup_tables$list_active_surveys)
  return(bind_years)
}