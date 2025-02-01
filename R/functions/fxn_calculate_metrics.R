#' Calculate Species Richness
#'
#' Calculates species richness for a specified subset of data
#'
#' @param index_data Processed data frame from fxn_prepare_input_data()
#' @param index_subset Type of subset ("frb", "nat", or "non")
#' @param index_native Native status filter ("n1" or "n0")
#' @param lookup_tables List of lookup tables from fxn_setup_lookup_tables()
#' @return A data frame of richness calculations
#' @import dplyr
#' @export
fxn_rich <- function(index_data, index_subset, index_native, lookup_tables) {
  # Validate packages
  required_packages <- c("dplyr")
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing")
  }
  
  # Input validation
  if(!index_subset %in% c("frb", "nat", "non")) {
    stop("index_subset must be 'frb', 'nat', or 'non'")
  }
  if(!index_native %in% c("n1", "n0")) {
    stop("index_native must be 'n1' or 'n0'")
  }
  
  # Create subset based on parameters
  subset <- if(index_subset == "frb") {
    index_data %>%
      dplyr::filter(f_native %in% "n1" & f_forb %in% "f1")
  } else {
    index_data %>%
      dplyr::filter(f_native %in% index_native)
  }
  
  # Calculate richness
  subset %>%
    dplyr::group_by(index_g) %>%
    dplyr::summarize(value = dplyr::n_distinct(genus_species)) %>%
    dplyr::ungroup() %>%
    dplyr::right_join(lookup_tables$lookup_annotation, "index_g") %>%
    dplyr::mutate(
      subset = index_subset,
      metric = "rich",
      met_sub = paste(metric, subset, sep = "_"),
      value = dplyr::if_else(is.na(value), 0, value),
      value_log = log(value + 1e-6), 
      value_std = as.vector(scale(value)),
      value_sqrt = sqrt(value)
    ) %>%
    dplyr::select(
      treatment,
      dplyr::starts_with("value"),
      plot_name,
      plot_type,
      year,
      grazer,
      f_year,
      f_break,
      f_new,
      f_one_yr,
      f_two_yr,
      met_sub
    )
}

#' Calculate Species Abundance
#'
#' Calculates species abundance for a specified subset of data
#'
#' @param index_data Processed data frame from fxn_prepare_input_data()
#' @param index_subset Type of subset ("frb", "nat", or "non")
#' @param index_native Native status filter ("n1" or "n0")
#' @param lookup_tables List of lookup tables from fxn_setup_lookup_tables()
#' @return A data frame of abundance calculations
#' @import dplyr
#' @export
fxn_abun <- function(index_data, index_subset, index_native, lookup_tables) {
  # Validate packages
  required_packages <- c("dplyr")
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing")
  }
  
  # Input validation
  if(!index_subset %in% c("frb", "nat", "non")) {
    stop("index_subset must be 'frb', 'nat', or 'non'")
  }
  if(!index_native %in% c("n1", "n0")) {
    stop("index_native must be 'n1' or 'n0'")
  }
  
  # Create subset based on parameters
  subset <- if(index_subset == "frb") {
    index_data %>%
      dplyr::filter(f_native %in% "n1" & f_forb %in% "f1")
  } else {
    index_data %>%
      dplyr::filter(f_native %in% index_native)
  }
  
  # Calculate abundance
  subset %>%
    dplyr::group_by(index_g) %>%
    dplyr::summarize(value = sum(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::right_join(lookup_tables$lookup_annotation, "index_g") %>%
    dplyr::mutate(
      subset = index_subset,
      metric = "abun",
      met_sub = paste(metric, subset, sep = "_"),
      value = dplyr::if_else(is.na(value), 0, value),
      value_log = log(value + 1e-6),
      value_std = as.vector(scale(value)),
      value_sqrt = sqrt(value)
    ) %>%
    dplyr::select(
      treatment,
      dplyr::starts_with("value"),
      plot_name,
      plot_type,
      year,
      grazer,
      f_year,
      f_break,
      f_new,
      f_one_yr,
      f_two_yr,
      met_sub
    )
}

#' Calculate and Save Richness and Abundance Metrics
#'
#' Calculates richness and abundance for different subsets and saves to CSV
#'
#' @param processed_data Processed data frame from fxn_create_input_data()
#' @param lookup_tables List of lookup tables from fxn_setup_lookup_tables()
#' @param project_paths List of file paths from fxn_setup_file_paths()
#' @return NULL (saves files to disk)
#' @import dplyr readr here
#' @export
fxn_calculate_metrics <- function(processed_data, lookup_tables, project_paths) {
  # Validate packages
  required_packages <- c("dplyr", "readr", "here")
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing")
  }
  
  # Calculate richness for all subsets
  rich_frb <- fxn_rich(
    index_data = processed_data,
    index_subset = "frb",
    index_native = "n1",
    lookup_tables = lookup_tables 
  )
  
  rich_nat <- fxn_rich(
    index_data = processed_data,
    index_subset = "nat",
    index_native = "n1",
    lookup_tables = lookup_tables 
  )
  
  rich_non <- fxn_rich(
    index_data = processed_data,
    index_subset = "non",
    index_native = "n0",
    lookup_tables = lookup_tables 
  )
  
  # Combine and save richness
  dplyr::bind_rows(rich_non, rich_nat, rich_frb) %>%
    readr::write_csv(here::here(
      project_paths$path_in_data_derived, 
      "richness.csv"
    ))
  
  # Calculate abundance for all subsets
  abun_frb <- fxn_abun(
    index_data = processed_data,
    index_subset = "frb",
    index_native = "n1",
    lookup_tables = lookup_tables 
  )
  
  abun_nat <- fxn_abun(
    index_data = processed_data,
    index_subset = "nat",
    index_native = "n1",
    lookup_tables = lookup_tables 
  )
  
  abun_non <- fxn_abun(
    index_data = processed_data,
    index_subset = "non",
    index_native = "n0",
    lookup_tables = lookup_tables 
  )
  
  # Combine and save abundance
  dplyr::bind_rows(abun_non, abun_nat, abun_frb) %>%
    readr::write_csv(here::here(
      project_paths$path_in_data_derived,
      "abundance.csv"
    ))
}

# Run script to calculate metrics
fxn_calculate_metrics(processed_data = processed_data, 
                      project_paths = project_paths, 
                      lookup_tables = lookup_tables)