#' Load and Subset Richness and Abundance Data
#' 
#' This function loads richness and abundance data from CSV files and processes
#' them into separate dataframes for native, non-native, and forb species.
#' 
#' @param project_paths List of file paths from fxn_setup_file_paths()
#' 
#' @details
#' This function requires project paths that must be generated first:
#' ```r
#' project_paths <- fxn_setup_file_paths()
#' ```
#' The function processes both richness and abundance data, creating separate
#' dataframes for different species categories (native, non-native, forb).
#' 
#' @return A list containing the following data frames:
#'   \item{richness}{List of richness data frames:
#'     \itemize{
#'       \item{rich_nat}{Native species richness}
#'       \item{rich_frb}{Native forb richness}
#'       \item{rich_non}{Non-native species richness}
#'     }
#'   }
#'   \item{abundance}{List of abundance data frames:
#'     \itemize{
#'       \item{abun_nat}{Native species abundance}
#'       \item{abun_frb}{Native forb abundance}
#'       \item{abun_non}{Non-native species abundance}
#'     }
#'   }
#' 
#' @import here
#' @import dplyr
#' @import readr
#' @importFrom stats sqrt
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # First generate required objects
#' project_paths <- fxn_setup_file_paths()
#' 
#' # Then create richness and abundance tables
#' rich_abun <- fxn_load_rich_abun(project_paths = project_paths)
#' }
fxn_load_rich_abun <- function(project_paths) {
  # Check for required packages
  required_packages <- c("here", "dplyr", "readr")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Load and process richness data
  rich <- readr::read_csv(
    here::here(
      project_paths$path_in_data_derived,
      "richness.csv"
    )
  ) %>%
    dplyr::arrange(year, f_break) %>%
    dplyr::mutate(
      value_sqrt = sqrt(value),
      dplyr::across(where(is.character), as.factor) %>%
      dplyr::mutate(treatment = factor(treatment, levels = c("Ungrazed", "Grazed")))
    )
  
  # Create richness subsets
  rich_nat <- dplyr::filter(rich, met_sub %in% "rich_nat")
  rich_frb <- dplyr::filter(rich, met_sub %in% "rich_frb")
  rich_non <- dplyr::filter(rich, met_sub %in% "rich_non")
  
  # Load and process abundance data
  abun <- readr::read_csv(
    here::here(
      project_paths$path_in_data_derived,
      "abundance.csv"
    )
  ) %>%
    dplyr::arrange(dplyr::desc(treatment), year, f_break) %>%
    dplyr::mutate(
      dplyr::across(where(is.character), as.factor) %>%
        dplyr::mutate(treatment = factor(treatment, levels = c("Ungrazed", "Grazed")))
    )
  
  # Create abundance subsets
  abun_nat <- dplyr::filter(abun, met_sub %in% "abun_nat")
  abun_frb <- dplyr::filter(abun, met_sub %in% "abun_frb")
  abun_non <- dplyr::filter(abun, met_sub %in% "abun_non")
  
  # Return results as a list
  return(list(
    richness = list(
      rich_nat = rich_nat,
      rich_frb = rich_frb,
      rich_non = rich_non
    ),
    abundance = list(
      abun_nat = abun_nat,
      abun_frb = abun_frb,
      abun_non = abun_non
    )
  ))
}

