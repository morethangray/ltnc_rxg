#' Install Missing Project Dependencies
#'
#' Checks and installs any missing packages required for the project
#'
#' @return NULL
#' @export
fxn_install_dependencies <- function() {
  packages_needed <- c(
    "here", 
    "readxl", 
    "dplyr", 
    "forcats", 
    "sessioninfo", 
    "lme4", 
    "MuMIn", 
    "ggplot2", 
    "DHARMa", 
    "performance", 
    "glmmTMB", 
    "broom.mixed",
    "knitr", 
    "kableExtra", 
    "car", 
    "MASS"
  )
  
  # Identify and install missing packages
  new_packages <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
  if(length(new_packages) > 0) {
    message("Installing missing packages: ", paste(new_packages, collapse=", "))
    install.packages(new_packages, dependencies = TRUE)
  }
  
  # Optionally, load packages silently
  suppressPackageStartupMessages({
    invisible(lapply(packages_needed, library, character.only = TRUE))
  })
}

# Optional: Run dependency check on script load
fxn_install_dependencies()

#' Check Package Availability
#'
#' Verifies if specified packages are installed
#'
#' @param packages Character vector of package names to check
#' @return Logical indicating whether all packages are available
#' @export
fxn_check_packages <- function(packages) {
  missing_packages <- setdiff(packages, installed.packages()[,"Package"])
  
  if(length(missing_packages) > 0) {
    warning(
      "Missing packages: ", 
      paste(missing_packages, collapse=", "), 
      "\nRun fxn_install_dependencies() to install them."
    )
    return(FALSE)
  }
  return(TRUE)
}
