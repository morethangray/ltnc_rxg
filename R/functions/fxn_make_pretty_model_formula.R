#' Make Model Formula Readable
#'
#' This function replaces shorthand variable names in a model formula with more descriptive terms and removes `"1 + "` when it appears inside parentheses.
#'
#' @param index_formula A character string representing a model formula.
#'
#' @return A character string with replaced terms and formatted structure.
#' @export
#'
#' @examples
#' index_formula <- "treatment + f_year + plot_type + (1 + treatment | plot_name)"
#' fxn_make_pretty_model_formula(index_formula)
fxn_make_pretty_model_formula <- function(index_formula){
  
  # Validate required packages
  required_packages <- c("stringr")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  replacements <- c(
    "treatment" = "Grazing treatment",
    "f_year" = "Survey year",
    "plot_type" = "Plot type",
    "plot_name" = "Plot pair",
    "f_break" = "Grazing interval",
    "f_new" = "New plot",
    "f_one_yr" = "One-year interval",
    "f_two_yr" = "Two-year interval",
    "grazer" = "Grazer"
    
  )
  # Apply replacements
  pretty_formula <- stringr::str_replace_all(index_formula, replacements)
  
  # Remove "1 + " when it appears inside parentheses
  pretty_formula <- stringr::str_replace_all(pretty_formula, "\\(1 \\+ ", "(")
  
  return(pretty_formula)
}
