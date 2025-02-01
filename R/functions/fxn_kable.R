#' Create Formatted Kable Table
#'
#' Generates a formatted HTML table using knitr and kableExtra
#'
#' @param df A data frame to be converted to a formatted table
#' @return A formatted HTML table
#' 
#' @import knitr
#' @import kableExtra
#' @export
fxn_kable <- function(data) {
  # Validate required packages
  required_packages <- c("knitr", "kableExtra")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Create formatted table
  data %>%
    knitr::kable() %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"), 
      full_width = TRUE,  
      position = "left", 
      fixed_thead = TRUE
    )
}

# Example usage
# formatted_table <- fxn_kable(your_dataframe)
