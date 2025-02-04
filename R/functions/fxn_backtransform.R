#' Data Transformation Utility
#'
#' @description 
#' Applies inverse transformations to index values based on specified transformation type.
#'
#' @param index_data A data frame containing the index values to transform
#' @param index_value Name of the column to be transformed
#' @param index_transform Type of transformation to apply ('sqrt', 'log', or default)
#'
#' @return Modified data frame with back-transformed values
#' @details 
#' Supports three transformation types:
#' - 'sqrt': Squares the original value
#' - 'log': Applies exponential transformation with a small offset
#' - Default: Returns original value unchanged
#'
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @export
fxn_backtransform <- function(index_data, index_value, index_transform) {
  
  # Validate required packages
  required_packages <- c("dplyr", "glue")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Generate a name for the back-transformed column
  bt_name <- glue::glue("bt_{index_value}")
  
  # Apply appropriate back-transformation
  index_data %>%
    dplyr::mutate(
      "{bt_name}" := dplyr::case_when(
        index_transform == "sqrt" ~ (!!rlang::sym(index_value))^2, 
        index_transform == "log" ~ exp(!!rlang::sym(index_value) - 1e-6), 
        TRUE ~ !!rlang::sym(index_value)
      )
    )
}
