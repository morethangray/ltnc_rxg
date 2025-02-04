#' Summarize Marginal Means for Statistical Models
#'
#' @description 
#' Computes and formats marginal means for a given statistical model
#'
#' @param index_model_name Name of the model object to analyze
#' @param lookup_tables A list of lookup tables for additional metadata
#'
#' @return A data frame with summarized marginal means and associated metadata
#' 
#' @details 
#' This function performs the following key steps:
#' 1. Retrieves the specified model object
#' 2. Determines response type (abundance or richness)
#' 3. Computes back-transformed marginal means using emmeans
#' 4. Applies delta method for standard error handling
#' 5. Joins with lookup tables for additional context
#'
#' @importFrom emmeans emmeans
#' @importFrom broom tidy
#' @importFrom dplyr mutate left_join select
#' @importFrom stringr str_detect str_sub str_remove
#' @export
fxn_summarize_marginal_means <- function(index_model_name, lookup_tables) {
  
  # Validate required packages
  required_packages <- c("broom", "dplyr", "emmeans", "janitor", "stringr")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Validate inputs
  stopifnot(
    !missing(index_model_name),
    !missing(lookup_tables),
    exists(index_model_name)
  )
  
  # Retrieve the model object
  index_model <- get(index_model_name)
  
  # Extract model metadata
  model_formula <- as.character(formula(index_model))[3]
  response <- if(stringr::str_detect(index_model_name, "abun")) "abundance" else "richness"
  subset <- stringr::str_sub(index_model_name, -3)
  transformation <- stringr::str_remove(as.character(formula(index_model))[2], "value_")
  
  # Compute marginal means
  marginal_means <- 
    emmeans::emmeans(index_model, ~ treatment) %>%
    broom::tidy() %>%
    janitor::clean_names() %>%
    fxn_backtransform(index_value = "estimate", index_transform = transformation) %>%
    fxn_backtransform(index_value = "std_error", index_transform = transformation) %>%
    dplyr::mutate(
      bt_std_error = std_error * bt_std_error, 
      response = response, 
      abbr_subset = subset,
      model_name = index_model_name, 
      model_formula = model_formula
    ) %>%
    dplyr::left_join(lookup_tables$lookup_model_subset, "abbr_subset") %>%
    dplyr::select(
      response, 
      subset, 
      treatment, 
      starts_with("bt"), 
      p_value, 
      statistic,
      df, 
      model_name, 
      model_formula, 
      starts_with("abbr_")
    )
  
  return(marginal_means)
}

