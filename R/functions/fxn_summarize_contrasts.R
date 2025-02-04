#' Compute Statistical Contrasts for Model Comparisons
#'
#' @description 
#' Performs pairwise comparisons and computes contrasts for statistical models
#'
#' @param index_model_name Name of the model object to analyze
#' @param lookup_tables A list of lookup tables for additional metadata
#'
#' @return A data frame with computed contrasts and associated statistical information
#' 
#' @details 
#' This function:
#' 1. Extracts model predictors
#' 2. Computes pairwise comparisons using marginaleffects
#' 3. Applies back-transformation
#' 4. Joins with lookup tables for comprehensive context
#'
#' @note Requires marginaleffects package version 0.5.0
#' @importFrom marginaleffects comparisons
#' @importFrom broom tidy
#' @importFrom dplyr group_by summarize mutate left_join select
#' @export
fxn_summarize_contrasts <- function(index_model_name, lookup_tables) {
  
  # Validate required packages
  required_packages <- c("broom", "dplyr", "janitor", "marginaleffects", "stringr")
  
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
  
  # Extract predictor variables
  predictors <- attr(terms(index_model), "term.labels")
  
  # Create pairwise comparison list
  variable_list <- setNames(rep(list("pairwise"), length(predictors)), predictors)
  
  # Compute and process contrasts
  contrasts <- 
    marginaleffects::comparisons(index_model, variables = variable_list) %>%
    broom::tidy() %>%
    janitor::clean_names() %>%
    group_by(term, contrast) %>%
    summarize(
      estimate = mean(estimate), 
      p_value = min(p_value), 
      statistic = mean(statistic), 
      s_value = mean(s_value), 
      conf_low = mean(conf_low), 
      conf_high = mean(conf_high)
    ) %>%
    fxn_backtransform(index_value = "estimate", index_transform = transformation) %>%
    fxn_backtransform(index_value = "conf_low", index_transform = transformation) %>%
    fxn_backtransform(index_value = "conf_high", index_transform = transformation) %>%
    dplyr::rename(
      abbr_term = term, 
      abbr_contrast = contrast
    ) %>%
    dplyr::mutate(
      response = response, 
      abbr_subset = subset,
      model_name = index_model_name, 
      model_formula = model_formula
    ) %>%
    dplyr::left_join(lookup_tables$lookup_model_subset, "abbr_subset") %>%
    dplyr::left_join(lookup_tables$lookup_model_term, "abbr_term") %>%
    dplyr::left_join(lookup_tables$lookup_model_contrast, "abbr_contrast") %>%
    dplyr::select(
      response, 
      subset,
      term,
      contrast,
      starts_with("bt"), 
      p_value, 
      statistic,
      estimate_raw = estimate, 
      model_name, 
      model_formula, 
      starts_with("abbr_")
    )
  
  return(contrasts)
}