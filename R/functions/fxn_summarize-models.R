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