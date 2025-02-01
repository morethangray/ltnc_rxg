#' Fit and Select Fixed Effects Models
#' 
#' This function fits various fixed effects models and selects the best model
#' based on AIC values.
#' 
#' @param index_model The initial model to build upon
#' 
#' @return A list containing:
#'   \item{best_model}{List with details of the best performing model:
#'     \itemize{
#'       \item{model}{The fitted model object}
#'       \item{model_name}{Name of the best model}
#'       \item{terms}{Terms included in the model}
#'       \item{terms_count}{Number of terms in the model}
#'       \item{formula}{Model formula}
#'       \item{aic}{AIC value of the best model}
#'     }
#'   }
#'   \item{summary_table}{Data frame containing AIC values for all fitted models}
#'
#' @import dplyr
#' @importFrom stats AIC formula update
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fixed_effects_results <- fxn_fixed_effects(initial_model)
#' }
fxn_lmer_fixed_effects <- function(index_model) {
  # Validate required packages
  required_packages <- c("dplyr")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }

  # Create empty data frame for summary of all models
  summary_df <- data.frame(
    model_name = character(),
    terms_count = integer(),
    terms = character(),
    aic = numeric(),
    model = I(list())
  )
  
  # One-term models
  model_terms_one <- c("plot_type", "f_break", "f_new", "f_one_yr", "f_two_yr", "f_year")
  fix_one <- lapply(model_terms_one, function(term) {
    model <- stats::update(index_model, stats::as.formula(paste(". ~ . +", term)))
    aic <- stats::AIC(model)
    data.frame(
      model_name = term,
      terms_count = 1,
      terms = term,
      aic = aic,
      model = I(list(model))
    )
  }) %>% dplyr::bind_rows()
  
  # Two-term models
  model_terms_two <- list(
    c("f_year", "plot_type"), c("f_year", "f_break"), c("f_year", "f_new"),
    c("f_year", "f_one_yr"), c("f_year", "f_two_yr")
  )
  
  fix_two <- lapply(model_terms_two, function(terms) {
    model <- stats::update(index_model, stats::as.formula(paste(". ~ . +", paste(terms, collapse = " + "))))
    aic <- stats::AIC(model)
    model_name <- paste(terms, collapse = " + ")
    data.frame(
      model_name = model_name,
      terms_count = 2,
      terms = model_name,
      aic = aic,
      model = I(list(model))
    )
  }) %>% dplyr::bind_rows()
  
  # Combine all models
  all_models <- dplyr::bind_rows(fix_one, fix_two)
  
  # Best overall model
  best_model_row <- all_models %>%
    dplyr::arrange(aic) %>%
    dplyr::slice(1)
  
  best_model <- best_model_row$model[[1]]
  best_model_name <- best_model_row$model_name
  best_aic <- best_model_row$aic
  best_terms <- best_model_row$terms
  best_terms_count <- best_model_row$terms_count
  
  # Return both the best model and the complete summary table
  return(list(
    best_model = list(
      model = best_model,
      model_name = best_model_name,
      terms = best_terms,
      terms_count = best_terms_count,
      formula = stats::formula(best_model),
      aic = best_aic
    ),
    summary_table = all_models %>% dplyr::arrange(aic)
  ))
}
