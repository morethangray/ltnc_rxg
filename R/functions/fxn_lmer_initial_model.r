#' Fit Initial Models for Data Analysis
#' 
#' This function fits various linear and mixed-effects models to the data and
#' selects the best model based on AIC values.
#' 
#' @param data A data frame containing the response variables and predictors
#' 
#' @return A list containing:
#'   \item{best_model}{List with details of the best performing model:
#'     \itemize{
#'       \item{response}{Name of the response variable}
#'       \item{model}{The fitted model object}
#'       \item{model_name}{Name of the best model}
#'       \item{formula}{Model formula}
#'       \item{aic}{AIC value of the best model}
#'     }
#'   }
#'   \item{summary_table}{Data frame containing AIC values for all fitted models}
#'
#' @import dplyr
#' @import lme4
#' @import MuMIn
#' @importFrom stats AIC formula lm terms
#'
#' @export
#'
#' @examples
#' \dontrun{
#' model_results <- fxn_initial_model(analysis_data)
#' }
fxn_lmer_initial_model <- function(data) {
  # Validate required packages
  required_packages <- c("dplyr", "lme4", "MuMIn")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }

  responses <- c("value_std", "value_log", "value_sqrt")
  best_models <- list()
  
  # Create empty data frame for summary of all models
  summary_df <- data.frame(
    response = character(),
    model_name = character(),
    aic = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (response in responses) {
    models <- list(
      lm = stats::lm(stats::as.formula(paste(response, "~ treatment")), data = data),
      lm_null = stats::lm(stats::as.formula(paste(response, "~ 1")), data = data),
      lmer_1 = lme4::lmer(stats::as.formula(paste(response, "~ treatment + (1 | plot_name)")), 
                         data = data, REML = FALSE),
      lmer_2 = lme4::lmer(stats::as.formula(paste(response, "~ treatment + (1 + treatment | plot_name)")), 
                         data = data, REML = FALSE),
      lmer_null_1 = lme4::lmer(stats::as.formula(paste(response, "~ 1 + (1 | plot_name)")), 
                              data = data, REML = FALSE),
      lmer_null_2 = lme4::lmer(stats::as.formula(paste(response, "~ 1 + (1 + treatment | plot_name)")), 
                              data = data, REML = FALSE)
    )
    
    # Calculate AIC values for all models
    aic_values <- stats::AIC(models$lm, models$lm_null)
    rownames(aic_values) <- c("lm", "lm_null")
    mixed_model_selection <- MuMIn::model.sel(
      lmer_1 = models$lmer_1, 
      lmer_2 = models$lmer_2,
      lmer_null_1 = models$lmer_null_1,
      lmer_null_2 = models$lmer_null_2
    ) 
    
    # Add linear models to summary
    lm_summary <- data.frame(
      response = response,
      model_name = rownames(aic_values),
      aic = aic_values$AIC,
      stringsAsFactors = FALSE
    )
    
    # Add mixed models to summary
    mixed_summary <- data.frame(
      response = response,
      model_name = rownames(mixed_model_selection),
      aic = mixed_model_selection$AICc,
      stringsAsFactors = FALSE
    )
    
    # Combine all models for this response
    response_summary <- dplyr::bind_rows(lm_summary, mixed_summary)
    summary_df <- dplyr::bind_rows(summary_df, response_summary)
    
    # Find best model for this response
    if (min(aic_values$AIC) < min(mixed_model_selection$AICc)) {
      best_model_name <- names(models)[which.min(aic_values$AIC)]
      best_model <- models[[best_model_name]]
      best_aic <- min(aic_values$AIC)
    } else {
      best_model_name <- rownames(mixed_model_selection)[which.min(mixed_model_selection$AICc)]
      best_model <- models[[best_model_name]]
      best_aic <- min(mixed_model_selection$AICc)
    }
    
    best_models[[response]] <- list(
      response = response,
      model = best_model,
      model_name = best_model_name,
      formula = stats::formula(best_model),
      aic = best_aic
    )
  }
  
  best_response <- names(best_models)[which.min(sapply(best_models, function(x) x$aic))]
  best_overall <- best_models[[best_response]]
  
  # Return both the best model and the complete summary table
  return(list(
    best_model = best_overall,
    summary_table = summary_df %>%
      dplyr::arrange(aic)
  ))
}
