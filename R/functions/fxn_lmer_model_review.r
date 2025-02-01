#' Review Model Performance
#' 
#' This function performs diagnostic checks on the fitted model and creates
#' residual plots for included terms.
#' 
#' @param index_model The fitted model object to review
#' @param data The data frame used to fit the model
#' 
#' @return No return value, called for side effects (plots and diagnostic tests)
#'
#' @import performance
#' @importFrom stats terms
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fxn_model_review(fitted_model, analysis_data)
#' }
fxn_lmer_model_review <- function(index_model, data) {
  # Validate required packages
  required_packages <- c("performance", "see")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }

  model_terms <- attr(stats::terms(index_model), "term.labels")
  
  # Diagnostic checks
  performance::check_model(index_model, check = "all")
  DHARMa::testDispersion(index_model)
  DHARMa::testUniformity(index_model)
  DHARMa::testOutliers(index_model)
  
  # Residual plots only for included terms
  for (term in model_terms) {
    if (term %in% names(data)) {
      DHARMa::plotResiduals(index_model, data[[term]])
    }
  }
}
