#' Fit Random Effects Models
#' 
#' This function fits various random effects models and evaluates their performance
#' using AIC values and singularity checks.
#' 
#' @param index_model The initial model to build upon
#' 
#' @return A tibble containing model comparison results and singularity checks
#'
#' @import dplyr
#' @import lme4
#' @import MuMIn
#' @import tibble
#'
#' @export
#'
#' @examples
#' \dontrun{
#' random_effects_results <- fxn_random_effects(initial_model)
#' }
fxn_lmer_random_effects <- function(index_model) {
  # Validate required packages
  required_packages <- c("dplyr", "lme4", "MuMIn", "tibble")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }

  # Add each random intercept one at a time
  plot_type <- stats::update(index_model, . ~ . + (1 | plot_type))
  f_year <- stats::update(index_model, . ~ . + (1 | f_year))
  f_year_grazer <- stats::update(index_model, . ~ . + (1 | f_year/grazer))
  grazer <- stats::update(index_model, . ~ . + (1 | grazer))
  f_break <- stats::update(index_model, . ~ . + (1 | f_break))
  f_one_yr <- stats::update(index_model, . ~ . + (1 | f_one_yr))
  f_two_yr <- stats::update(index_model, . ~ . + (1 | f_two_yr))
  f_new <- stats::update(index_model, . ~ . + (1 | f_new))
  
  # Add each random slope one at a time
  f_year_plot_type <- stats::update(index_model, . ~ . + (1 + f_year | plot_type))
  treatment_f_year <- stats::update(index_model, . ~ . + (1 + treatment | f_year))
  treatment_f_year_grazer <- stats::update(index_model, . ~ . + (1 + treatment | f_year/grazer))
  treatment_grazer <- stats::update(index_model, . ~ . + (1 + treatment | grazer))
  treatment_f_break <- stats::update(index_model, . ~ . + (1 + treatment | f_break))
  treatment_f_one_yr <- stats::update(index_model, . ~ . + (1 + treatment | f_one_yr))
  treatment_f_two_yr <- stats::update(index_model, . ~ . + (1 + treatment | f_two_yr))
  treatment_f_new <- stats::update(index_model, . ~ . + (1 + treatment | f_new))
  
  # Check for singularity using a robust method 
  list_model_names <- c("index_model", 
                        "plot_type", 
                        "f_year", 
                        "f_year_grazer", 
                        "grazer", 
                        "f_break", 
                        "f_one_yr",
                        "f_two_yr",
                        "f_new",
                        "f_year_plot_type",
                        "treatment_f_year",
                        "treatment_f_year_grazer",
                        "treatment_grazer",
                        "treatment_f_break",
                        "treatment_f_one_yr",
                        "treatment_f_two_yr",
                        "treatment_f_new")
  
  datalist <- list()
  for(n in list_model_names){
    
    model <- get(n)
    
    datalist[[n]] <- tibble::tibble(
      model_name = n, 
      is_singular = lme4::isSingular(model)
    )
  }
  
  bind_datalist <- do.call(bind_rows, datalist)
  
  all_models <- MuMIn::model.sel(index_model, 
                                 plot_type, 
                                 f_year, 
                                 f_year_grazer, 
                                 grazer, 
                                 f_break, 
                                 f_one_yr,
                                 f_two_yr,
                                 f_new,
                                 f_year_plot_type,
                                 treatment_f_year,
                                 treatment_f_year_grazer,
                                 treatment_grazer,
                                 treatment_f_break,
                                 treatment_f_one_yr,
                                 treatment_f_two_yr,
                                 treatment_f_new)  
  
  all_models_tbl <-   all_models %>%
    tibble::as_tibble() %>%
    # Add rownames as a column since they contain model names
    tibble::rownames_to_column("model")
  
  # To see what columns are being dropped:
  dropped_cols <- setdiff(names(all_models), names(all_models_tbl))
  
  summary_table <- tibble::tibble(
    model_name = rownames(all_models),
    AICc = all_models$AICc,
    delta_AICc = round(all_models$delta, 2),
    weight = round(all_models$weight, 4),
    df = all_models$df
  ) %>%
    dplyr::arrange(AICc) %>%
    dplyr::left_join(bind_datalist, "model_name")
}