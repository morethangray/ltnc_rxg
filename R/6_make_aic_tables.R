# Run setup scripts -----
source(here::here("R/1_setup.R"))

# Define functions for model selection ----
# Fit the initial models
source(here::here("R/functions/fxn_lmer_initial_model.R"))

# Review model diagnostics 
source(here::here("R/functions/fxn_lmer_model_review.R"))

# Evaluate fixed effects
source(here::here("R/functions/fxn_lmer_fixed_effects.R"))

# Evaluate random effects
source(here::here("R/functions/fxn_lmer_random_effects.R"))

# Write functions ----
get_aic_lmer_initial_model <- function(index_data_name) {
  
  index_data <- get(index_data_name)
  
  # Define model responses
  responses <- c("value_std", "value_log", "value_sqrt")
  all_results <- list()
  
  for (response in responses) {
    models <- list(
      lm = stats::lm(stats::as.formula(paste(response, "~ treatment")), data = index_data),
      lm_null = stats::lm(stats::as.formula(paste(response, "~ 1")), data = index_data),
      lmer_1 = lme4::lmer(stats::as.formula(paste(response, "~ treatment + (1 | plot_name)")), 
                          data = index_data, REML = FALSE),
      lmer_2 = lme4::lmer(stats::as.formula(paste(response, "~ treatment + (1 + treatment | plot_name)")), 
                          data = index_data, REML = FALSE),
      lmer_null_1 = lme4::lmer(stats::as.formula(paste(response, "~ 1 + (1 | plot_name)")), 
                               data = index_data, REML = FALSE),
      lmer_null_2 = lme4::lmer(stats::as.formula(paste(response, "~ 1 + (1 + treatment | plot_name)")), 
                               data = index_data, REML = FALSE)
    )
    
    # Calculate AIC values for all models
    aic_values <- stats::AIC(models$lm, models$lm_null)
    rownames(aic_values) <- c("lm", "lm_null")
    
    # Add linear models to summary
    lm_summary <- tibble::tibble(
      model_id = index_data_name,
      response = response,
      model_name = rownames(aic_values),
      aic = aic_values$AIC, 
      is_singular = NA  # Not applicable to lm
    )
    
    # Mixed model selection
    mixed_model_selection <- MuMIn::model.sel(
      lmer_1 = models$lmer_1, 
      lmer_2 = models$lmer_2,
      lmer_null_1 = models$lmer_null_1,
      lmer_null_2 = models$lmer_null_2
    )
    
    # Check singularity for lmer models
    singular_flags <- c(
      lme4::isSingular(models$lmer_1),
      lme4::isSingular(models$lmer_2),
      lme4::isSingular(models$lmer_null_1),
      lme4::isSingular(models$lmer_null_2)
    )
    
    # Add mixed models to summary
    mixed_summary <- tibble::tibble(
      model_id = index_data_name,
      response = response,
      model_name = rownames(mixed_model_selection),
      aic = mixed_model_selection$AICc,
      is_singular = singular_flags 
    )
    
    # Combine all models for this response
    all_results[[response]] <- dplyr::bind_rows(lm_summary, mixed_summary)
  }
  # Combine results across all responses
  return(dplyr::bind_rows(all_results))
}
aic_init_abun_nat <- get_aic_lmer_initial_model(index_data_name = "abun_nat") 
aic_init_abun_frb <- get_aic_lmer_initial_model(index_data_name = "abun_frb") 
aic_init_abun_non <- get_aic_lmer_initial_model(index_data_name = "abun_non") 

get_aic_lmer_fixed_effects <- function(index_model_name) {
  
  index_model <- get(index_model_name)
  
  # Define model terms
  model_terms_one <- model_terms_one <- c("plot_type", "f_break", "f_new", "f_one_yr", "f_two_yr", "f_year")
  model_terms_two <- list(
    c("f_year", "plot_type"),
    c("f_year", "f_break"),
    c("f_year", "f_new"),
    c("f_year", "f_one_yr"),
    c("f_year", "f_two_yr")
  )
  
  # Helper function to fit model and summarize results
  fit_model <- function(terms) {
    model_formula <- stats::as.formula(paste(". ~ . +", paste(terms, collapse = " + ")))
    model <- stats::update(index_model, model_formula)
    tibble::tibble(
      model_id = index_model_name,
      model_terms = paste(terms, collapse = " + "),
      aic = stats::AIC(model),
      is_singular = lme4::isSingular(model),
      n_terms = length(terms),
      model_call = paste(deparse(model@call), collapse = " ")
    ) 
  }

  # Apply to one-term and two-term models
  results_one <- purrr::map_dfr(model_terms_one, ~ fit_model(.x))
  results_two <- purrr::map_dfr(model_terms_two, ~ fit_model(.x))
  
  # Combine all models
  all_models <- dplyr::bind_rows(results_one, results_two) |>
    mutate(response = 
             case_when(
               str_detect(model_call, "value_sqrt") ~ "value_sqrt", 
               str_detect(model_call,"value_log") ~ "value_log", 
               str_detect(model_call,"value_std") ~ "value_std")) |>
    relocate(response, .after = model_id)
  
  return(all_models)
  
  
}
aic_fixed_abun_nat <- get_aic_lmer_fixed_effects("mod_abun_nat") 
aic_fixed_abun_frb <- get_aic_lmer_fixed_effects("mod_abun_frb") 
aic_fixed_abun_non <- get_aic_lmer_fixed_effects("mod_abun_non") 

get_aic_lmer_random_effects <- function(index_model_name) {
  
  index_model <- get(index_model_name)
  
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
  
  random_terms <- attr(all_models, "random.terms") 
  enframe_random_terms <- tibble::enframe(random_terms, name = "model_name", value = "random_terms")  
  
  all_models_tbl <-   all_models %>%
    tibble::as_tibble() %>%
    # Add rownames as a column since they contain model names
    tibble::rownames_to_column("model") |>
    bind_cols(enframe_random_terms)
  
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
  
  # Define model terms
  model_terms_one <- model_terms_one <- c("plot_type", "f_break", "f_new", "f_one_yr", "f_two_yr", "f_year")
  model_terms_two <- list(
    c("f_year", "plot_type"),
    c("f_year", "f_break"),
    c("f_year", "f_new"),
    c("f_year", "f_one_yr"),
    c("f_year", "f_two_yr")
  )
  
  # Helper function to fit model and summarize results
  fit_model <- function(terms) {
    model_formula <- stats::as.formula(paste(". ~ . +", paste(terms, collapse = " + ")))
    model <- stats::update(index_model, model_formula)
    tibble::tibble(
      model_id = index_model_name,
      model_terms = paste(terms, collapse = " + "),
      aic = stats::AIC(model),
      is_singular = lme4::isSingular(model),
      n_terms = length(terms),
      model_call = paste(deparse(model@call), collapse = " ")
    ) 
  }
  
  # Apply to one-term and two-term models
  results_one <- purrr::map_dfr(model_terms_one, ~ fit_model(.x))
  results_two <- purrr::map_dfr(model_terms_two, ~ fit_model(.x))
  
  # Combine all models
  all_models <- dplyr::bind_rows(results_one, results_two) |>
    mutate(response = 
             case_when(
               str_detect(model_call, "value_sqrt") ~ "value_sqrt", 
               str_detect(model_call,"value_log") ~ "value_log", 
               str_detect(model_call,"value_std") ~ "value_std")) |>
    relocate(response, .after = model_id)
  
  return(all_models)
  
  
}
# ========================================================== -----
# ========================================================== -----
# Create input data subsets & fit models ----
#   Abundance ----
abun <- rich_abun$abundance  
abun_nat <- abun$abun_nat 
abun_frb <- abun$abun_frb
abun_non <- abun$abun_non

mod_abun_nat <- 
  lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), 
             data = abun_nat, REML = FALSE)
mod_abun_frb <- 
  lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), 
             data = abun_frb, REML = FALSE)

mod_abun_non <- 
  lme4::lmer(value_sqrt ~ treatment + f_year + f_two_yr + (1 | plot_name), 
             data = abun_non, REML = FALSE)

# ---------------------------------------------------------- -----
#   Richness ----
rich <- rich_abun$richness
rich_nat <- rich$rich_nat
rich_frb <- rich$rich_frb
rich_non <- rich$rich_non

mod_rich_nat <- 
  glmmTMB::glmmTMB(
    value ~ treatment + f_year + (1 | plot_type) + (1 + treatment | plot_name),
    data = rich_nat, family = nbinom2,
    control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

mod_rich_frb <- 
  glmmTMB::glmmTMB(
    value ~ treatment + f_year + (1 | plot_type) + (1 + treatment | plot_name),
    data = rich_frb, family = nbinom2,
    control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

mod_rich_non <- 
  lme4::lmer(
    value_sqrt ~ treatment + f_year + (1 + treatment | plot_name), 
    data = rich_non, REML = FALSE)

# ========================================================== -----
