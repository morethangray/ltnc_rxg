source(here::here("R/functions/fxn_install_dependencies.R"))
# fxn_initial_model: Function to fit initial models ----
fxn_initial_model <- function(data) {
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
      lm = lm(as.formula(paste(response, "~ treatment")), data = data),
      lm_null = lm(as.formula(paste(response, "~ 1")), data = data),
      lmer_1 = lmer(as.formula(paste(response, "~ treatment + (1 | plot_name)")), data = data, REML = FALSE),
      lmer_2 = lmer(as.formula(paste(response, "~ treatment + (1 + treatment | plot_name)")), data = data, REML = FALSE),
      lmer_null_1 = lmer(as.formula(paste(response, "~ 1 + (1 | plot_name)")), data = data, REML = FALSE),
      lmer_null_2 = lmer(as.formula(paste(response, "~ 1 + (1 + treatment | plot_name)")), data = data, REML = FALSE)
    )
    
    # Calculate AIC values for all models
    aic_values <- AIC(models$lm, models$lm_null)
    rownames(aic_values) <- c("lm", "lm_null")
    mixed_model_selection <- model.sel(lmer_1 = models$lmer_1, 
                                       lmer_2 = models$lmer_2,
                                       lmer_null_1 = models$lmer_null_1,
                                       lmer_null_2 = models$lmer_null_2) 
    
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
    response_summary <- rbind(lm_summary, mixed_summary)
    summary_df <- rbind(summary_df, response_summary)
    
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
      formula = formula(best_model),
      aic = best_aic
    )
  }
  
  best_response <- names(best_models)[which.min(sapply(best_models, function(x) x$aic))]
  best_overall <- best_models[[best_response]]
  
  # Return both the best model and the complete summary table
  return(list(
    best_model = best_overall,
    summary_table = summary_df %>%
      arrange(aic)
  ))
}

# fxn_model_review: Function to review model performance ----
fxn_model_review <- function(model, data) {
  model_terms <- attr(terms(model), "term.labels")
  
  # Diagnostic checks
  testDispersion(model)
  check_singularity(model)
  testUniformity(model)
  testOutliers(model)
  
  # Residual plots only for included terms
  for (term in model_terms) {
    if (term %in% names(data)) {
      plotResiduals(model, data[[term]])
    }
  }
}

# fxn_fixed_effects: Function to fit and select fixed effects models ----
fxn_fixed_effects <- function(input_model) {
  # Create empty data frame for summary of all models
  summary_df <- data.frame(
    model_name = character(),
    terms_count = integer(), # Number of terms
    terms = character(),
    aic = numeric(),
    model = I(list()) # Store model as a list column
  )
  
  # One-term models
  model_terms_one <- c("plot_type", "f_break", "f_new", "f_one_yr", "f_two_yr", "f_year")
  fix_one <- lapply(model_terms_one, function(term) {
    model <- update(input_model, as.formula(paste(". ~ . +", term)))
    aic <- AIC(model)
    data.frame(
      model_name = term,
      terms_count = 1,
      terms = term,
      aic = aic,
      model = I(list(model))
    )
  }) %>% bind_rows()
  
  # Two-term models
  model_terms_two <- list(
    c("f_year", "plot_type"), c("f_year", "f_break"), c("f_year", "f_new"),
    c("f_year", "f_one_yr"), c("f_year", "f_two_yr")
  )
  
  fix_two <- lapply(model_terms_two, function(terms) {
    model <- update(input_model, as.formula(paste(". ~ . +", paste(terms, collapse = " + "))))
    aic <- AIC(model)
    model_name <- paste(terms, collapse = " + ")
    data.frame(
      model_name = model_name,
      terms_count = 2,
      terms = model_name, # Store combined term name
      aic = aic,
      model = I(list(model))
    )
  }) %>% bind_rows()
  
  # Combine all models
  all_models <- bind_rows(fix_one, fix_two)
  
  # Best overall model
  best_model_row <- all_models %>%
    arrange(aic) %>%
    slice(1)
  
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
      formula = formula(best_model),
      aic = best_aic
    ),
    summary_table = all_models %>% arrange(aic) # Sort summary table by AIC
  ))
}

# fxn_random_effects: Function to fit random effects  ----
fxn_random_effects <- function(best_fixed_model) {
  
  # Add each random intercept one at a time
  plot_type <- update(best_fixed_model, . ~ . + (1 | plot_type))
  f_year <- update(best_fixed_model, . ~ . + (1 | f_year))
  f_year_grazer <- update(best_fixed_model, . ~ . + (1 | f_year/grazer))
  grazer <- update(best_fixed_model, . ~ . + (1 | grazer))
  f_break <- update(best_fixed_model, . ~ . + (1 | f_break))
  f_one_yr <- update(best_fixed_model, . ~ . + (1 | f_one_yr))
  f_two_yr <- update(best_fixed_model, . ~ . + (1 | f_two_yr))
  f_new <- update(best_fixed_model, . ~ . + (1 | f_new))
  # 
  # # Add each random slope one at a time
  f_year_plot_type <- update(best_fixed_model, . ~ . + (1 + f_year | plot_type))
  treatment_f_year <- update(best_fixed_model, . ~ . + (1 + treatment | f_year))
  treatment_f_year_grazer <- update(best_fixed_model, . ~ . + (1 + treatment | f_year/grazer))
  treatment_grazer <- update(best_fixed_model, . ~ . + (1 + treatment | grazer))
  treatment_f_break <- update(best_fixed_model, . ~ . + (1 + treatment | f_break))
  treatment_f_one_yr <- update(best_fixed_model, . ~ . + (1 + treatment | f_one_yr))
  treatment_f_two_yr  <- update(best_fixed_model, . ~ . + (1 + treatment | f_two_yr))
  treatment_f_new <- update(best_fixed_model, . ~ . + (1 + treatment | f_new))
  
  
  # Check for singularity using a robust method 
  list_model_names <- c("best_fixed_model", 
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
    
    datalist[[n]] <- tibble(
      model_name = n, 
      is_singular = isSingular(model)
    )
  }
  
  bind_datalist <- do.call(bind_rows, datalist)
  
  # not_singular <-  bind_datalist %>%
  #   filter(is_singular == FALSE) %>%
  #   filter(model_name != "best_fixed_model") 
  # 
  # if(nrow(not_singular) == 0){
  #   print("All random effect models are singular")
  # }else{
  #   print(not_singular)
  # }
  
  all_models <- model.sel(best_fixed_model, 
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
  
  names(all_models)
  
  all_models_tbl <-   all_models %>%
    as_tibble() %>%
    # Add rownames as a column since they often contain model names
    tibble::rownames_to_column("model")
  
  # To see what columns are being dropped:
  dropped_cols <- setdiff(names(all_models), names(all_models_tbl))
  
  summary_table <- tibble(
    model_name = rownames(all_models),
    AICc = all_models$AICc,
    delta_AICc = round(all_models$delta, 2),
    weight = round(all_models$weight, 4),
    df = all_models$df
  ) %>%
    arrange(AICc) %>%
    left_join(bind_datalist, "model_name")
}