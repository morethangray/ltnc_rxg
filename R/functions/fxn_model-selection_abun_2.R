# Function to fit initial models ----
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

# Function to review model performance ----
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


# Function to review model performance with combined plots [NOT WORKING]----
fxn_model_review_grid <- function(model, data) {
  # Simulated residuals
  sim_res <- simulateResiduals(model)
  
  # Generate diagnostic plots (using as.ggplot to convert base R plots)
  p1 <- ggplotify::as.ggplot(function() plot(sim_res) + ggtitle("Simulated Residuals"))
  p2 <- ggplotify::as.ggplot(function() plotResiduals(sim_res, rank = TRUE)) + ggtitle("Ranked Residuals") # Convert to ggplot
  p3 <- ggplotify::as.ggplot(function() plotResiduals(sim_res)) + ggtitle("Residuals vs Predicted") # Convert to ggplot
  p4 <- ggplotify::as.ggplot(function() qqnorm(residuals(model)) + qqline()) + ggtitle("Normal Q-Q") # Convert to ggplot
  
  # Combine main diagnostic plots (2x2 layout)
  diagnostic_plot <- p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 2)
  
  # Print combined diagnostic plots
  print(diagnostic_plot)
  
  # Residual plots for included terms (using as.ggplot here too)
  model_terms <- attr(terms(model), "term.labels")
  residual_plots <- list()
  for (term in model_terms) {
    if (term %in% names(data)) {
      residual_plots[[term]] <- ggplotify::as.ggplot(function() plotResiduals(sim_res, data[[term]])) + 
        ggtitle(paste("Residuals vs", term)) # Convert to ggplot
    }
  }
  
  # Print residual plots if available, also in a 2x2 grid
  if (length(residual_plots) > 0) {
    print(patchwork::wrap_plots(residual_plots, ncol = 2))
  }
}



# Function to fit and select fixed effects models ----
input_model <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)

fxn_fixed_effects <- function(input_model) {
  # Create empty data frame for summary of all models
  summary_df <- data.frame(
    model_name = character(),
    terms_count = integer(), # Number of terms
    terms = character(),
    aic = numeric(),
    model = I(list()) # Store model as a list column
  )
  # Create an empty list to store summary information
  model_summaries <- list()
  
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


# Function to fit and select random effects models ----
fxn_random_effects <- function(best_fixed_model) {
  random_terms <- c(
    "(1 | plot_type)", "(1 | f_year)", "(1 | f_year / grazer)", "(1 | grazer)",
    "(1 | f_break)", "(1 | f_one_yr)", "(1 | f_two_yr)", "(1 | f_new)",
    "(1 + f_year | plot_type)", "(1 + treatment | f_year)", "(1 + treatment | f_year / grazer)",
    "(1 + treatment | grazer)", "(1 + treatment | f_break)", "(1 + treatment | f_one_yr)",
    "(1 + treatment | f_two_yr)", "(1 + treatment | f_new)"
  )
  
  models <- lapply(random_terms, function(term) update(best_fixed_model, as.formula(paste(". ~ . +", term))))
  model.sel(best_fixed_model, models)
}

# Function to fit final model ----
fit_final_model <- function(data, response) {
  lmer(
    as.formula(paste(response, "~ treatment + (1 + treatment | plot_name) + f_year + plot_type")), 
    data = data, REML = FALSE
  )
}


# ========================================================== -----

# Function to fit and assess random effects  ----
best_fixed_model = best_model_fixed_abun_nat
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
  
  not_singular <-  bind_datalist %>%
    filter(is_singular == FALSE) %>%
    filter(model_name != "best_fixed_model") 
  
  if(nrow(not_singular) == 0){
    print("All random effect models are singular")
  }else{
    print(not_singular)
  }
  
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
    delta_AICc = all_models$delta,
    weight = round(all_models$weight, 4),
    df = all_models$df
  ) %>%
    arrange(AICc) %>%
    left_join(bind_datalist, "model_name")
}



fxn_select_random_effects <- function(model, data) {
  # Fit models with random effects
  ran_init <- fxn_random(index_model = model)
  ran_fix <- fxn_random(index_model = model)
  
  # Fit the best random effects model (assumed same formula as input)
  model_ran <- update(model, REML = FALSE)
  
  # Assess model diagnostics
  testDispersion(model_ran)
  check_singularity(model_ran)
  testUniformity(model_ran)
  testOutliers(model_ran)
  
  # Residual analysis for included terms
  included_terms <- names(fixef(model_ran))
  lapply(included_terms, function(term) plotResiduals(model_ran, data[[term]]))
  
  return(model_ran)
}

# Function to finalize model selection and review residuals  ----
fxn_final_model_review <- function(model, data) {
  # Fit the final model (assumed same formula as input)
  final_model <- update(model, REML = FALSE)
  
  # Residual analysis for included terms
  included_terms <- names(fixef(final_model))
  lapply(included_terms, function(term) plotResiduals(final_model, data[[term]]))
  
  return(final_model)
}

# Function to fit final model ----
fit_final_model <- function(data, response) {
  lmer(as.formula(paste(response, "~ treatment + (1 + treatment | plot_name) + f_year + plot_type")), 
       data = data, REML = FALSE)
}
# --- Running the workflow ---- 
# Apply functions to native species ----
# Fit the models with transformations
model_init_abun_nat <- fxn_initial_model(abun_nat)

# lm: lm(abundance ~ treatment)
# lm_null: lm(abundance ~ 1)
# lmer_1: lmer(abundance ~ treatment + (1 | plot_name))
# lmer_2: lmer(abundance ~ treatment + (1 + treatment | plot_name))
# lmer_null_1: lmer(abundance ~ 1 + (1 | plot_name))
# lmer_null_2: lmer(abundance ~ 1 + (1 + treatment | plot_name))

# View the summary table
print(model_init_abun_nat$summary_table)

# Review model diagnostics
# Based on AIC the standardized values had a better fit
# However, model diagnostics for this response didn't look good
fxn_model_review(model_init_abun_nat$best_model$model, abun_nat)

# Next best AIC set was value_log, These model diagnostics looked better but not perfect
best_model_init_abun_nat_log_2 <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
best_model_init_abun_nat_log_1 <- lmer(
  value_log ~ treatment + (1 | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
fxn_model_review(best_model_init_abun_nat_log_1, abun_nat)
fxn_model_review(best_model_init_abun_nat_log_2, abun_nat)

# Also looked at square root transformed diagnostics; super skewed

best_model_init_abun_nat_sqrt_2 <- lmer(
  value_sqrt ~ treatment + (1 + treatment | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
best_model_init_abun_nat_sqrt_1 <- lmer(
  value_sqrt ~ treatment + (1 | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
fxn_model_review(best_model_init_abun_nat_sqrt_1, abun_nat)
fxn_model_review(best_model_init_abun_nat_sqrt_2, abun_nat)

# After comparing model diagnostics for the three transformations I used log-transformed values for subsequent model selection
# Get the best overall model
best_model_init_abun_nat <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)




# Fit and select fixed effects model
best_fixed_abun_nat_dredge <- fxn_fixed_effects("best_model_abun_nat", method = "dredge")  # Choose "one" or "two" as needed
best_fixed_abun_nat_one <- fxn_fixed_effects("best_model_abun_nat", method = "one")  # Choose "one" or "two" as needed
best_fixed_abun_nat_two <- fxn_fixed_effects("best_model_abun_nat", method = "two")  # Choose "one" or "two" as needed

best_fixed_abun_nat

# Review model diagnostics
fxn_model_review(best_fixed_abun_nat$mixed_model_selection, abun_nat)

# Fit and select random effects model
best_random_abun_nat <- fxn_random_effects(best_fixed_abun_nat)
# Review model diagnostics
fxn_model_review(best_random_abun_nat$mixed_model_selection, abun_nat)

# Fit final model
final_model_abun_nat <- fit_final_model(abun_nat, response = "abun_nat")
# Review model diagnostics
fxn_model_review(final_model_abun_nat$mixed_model_selection, abun_nat)



nat_models <- define_models(abun_nat, "value_log")

# Initial model selection
nat_selection <- select_best_model(nat_models)
best_model_nat <- nat_models$lmer_2   
nat_assessment <- assess_model(best_model_nat, abun_nat)

model_abun_nat_fix <- fxn_select_fixed_effects(index_model = "model_abun_nat_init", data = abun)
model_abun_nat_ran <- fxn_select_random_effects(model = model_abun_nat_fix, data = abun)
model_abun_nat_final <- fxn_final_model_review(model = model_abun_nat_ran, data = abun)

final_model_nat <- fit_final_model(data_nat, "value_log")


# Apply functions to native forb species ----
frb_models <- define_models(abun_frb, "value_log")
frb_selection <- select_best_model(frb_models)
best_model_frb <- frb_models$lmer_2   
frb_assessment <- assess_model(best_model_frb, data_frb)
final_model_frb <- fit_final_model(data_frb, "value_log")
 
# Apply functions to non-native species ----
non_models <- define_models(abun_non, "value_log")
non_selection <- select_best_model(non_models)
best_model_non <- non_models$lmer_2   
non_assessment <- assess_model(best_model_non, data_non)
final_model_non <- fit_final_model(data_non, "value_log")
