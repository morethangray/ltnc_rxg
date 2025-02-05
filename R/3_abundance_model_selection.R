# Run setup scripts -----
source(here::here("R/1_setup.R"))

# Define functions for model selection ----
# Fit the initial models
source(here::here("R/functions/fxn_lmer_initial_model.R"))

# View the summary table
source(here::here("R/functions/fxn_kable.R"))

# Review model diagnostics 
source(here::here("R/functions/fxn_lmer_model_review.R"))

# Evaluate fixed effects
source(here::here("R/functions/fxn_lmer_fixed_effects.R"))

# Evaluate random effects
source(here::here("R/functions/fxn_lmer_random_effects.R"))

# Create input data subsets -----
abun <- rich_abun$abundance
abun_nat <- abun$abun_nat
abun_frb <- abun$abun_frb
abun_non <- abun$abun_non

# ========================================================== -----
# Native species abundance ----
#   Fit the initial models ----
models_init_abun_nat <- fxn_lmer_initial_model(abun_nat)

# View the summary table
models_init_abun_nat$summary_table

# Review model diagnostics 
fxn_lmer_model_review(models_init_abun_nat$best_model$model, abun_nat)

# Identify the model with the lowest AIC for each transformation
models_init_abun_nat$summary_table %>%
  group_by(response) %>%
  slice(1) %>%
  arrange(aic) %>%
  ungroup()

# Review alternate model diagnostics 
model_init_abun_nat_log <- lmer(value_log ~ treatment + (1 + treatment | plot_name),
  data = abun_nat, REML = FALSE)

fxn_lmer_model_review(model_init_abun_nat_log, abun_nat)

# Define the best initial model
best_model_init_abun_nat <- lmer(value_log ~ treatment + (1 + treatment | plot_name),
                                 data = abun_nat, REML = FALSE)

#   Fit the fixed effects models ----
models_fixed_abun_nat <- fxn_lmer_fixed_effects(best_model_init_abun_nat)

# View the summary table
fixed_summary_table_abun_nat <- models_fixed_abun_nat$summary_table %>%
  as_tibble() %>%
  dplyr::select(model_name, terms_count, aic) %>%
  dplyr::mutate(delta = round(aic - min(aic), 1))

fixed_summary_table_abun_nat 

# View the best model 
best_fixed_model_abun_nat <- models_fixed_abun_nat$best_model_abun_nat$formula

# Review the model diagnostics 
model_fixed_abun_nat_5 <- update(best_model_init_abun_nat, . ~ . + f_year)
fxn_lmer_model_review(model_fixed_abun_nat_5, abun_nat)

#   Evaluate random effects ----
# Define the best fixed effects model
best_model_fixed_abun_nat <- update(best_model_init_abun_nat, . ~ . + f_year)
models_random_abun_nat <- fxn_lmer_random_effects(best_model_fixed_abun_nat)

models_random_abun_nat

#   Define the final model ----
mod_abun_nat <- 
  lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), 
             data = abun_nat, REML = FALSE)

# ---------------------------------------------------------- -----
# Native forb species abundance ----
#   Fit the initial models ----
models_init_abun_frb <- fxn_lmer_initial_model(abun_frb)
# View the summary table
models_init_abun_frb$summary_table %>%
  fxn_kable()
# Review model diagnostics 
fxn_lmer_model_review(models_init_abun_frb$best_model$model, abun_frb)
# Identify the model with the lowest AIC for each transformation
models_init_abun_frb$summary_table %>%
  group_by(response) %>%
  slice(1) %>%
  arrange(aic) %>%
  ungroup()
# Review alternate model diagnostics 
model_init_abun_frb_log <- lmer(value_log ~ treatment + (1 + treatment | plot_name),
                                data = abun_frb, REML = FALSE)
fxn_lmer_model_review(model_init_abun_frb_log, abun_frb)
# Define the best initial model
best_model_init_abun_frb <- lmer(value_log ~ treatment + (1 + treatment | plot_name),
                                 data = abun_frb, REML = FALSE)

#   Fit the fixed effects models ----
models_fixed_abun_frb <- fxn_lmer_fixed_effects(best_model_init_abun_frb)
# View the summary table
fixed_summary_table_abun_frb <- models_fixed_abun_frb$summary_table %>%
  as_tibble() %>%
  dplyr::select(model_name, terms_count, aic) %>%
  dplyr::mutate(delta = round(aic - min(aic), 1))
fixed_summary_table_abun_frb %>%
  fxn_kable()
# View the best model 
best_fixed_model_abun_frb <- models_fixed_abun_frb$best_model_abun_frb$formula

# Review the model diagnostics 
model_fixed_abun_frb_5 <- update(best_model_init_abun_frb, . ~ . + f_year)
fxn_lmer_model_review(model_fixed_abun_frb_5, abun_frb)

#   Evaluate random effects ----
# Define the best fixed effects model
best_model_fixed_abun_frb <- update(best_model_init_abun_frb, . ~ . + f_year)
models_random_abun_frb <- fxn_lmer_random_effects(best_model_fixed_abun_frb)
models_random_abun_frb %>%
  fxn_kable()
#   Define the final model ----
mod_abun_frb <- 
  lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), 
             data = abun_frb, REML = FALSE)

# ---------------------------------------------------------- -----
# Non-native species abundance ----
#   Fit the initial models ----
models_init_abun_non <- fxn_lmer_initial_model(abun_non)
# View the summary table
models_init_abun_non$summary_table %>%
  fxn_kable()
# Review model diagnostics 
fxn_lmer_model_review(models_init_abun_non$best_model$model, abun_non)
# Identify the model with the lowest AIC for each transformation
models_init_abun_non$summary_table %>%
  group_by(response) %>%
  slice(1) %>%
  arrange(aic) %>%
  ungroup()
# Review alternate model diagnostics 
model_init_abun_non_log <- lmer(value_log ~ treatment + (1 + treatment | plot_name),
                                data = abun_non, REML = FALSE)
fxn_lmer_model_review(model_init_abun_non_log, abun_non)
# Define the best initial model
best_model_init_abun_non <- lmer(value_log ~ treatment + (1 + treatment | plot_name),
                                 data = abun_non, REML = FALSE)

#   Fit the fixed effects models ----
models_fixed_abun_non <- fxn_lmer_fixed_effects(best_model_init_abun_non)
# View the summary table
fixed_summary_table_abun_non <- models_fixed_abun_non$summary_table %>%
  as_tibble() %>%
  dplyr::select(model_name, terms_count, aic) %>%
  dplyr::mutate(delta = round(aic - min(aic), 1))
fixed_summary_table_abun_non %>%
  fxn_kable()
# View the best model 
best_fixed_model_abun_non <- models_fixed_abun_non$best_model_abun_non$formula

# Review the model diagnostics 
model_fixed_abun_non_5 <- update(best_model_init_abun_non, . ~ . + f_year)
fxn_lmer_model_review(model_fixed_abun_non_5, abun_non)

#   Evaluate random effects ----
# Define the best fixed effects model
best_model_fixed_abun_non <- update(best_model_init_abun_non, . ~ . + f_year)
models_random_abun_non <- fxn_lmer_random_effects(best_model_fixed_abun_non)
models_random_abun_non %>%
  fxn_kable()
#   Define the final model ----
mod_abun_non <- 
  lme4::lmer(value_sqrt ~ treatment + f_year + f_two_yr + (1 | plot_name), 
             data = abun_non, REML = FALSE)


# ---------------------------------------------------------- -----