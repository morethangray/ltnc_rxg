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
rich <- rich_rich$richdance
rich_nat <- rich$rich_nat
rich_frb <- rich$rich_frb
rich_non <- rich$rich_non

# ========================================================== -----
# Native species richness ----
#   Fit the initial models ----
# View the summary table

# Review model diagnostics 

# Identify the model with the lowest AIC for each transformation

# Review alternate model diagnostics 

# Define the best initial model

#   Fit the fixed effects models ----

# View the summary table

# View the best model 

# Review the model diagnostics 

#   Evaluate random effects ----
# Define the best fixed effects model

# View the summary table

#   Define the final model ----
# ---------------------------------------------------------- -----
# Native forb species richness ----
#   Fit the initial models ----
# View the summary table

# Review model diagnostics 

# Identify the model with the lowest AIC for each transformation

# Review alternate model diagnostics 

# Define the best initial model

#   Fit the fixed effects models ----

# View the summary table

# View the best model 

# Review the model diagnostics 

#   Evaluate random effects ----
# Define the best fixed effects model

# View the summary table

#   Define the final model ----
# ---------------------------------------------------------- -----
# Non-native species richness ----
#   Fit the initial models ----
# View the summary table

# Review model diagnostics 

# Identify the model with the lowest AIC for each transformation

# Review alternate model diagnostics 

# Define the best initial model

#   Fit the fixed effects models ----

# View the summary table

# View the best model 

# Review the model diagnostics 

#   Evaluate random effects ----
# Define the best fixed effects model

# View the summary table

#   Define the final model ----
# ---------------------------------------------------------- -----
