# LOAD LIBRARIES ----
library(lme4) ## Model fitting
library(MuMIn) ## Model selection and multi-model inference 
options(na.action = "na.fail")
library(performance)
library(DHARMa) ## Model diagnostics and residual checks 

# For richness models
library(glmmTMB)
library(broom.mixed)
library(car)
library(MASS)
#
# ========================================================== -----
# WRITE FUNCTIONS ----

#   fxn_fixed_dredge_by_model ----
fxn_fixed_dredge_by_model <- function(index_model) {
  input_model <- get(index_model)
  
  f0 <- update(input_model, . ~ . +
                 plot_type +
                 f_year +
                 # grazer +
                 f_break +
                 f_new +
                 f_one_yr +
                 f_two_yr)
  
  dd <- dredge(f0)
  
  aic_fixed <-
    subset(dd, delta < 6) %>%
    as_tibble(rownames = "model_name") %>%
    mutate(
      rank = 1:n(),
      input_model = index_model
    ) %>%
    clean_names() %>%
    dplyr::relocate(input_model,
                    rank,
                    treatment,
                    f_one_yr,
                    f_two_yr,
                    f_break,
                    f_new,
                    plot_type,
                    f_year,
                    # grazer,
                    delta,
                    weight,
                    aicc = ai_cc,
                    df,
                    model_name
    )
}

#   fxn_fixed_one ----
# index_subset <- "n_val_i"
fxn_fixed_one <- function(index_subset) {
  input_model <- get(paste0(index_subset, "_red"))
  
  f01 <- update(input_model, . ~ . + plot_type)
  # f02 <- update(input_model, . ~ . + grazer)
  f03 <- update(input_model, . ~ . + f_break)
  f04 <- update(input_model, . ~ . + f_new)
  f05 <- update(input_model, . ~ . + f_one_yr)
  f06 <- update(input_model, . ~ . + f_two_yr)
  f07 <- update(input_model, . ~ . + f_year)
  
  model.sel(
    f01,
    # f02,
    f03,
    f04,
    f05, f06, f07
  ) %>%
    as_tibble(rownames = "model_name") %>%
    clean_names() %>%
    mutate(
      input_model = index_subset,
      rank = 1:n()
    ) %>%
    dplyr::relocate(input_model,
                    rank,
                    treatment,
                    f_one_yr,
                    f_two_yr,
                    f_break,
                    f_new,
                    plot_type,
                    f_year,
                    # grazer,
                    delta,
                    weight,
                    aicc = ai_cc,
                    df,
                    model_name
    )
}
#   fxn_fixed_two ----
fxn_fixed_two <- function(index_subset) {
  input_model <- get(paste0(index_subset, "_red"))
  
  f08 <- update(input_model, . ~ . + f_year + plot_type)
  # f09 <- update(input_model, . ~ . + f_year + grazer)
  f10 <- update(input_model, . ~ . + f_year + f_break)
  f11 <- update(input_model, . ~ . + f_year + f_new)
  f12 <- update(input_model, . ~ . + f_year + f_one_yr)
  f13 <- update(input_model, . ~ . + f_year + f_two_yr)
  
  model.sel(
    f08,
    # f09,
    f10, f11, f12, f13
  ) %>%
    as_tibble(rownames = "model_name") %>%
    clean_names() %>%
    mutate(
      input_model = index_subset,
      rank = 1:n()
    ) %>%
    dplyr::relocate(input_model,
                    rank,
                    treatment,
                    f_one_yr,
                    f_two_yr,
                    f_break,
                    f_new,
                    plot_type,
                    f_year,
                    # grazer,
                    delta,
                    weight,
                    aicc = ai_cc,
                    df,
                    model_name
    )
}

#   fxn_fixed_one_glmmTMB ----
fxn_fixed_one_glmmTMB <- function(index_data) {
  f01 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      plot_type,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f02 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      grazer,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f03 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_break,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f04 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_new,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f05 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_one_yr,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f06 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_two_yr,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f07 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  model.sel(
    f01, f02, f03,
    f04,
    f05, f06, f07
  )
}
# Richness: Non-native ---


#   fxn_fixed_two_glmmTMB ----
fxn_fixed_two_glmmTMB <- function(index_data) {
  # Add a second fixed effect
  f08 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + plot_type,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f09 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + f_break,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f10 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + f_new,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  f11 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + f_one_yr,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f12 <- glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + f_two_yr,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  model.sel(f08, f09, f10, f11, f12)
}
#   fxn_reduced ----
fxn_reduced <- function(index_data) {
  g0 <- MASS::glm.nb(value ~ treatment, data = index_data)
  r0 <- glmmTMB(value ~ treatment + (1 | plot_name),
                data = index_data, family = nbinom2,
                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  r1 <- glmmTMB(value ~ treatment + (1 + treatment | plot_name),
                data = index_data, family = nbinom2,
                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  aic_table <- model.sel(g0, r0, r1)
}
#   fxn_random ----
fxn_random <- function(index_model) {
  best_fixed <- index_model
  
  # Add each random intercept one at a time
  m11 <- update(best_fixed, . ~ . + (1 | plot_type))
  m12 <- update(best_fixed, . ~ . + (1 | f_year))
  m13 <- update(best_fixed, . ~ . + (1 | f_year / grazer))
  m14 <- update(best_fixed, . ~ . + (1 | grazer))
  m15 <- update(best_fixed, . ~ . + (1 | f_break))
  m16 <- update(best_fixed, . ~ . + (1 | f_one_yr))
  m17 <- update(best_fixed, . ~ . + (1 | f_two_yr))
  m18 <- update(best_fixed, . ~ . + (1 | f_new))
  #
  # # Add each random slope one at a time
  m21 <- update(best_fixed, . ~ . + (1 + f_year | plot_type))
  m22 <- update(best_fixed, . ~ . + (1 + treatment | f_year))
  m23 <- update(best_fixed, . ~ . + (1 + treatment | f_year / grazer))
  m24 <- update(best_fixed, . ~ . + (1 + treatment | grazer))
  m25 <- update(best_fixed, . ~ . + (1 + treatment | f_break))
  m26 <- update(best_fixed, . ~ . + (1 + treatment | f_one_yr))
  m27 <- update(best_fixed, . ~ . + (1 + treatment | f_two_yr))
  m28 <- update(best_fixed, . ~ . + (1 + treatment | f_new))
  
  model.sel(
    best_fixed,
    m11,
    m12, m13, m14,
    m15, m16, m17, m18,
    m21, m22, m23, m24,
    m25, m26, m27, m28
  )
}

# ========================================================== -----