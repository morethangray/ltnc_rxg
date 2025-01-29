# LOAD LIBRARIES ----
# For richness models
library(glmmTMB)
library(broom.mixed)
library(car)
library(MASS)
#
source(here("R/functions/fxn_model-selection.R"))
# ========================================================== -----
# WRITE FUNCTIONS ----

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
# ========================================================== -----