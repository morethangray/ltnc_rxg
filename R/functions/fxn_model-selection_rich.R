source(here::here("R/functions/fxn_install_dependencies.R"))
# ========================================================== -----
# WRITE FUNCTIONS ----

#   fxn_reduced_glmmTMB ----
fxn_reduced_glmmTMB <- function(index_data){
  
  g0 <- MASS::glm.nb(value ~ treatment, data = index_data)
  r0 <- glmmTMB::glmmTMB(value ~ treatment + (1 | plot_name), data = index_data, family = nbinom2, 
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
  r1 <- glmmTMB::glmmTMB(value ~ treatment + (1 + treatment | plot_name), data = index_data, family = nbinom2, 
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
  
  aic_table <- model.sel(g0, r0, r1)
  
}
#   fxn_fixed_one_glmmTMB ----
fxn_fixed_one_glmmTMB <- function(index_data) {
  plot_type <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      plot_type,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  grazer <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      grazer,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f_break <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_break,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f_new <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_new,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f_one_yr <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_one_yr,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f_two_yr <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_two_yr,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  f_year <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  model.sel(
    plot_type, 
    grazer, 
    f_break, 
    f_new,
    f_one_yr,
    f_two_yr,
    f_year
  )
}



#   fxn_fixed_two_glmmTMB ----
fxn_fixed_two_glmmTMB <- function(index_data) {
  # Add a second fixed effect
  year_plot_type <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + plot_type,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  year_f_break <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + f_break,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  year_f_new <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + f_new,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  year_f_one_yr <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + f_one_yr,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  year_f_two_yr <- glmmTMB::glmmTMB(
    value ~ treatment + (1 + treatment | plot_name) +
      f_year + f_two_yr,
    data = index_data, family = nbinom2,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  model.sel(year_plot_type, 
            year_f_break, 
            year_f_new,
            year_f_one_yr,
            year_f_two_yr)
}
# ========================================================== -----