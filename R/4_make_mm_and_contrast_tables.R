# ---------------------------------------------------------- -----
#   !!!   VERY IMPORTANT NOTE   !!!   ----
# You must run the scripts in fxn_summarize_models with marginaleffects_0.5.0  
# If you use a later version you could (and probably will) get different results
#
# Get marginaleffects_0.5.0 from cran by running these 2 lines:
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/marginaleffects/marginaleffects_0.5.0.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
#
# Then confirm the package version:
# sessionInfo()
#
# ---------------------------------------------------------- -----
# Run setup scripts ----
source(here::here("R/1_setup.R"))

# Load necessary functions and helpers ----
# Create functions to summarize the models
source(here::here("R/functions/fxn_summarize_models.R"))
#
# ========================================================== -----

# Create input data subsets & fit models ----
#   Abundance ----
abun <- rich_abun$abundance  
abun_nat <- abun$abun_nat 
abun_frb <- abun$abun_frb
abun_non <- abun$abun_non

mod_abun_nat <- lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)
mod_abun_frb <- lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), data = abun_frb, REML = FALSE)
mod_abun_non <- lme4::lmer(value_sqrt ~ treatment + f_year + f_two_yr + (1 | plot_name), data = abun_non, REML = FALSE)
# ---------------------------------------------------------- -----
#   Richness ----
rich <- rich_abun$richness
rich_nat <- rich$rich_nat
rich_frb <- rich$rich_frb
rich_non <- rich$rich_non

mod_rich_nat <- glmmTMB::glmmTMB(value ~ treatment + f_year + (1 | plot_type) + (1 + treatment | plot_name),
                        data = rich_nat, family = nbinom2,
                        control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

mod_rich_frb <- glmmTMB::glmmTMB(value ~ treatment + f_year + (1 | plot_type) + (1 + treatment | plot_name),
                                 data = rich_frb, family = nbinom2,
                                 control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

mod_rich_non <- lme4::lmer(value_sqrt ~ treatment + f_year + (1 + treatment | plot_name), 
                           data = rich_non, REML = FALSE)


# List the models ----
list_model_names <- c("mod_abun_nat", 
                      "mod_abun_frb", 
                      "mod_abun_non", 
                      "mod_rich_nat", 
                      "mod_rich_frb", 
                      "mod_rich_non")
# Summarize marginal means  ----
marginal_means <- list_model_names %>%
  purrr::map(~ fxn_summarize_marginal_means(.x, lookup_tables)) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(response, subset, treatment)

readr::write_csv(marginal_means,
                 here(project_paths$path_out,
                        "marginal-means.csv"))

# Summarize contrasts ----
contrasts <- list_model_names %>%
  purrr::map(~ fxn_summarize_contrasts(.x, lookup_tables)) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(response, subset, term, contrast)

readr::write_csv(contrasts,
                 here(project_paths$path_out,
                      "contrasts.csv"))
