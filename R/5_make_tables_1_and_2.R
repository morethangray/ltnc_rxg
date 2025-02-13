# Run setup scripts ----
source(here::here("R/1_setup.R"))

# Load necessary functions and helpers ----
# Create functions to summarize the models
source(here::here("R/functions/fxn_summarize_models.R"))

# Create readable model formula 
source(here::here("R/functions/fxn_make_pretty_model_formula.R"))

# Create table 1: summary of model attributes 
source(here::here("R/functions/fxn_make_table_1.R"))

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

# ========================================================== -----
# List the models for iteration ----
list_model_names <- c("mod_abun_nat", 
                      "mod_abun_frb", 
                      "mod_abun_non", 
                      "mod_rich_nat", 
                      "mod_rich_frb", 
                      "mod_rich_non")
# ========================================================== -----
# Make Table 1: Summary of model attributes ----
# Response variable, transformation, formula, distribution
table_1 <- list_model_names %>%
  purrr::map_dfr(fxn_make_table_1)

readr::write_csv(table_1,
                 here(project_paths$path_out,
                      "table_1.csv"))

# Make Table 2: Summary of marginal means and effect size ----
# Response variable, marginal mean (ungrazed, grazed, % change), effect size (effect, p-value)


table_2_marginal_means <- list_model_names %>%
  purrr::map(~ fxn_summarize_marginal_means(.x, lookup_tables)) %>%
  dplyr::bind_rows()  %>%
  dplyr::mutate(response_variable = glue::glue("{subset} {response}")) %>%
  dplyr::select(response_variable, treatment, bt_estimate) %>%
  tidyr::spread(treatment, bt_estimate) %>%
  dplyr::mutate(
    Ungrazed = round(Ungrazed, 2), 
    Grazed = round(Grazed, 2), 
    percent_change = round(((Grazed - Ungrazed)/Ungrazed)*100, 0))

table_2_effects <- list_model_names %>%
  purrr::map(~ fxn_summarize_contrasts(.x, lookup_tables)) %>%
  dplyr::bind_rows()%>% 
  dplyr::filter(term == "Grazing treatment") %>%
  dplyr::mutate(response_variable = glue::glue("{subset} {response}")) %>%
  
  dplyr::rename(effect_size = bt_estimate) %>%
  dplyr::mutate(
    effect_size = round(effect_size, 2), 
    signif = 
      dplyr::case_when(
        p_value < 0.05 & p_value > 0.01 ~ "<0.05", 
        p_value < 0.01 & p_value > 0.001 ~ "<0.01", 
        p_value < 0.001  ~ "<0.001", 
        TRUE ~ "n.s."
        )
  ) %>%
  dplyr::select(response_variable, effect_size, signif) 

table_2 <- table_2_marginal_means %>%
  left_join(table_2_effects, "response_variable")

readr::write_csv(table_2,
                 here(project_paths$path_out,
                      "table_2.csv"))
