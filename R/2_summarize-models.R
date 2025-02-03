rm(list = ls())
# Set up for analysis ----
# Set seed for reproducible simulations during model selection
set.seed(912)
# Manage dependencies
source(here::here("R/functions/setup.R"))
# Define file paths
source(here::here("R/functions/fxn_setup_file_paths.R"))
# Create lookup tables
source(here::here("R/functions/fxn_setup_lookup_tables.R"))
# Prepare data from xlsx 
source(here::here("R/functions/fxn_prepare_input_data.R"))
# Create the processed data 
processed_data <- fxn_prepare_input_data(list_years = 2019:2022, 
                                         project_paths = project_paths, 
                                         lookup_tables = lookup_tables)
# Calculate metrics for richness and abundance
source(here::here("R/functions/fxn_calculate_metrics.R"))
# Load richness and abundance
source(here::here("R/functions/fxn_load_rich_abun.R"))
# Run function to load richness and abundance data
rich_abun <- fxn_load_rich_abun(project_paths = project_paths)

# Create input data subsets
abun <- rich_abun$abundance
abun_nat <- abun$abun_nat
abun_frb <- abun$abun_frb
abun_non <- abun$abun_non

rich <- rich_abun$richness
rich_nat <- rich$rich_nat
rich_frb <- rich$rich_frb
rich_non <- rich$rich_non

# sessionInfo()
# ========================================================== -----
# Define models ----
# List model names ----
list_model_names <- c("mod_abun_nat", 
                      "mod_abun_frb", 
                      "mod_abun_non", 
                      "mod_rich_nat", 
                      "mod_rich_frb", 
                      "mod_rich_non")

# Richness ----
# Native 
# mod_rich_nat <- glmmTMB::glmmTMB(value ~ treatment + f_year + (1 + treatment | plot_name),
#                                  data = rich_nat, family = nbinom2,
#                                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
# mod_rich_nat <- glmmTMB::glmmTMB(value ~ treatment + f_year + (1 + treatment | plot_name) + (1 | plot_type),
#                         data = rich_nat, family = nbinom2,
#                         control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

mod_rich_nat <- glmmTMB::glmmTMB(value ~ treatment + f_year + plot_type + (1 + treatment | plot_name),
                        data = rich_nat, family = nbinom2,
                        control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Native forb
# mod_rich_frb <- glmmTMB::glmmTMB(value ~ treatment + f_year + (1 + treatment | plot_name) + (1 | plot_type),
#                         data = rich_frb, family = nbinom2,
#                         control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

mod_rich_frb <- glmmTMB::glmmTMB(value ~ treatment + f_year + plot_type + (1 + treatment | plot_name),
                        data = rich_frb, family = nbinom2,
                        control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Non-native
mod_rich_non <- lme4::lmer(value_sqrt ~ treatment + (1 + treatment | plot_name) + f_year,
                 data = rich_non, REML = FALSE)
# Abundance ----
mod_abun_nat <- lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)
mod_abun_frb <- lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), data = abun_frb, REML = FALSE)
mod_abun_non <- lme4::lmer(value_sqrt ~ treatment + f_year + f_two_yr + (1 | plot_name), data = abun_non, REML = FALSE)

# ========================================================== -----
# DEFINE FUNCTIONS ----
#   fxn_backtransform ----
fxn_backtransform <- function(index_data, index_value, index_transform){
  
  bt_name <- glue::glue("bt_{index_value}")
  
 index_data %>%
    dplyr::mutate(
      "{bt_name}" :=  
        dplyr::case_when(
          index_transform == "sqrt" ~ (!!sym(index_value))^2, 
          index_transform == "log" ~ exp(!!sym(index_value) - 1e-6), 
          TRUE ~ !!sym(index_value)
        )
    )
  
}

# Summarize marginal means  ----

#   fxn_summarize_marginal_means ----
fxn_summarize_marginal_means <- function(index_model_name){
  
  # Retrieve the model object
  index_model <- get(index_model_name)
  
  # Extract model formula and attributes for reference
  model_formula <- as.character(formula(index_model))[3]
  response <- if(stringr::str_detect(index_model_name, "abun")){
    "abundance"
  }else{
    "richness"
  }
  subset <- stringr::str_sub(index_model_name, -3)
  
  # Extract response variable and determine transformation
  transformation <- stringr::str_remove(as.character(formula(index_model))[2], "value_")
  
  # Compute marginal means dynamically
  marginal_means <- 
    emmeans::emmeans(index_model, ~ treatment) %>%
    broom::tidy() %>%
    janitor::clean_names() %>%
    fxn_backtransform(index_value = "estimate", index_transform = transformation) %>%
    # Use the delta method for handling standard errors
    fxn_backtransform(index_value = "std_error", index_transform = transformation) %>%
    dplyr::mutate(bt_std_error = std_error * bt_std_error, 
                  response = response, 
                  subset = subset,
                  model_name = index_model_name, 
                  model_formula = model_formula) %>%
    dplyr::select(response, 
                  treatment, 
                  starts_with("bt"), 
                  p_value, 
                  statistic,
                  df, 
                  subset, 
                  model_name, 
                  model_formula) 
  
  return(marginal_means)

}
#
# Apply function to the models ----
bind_results_mm <- list_model_names %>%
  purrr::map(fxn_summarize_marginal_means) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(lookup_tables$lookup_model_subset, "subset") %>%
  dplyr::rename(abbr_subset = subset, 
         subset = lab_subset) %>%
  dplyr::relocate(response, subset) %>%
  dplyr::relocate(starts_with("abbr_"), .after = last_col())  %>%
  dplyr::arrange(response, subset, treatment)

readr::write_csv(bind_results_mm, 
                 here(project_paths$path_out_summary, 
                        "marginal-means.csv"))

# Summarize contrasts ----
#   fxn_summarize_contrasts ----
# use marginaleffects_0.5.0  
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/marginaleffects/marginaleffects_0.5.0.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# sessionInfo()

fxn_summarize_contrasts <- function(index_model_name) {
  
  # Retrieve the model object
  index_model <- get(index_model_name)
  
  # Extract model formula and attributes for reference
  model_formula <- as.character(formula(index_model))[3]
  response <- if(stringr::str_detect(index_model_name, "abun")){
    "abundance"
  }else{
    "richness"
  }
  subset <- stringr::str_sub(index_model_name, -3)
  
  # Extract response variable and determine transformation
  transformation <- stringr::str_remove(as.character(formula(index_model))[2], "value_")
  
  # Extract predictor variables from model formula
  predictors <- attr(terms(index_model), "term.labels")
  
  # Create a named list for pairwise comparisons
  variable_list <- setNames(rep(list("pairwise"), length(predictors)), predictors)
  
  # Compute contrasts dynamically
  contrasts <- 
    marginaleffects::comparisons(index_model, variables = variable_list)  %>%
    broom::tidy() %>%
    janitor::clean_names() %>%
    group_by(term, contrast) %>%
    summarize(estimate = mean(estimate), 
              p_value = min(p_value), 
              statistic = mean(statistic), 
              s_value = mean(s_value), 
              conf_low = mean(conf_low), 
              conf_high = mean(conf_high)) %>%
      fxn_backtransform(index_value = "estimate", index_transform = transformation) %>%
      fxn_backtransform(index_value = "conf_low", index_transform = transformation) %>%
      fxn_backtransform(index_value = "conf_high", index_transform = transformation) %>%
    dplyr::mutate(response = response, 
                  subset = subset,
                  model_name = index_model_name, 
                  model_formula = model_formula) %>%
    dplyr::select(response, 
                  term,
                  contrast,
                  starts_with("bt"), 
                  p_value, 
                  statistic,
                  subset, 
                  model_name, 
                  model_formula) 
  
  return(contrasts)
}
# Apply function to the models ----

 <- list_model_names %>%
  purrr::map(fxn_summarize_contrasts) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(lookup_tables$lookup_model_subset, "subset") %>%
  dplyr::left_join(lookup_tables$lookup_model_subset, "term") %>%
  dplyr::left_join(lookup_tables$lookup_model_subset, "contrast") %>%
  dplyr::rename(abbr_subset = subset, 
                abbr_term = term, 
                abbr_contrast = contrast, 
                subset = lab_subset, 
                term = lab_term, 
                contrast = lab_contrast) %>%
  dplyr::relocate(response, subset, term, contrast) %>%
  dplyr::relocate(starts_with("abbr_"), .after = last_col())  %>%
  dplyr::arrange(response, subset, treatment)

readr::write_csv(bind_results_contrasts, 
                 here(project_paths$path_out_summary, 
                      "contrasts.csv"))

# ========================================================== -----

# rich_non: Create model summaries ----
mod_rich_non_sqrt <- lmer(value_sqrt ~ treatment + f_year + (1 + treatment | plot_name), 
                     data = rich_non, REML = FALSE)
mod_rich_non_sqrt_1 <- lmer(value_sqrt ~ treatment + f_year + (1 | plot_name), 
                          data = rich_non, REML = FALSE)
# mod_rich_non_sqrt_2 <- lmer(value_sqrt ~ treatment + f_year, 
#                             data = rich_non, REML = FALSE)

# mod_rich_non_log_1 <- lmer(value_log ~ treatment + f_year, 
#                 data = r_non, REML = FALSE)

mod_rich_non_log_2 <- lmer(value_log ~ treatment + f_year + (1 | plot_name), 
                       data = rich_non, REML = FALSE)

mod_rich_non_log_3 <- lmer(value_log ~ treatment + f_year + f_break + (1 | plot_name), 
                       data = rich_non, REML = FALSE)


mod_rich_non_log_4 <- lmer(value_log ~ treatment + f_year + (1 + treatment | plot_name),
                       data = rich_non, REML = FALSE)



mod_rich_non_log_5 <- lmer(value_log ~ treatment + f_year + f_two_yr + (1 + treatment | plot_name),
                  data = rich_non)


# 
model.sel(mod_rich_non_sqrt,
          mod_rich_non_sqrt_1,
          mod_rich_non_log_2, 
          mod_rich_non_log_3, 
          mod_rich_non_log_4, 
          mod_rich_non_log_5) %>%
  arrange(AICc)
# 
# summary(mod_rich_non)
# performance::performance(mod_rich_non_log_4)


# List effects for marginal summaries 
list_effects_rich_non_x <- c("treatment", "f_year")

mfx <- 
  marginaleffects::slopes(mod_rich_non_log_4, 
                          re.form = NULL) %>%
  tidy() %>%
  clean_names() %>%
  mutate(method = "effect")  %>%
  rename(value = contrast) %>%
  relocate(term, value, starts_with("p_"))

mm <- 
  marginalmeans(mod_rich_non_log_4, 
                variables = all_of(list_effects_rich_non_x)) %>%
  tidy() %>%
  mutate(method = "mean")  

bind_summaries <- 
  bind_rows(mfx, 
            mm) %>%
  dplyr::select(method, term, value, estimate, conf.low, conf.high, p.value)  %>%
  write_csv(here(path_out, "richness_nonnative_glmm-summary-tables.csv"))

# predictions(mod_rich_non_log_4, 
#             re.form = NULL) %>%
#   as_tibble() %>%
#   mutate(treatment = as.character(treatment)) %>%
#   arrange(treatment) %>%
#   mutate(treatment = as_factor(treatment)) %>%
#   ggplot(aes(x = treatment, 
#              y = estimate, 
#              color = treatment)) + 
#   geom_quasirandom(aes(color = treatment),
#                    # size = 1,
#                    cex = 1,
#                    alpha = 0.8, 
#                    # dodge.width=0.5, 
#                    method = "pseudorandom") +
#   # scale_color_manual(values = palette_treatment) + 
#   ylab("Richness") +
#   theme_minimal() + 
#   theme(legend.position = "none",
#         axis.title.x = element_blank())  
# 
# ggsave(filename =
#          here(path_out, "richness_nonnative_best-model.png"),
#        width = 3, height = 5)
