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
# Create functions to summarize the models
source(here::here("R/functions/fxn_backtransform.R"))

# ========================================================== -----
# Abundance ----
#   Create input data subsets & fit models ----
abun <- rich_abun$abundance  
abun_nat <- abun$abun_nat %>%
  mutate(treatment = factor(treatment, levels = c("Ungrazed", "Grazed")))
abun_frb <- abun$abun_frb
abun_non <- abun$abun_non

levels(abun_nat$treatment)


mod_abun_nat <- lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 | plot_name), data = abun_nat, REML = FALSE)
mod_abun_frb <- lme4::lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), data = abun_frb, REML = FALSE)
mod_abun_non <- lme4::lmer(value_sqrt ~ treatment + f_year + f_two_yr + (1 | plot_name), data = abun_non, REML = FALSE)


list_models_abun <- c("mod_abun_nat", 
                      "mod_abun_frb", 
                      "mod_abun_non")

#   Marginal means ----
marginal_means_abun <- list_models_abun %>%
  purrr::map(~ fxn_summarize_marginal_means(.x, lookup_tables)) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(response, subset, treatment) 

readr::write_csv(marginal_means_abun,
                 here(project_paths$path_out_summary,
                        "marginal-means_abun.csv"))

#   Contrasts ----
# use marginaleffects_0.5.0  
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/marginaleffects/marginaleffects_0.10.0.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# sessionInfo()

# Retrieve the model object
index_model_name = "mod_abun_nat"
index_model <- get(index_model_name)

# Extract model metadata
model_formula <- as.character(formula(index_model))[3]
response <- if(stringr::str_detect(index_model_name, "abun")) "abundance" else "richness"
subset <- stringr::str_sub(index_model_name, -3)
transformation <- stringr::str_remove(as.character(formula(index_model))[2], "value_")

# Extract predictor variables
predictors <- attr(terms(index_model), "term.labels")

# Create pairwise comparison list
variable_list <- setNames(rep(list("pairwise"), length(predictors)), predictors)

# Compute and process contrasts
# contrasts <- 
  marginaleffects::comparisons(index_model, variables = variable_list) %>%
  broom::tidy() %>%
  # janitor::clean_names() %>%
  #   tibble() %>%
  # group_by(term, contrast) %>%
  # summarize(
  #   # estimate = mean(estimate), 
  #   estimate = mean(comparison), 
  #   p_value = min(p_value), 
  #   statistic = mean(statistic), 
  #   # s_value = mean(s_value), 
  #   conf_low = mean(conf_low), 
  #   conf_high = mean(conf_high)
  # ) %>%
  fxn_backtransform(index_value = "estimate", index_transform = transformation) %>%
    relocate(term, contrast, bt_estimate)
  
#   dplyr::rename(
#     abbr_term = term, 
#     abbr_contrast = contrast
#   ) %>%
#   dplyr::mutate(
#     response = response, 
#     abbr_subset = subset,
#     model_name = index_model_name, 
#     model_formula = model_formula
#   ) %>%
#   dplyr::left_join(lookup_tables$lookup_model_subset, "abbr_subset") %>%
#   dplyr::left_join(lookup_tables$lookup_model_term, "abbr_term") %>%
#   dplyr::left_join(lookup_tables$lookup_model_contrast, "abbr_contrast") %>%
#   dplyr::select(
#     response, 
#     subset,
#     term,
#     contrast,
#     starts_with("bt"), 
#     p_value, 
#     statistic,
#     estimate_raw = estimate, 
#     model_name, 
#     model_formula, 
#     starts_with("abbr_")
#   )
# 
# contrasts_abun <- list_models_abun %>%
#   purrr::map(~ fxn_summarize_contrasts(.x, lookup_tables)) %>%
#   dplyr::bind_rows() %>%
#   dplyr::arrange(response, subset, term, contrast)
# 
# readr::write_csv(contrasts_abun,
#                  here(project_paths$path_out_summary,
#                       "contrasts_abun.csv"))
# ---------------------------------------------------------- -----
# Richness ----
# Create input data subsets
rich <- rich_abun$richness
rich_nat <- rich$rich_nat
rich_frb <- rich$rich_frb
rich_non <- rich$rich_non

# Native 
# CLOSEST SO FAR
# mod_rich_nat <- glmmTMB::glmmTMB(value ~ treatment + f_year + plot_type + (1 + treatment | plot_name),
#                         data = rich_nat, family = nbinom2,
#                         control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
mod_rich_nat <- glmmTMB::glmmTMB(value ~ treatment + f_year,
                                 data = rich_nat, family = nbinom2,
                                 control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

mod_rich_nat <- lme4::lmer(value_sqrt ~ treatment + f_year + plot_type + (1 | plot_name),
                           data = rich_nat, REML = FALSE)

fxn_summarize_marginal_means("mod_rich_nat")  %>%
  mutate(bt_log = exp(bt_estimate)) %>%
  mutate(bt_sqrt = estimate_raw^2) %>%
  mutate(bt_sqrt1 = exp(estimate_raw)) %>%
  relocate(starts_with("bt"))

# Native forb
mod_rich_frb <- glmmTMB::glmmTMB(value ~ treatment + f_year + (1 + treatment | plot_name) + (1 | plot_type),
                                 data = rich_frb, family = nbinom2,
                                 control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

fxn_summarize_marginal_means("mod_rich_frb")  %>%
  mutate(bt_log = exp(bt_estimate)) %>%
  relocate(starts_with("bt"))


# mod_rich_frb <- glmmTMB::glmmTMB(value ~ treatment + f_year + plot_type + (1 + treatment | plot_name),
#                         data = rich_frb, family = nbinom2,
#                         control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Non-native
mod_rich_non <- lme4::lmer(value_sqrt ~ treatment + (1 + treatment | plot_name) + f_year,
                           data = rich_non, REML = FALSE)

list_models_rich <- c("mod_rich_nat", 
                      "mod_rich_frb", 
                      "mod_rich_non")

marginal_means_rich <- list_models_rich %>%
  purrr::map(~ fxn_summarize_marginal_means(.x, lookup_tables)) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(response, subset, treatment) 

readr::write_csv(marginal_means_rich,
                 here(project_paths$path_out_summary,
                      "marginal-means_rich.csv"))

contrasts_rich <- list_models_rich %>%
  purrr::map(~ fxn_summarize_contrasts(.x, lookup_tables)) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(response, subset, term, contrast)

readr::write_csv(contrasts_rich,
                 here(project_paths$path_out_summary,
                      "contrasts_rich.csv"))
# ========================================================== -----
# All models ----
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

# readr::write_csv(marginal_means, 
#                  here(project_paths$path_out_summary, 
#                         "marginal-means.csv"))

# Summarize contrasts ----
contrasts <- list_model_names %>%
  purrr::map(~ fxn_summarize_contrasts(.x, lookup_tables)) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(response, subset, term, contrast)

# readr::write_csv(contrasts, 
#                  here(project_paths$path_out_summary, 
#                       "contrasts.csv"))

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
