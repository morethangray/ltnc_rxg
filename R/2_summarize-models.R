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

abun <- rich_abun$abundance
abun_nat <- abun$abun_nat
abun_frb <- abun$abun_frb
abun_non <- abun$abun_non

rich <- rich_abun$richness
rich_nat <- rich$rich_nat
rich_frb <- rich$rich_frb
rich_non <- rich$rich_non
# ========================================================== -----
# DEFINE MODELS ----
# Richness ----
# Native 
# mod_rich_nat <- glmmTMB(value ~ treatment + f_year + (1 + treatment | plot_name) + (1 | plot_type),
#                         data = rich_nat, family = nbinom2,
#                         control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

mod_rich_nat <- glmmTMB(value ~ treatment + f_year + plot_type + (1 + treatment | plot_name),
                        data = rich_nat, family = nbinom2,
                        control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Native forb
mod_rich_frb <- glmmTMB(value ~ treatment + f_year + plot_type + (1 + treatment | plot_name),
                        data = rich_frb, family = nbinom2,
                        control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# mod_rich_frb <- glmmTMB(value ~ treatment + f_year + (1 + treatment | plot_name) + (1 | plot_type),
#                         data = rich_frb, family = nbinom2,
#                         control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Non-native
# mod_rich_non <- lmer(value_sqrt ~ treatment + (1 + treatment | plot_name) + f_year, 
#                  data = rich_non, REML = FALSE)
mod_rich_non  <- lmer(value_log ~ treatment + f_year + (1 + treatment | plot_name) ,
                     data = rich_non, REML = FALSE)

# mod_rich_non_final <- lmer(value_sqrt ~ treatment + (1 + treatment | plot_name) + f_year, 
#                            data = rich_non, REML = FALSE)
# 
# model.sel(mod_rich_non, 
#           mod_rich_non_1) %>%
#   arrange(AICc)
# 
# summary(mod_rich_non)
# performance::performance(mod_rich_non)

# List effects for marginal summaries 
list_effects_rich_nat <- c("treatment", "f_year", "plot_type")
list_effects_rich_frb <- list_effects_rich_nat
list_effects_rich_non <- c("treatment", "f_year")

# Abundance ----
# Native
mod_abun_nat <- lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)

# Native forb
mod_abun_frb <- lmer(value_log ~ treatment + f_year + plot_type + (1 + treatment | plot_name), data = abun_frb, REML = FALSE)

# Non-native
mod_abun_non <- lmer(value_sqrt ~ treatment + f_year + f_two_yr + (1 | plot_name), data = abun_non, REML = FALSE)

# List effects for marginal summaries 
list_effects_abun_nat <- c("treatment", "plot_type", "f_year")
list_effects_abun_frb <- list_effects_abun_nat
list_effects_abun_non <- c("treatment", "f_year", "f_two_yr")
# ========================================================== -----
# SUMMARIZE MODEL ----
# index_subset <- "abun_non"
# index_effects <- list_effects_abun_non
# index_bt <- "sqrt"

#   fxn_annotate_summary ----
fxn_annotate_summary <- function(df, index_subset){
  df %>%
    tidy() %>%
    janitor::clean_names() %>%
    mutate(subset = index_subset) %>%
    dplyr::relocate(subset, 
                    term,
                    estimate, 
                    conf_low, 
                    conf_high, 
                    starts_with("p_")) 
}
#   fxn_bt_log ----
fxn_bt_log <- function(df, index_bt){
  
     df %>%
      mutate(bt_estimate = exp(estimate), 
             bt_conf_low = exp(conf_low),
             bt_conf_high = exp(conf_high))
}
#   fxn_bt_log1 ----
fxn_bt_log1 <- function(df, index_bt){
  
  df %>%
    mutate(bt_estimate = exp(estimate) - 1e-6, 
           bt_conf_low = exp(conf_low) - 1e-6, 
           bt_conf_high = exp(conf_high) - - 1e-6)
}
#   fxn_bt_sqrt ----
fxn_bt_sqrt <- function(df, index_bt){
  
  df %>%
    mutate(bt_estimate = (estimate)^2, 
           bt_conf_low = (conf_low)^2,
           bt_conf_high = (conf_high)^2)
}
#   fxn_summarize_contrasts ----
fxn_summarize_contrasts <- function(index_subset, index_bt){
  
  best_fit <- get(glue::glue("mod_{index_subset}"))
  
  contrasts <- 
    marginaleffects::comparisons(
      best_fit, 
      variables = list(f_year = "pairwise",
                       f_two_yr = "pairwise",
                       plot_type = "pairwise",
                       treatment = "pairwise")) %>%
    fxn_annotate_summary(index_subset = index_subset) %>%
    dplyr::mutate(summary_type = "contrast")  
  
  if(index_bt == "none"){
    bt_data <- contrasts  %>%
      dplyr::mutate(bt_estimate = estimate,
                    bt_conf_low = conf_low,
                    bt_conf_high = conf_high)
  }
  
  if(index_bt == "log1"){
    bt_data <- fxn_bt_log1(contrasts)
  }
  
  if(index_bt == "log"){
    bt_data <- fxn_bt_log(contrasts)
  }
  
  if(index_bt == "sqrt"){
    bt_data <- fxn_bt_sqrt(contrasts)
  }
  
  bt_data %>%
    dplyr::relocate(subset, 
                    summary_type, 
                    term,
                    contrast, 
                    starts_with("bt"), 
                    starts_with("p_")) %>%
    readr::write_csv(here(project_paths$path_out_summary, 
                          glue::glue("{index_subset}_contrasts.csv")))
  
}

#   fxn_summarize_marginal_mean ----
fxn_summarize_marginal_mean <- function(index_subset, index_effects, index_bt){
  
  best_fit <- get(glue::glue("mod_{index_subset}"))

  mm <- 
    marginaleffects::marginalmeans(
      best_fit, 
      variables = all_of(index_effects)) %>%
    fxn_annotate_summary(index_subset = index_subset) %>%
    dplyr::mutate(summary_type = "marginal_mean")  
  
  if(index_bt == "none"){
    bt_data <- mm  %>%
      dplyr::mutate(bt_estimate = estimate,
                    bt_conf_low = conf_low,
                    bt_conf_high = conf_high)
  }
  
  if(index_bt == "log1"){
    bt_data <- fxn_bt_log1(mm)
  }
  
  if(index_bt == "log"){
    bt_data <- fxn_bt_log(mm)
  }
  
  if(index_bt == "sqrt"){
    bt_data <- fxn_bt_sqrt(mm)
  }
  
  bt_data %>%
    dplyr::relocate(subset, 
                    summary_type, 
                    term,
                    value, 
                    starts_with("bt"), 
                    starts_with("p_")) %>%
    readr::write_csv(here(project_paths$path_out_summary, 
                          glue::glue("{index_subset}_marginal-mean.csv")))
  
}


#   fxn_summarize_all ----
fxn_summarize_all <- function(index_subset, 
                          index_effects, 
                          index_bt){
  
  fxn_summarize_contrasts(index_subset, 
                          index_bt)
  
  fxn_summarize_marginal_mean(index_subset, 
                              index_effects, 
                              index_bt)
  
  # fxn_create_plots(index_subset, 
  #                  index_effects,
  #                  index_bt)
}

# Create model summaries ----
# fxn_summarize_contrasts(index_subset = "abun_non", index_bt = "sqrt")
# fxn_summarize_marginal_mean(index_subset = "abun_non",
#                   index_effects = list_effects_abun_non,
#                   index_bt = "sqrt")

fxn_summarize_all(index_subset = "rich_nat",
                  index_effects = list_effects_rich_nat,
                  index_bt = "none")

fxn_summarize_all(index_subset = "rich_frb",
                  index_effects = list_effects_rich_frb,
                  index_bt = "none")

fxn_summarize_all(index_subset = "rich_non",
                  index_effects = list_effects_rich_non,
                  index_bt = "log1")

fxn_summarize_all(index_subset = "abun_nat",
                  index_effects = list_effects_abun_nat,
                  index_bt = "log1")

fxn_summarize_all(index_subset = "abun_frb",
                  index_effects = list_effects_abun_frb,
                  index_bt = "log1")
#
fxn_summarize_all(index_subset = "abun_non",
                  index_effects = list_effects_abun_non,
                  index_bt = "sqrt")

# Combine model summaries ----
list_csv_contrast <- 
  list.files(here(project_paths$path_out_summary),
             pattern = "*contrasts.csv")

list_csv_mm <- 
  list.files(here(project_paths$path_out_summary),
             pattern = "*marginal-mean.csv")

bind_contrasts <- 
  sapply(here(project_paths$path_out_summary,
              list_csv_contrast),
         readr::read_csv, 
         show_col_types = FALSE,
         simplify = FALSE) %>%
  bind_rows() %>%
  rename(met_sub = subset) %>%
  mutate(metric = stringr::str_sub(met_sub, 1, 4), 
         subset = stringr::str_sub(met_sub, 6, 8)) %>%
  left_join(lookup_tables$lookup_model_subset, "subset") %>%
  left_join(lookup_tables$lookup_model_term, "term") %>%
  left_join(lookup_tables$lookup_model_contrast, "contrast")  %>%
  dplyr::select(metric, 
                starts_with("lab"),
                contrast, 
                term, 
                treatment, 
                starts_with("f_"), 
                starts_with("plot"),
                bt_estimate, 
                starts_with("bt_"), 
                p_value, 
                s_value, 
                predicted,
                predicted_lo, 
                predicted_hi, 
                raw_estimate = estimate, 
                raw_conf_low = conf_low, 
                raw_conf_high = conf_high, 
                raw_std_error = std_error, 
                raw_statistic = statistic, 
                value,
                starts_with("value"),
                summary_type) %>%
  distinct() %>%
  arrange(metric, lab_subset, lab_term) %>%
  readr::write_csv(here(project_paths$path_out_summary, 
                        "contrasts.csv"))

bind_mm <- 
  sapply(here(project_paths$path_out_summary,
              list_csv_mm),
         readr::read_csv, 
         show_col_types = FALSE,
         simplify = FALSE) %>%
  bind_rows() %>%
  rename(met_sub = subset) %>%
  mutate(metric = stringr::str_sub(met_sub, 1, 4), 
         subset = stringr::str_sub(met_sub, 6, 8)) %>%
  left_join(lookup_tables$lookup_model_subset, "subset") %>%
  left_join(lookup_tables$lookup_model_term, "term") %>%
  left_join(lookup_tables$lookup_model_mean, "value") %>%
  dplyr::select(metric, 
                starts_with("lab"),
                bt_estimate, 
                starts_with("bt_"), 
                p_value, 
                raw_estimate = estimate, 
                raw_conf_low = conf_low, 
                raw_conf_high = conf_high, 
                raw_std_error = std_error, 
                raw_statistic = statistic, 
                s_value, 
                term, 
                term_id = value, 
                summary_type)  %>%
  arrange(metric, lab_subset, lab_term)  %>%
  readr::write_csv(here(project_paths$path_out_summary, 
                        "marginal-means.csv"))

# Review marginal means ----
# Summary table of back-transformed values 
# bind_mm_wide <-
  bind_mm  %>%
      dplyr::mutate(
        signif =
          case_when(
            p_value < 0.001 ~ 0.001,
            p_value > 0.001 & p_value < 0.05 ~ 0.05,
            p_value > 0.05 & p_value < 0.01 ~ 0.01,
            TRUE ~ NA
          )
      ) %>%
  dplyr::select(metric, 
                starts_with("lab"),
                bt_estimate) %>%
  
  filter(lab_term == "Grazing treatment") %>%
  tidyr::spread(lab_value, bt_estimate) %>%
  dplyr::select(metric, lab_subset, Ungrazed, Grazed) 
# ========================================================== -----
# CREATE PLOTS ----
#   fxn_create_plots ----
# Define plot settings
library(ggbeeswarm)
source(here("R/functions/fxn_plot_settings.R"))
fxn_create_plots <- function(index_subset, index_effects, index_bt){
  
  best_fit <- get(paste0("mod_", index_subset))
  index_metric <- str_sub(index_subset, 1, 1)
  
  input_predictions <-
    predictions(best_fit) %>%
    as_tibble() %>%
    clean_names()
  
  mm <-
    marginalmeans(best_fit,
                  variables = all_of(index_effects)) %>%
    fxn_annotate_summary(index_subset = index_subset)
  
  # Backtransform values -----
  if(index_bt == "none"){
    predictions <- input_predictions %>%
      mutate(bt_predicted = predicted,
             bt_value = predicted)
    
    bt_mm <- mm  %>%
      mutate(bt_estimate = estimate,
             bt_conf_low = conf_low,
             bt_conf_high = conf_high)
  }
  if(index_bt == "log1"){
    predictions <-
      input_predictions %>%
      mutate(bt_predicted = exp(predicted) - 1e-06,
             bt_value = exp(value_log) - 1e-06)
    
    bt_mm <- fxn_bt_log1(mm)
  }
  if(index_bt == "log"){
    predictions <-
      input_predictions %>%
      mutate(bt_predicted = exp(predicted),
             bt_value = exp(value_log))
    
    bt_mm <- fxn_bt_log(mm)
  }
  if(index_bt == "sqrt"){
    predictions <-
      input_predictions %>%
      mutate(bt_predicted = (predicted)^2,
             bt_value = value_sqrt^2)
    
    bt_mm <- fxn_bt_sqrt(mm)
  }
  
  
  # Calculate mean by treatment ----
  treatment_mean <-
    bt_mm %>%
    filter(term %in% "treatment") %>%
    dplyr::select(treatment = value,
                  bt_mean = bt_estimate)
  
  # Create plot ----
  plot <-
    predictions %>%
    ggplot(aes(x = treatment,
               y = bt_predicted,
               color = treatment)) +
    geom_quasirandom(aes(color = treatment),
                     cex = 1,
                     alpha = 0.6,
                     method = "pseudorandom") +
    stat_summary(fun = mean,
                 fun.min = mean,
                 fun.max = mean,
                 geom = "crossbar",
                 width = 0.8,
                 size = 0.4,
                 color = "gray30") +
    scale_color_manual(values = palette_treatment_line_light) +
    ylab(index_subset) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank())
  
  if(index_metric == "r"){
    plot + ylim(0, 16)
  }else{
    plot + ylim(0, 130)
  }
  
  ggsave(filename =
           here(path_out_results, "plots",
                paste0(index_subset,
                       "_best-model_",
                       Sys.Date(),
                       ".png")),
         width = 2, height = 5)
}
#   fxn_plot_subset ----
fxn_plot_subset <- function(index_factor, index_data){
  
  # Summarize by treatment ----
  by_trmt <- 
    index_data %>%
    group_by(treatment) %>%
    summarize(trmt_mean = mean(value, na.rm = TRUE), 
              n = n(), 
              df = n-1,
              sd = sd(value, na.rm = TRUE),
              se = sd/sqrt(n), 
              t = qt(p = 0.025, 
                     df = df,
                     lower.tail = FALSE), 
              margin = t * se,
              lower = trmt_mean - margin, 
              upper = trmt_mean + margin) 
  
  # Summarize by treatment x factor ----
  df <-   
    index_data %>%
    arrange(year) %>%
    mutate(year = as_factor(year)) %>%
    select(treatment, 
           subset_name = all_of(index_factor), 
           value) 
  
  
  # Plot by treatment x factor ----
  df %>%
    group_by(treatment, subset_name) %>%
    summarize(subset_mean = mean(value, na.rm = TRUE), 
              n = n(), 
              df = n-1,
              sd = sd(value, na.rm = TRUE),
              se = sd/sqrt(n), 
              t = qt(p = 0.025, 
                     df = df,
                     lower.tail = FALSE), 
              margin = t * se,
              lower = subset_mean - margin, 
              upper = subset_mean + margin) %>%
    ungroup() %>%
    ggplot(aes(x = treatment,
               y = subset_mean)) + 
    geom_quasirandom(data = df, 
                     aes(x = treatment, 
                         y = value, 
                         color = subset_name),
                     size = 0.6,
                     alpha = 0.5, 
                     dodge.width=0.5, 
                     method = "pseudorandom") +
    geom_line(aes(group = subset_name, 
                  color = subset_name),
              position = position_dodge(width = 0.5), 
              color = "gray70") +
    geom_point(aes(group = subset_name, 
                   color = subset_name),
               size = 2.5,
               position = position_dodge(width = 0.5)) + 
    geom_point(data = by_trmt, 
               aes(x = treatment, 
                   y = trmt_mean), 
               color = "gray10",
               size = 4) + 
    geom_errorbar(data = by_trmt, 
                  aes(y = trmt_mean,
                      ymin = lower,
                      ymax = upper),
                  color = "gray10",
                  width = .05) +
    # ylab("Richness") +
    theme_minimal() +
    theme(axis.title.x = element_blank(), 
          legend.position = "none",
          axis.title.y = element_blank()) + 
    ylim(0, 20)
  # 
  # df %>%
  #   group_by(treatment, subset_name) %>%
  #   summarize(subset_mean = mean(value, na.rm = TRUE), 
  #             n = n(), 
  #             df = n-1,
  #             sd = sd(value, na.rm = TRUE),
  #             se = sd/sqrt(n), 
  #             t = qt(p = 0.025, 
  #                    df = df,
  #                    lower.tail = FALSE), 
  #             margin = t * se,
  #             lower = subset_mean - margin, 
  #             upper = subset_mean + margin) %>%
  #   ungroup() %>%
  #   ggplot(aes(x = treatment,
  #              y = subset_mean)) + 
  # geom_quasirandom(data = df, 
  #               aes(x = treatment, 
  #                   y = value, 
  #                   color = subset_name),
  #               size = .8,
  #               alpha = 0.5, 
  #               dodge.width=0.5, 
  #               method = "pseudorandom") +
  #   geom_line(aes(group = subset_name, 
  #                 color = subset_name),
  #             position = position_dodge(width = 0.5), 
  #             color = "gray80") +
  #   geom_point(aes(group = subset_name, 
  #                  color = subset_name),
  #              size = 2.5,
  #              position = position_dodge(width = 0.5)) + 
  #   geom_point(data = by_trmt, 
  #              aes(x = treatment, 
  #                  y = trmt_mean), 
  #              color = "gray10",
  #              size = 4) + 
  #   geom_errorbar(data = by_trmt, 
  #                 aes(y = trmt_mean,
  #                     ymin = lower,
  #                     ymax = upper),
  #                 color = "gray10",
  #                 width = .05) +
  #   ylab("Richness") +
  #   theme_minimal() +
  #   theme(axis.title.x = element_blank())
  
  # Potential beeswarm methods:
  # “quasirandom”, “pseudorandom”, “smiley”, “maxout”, “frowney”, “minout”, “tukey”, “tukeyDense”
}

#   fxn_plot_variables ----
fxn_plot_variables <- function(index_model, index_subset, index_data){
  
  best_fit <- index_model
  file_prefix <- paste0("r_", index_subset, "_")
  
  # Summarize by treatment ----
  by_trmt <- 
    index_data %>%
    group_by(treatment) %>%
    summarize(trmt_mean = mean(value, na.rm = TRUE), 
              n = n(), 
              df = n-1,
              sd = sd(value, na.rm = TRUE),
              se = sd/sqrt(n), 
              t = qt(p = 0.025, 
                     df = df,
                     lower.tail = FALSE), 
              margin = t * se,
              lower = trmt_mean - margin, 
              upper = trmt_mean + margin) 
  
  
  # Create plots ----
  # index_variable_subset <- list_variables_subset[1]

  for(index_variable_subset in list_variables_subset){
    
    subset <- tbl_variables %>%
      filter(variable %in% index_variable_subset)
    palette <- get(subset$palette)
    
    fxn_plot_subset(index_factor = subset$variable, 
                    index_data = index_data) + 
      scale_color_manual(values = palette)
    
    ggsave(filename =
             here(path_out_results, paste0(file_prefix, "treatment-by-", 
                                   str_to_lower(subset$label), 
                                   "_points.png")),
           width = 3, height = 4.5)
  }
}

# Create plots by variable ----
fxn_plot_variables(index_model = mod_rich_nat, 
                   index_subset = "native", 
                   index_data = rich_nat)

fxn_plot_variables(index_model = mod_r_for, 
                   index_subset = "native-forb", 
                   index_data = rich_frb)
# NONNATIVE SPECIES ----
#   Plot the residuals for each variable ----
for(index_variable in list_variables){
  
  variable <- rich_non %>%
    dplyr::select(variable = all_of(index_variable)) %>%
    pull(variable)
  
  png(file = here(path_out_results, paste0("rich_nonnative_residual-plot_model-output_", 
                                   index_variable, 
                                   ".png")), 
      width = 300, height = 600)
  plotResiduals(lmm_non, variable)  
  dev.off()
}


#   fxn_plot_subset ----
fxn_plot_subset <- function(index_factor){
  
  # Summarize by treatment ----
  by_trmt <- 
    rich_non %>%
    group_by(treatment) %>%
    summarize(trmt_mean = mean(value, na.rm = TRUE), 
              n = n(), 
              df = n-1,
              sd = sd(value, na.rm = TRUE),
              se = sd/sqrt(n), 
              t = qt(p = 0.025, 
                     df = df,
                     lower.tail = FALSE), 
              margin = t * se,
              lower = trmt_mean - margin, 
              upper = trmt_mean + margin) 
  
  # Summarize by treatment x factor ----
  df <-   
    rich_non %>%
    arrange(year) %>%
    mutate(year = as_factor(year)) %>%
    dplyr::select(treatment, 
                  subset_name = all_of(index_factor), 
                  value) 
  
  # Plot by treatment x factor ----
  df %>%
    group_by(treatment, subset_name) %>%
    summarize(subset_mean = mean(value, na.rm = TRUE), 
              n = n(), 
              df = n-1,
              sd = sd(value, na.rm = TRUE),
              se = sd/sqrt(n), 
              t = qt(p = 0.025, 
                     df = df,
                     lower.tail = FALSE), 
              margin = t * se,
              lower = subset_mean - margin, 
              upper = subset_mean + margin) %>%
    ungroup() %>%
    ggplot(aes(x = treatment,
               y = subset_mean)) + 
    geom_quasirandom(data = df, 
                     aes(x = treatment, 
                         y = value, 
                         color = subset_name),
                     size = 0.6,
                     alpha = 0.5, 
                     dodge.width=0.5, 
                     method = "pseudorandom") +
    geom_line(aes(group = subset_name, 
                  color = subset_name),
              position = position_dodge(width = 0.5), 
              color = "gray70") +
    geom_point(aes(group = subset_name, 
                   color = subset_name),
               size = 2.5,
               position = position_dodge(width = 0.5)) + 
    geom_point(data = by_trmt, 
               aes(x = treatment, 
                   y = trmt_mean), 
               color = "gray10",
               size = 4) + 
    geom_errorbar(data = by_trmt, 
                  aes(y = trmt_mean,
                      ymin = lower,
                      ymax = upper),
                  color = "gray10",
                  width = .05) +
    # ylab("Richness") +
    theme_minimal() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          legend.position = "none") + 
    ylim(0, 20)
  
  # Potential beeswarm methods:
  # “quasirandom”, “pseudorandom”, “smiley”, “maxout”, “frowney”, “minout”, “tukey”, “tukeyDense”
}

# Summarize by treatment 
by_trmt <- 
  rich_non %>%
  group_by(treatment) %>%
  summarize(trmt_mean = mean(value, na.rm = TRUE), 
            n = n(), 
            df = n-1,
            sd = sd(value, na.rm = TRUE),
            se = sd/sqrt(n), 
            t = qt(p = 0.025, 
                   df = df,
                   lower.tail = FALSE), 
            margin = t * se,
            lower = trmt_mean - margin, 
            upper = trmt_mean + margin) 


# Create plots 
tbl_variables <- 
  tibble(variable = c("year", 
                      "grazer", 
                      "plot_type",
                      "f_break"), 
         label = c("Year", 
                   "Grazer", 
                   "Plot type", 
                   "Grazing interval"), 
         palette = c("palette_yr", 
                     "palette_grazer",
                     "palette_plots",
                     "palette_delta"))  

list_variables_subset <- unique(tbl_variables$variable)

# index_variable_subset <- list_variables_subset[1]
for(index_variable_subset in list_variables_subset){
  
  subset <- tbl_variables %>%
    filter(variable %in% index_variable_subset)
  palette <- get(subset$palette)
  
  fxn_plot_subset(index_factor = subset$variable) + 
    scale_color_manual(values = palette)
  
  ggsave(filename =
           here(path_out_results, paste0("rich_nonnative_treatment-by-", 
                                 str_to_lower(subset$label), 
                                 "_points.png")),
         width = 3, height = 4.5)
}


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
