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


