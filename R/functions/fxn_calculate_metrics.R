#   fxn_rich ----
fxn_rich <- function(index_subset, index_native) {
  if (index_subset == "frb") {
    subset <-
      input_surveys %>%
      filter(f_native %in% "n1" & f_forb %in% "f1")
  } else {
    (
      subset <-
        input_surveys %>%
        filter(f_native %in% index_native)
    )
  }

  rich <-
    subset %>%
    group_by(index_g) %>%
    summarize(value = n_distinct(genus_species)) %>%
    ungroup() %>%
    right_join(lookup_annotation, "index_g") %>%
    mutate(
      subset = index_subset,
      metric = "rich",
      met_sub = paste(metric, subset, sep = "_"),
      value = ifelse(is.na(value), 0, value),
      value_log = log(value + 1e-6)
    ) %>%
    select(
      treatment,
      starts_with("value"),
      plot_name,
      plot_type,
      year,
      interval,
      grazer,
      f_year,
      f_break,
      f_new,
      f_one_yr,
      f_two_yr,
      met_sub,
      subset
    )
}
#   fxn_abun ----
fxn_abun <- function(index_subset, index_native) {
  if (index_subset == "frb") {
    subset <-
      input_surveys %>%
      filter(f_native %in% "n1" & f_forb %in% "f1")
  } else {
    (
      subset <-
        input_surveys %>%
        filter(f_native %in% index_native)
    )
  }

  abun <-
    subset %>%
    group_by(index_g) %>%
    summarize(value = sum(value, na.rm = TRUE)) %>%
    ungroup() %>%
    right_join(lookup_annotation, "index_g") %>%
    mutate(
      subset = index_subset,
      metric = "abun",
      met_sub = paste(metric, subset, sep = "_"),
      value = ifelse(is.na(value), 0, value),
      value_log = log(value + 1e-6),
      value_std = as.vector(scale(value)),
      value_sqrt = sqrt(value)
    ) %>%
    select(
      treatment,
      starts_with("value"),
      plot_name,
      plot_type,
      year,
      interval,
      grazer,
      f_year,
      f_break,
      f_new,
      f_one_yr,
      f_two_yr,
      met_sub,
      subset
    )
}
