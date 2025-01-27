# Richness ----
rich <-
  read_csv(here(
    path_in_data_derived,
    "rich_nat-frb-non_2019-2022.csv"
  )) %>%
  # arrange(desc(treatment), year, f_break) %>%
  arrange(year, f_break) %>%
  mutate(value_sqrt = sqrt(value)) %>%
  mutate_if(is.character, as_factor) %>%
  dplyr::select(
    treatment,
    starts_with("value"),
    plot_name,
    plot_type,
    year,
    grazer,
    f_year,
    f_break,
    f_new,
    f_one_yr,
    f_two_yr,
    met_sub
  )

# Native
rich_nat <- rich %>%
  filter(met_sub %in% "rich_nat")

# Native forb
rich_frb <- rich %>%
  filter(met_sub %in% "rich_frb")

# Non-native
rich_non <- rich %>%
  filter(met_sub %in% "rich_non")

rich_non_s <- as.vector(scale(rich_non$value))
rich_non$value_std <- rich_non_s
rich_non$value_sqrt <- sqrt(rich_non$value)

# Abundance ----
abun <-
  read_csv(here(
    path_in_data_derived,
    "abun_nat-frb-non_2019-2022.csv"
  )) %>%
  arrange(desc(treatment), year, f_break) %>%
  mutate_if(is.character, as_factor) %>%
  dplyr::select(
    treatment,
    starts_with("value"),
    plot_name,
    plot_type,
    year,
    grazer,
    f_year,
    f_break,
    f_new,
    f_one_yr,
    f_two_yr,
    met_sub
  )


# Native
abun_nat <- abun %>%
  filter(met_sub %in% "abun_nat")

# Native forb
abun_frb <- abun %>%
  filter(met_sub %in% "abun_frb")

# Non-native
abun_non <- abun %>%
  filter(met_sub %in% "abun_non")
