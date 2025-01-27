# Convert xlsx to csv
# list_years <- 2019:2022
# index_year <- list_years[2]
fxn_xlsx_to_csv <- function(list_years) {
  datalist_year <- list()
  
  for (index_year in list_years) {
    # Read xlsx ----
    path <- here(
      path_in_data_derived,
      glue(
        "wantrup_{index_year}.xlsx"
      )
    )
    input_data <-
      path %>%
      excel_sheets() %>%
      set_names() %>%
      map_dfr(read_excel,
        path = path,
        .id = "sheet"
      )

    # Reshape and annotate survey table ----
    column_count <- ncol(input_data)

    # Gather into long format
    long_data <-
      input_data %>%
      gather(orig_plot, value, 3:all_of(column_count)) %>%
      drop_na(value) %>%
      mutate(
        year = index_year,
        value = as.numeric(value)
      ) %>%
      # Join plant attributes
      left_join(lookup_plants, "orig_name") %>%
      filter(include_plant == TRUE) %>%
      mutate(
        native = ifelse(is_native == TRUE, "Native", "Non-native"),
        forb = ifelse(is_forb == TRUE, "Forb", "Non-forb"),
        f_native = ifelse(is_native == TRUE, "n1", "n0"),
        f_forb = ifelse(is_forb == TRUE, "f1", "f0")
      ) %>%
      # Join plot attributes
      left_join(lookup_plot_names, "orig_plot") %>%
      mutate(
        index = paste(plot_type, plot_name, year, sep = "_"),
        index_g = paste(index, f_grazed, sep = "_"),
        index_g_n_f = paste(index_g, f_native, f_forb, sep = "_")
      ) %>%
      # Join temporal attributes
      left_join(lookup_grazing_interval, "index") %>%
      # Reorder columns
      select(
        plot_type,
        plot_name,
        year,
        month,
        treatment,
        grazer,
        native,
        forb,
        interval,
        genus_species,
        value,
        starts_with("f_"),
        starts_with("index"),
        taxonomic_name
      ) %>%
      distinct()

    # Combine counts for taxa with ssp. or var. ----
    list_multi_gs <-
      long_data %>%
      distinct(genus_species, taxonomic_name) %>%
      group_by(genus_species) %>%
      count() %>%
      ungroup() %>%
      filter(n > 1) %>%
      pull(genus_species)

    one_gs <-
      long_data %>%
      filter(genus_species %nin% list_multi_gs) %>%
      select(-taxonomic_name)

    # Calculate the maximum count for Apr & May surveys ----
    percent_cover <-
      long_data %>%
      filter(genus_species %in% list_multi_gs) %>%
      group_by(across(c(-value, -taxonomic_name))) %>%
      summarize(value = sum(value)) %>%
      ungroup() %>%
      bind_rows(one_gs) %>%
      distinct() %>%
      group_by(across(c(-month, -value))) %>%
      summarize(value = max(value)) %>%
      ungroup() %>%
      relocate(value, .after = genus_species)

    datalist_year[[index_year]] <- percent_cover
  }

  bind_years <-
    do.call(
      bind_rows,
      datalist_year
    )

  return(bind_years)
}
