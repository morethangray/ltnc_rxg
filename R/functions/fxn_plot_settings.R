# Define plot layouts for percent cover ----
#   fxn_boxplot ----
fxn_boxplot <- function(index_x, index_fill) {
  percent_cover %>%
    select(value,
      index_x = all_of(index_x),
      index_fill = all_of(index_fill),
      nativity,
      forb_status
    ) %>%
    ggplot(aes(
      x = index_x,
      y = value,
      fill = index_fill
    )) +
    geom_boxplot(
      width = 0.8,
      outlier.shape = NA,
      alpha = 0.8
    ) +
    geom_point(
      alpha = 0.3,
      size = 0.4,
      position = position_jitterdodge(
        jitter.width = 0.2,
        jitter.height = 0.3,
        dodge.width = 0.8
      )
    ) +
    theme_minimal() +
    scale_fill_manual(values = palette_plots) +
    ylab("Percent cover") +
    theme(
      axis.title.x = element_blank(),
      panel.spacing = unit(2, "lines")
    )
}
#   fxn_violin ----
fxn_violin <- function(index_x, index_fill) {
  percent_cover %>%
    select(value,
      index_x = all_of(index_x),
      index_fill = all_of(index_fill),
      nativity,
      forb_status
    ) %>%
    ggplot(aes(
      x = index_x,
      y = value,
      fill = index_fill
    )) +
    geom_violin(
      width = 0.8,
      outlier.shape = NA,
      alpha = 0.8
    ) +
    theme_minimal() +
    scale_fill_manual(values = palette_plots) +
    ylab("Percent cover") +
    theme(
      axis.title.x = element_blank(),
      panel.spacing = unit(2, "lines")
    )
}

# ---------------------------------------------------------- -----
# Define palette colors ----
#   plot_colors ----
plot_colors <-
  read_csv(here(
    path_in_lookup,
    "attributes_palette.csv"
  )) %>%
  arrange(palette, palette_subset, levels)

#   palette_treatment ----
palette_treatment_line_light <-
  plot_colors %>%
  filter(
    palette %in% "treatment",
    palette_subset %in% "line_light"
  ) %>%
  arrange(desc(levels)) %>%
  pull(hex_code)

#   palette_treatment_dark ----
palette_treatment <-
  plot_colors %>%
  filter(palette %in% "treatment") %>%
  pull(hex_code)

#   palette_plots ----
palette_plots <-
  plot_colors %>%
  filter(
    palette %in% "plot_type",
    palette_subset %in% "teal"
  ) %>%
  pull(hex_code)

#   palette_plots_dark ----
palette_plots <-
  plot_colors %>%
  filter(
    palette %in% "plot_type",
    palette_subset %in% "dark"
  ) %>%
  pull(hex_code)

palette_yr <- pals::parula(5)[1:4]

#   palette_grazer ----
palette_grazer <-
  plot_colors %>%
  filter(
    palette %in% "grazer",
    palette_subset %in% "points"
  ) %>%
  pull(hex_code)

#   palette_delta ----
palette_delta <-
  plot_colors %>%
  filter(
    palette %in% "delta_t",
    palette_subset %in% "points"
  ) %>%
  pull(hex_code)

# ---------------------------------------------------------- -----

