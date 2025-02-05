# Run setup scripts
source(here::here("R/1_setup.R"))

# Create functions to summarize the metrics by variable group
source(here::here("R/functions/fxn_summarize_models.R"))

# Create function for pretty tables 
source(here::here("R/functions/fxn_kable.R"))

# # Create functions to tidy tables 
# fxn_tidy_table <- function(index_table){
#   
#   index_table %>%
#     dplyr::mutate(estimate = round(bt_estimate, 2), 
#                   conf_low = round(bt_conf_low, 2),
#                   conf_high = round(bt_conf_high, 2), 
#                   statistic = round(statistic, 2)) %>%
#     dplyr::select(dplyr::any_of(c(
#       # "subset", 
#       "treatment",
#       "term",
#       "contrast",
#       "estimate", 
#       "conf_low", 
#       "conf_high",
#       "p_value", 
#       "statistic")))
#   
# }

# Create input data subsets 
abun <- rich_abun$abundance  
abun_nat <- abun$abun_nat 
abun_frb <- abun$abun_frb
abun_non <- abun$abun_non
abun_all <- do.call(bind_rows, abun)

rich <- rich_abun$richness
rich_nat <- rich$rich_nat
rich_frb <- rich$rich_frb
rich_non <- rich$rich_non
rich_all <- do.call(bind_rows, rich)

# Abundance ----
 