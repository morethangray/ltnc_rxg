fxn_make_table_1 <- function(index_model_name){
  
  # Validate required packages
  required_packages <- c("dplyr", "glue", "stringr")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Retrieve the model object
  index_model <- get(index_model_name)
  
  # Extract and collate model metadata
  model_formula <- as.character(formula(index_model))[3]
  model_formula_pretty <- fxn_make_pretty_model_formula(index_formula = model_formula)
  
  subset_name <- lookup_tables$lookup_species_subset %>%
    dplyr::filter(abbr_subset == stringr::str_sub(index_model_name, -3)) %>%
    dplyr::pull(subset)
  
  response <- if(stringr::str_detect(index_model_name, "abun")) "abundance" else "richness"
  response_upper <- stringr::str_to_title(response)
  
  transformation_init <- stringr::str_remove(as.character(formula(index_model))[2], "value_")
  transformation <- if(transformation_init == "value") "none" else transformation_init
  
  model_family <- summary(index_model)$objClass[1]
  distribution_name <- if(is.null(model_family)) "Negative binomial" else "Normal"
  
  # Populate table  
  table_subset <- tibble::tibble(
    response_variable = glue::glue("{subset_name} {response}"), 
    formula = dplyr::case_when(
      transformation == "none" ~ glue::glue("{response_upper} ~ {model_formula_pretty}"),
      TRUE ~ glue::glue("{transformation}({response_upper}) ~ {model_formula_pretty}")),
    distribution = glue::glue("{distribution_name} distribution")
  )
  
  return(table_subset)
  
}
