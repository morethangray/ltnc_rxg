# Function for formatted Quarto tables 
fxn_kable <- function(df){
  
  df  %>%
    knitr::kable() %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"), 
      full_width = TRUE,  
      position = "left", 
      fixed_thead = TRUE)
}
