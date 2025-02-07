fxn_get_citation <- function(index_package){
  
  result <- 
    as_tibble(capture.output(citation(index_package))) %>%
    dplyr::mutate(is_blank = ifelse(value == "", 1, 0), 
                  cum_sum = cumsum(is_blank),
                  is_ref = is_blank == 0 & cum_sum == 1, 
                  row_n = 1:n(), 
                  package = index_package, 
                  value = stringr::str_remove_all(value, "  ")) %>%
    dplyr::filter(is_ref == TRUE) %>%
    dplyr::select(package, ref = value, row_n)
  
  return(result)
}
# 
# datalist <- list()
# for(p in packages_needed){
#   datalist[[p]] <- fxn_get_citation(p)
# }
# 
# all_citations <- do.call(bind_rows, datalist) 
# 
# print(all_citations, n = 100)
# 
# all_citations %>%
#   write_csv("r_package_citations.csv")
