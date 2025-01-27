# ========================================================== -----
# Load libraries ----
# File management ---
library(fs)   ## To manage directories
library(hms)  ## For working with times
library(janitor)   ## To clean data tables
library(lubridate)   ## To work with dates and times
library(stringr)   ## To wrangle character variables
library(openxlsx)   ## To write xlsx files
library(readxl)   ## To read xlsx files
library(testthat)  ## For data checks
# library(exiftoolr)   ## To work with exif details from jpg files
# library(plotrix)   ## To get color names from hex codes
#
# Define basic functions ----
#   %nin%: Negate %in% 
"%nin%" <- Negate("%in%")
#   substrRight: Subset character string from the right 
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}
# Define file paths ----
fxn_define_file_paths <- function(){
  #
  # For R scripts and work
  path_r <- here()
  path_in <- here("input")
  #
  # Input files
  path_in_data_raw <<- here(path_in, "data_raw")
  path_in_data_derived <<- here(path_in, "data_derived")
  path_in_lookup <<- here(path_in, "lookup-tables")
  #
  # Output files
  path_out <<- here("output")
  path_out_explore <<- here(path_out, "1_explore")
  path_out_prepare <<- here(path_out, "2_prepare")
  path_out_summary <<- here(path_out, "3_summary")
  
  path_out_prepare <<- here(path_out, "0_prepare")
  path_out_clean <<- here(path_out, "1_clean")
  path_out_explore <<- here(path_out, "2_explore")
  path_out_analysis1 <<- here(path_out, "4_analysis-1")
  path_out_analysis2 <<- here(path_out, "5_analysis-2")
  #
}
#
# Execute function to define file paths
fxn_define_file_paths()
# ========================================================== -----
