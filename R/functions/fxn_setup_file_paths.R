#' Define and Create Project File Paths
#'
#' Sets up input and output directory paths for the project
#' Creates directories if they do not exist
#'
#' @return List of file paths
#' @import here
#' @export
fxn_setup_file_paths <- function() {
  # Validate required packages
  required_packages <- c("here")
  
  if(!fxn_check_packages(required_packages)) {
    stop("Required packages are missing. Please install them.")
  }
  
  # Base project paths
  path_r <- here::here()
  path_in <- here::here("input")
  path_out <- here::here("output")
  
  # Create directories if they don't exist
  dirs_to_create <- c(
    path_in,
    here::here(path_in, "data"),
    here::here(path_in, "lookup_tables"),
    path_out
  )
  
  # Create directories
  lapply(dirs_to_create, function(dir_path) {
    if (!dir.exists(dir_path)) {
      tryCatch(
        {
          dir.create(dir_path, recursive = TRUE)
          message("Created directory: ", dir_path)
        },
        error = function(e) {
          warning("Could not create directory: ", dir_path, "\n", e$message)
        }
      )
    }
  })
  
  # Return a list of paths  
  return(list(
    path_r = path_r,
    path_in = path_in,
    path_in_data = here::here(path_in, "data"),
    path_in_lookup = here::here(path_in, "lookup_tables"),
    path_out = path_out
  ))
}

# Store paths in a global variable
project_paths <- fxn_setup_file_paths()

# If you need global variables (not recommended, but possible)
# list2env(project_paths, envir = .GlobalEnv)