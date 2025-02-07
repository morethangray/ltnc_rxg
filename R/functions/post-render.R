#' Post-render Script for Quarto PDF Documents
#' 
#' This script is designed to be used as a post-render hook for Quarto PDF documents.
#' It renames the rendered PDF file to include a timestamp and moves it to an output
#' directory. The script is intended to help maintain versioned copies of rendered
#' documents with clear timestamps. The script works from the project root directory.

# Get the list of output files from Quarto
output_files <- Sys.getenv("QUARTO_PROJECT_OUTPUT_FILES")
if (output_files == "") {
  stop("No output files found in QUARTO_PROJECT_OUTPUT_FILES")
}

# Split the newline-separated list into a vector
files <- strsplit(output_files, "\n")[[1]]

# Filter to only get PDF files
pdf_files <- files[grep("\\.pdf$", files)]

# Process each PDF file
for (pdf_path in pdf_files) {
  # Extract the base filename without extension
  base_name <- tools::file_path_sans_ext(basename(pdf_path))

  # Create timestamp string
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M")

  # Construct new filename with timestamp
  new_filename <- paste0(base_name, "_", timestamp, ".pdf")

  # Define output directory (from project root)
  output_dir <- "output"

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Construct new path
  new_path <- file.path(output_dir, new_filename)

  # Move and rename the file
  if (file.exists(pdf_path)) {
    file.rename(pdf_path, new_path)
    cat(sprintf("Moved rendered PDF to: %s\n", new_path))
  } else {
    warning(sprintf("Could not find PDF file: %s", pdf_path))
  }
}
