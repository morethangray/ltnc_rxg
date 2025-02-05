#' Plot histograms for variable transformations
#'
#' Creates a series of histograms to visualize the distribution of data values
#' and their transformations (standardized, log-transformed, and square root-transformed).
#' This helps in assessing the effectiveness of different data transformations for
#' achieving normality.
#'
#' @param data A data frame containing the response variable and its transformations.
#'        Must be processed using fxn_prepare_input_data().
#' @param data_name Character string specifying the name of response and species group
#'        pair (e.g., "native richness"). Used for plot titles.
#'
#' @return Plots multiple histograms showing the distribution of original and 
#'         transformed data.
#'
#' @examples
#' \dontrun{
#' data <- fxn_prepare_input_data(your_data)
#' fxn_histogram(data, "native richness")
#' }
#'
#' @export
fxn_histogram <- function(data, data_name) {
  
  # Create base histogram
  hist(data$value, 
       breaks = 15, 
       col = "gray", 
       main = paste("Histogram of", data_name), 
       xlab = "Value", 
       ylab = "Count")
  
  # Create transformed histograms if available
  if ("value_std" %in% names(data)) {
    hist(data$value_std, 
         breaks = 25, 
         col = "gray", 
         main = paste("Histogram of standardized", data_name), 
         xlab = "Standardized Value", 
         ylab = "Count")
  }
  
  if ("value_log" %in% names(data)) {
    hist(data$value_log, 
         breaks = 25, 
         col = "gray", 
         main = paste("Histogram of log-transformed", data_name), 
         xlab = "Log Value", 
         ylab = "Count")
  }
  
  if ("value_sqrt" %in% names(data)) {
    hist(data$value_sqrt, 
         breaks = 25, 
         col = "gray", 
         main = paste("Histogram of square root-transformed", data_name), 
         xlab = "Sqrt Value", 
         ylab = "Count")
  }
}

#' Conduct goodness of fit tests for species richness data
#'
#' Performs goodness of fit tests for Poisson and negative binomial distributions
#' on species richness data. Creates diagnostic plots to visualize the fit of each
#' distribution to the observed data.
#'
#' @param data A data frame containing the response variable.
#'        Must be processed using fxn_prepare_input_data().
#' @param data_name Character string specifying the name of response and species group
#'        pair (e.g., "native richness"). Used for documentation.
#'
#' @return Prints goodness of fit test results and creates diagnostic plots for
#'         Poisson and negative binomial distributions.
#'
#' @import vcd
#'
#' @examples
#' \dontrun{
#' data <- fxn_prepare_input_data(your_data)
#' fxn_gof_test_rich(data, "native richness")
#' }
#'
#' @export
fxn_gof_test_rich <- function(data, data_name) {
  # Validate required packages
  required_packages <- c("vcd")
  if (!fxn_check_packages(required_packages)) {
    stop("Required package 'vcd' is missing. Please install it.")
  }
  
  # Fit and plot Poisson distribution
  tryCatch({
    gf_pois <- vcd::goodfit(data$value, type = "poisson")
    pois_first_line <- summary(gf_pois)[1, ]
    plot(gf_pois, main = "") 
  }, error = function(e) {
    warning("Error fitting Poisson distribution: ", e$message)
  })
  
  # Fit and plot negative binomial distribution
  tryCatch({
    gf_nbin <- vcd::goodfit(data$value, type = "nbinomial")
    nbin_first_line <- summary(gf_nbin)[1, ]
    plot(gf_nbin, main = "") 
  }, error = function(e) {
    warning("Error fitting negative binomial distribution: ", e$message)
  })
}

#' Conduct goodness of fit tests for abundance data
#'
#' Performs goodness of fit tests for normal distribution on original and
#' transformed abundance data. Creates diagnostic plots and conducts Shapiro-Wilk
#' tests for each transformation.
#'
#' @param data A data frame containing the response variable and its transformations.
#'        Must be processed using fxn_prepare_input_data().
#' @param data_name Character string specifying the name of response and species group
#'        pair (e.g., "native abundance"). Used for documentation.
#'
#' @return Prints Shapiro-Wilk test results and creates diagnostic plots for
#'         normal distribution fits on original and transformed data.
#'
#' @import fitdistrplus
#' @import stats
#'
#' @examples
#' \dontrun{
#' data <- fxn_prepare_input_data(your_data)
#' fxn_gof_test_abun(data, "native abundance")
#' }
#'
#' @export
fxn_gof_test_abun <- function(data, data_name) {
  # Validate required packages
  required_packages <- c("fitdistrplus", "stats")
  if (!fxn_check_packages(required_packages)) {
    stop("Required packages 'fitdistrplus' and/or 'stats' are missing. Please install them.")
  }
  
  # Test and plot original values
  print(paste("Normal Distribution Fit for", data_name))
  print(stats::shapiro.test(data$value))
  tryCatch({
    fitdistrplus::fitdist(data$value, dist = "norm") %>%
      plot(histo = TRUE, demp = TRUE)
  }, error = function(e) {
    warning("Error fitting normal distribution to original values: ", e$message)
  })
  cat("\n")
  
  # Test and plot standardized values
  print(paste("Normal Distribution Fit for Standardized", data_name))
  print(stats::shapiro.test(data$value_std))
  tryCatch({
    fitdistrplus::fitdist(data$value_std, dist = "norm") %>%
      plot(histo = TRUE, demp = TRUE)
  }, error = function(e) {
    warning("Error fitting normal distribution to standardized values: ", e$message)
  })
  cat("\n")
  
  # Test and plot log-transformed values
  print(paste("Normal Distribution Fit for Log-transformed", data_name))
  print(stats::shapiro.test(data$value_log))
  tryCatch({
    fitdistrplus::fitdist(data$value_log, dist = "norm") %>%
      plot(histo = TRUE, demp = TRUE)
  }, error = function(e) {
    warning("Error fitting normal distribution to log-transformed values: ", e$message)
  })
  cat("\n")
  
  # Test and plot square root-transformed values
  print(paste("Normal Distribution Fit for Square Root-transformed", data_name))
  print(stats::shapiro.test(data$value_sqrt))
  tryCatch({
    fitdistrplus::fitdist(data$value_sqrt, dist = "norm") %>%
      plot(histo = TRUE, demp = TRUE)
  }, error = function(e) {
    warning("Error fitting normal distribution to square root-transformed values: ", e$message)
  })
}
  
#' Conduct Vuong tests for comparing model fits
#'
#' Performs Vuong tests to compare the fit of different count data models
#' (negative binomial, Poisson, and zero-inflated models if applicable).
#' Tests are only performed when a treatment variable is present in the data.
#'
#' @param data A data frame containing the response variable and treatment.
#'        Must be processed using fxn_prepare_input_data().
#' @param data_name Character string specifying the name of response and species group
#'        pair (e.g., "native richness"). Used for documentation.
#'
#' @return Prints results of Vuong tests comparing different model fits.
#'
#' @import MASS
#' @import pscl
#' @import stats
#'
#' @examples
#' \dontrun{
#' data <- fxn_prepare_input_data(your_data)
#' fxn_vuong_test(data, "native richness")
#' }
#'
#' @export
fxn_vuong_test <- function(data, data_name) {
  # Validate required packages
  required_packages <- c("MASS", "pscl", "stats")
  if (!fxn_check_packages(required_packages)) {
    stop("Required packages 'MASS', 'pscl', and/or 'stats' are missing. Please install them.")
  }
  
  # Only proceed if there are multiple treatment levels
  if (length(unique(data$treatment)) > 1) {
    tryCatch({
      nb <- MASS::glm.nb(value ~ treatment, data = data)
      p <- stats::glm(value ~ treatment, data = data, family = "poisson")
      
      # Compare negative binomial vs. Poisson
      print(paste("Vuong test: Negative Binomial vs. Poisson for", data_name))
      print(pscl::vuong(nb, p))
      cat("\n\n")
      
      # If zero values present, also fit and compare zero-inflated model
      if (min(data$value) == 0) {
        z <- pscl::zeroinfl(value ~ treatment, data = data)
        print(paste("Vuong test: Poisson vs. Zero-Inflated for", data_name))
        print(pscl::vuong(p, z))
        cat("\n\n")
        
        print(paste("Vuong test: Negative Binomial vs. Zero-Inflated for", data_name))
        print(pscl::vuong(nb, z))
      }
    }, error = function(e) {
      message(paste("Error in Vuong tests for", data_name, ":", e$message))
    })
  } else {
    message("Skipping Vuong tests as data contains only one treatment level.")
  }
}

#' Create quantile-quantile (QQ) plots for distribution fitting
#'
#' Creates QQ plots to assess the fit of various theoretical distributions
#' to the observed data. For non-native richness data, additional distributions
#' (normal, lognormal, and gamma) are tested.
#'
#' @param data A data frame containing the response variable.
#'        Must be processed using fxn_prepare_input_data().
#' @param data_name Character string specifying the name of response and species group
#'        pair (e.g., "native richness"). Used for documentation.
#'
#' @return Creates QQ plots for various distributions.
#'
#' @import car
#' @import MASS
#'
#' @examples
#' \dontrun{
#' data <- fxn_prepare_input_data(your_data)
#' fxn_qq_plot(data, "native richness")
#' }
#'
#' @export
fxn_qq_plot <- function(data, data_name) {
  # Validate required packages
  required_packages <- c("car", "MASS")
  if (!fxn_check_packages(required_packages)) {
    stop("Required packages 'car' and/or 'MASS' are missing. Please install them.")
  }
  
  # Add 1 to avoid issues with zero values
  data$value_1 <- data$value + 1
  
  # Fit and plot Poisson QQ plot
  tryCatch({
    poisson <- MASS::fitdistr(data$value_1, "Poisson")
    car::qqp(data$value_1, "pois", lambda = poisson$estimate, 
             main = paste("Poisson QQ plot for", data_name))
  }, error = function(e) {
    warning("Error creating Poisson QQ plot: ", e$message)
  })
  
  # Additional distributions for non-native richness
  if (data_name == "nonnative richness") {
    # Normal QQ plot
    car::qqp(data$value_1, "norm", 
             main = paste("Normal QQ plot for", data_name))
    
    # Lognormal QQ plot
    car::qqp(data$value_1, "lnorm", 
             main = paste("Lognormal QQ plot for", data_name))
    
    # Gamma QQ plot
    tryCatch({
      gamma <- MASS::fitdistr(data$value_1, "gamma")
      car::qqp(data$value_1, "gamma", 
               shape = gamma$estimate[[1]], 
               rate = gamma$estimate[[2]], 
               main = paste("Gamma QQ plot for", data_name))
    }, error = function(e) {
      warning("Error creating Gamma QQ plot: ", e$message)
    })
  }
}