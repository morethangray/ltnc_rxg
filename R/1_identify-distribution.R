# Load libraries, functions, workflows -----
rm(list = ls())
#
library(tidyverse) ## To manipulate data frames
library(here) ## To manage directories
#
source(here("R/functions/fxn_utilities.R"))
#
# ========================================================== -----
# PREPARE DATA  ----
source(here("R/functions/fxn_load_rich_abun.R"))

# Load libraries ----
library(vcd)
library(pscl)
library(MASS)
library(car)
library(fitdistrplus)
library(tidyr) # For pivot_longer

# fig : Simple theme for ggplot figures  ----
fig <-
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), panel.background = element_blank(),
    strip.background = element_blank(), strip.text.y = element_text(),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  )
# ========================================================== -----
# Gemini ----
# Function to perform distribution checks for count data (Richness)
data = datasets_count[[1]]
# data_name = "rich_nat"
check_count_distribution <- function(data, data_name) {
  message(paste("Checking distribution for:", data_name))
  
  # Basic histogram
  hist(data$value, breaks = 15, col = "gray", main = paste("Histogram of", data_name), xlab = "Value", ylab = "Count")
  
  # Goodness-of-fit tests
  gf_pois <- goodfit(data$value, type = "poisson")
  print(paste("Poisson Goodness of Fit for", data_name))
  print(summary(gf_pois))
  plot_gf_pois <- plot(gf_pois, main = paste("Poisson Goodness of Fit Plot for", data_name))
  
  gf_nbin <- goodfit(data$value, type = "nbinomial")
  print(paste("Negative Binomial Goodness of Fit for", data_name))
  print(summary(gf_nbin))
  plot_gf_nbin <- plot(gf_nbin, main = paste("Negative Binomial Goodness of Fit Plot for", data_name))
  
  
  # Vuong tests (if applicable)
  if (length(unique(data$treatment)) > 1) { # Only run if there's a treatment variable
    z <- zeroinfl(value ~ treatment, data = data)
    nb <- glm.nb(value ~ treatment, data = data)
    p <- glm(value ~ treatment, data = data, family = "poisson")
    print(paste("Vuong test: Poisson vs. Zero-Inflated for", data_name))
    print(vuong(p, z))
    print(paste("Vuong test: Negative Binomial vs. Zero-Inflated for", data_name))
    print(vuong(nb, z))
    print(paste("Vuong test: Negative Binomial vs. Poisson for", data_name))
    print(vuong(nb, p))
  } else {
    message("Skipping Vuong tests as no treatment variable is present.")
  }
  
  # Q-Q plots with shifted data
  data$value_1 <- data$value + 1
  plot_qq_norm <- qqp(data$value_1, "norm", main = paste("Normal Q-Q Plot for", data_name))
  plot_qq_lognorm <- qqp(data$value_1, "lnorm", main = paste("Lognormal Q-Q Plot for", data_name))
  
  # Fit distributions and Q-Q plots
  nbinom <- fitdistr(data$value_1, "Negative Binomial")
  plot_qq_nbinom <- qqp(data$value_1, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]], main = paste("Negative Binomial Q-Q Plot for", data_name))
  
  poisson <- fitdistr(data$value_1, "Poisson")
  plot_qq_poisson <- qqp(data$value_1, "pois", lambda = poisson$estimate, main = paste("Poisson Q-Q Plot for", data_name))
  
  gamma <- fitdistr(data$value_1, "gamma")
  plot_qq_gamma <- qqp(data$value_1, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]], main = paste("Gamma Q-Q Plot for", data_name))
  
  invisible(data) # Return the data invisibly for chaining

}

# Function to perform distribution checks for continuous data (Abundance)
check_continuous_distribution <- function(data, data_name) {
  message(paste("Checking distribution for:", data_name))
  # Histograms
  hist(data$value, breaks = 25, col = "gray", main = paste("Histogram of", data_name), xlab = "Value", ylab = "Count")
  hist(data$value_std, breaks = 25, col = "gray", main = paste("Histogram of Standardized", data_name), xlab = "Standardized Value", ylab = "Count")
  hist(data$value_log, breaks = 25, col = "gray", main = paste("Histogram of Log-transformed", data_name), xlab = "Log Value", ylab = "Count")
  hist(data$value_sqrt, breaks = 25, col = "gray", main = paste("Histogram of Square Root-transformed", data_name), xlab = "Sqrt Value", ylab = "Count")
  
  # Q-Q Plots
  qqp(data$value_std, "norm", main = paste("Normal Q-Q Plot of Standardized", data_name))
  qqp(data$value_sqrt, "norm", main = paste("Normal Q-Q Plot of Square Root-transformed", data_name))
  qqp(data$value_log, "norm", main = paste("Normal Q-Q Plot of Log-transformed", data_name))
  
  # fitdistrplus
  fitdist(data$value, dist = "norm") %>% plot(histo = FALSE, demp = TRUE, main = paste("Distribution Fit for", data_name))
  fitdist(data$value_log, dist = "norm") %>% plot(histo = FALSE, demp = TRUE, main = paste("Distribution Fit for Log-transformed", data_name))
  fitdist(data$value_sqrt, dist = "norm") %>% plot(histo = FALSE, demp = TRUE, main = paste("Distribution Fit for Square Root-transformed", data_name))
  fitdist(data$value, dist = "logis") %>% plot(histo = FALSE, demp = TRUE, main = paste("Logistic Distribution Fit for", data_name))
  fitdist(data$value_log, dist = "logis") %>% plot(histo = FALSE, demp = TRUE, main = paste("Logistic Distribution Fit for Log-transformed", data_name))
  fitdist(data$value_sqrt, dist = "logis") %>% plot(histo = FALSE, demp = TRUE, main = paste("Logistic Distribution Fit for Square Root-transformed", data_name))
  invisible(data)
}


# Prepare data for plotting effect variables (using tidyr::pivot_longer)
prepare_data_for_boxplot <- function(data, data_name) {
  data %>%
    mutate(value_s = as.vector(scale(value))) %>% # Scale within the function
    pivot_longer(cols = c(treatment, f_year, f_break, f_one_yr, f_two_yr, grazer),
                 names_to = "predictor", values_to = "predictor_value")
}

# Function to create boxplots
create_effect_boxplots <- function(data, data_name) {
  plot_data <- prepare_data_for_boxplot(data, data_name)
  ggplot(plot_data, aes(y = value_s, x = predictor_value)) +
    geom_boxplot() +
    fig + # Assuming 'fig' is your ggplot theme
    facet_wrap(~ predictor, scales = "free_x") + # Use facet_wrap for separate plots
    geom_hline(yintercept = c(0, 0.75, -0.75), linetype = c("solid", "dashed", "dashed")) +
    labs(title = paste("Effect of Predictors on Scaled Value for", data_name),
         x = "Predictor Value", y = "Scaled Value")
}



# Apply functions to your data
datasets_count <- list(rich_nat = rich_nat, rich_frb = rich_frb, rich_non = rich_non)
datasets_continuous <- list(abun_nat = abun_nat, abun_frb = abun_frb, abun_non = abun_non)

purrr::iwalk(datasets_count, check_count_distribution)
purrr::iwalk(datasets_continuous, check_continuous_distribution)

# Create and display boxplots
purrr::iwalk(datasets_count, create_effect_boxplots)
purrr::iwalk(datasets_continuous, create_effect_boxplots)
# Claude ----
# Function to perform comprehensive distribution analysis
# data = rich_non
# var_type = "count"

analyze_distribution <- function(data, response_var = "value", var_type = c("count", "continuous")) {
  var_type <- match.arg(var_type)
  
  # Basic histogram
  hist(data[[response_var]], 
       breaks = ifelse(var_type == "count", 15, 25),
       col = "gray", main = "",
       xlab = ifelse(var_type == "count", "Richness", "Abundance"),
       ylab = "Count")
  
  # For count data, perform discrete distribution tests
  if (var_type == "count") {
    # Poisson and negative binomial tests
    tests <- list(
      poisson = summary(vcd::goodfit(data[[response_var]], type = "poisson")),
      nbinom = summary(vcd::goodfit(data[[response_var]], type = "nbinomial"))
    )
    
    # Vuong tests for zero-inflation
    z <- pscl::zeroinfl(as.formula(paste(response_var, "~ treatment")), data = data)
    nb <- MASS::glm.nb(as.formula(paste(response_var, "~ treatment")), data = data)
    p <- glm(as.formula(paste(response_var, "~ treatment")), data = data, family = "poisson")
    
    vuong_tests <- list(
      poisson_vs_zeroinfl = pscl::vuong(p, z),
      nb_vs_zeroinfl = pscl::vuong(nb, z),
      nb_vs_poisson = pscl::vuong(nb, p)
    )
    
    return(list(goodfit_tests = tests, vuong_tests = vuong_tests))
  }
  
  # For continuous data, test different distributions using fitdistrplus
  if (var_type == "continuous") {
    distributions <- c("norm", "gamma", "lnorm", "logis")
    fits <- lapply(distributions, function(dist) {
      tryCatch(
        fitdistrplus::fitdist(data[[response_var]], dist),
        error = function(e) NULL
      )
    })
    names(fits) <- distributions
    return(fits)
  }
}

# Function to create standardized boxplots for all predictors
plot_predictors <- function(data, response_var = "value_std", predictors = c("treatment", "f_year", "f_break", "f_one_yr", "f_two_yr", "grazer")) {
  lapply(predictors, function(pred) {
    ggplot(data, aes(y = .data[[response_var]], x = .data[[pred]])) +
      geom_boxplot() +
      fig +  
      geom_hline(yintercept = c(0, 0.75, -0.75), 
                 linetype = c("solid", "dashed", "dashed"))
  })
}

# Usage example:
# For count data (richness):
rich_nat_dist <- analyze_distribution(rich_nat, var_type = "count")
rich_non_plots <- plot_predictors(rich_non)

# For continuous data (abundance):
abun_nat_dist <- analyze_distribution(abun_nat, var_type = "continuous")
abun_nat_plots <- plot_predictors(abun_nat)

# GPT ----
# Define a function to check and compare distributions 
check_distribution <- function(data, response_var, response_var_adjusted = NULL, predictors) {
  response <- data[[response_var]]
  
  # Plot histogram
  hist(response,
       breaks = 15, col = "gray", main = "",
       xlab = "Richness", ylab = "Count")
  
  # Perform Goodness-of-Fit Tests
  goodfit_pois <- vcd::goodfit(response, type = "poisson")
  goodfit_nbin <- vcd::goodfit(response, type = "nbinomial")
  
  # Summary and plots for both fits
  list(
    pois = list(summary = summary(goodfit_pois), plot = plot(goodfit_pois)),
    nbin = list(summary = summary(goodfit_nbin), plot = plot(goodfit_nbin))
  )
  
  # Model comparisons
  z <- pscl::zeroinfl(as.formula(paste(response_var, "~", predictors)), data = data)
  nb <- MASS::glm.nb(as.formula(paste(response_var, "~", predictors)), data = data)
  p <- glm(as.formula(paste(response_var, "~", predictors)), data = data, family = "poisson")
  
  # Vuong tests
  vuong_results <- list(
    poisson_vs_zeroinfl = pscl::vuong(p, z),
    nb_vs_zeroinfl = pscl::vuong(nb, z),
    nb_vs_poisson = pscl::vuong(nb, p)
  )
  
  # Adjust response for zero-inflation (if applicable)
  if (!is.null(response_var_adjusted)) {
    data[[response_var_adjusted]] <- response + 1
    
    # QQ Plots
    fit_nbinom <- MASS::fitdistr(data[[response_var_adjusted]], "Negative Binomial")
    fit_poisson <- MASS::fitdistr(data[[response_var_adjusted]], "Poisson")
    fit_gamma <- MASS::fitdistr(data[[response_var_adjusted]], "gamma")
    
    qq_results <- list(
      nbinom = qqp(data[[response_var_adjusted]], "nbinom", size = fit_nbinom$estimate[[1]], mu = fit_nbinom$estimate[[2]]),
      poisson = qqp(data[[response_var_adjusted]], "pois", lambda = fit_poisson$estimate),
      gamma = qqp(data[[response_var_adjusted]], "gamma", shape = fit_gamma$estimate[[1]], rate = fit_gamma$estimate[[2]])
    )
    
    return(list(vuong_results = vuong_results, qq_results = qq_results))
  }
  
  return(list(vuong_results = vuong_results))
}

# Function for boxplots 
create_boxplots <- function(data, response_var, factor_vars) {
  data[[paste0(response_var, "_scaled")]] <- scale(data[[response_var]])
  
  plots <- map(factor_vars, ~ ggplot(data, aes_string(y = paste0(response_var, "_scaled"), x = .x)) +
                 geom_boxplot() +
                 geom_hline(yintercept = c(0, 0.75, -0.75), linetype = c("solid", "dashed", "dashed")) +
                 labs(title = paste("Effect of", .x, "on", response_var)))
  
  return(plots)
}

# Apply the functions to datasets 
factor_vars <- c("treatment", "f_year", "f_break", "f_one_yr", "f_two_yr", "grazer")

# For each subset, perform the analysis
datasets <- list(
  rich_nat = list(data = rich_nat, response_var = "value", predictors = "treatment"),
  rich_frb = list(data = rich_frb, response_var = "value", predictors = "treatment"),
  rich_non = list(data = rich_non, response_var = "value", predictors = "treatment"),
  abun_nat = list(data = abun_nat, response_var = "value", predictors = "treatment"),
  abun_frb = list(data = abun_frb, response_var = "value", predictors = "treatment"),
  abun_non = list(data = abun_non, response_var = "value", predictors = "treatment")
)

results <- map(datasets, function(d) {
  list(
    distribution_check = check_distribution(d$data, d$response_var, response_var_adjusted = "value_1", d$predictors),
    boxplots = create_boxplots(d$data, d$response_var, factor_vars)
  )
})
