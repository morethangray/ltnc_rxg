# ========================================================== -----
fxn_gof_test_rich <- function(data, data_name){
  
  gf_pois <- vcd::goodfit(data$value, type = "poisson")
  # Extract the first line (Likelihood Ratio goodness-of-fit test results)
  pois_first_line <- summary(gf_pois)[1, ]
  plot(gf_pois, main = "") 
  
  gf_nbin <- vcd::goodfit(data$value, type = "nbinomial")
  nbin_first_line <- summary(gf_nbin)[1, ]
  plot(gf_nbin,  main = "") 
  
  if(data_name == "non-native richness"){
    
    print(paste("Square root-transformed", data_name)) 
    gf_pois_sqrt <- vcd::goodfit(data$value_sqrt, type = "poisson")
    # Extract the first line (Likelihood Ratio goodness-of-fit test results)
    pois_sqrt_first_line <- summary(gf_pois_sqrt)[1, ]
    plot(pois_sqrt_first_line, main = "") 
    
    print(paste("Log-transformed", data_name)) 
    
    gf_pois_log <- vcd::goodfit(data$value_log, type = "poisson")
    # Extract the first line (Likelihood Ratio goodness-of-fit test results)
    pois_log_first_line <- summary(gf_pois_log)[1, ]
    plot(pois_log_first_line, main = "") 
    
  }
  
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
  
  # fitdistrplus - Corrected plotting
  
  print(paste("Normal Distribution Fit for", data_name)) 
  fit_norm <- fitdist(data$value, dist = "norm")
  print(plot(fit_norm, histo = FALSE, demp = TRUE))
  
  print(paste("Normal Distribution Fit for Log-transformed", data_name))
  fit_lognorm <- fitdist(data$value_log, dist = "norm")
  plot(fit_lognorm, histo = FALSE, demp = TRUE)
  
  print(paste("Normal Distribution Fit for Square Root-transformed", data_name)) 
  fit_sqrtnorm <- fitdist(data$value_sqrt, dist = "norm")
  plot(fit_sqrtnorm, histo = FALSE, demp = TRUE)
  
  print(paste("Logistic Distribution Fit for", data_name))
  fit_logis <- fitdist(data$value, dist = "logis")
  plot(fit_logis, histo = FALSE, demp = TRUE)
  
  print(paste("Logistic Distribution Fit for Log-transformed", data_name))
  fit_loglogis <- fitdist(data$value_log, dist = "logis")
  plot(fit_loglogis, histo = FALSE, demp = TRUE)
  
  print(paste("Logistic Distribution Fit for Square Root-transformed", data_name))
  fit_sqrtlogis <- fitdist(data$value_sqrt, dist = "logis")
  plot(fit_sqrtlogis, histo = FALSE, demp = TRUE)
  
  invisible(data)
}

# Negative Binomial fit
nbinom <- fitdistr(data$value_1, "Negative Binomial")
qqp(data$value_1, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]], main = paste("Negative binomial QQ plot for", data_name))

tryCatch({
  nbinom_mle <- mle2(value_1 ~ dnegbin(mu = mu, theta = theta), 
                     start = list(mu = mean(data$value_1), theta = 1), 
                     data = data)
  qqp(data$value_1, 
      "nbinom", 
      size = coef(nbinom_mle)["theta"], 
      mu = coef(nbinom_mle)["mu"], 
      main = paste("Negative Binomial Q-Q Plot for", data_name))
}, 
error = function(e) {
  message(paste("Error fitting Negative Binomial with mle2 for", data_name, ":", e$message))
})
# Function to perform distribution checks for count data (Richness)
check_count_distribution <- function(data, data_name) {
  message(paste("Checking distribution for:", data_name))
  
  # Basic histogram
  hist(data$value, breaks = 15, col = "gray", main = paste("Histogram of", data_name), xlab = "Value", ylab = "Count")
  
  # Goodness-of-fit tests
  gf_pois <- goodfit(data$value, type = "poisson")
  print(paste("Poisson Goodness of Fit for", data_name))
  print(summary(gf_pois))
  plot(gf_pois, main = paste("Poisson Goodness of Fit Plot for", data_name))
  
  gf_nbin <- goodfit(data$value, type = "nbinomial")
  print(paste("Negative Binomial Goodness of Fit for", data_name))
  print(summary(gf_nbin))
  plot(gf_nbin, main = paste("Negative Binomial Goodness of Fit Plot for", data_name))
  
  # Vuong tests (if applicable and if models converge)
  if (length(unique(data$treatment)) > 1) {
    tryCatch({ # Wrap in tryCatch to handle potential errors
      nb <- glm.nb(value ~ treatment, data = data)
      p <- glm(value ~ treatment, data = data, family = "poisson")
      print(paste("Vuong test: Negative Binomial vs. Poisson for", data_name))
      print(vuong(nb, p))
      
      if (min(data$value) == 0) {
        z <- zeroinfl(value ~ treatment, data = data)
        print(paste("Vuong test: Poisson vs. Zero-Inflated for", data_name))
        print(vuong(p, z))
        print(paste("Vuong test: Negative Binomial vs. Zero-Inflated for", data_name))
        print(vuong(nb, z))
      }
    }, error = function(e) {
      message(paste("Error in Vuong tests for", data_name, ":", e$message))
    })
  } else {
    message("Skipping Vuong tests as no treatment variable is present.")
  }
  
  # Q-Q plots - Use mle2 for more robust Negative Binomial fitting
  data$value_1 <- data$value + 1 # Shifting values for better visualization
  tryCatch({
    nbinom_mle <- mle2(value_1 ~ dnegbin(mu = mu, theta = theta), start = list(mu = mean(data$value_1), theta = 1), data = data)
    qqp(data$value_1, "nbinom", size = coef(nbinom_mle)["theta"], mu = coef(nbinom_mle)["mu"], main = paste("Negative Binomial Q-Q Plot for", data_name))
  }, error = function(e) {
    message(paste("Error fitting Negative Binomial with mle2 for", data_name, ":", e$message))
  })
  
  qqp(data$value_1, "norm", main = paste("Normal Q-Q Plot for", data_name))
  qqp(data$value_1, "lnorm", main = paste("Lognormal Q-Q Plot for", data_name))
  
  poisson <- fitdistr(data$value_1, "Poisson")
  qqp(data$value_1, "pois", lambda = poisson$estimate, main = paste("Poisson Q-Q Plot for", data_name))
  
  gamma <- fitdistr(data$value_1, "gamma")
  qqp(data$value_1, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]], main = paste("Gamma Q-Q Plot for", data_name))
  
  invisible(data)
}
# ---------------------------------------------------------- -----
# ========================================================== -----
