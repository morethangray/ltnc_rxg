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

# ----
fxn_histogram <- function(data, data_name){
    
    hist(data$value, breaks = 15, col = "gray", main = paste("Histogram of", data_name), xlab = "Value", ylab = "Count")
    
    if("value_std" %in% names(data)){
      hist(data$value_std, breaks = 25, col = "gray", main = paste("Histogram of Standardized", data_name), xlab = "Standardized Value", ylab = "Count")
    }
    
    if("value_log" %in% names(data)){
      hist(data$value_log, breaks = 25, col = "gray", main = paste("Histogram of Log-transformed", data_name), xlab = "Log Value", ylab = "Count")
    }
    
    if("value_sqrt" %in% names(data)){
      hist(data$value_sqrt, breaks = 25, col = "gray", main = paste("Histogram of Square Root-transformed", data_name), xlab = "Sqrt Value", ylab = "Count")
    }
    
  }
fxn_gof_test_rich <- function(data, data_name){
  
  gf_pois <- vcd::goodfit(data$value, type = "poisson")
  # Extract the first line (Likelihood Ratio goodness-of-fit test results)
  pois_first_line <- summary(gf_pois)[1, ]
  plot(gf_pois, main = "") 
  
  gf_nbin <- vcd::goodfit(data$value, type = "nbinomial")
  nbin_first_line <- summary(gf_nbin)[1, ]
  plot(gf_nbin,  main = "") 
  
}
fxn_gof_test_abun <- function(data, data_name){
  
  print(paste("Normal Distribution Fit for", data_name)) 
  shapiro.test(data$value)
  fitdist(data$value, dist = "norm") %>%
    plot(histo = TRUE, demp = TRUE)
  cat("\n")
  
  print(paste("Normal Distribution Fit for Standardized", data_name)) 
  shapiro.test(data$value_std)
  fitdist(data$value_std, dist = "norm") %>%
    plot(histo = TRUE, demp = TRUE)
  cat("\n")
  
  print(paste("Normal Distribution Fit for Log-transformed", data_name))
  shapiro.test(data$value_log)
  fitdist(data$value_log, dist = "norm") %>%
    plot(histo = TRUE, demp = TRUE)
  cat("\n")
  
  print(paste("Normal Distribution Fit for Square Root-transformed", data_name)) 
  shapiro.test(data$value_sqrt)
  fitdist(data$value_sqrt, dist = "norm") %>%
    plot(histo = TRUE, demp = TRUE)
  
  
}
fxn_vuong_test <- function(data, data_name){
  # Vuong tests (if applicable and if models converge)
  if (length(unique(data$treatment)) > 1) {
    tryCatch({ # Wrap in tryCatch to handle potential errors
      nb <- MASS::glm.nb(value ~ treatment, data = data)
      p <- glm(value ~ treatment, data = data, family = "poisson")
      print(paste("Vuong test: Negative Binomial vs. Poisson for", data_name))
      print(vuong(nb, p))
      cat("\n\n")
      
      if (min(data$value) == 0) {
        z <- pscl::zeroinfl(value ~ treatment, data = data)
        print(paste("Vuong test: Poisson vs. Zero-Inflated for", data_name))
        print(vuong(p, z))
        cat("\n\n")
        print(paste("Vuong test: Negative Binomial vs. Zero-Inflated for", data_name))
        print(vuong(nb, z))
      }
    }, error = function(e) {
      message(paste("Error in Vuong tests for", data_name, ":", e$message))
    })
  } else {
    message("Skipping Vuong tests as no treatment variable is present.")
  }
}
fxn_qq_plot <- function(data, data_name) {
  # Q-Q plots - Use mle2 for more robust Negative Binomial fitting
  data$value_1 <- data$value + 1 # Shifting values for better visualization
  
  # Poisson fit
  poisson <- fitdistr(data$value_1, "Poisson")
  qqp(data$value_1, "pois", lambda = poisson$estimate, main = paste("Poisson QQ plot for", data_name))
  
  tryCatch({
    nbinom_mle <- mle2(value_1 ~ dnegbin(mu = mu, theta = theta), 
                       start = list(mu = mean(data$value_1), theta = 1), 
                       data = data)
    qqp(data$value_1, 
        "nbinom", 
        size = coef(nbinom_mle)["theta"], 
        mu = coef(nbinom_mle)["mu"], 
        main = paste("Negative binomial QQ plot for", data_name))
  }, 
  error = function(e) {
    message(paste("Error fitting Negative Binomial with mle2 for", data_name, ":", e$message))
  })
  
  # Conditional additional fits for specific data_name
  if (data_name == "nonnative richness") {
    qqp(data$value_1, "norm", main = paste("Normal QQ plot for", data_name))
    qqp(data$value_1, "lnorm", main = paste("Lognormal QQ plot for", data_name))
    
    gamma <- fitdistr(data$value_1, "gamma")
    qqp(data$value_1, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]], main = paste("Gamma QQ plot for", data_name))
  }
}
# ========================================================== -----
# ========================================================== -----
# abundance model selection ----
# ========================================================== -----
# ========================================================== -----
# NATIVE SPECIES  ----
# ---------------------------------------------------------- -----
# Fit initial models  -----
# Create a bunch of potential models: glm, glmer, glm.nb, glmmTMB, lmer

#   value_std: normal (no random effect)
n0 <- lm(value_std ~ treatment, data = abun_nat) #lm
n0n <- lm(value_std ~ 1, data = abun_nat) #lm_null
#   value: normal
n1 <- lmer(value_std ~ treatment + (1 | plot_name), data = abun_nat, REML = FALSE) #lmer_1
n2 <- lmer(value_std ~ treatment + (1 + treatment | plot_name), data = abun_nat, REML = FALSE) #lmer_2
n1n <- lmer(value_std ~ 1 + (1 | plot_name), data = abun_nat, REML = FALSE) #lmer_null_1
n2n <- lmer(value_std ~ 1 + (1 + treatment | plot_name), data = abun_nat, REML = FALSE) #lmer_null_2

#   value_log: normal
n3 <- lmer(value_log ~ treatment + (1 | plot_name), data = abun_nat, REML = FALSE)
n4 <- lmer(value_log ~ treatment + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)
n3n <- lmer(value_log ~ 1 + (1 | plot_name), data = abun_nat, REML = FALSE)
n4n <- lmer(value_log ~ 1 + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)

#   value_sqrt: normal
n5 <- lmer(value_sqrt ~ treatment + (1 | plot_name), data = abun_nat, REML = FALSE)
n6 <- lmer(value_sqrt ~ treatment + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)
n5n <- lmer(value_sqrt ~ 1 + (1 | plot_name), data = abun_nat, REML = FALSE)
n6n <- lmer(value_sqrt ~ 1 + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)

# Review initial model performance ----
AIC(n0, n0n) # n0
model.sel(n1, n2, n1n, n2n) # n2, n1
model.sel(n3, n4, n3n, n4n) # n4, n4n
model.sel(n5, n6, n5n, n6n) # n6, n5

plotQQunif(n0) # significant deviation for KS test
plotQQunif(n2) # significant deviation for KS test
plotQQunif(n4) # looks ok
plotQQunif(n6) # significant deviation for KS test

best_initial <- n4

testDispersion(best_initial)
check_singularity(best_initial)
testUniformity(best_initial)
plotQQunif(best_initial)
testOutliers(best_initial)

plotResiduals(best_initial)
plotResiduals(best_initial, abun_nat$treatment)
plotResiduals(best_initial, abun_nat$f_year) # y4
plotResiduals(best_initial, abun_nat$f_break)
plotResiduals(best_initial, abun_nat$f_one_yr)
plotResiduals(best_initial, abun_nat$f_two_yr)
plotResiduals(best_initial, abun_nat$f_new)
plotResiduals(best_initial, abun_nat$grazer) # grazer

# Best initial: n4 (log-transformed)----
lmer(value_log ~ treatment + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)

# Define subset and model ----
model_abun_nat_init <-
  lmer(value_log ~ treatment + (1 + treatment | plot_name),
       data = abun, subset = met_sub == "abun_nat", REML = FALSE
  )

# Fit fixed effects  ----
abun_nat_tbl_fix_dredge <- fxn_fixed_dredge_by_model(index_model = "model_abun_nat_init")
abun_nat_tbl_fix_one <- fxn_fixed_one(index_model = "model_abun_nat_init")
abun_nat_tbl_fix_two <- fxn_fixed_two(index_model = "model_abun_nat_init")

bind_fixed <- bind_rows(abun_nat_tbl_fix_dredge, 
                        abun_nat_tbl_fix_one, 
                        abun_nat_tbl_fix_two) %>%
  arrange(aicc) 


#   Identify best fixed model ----
model_abun_nat_fix <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name) + f_year,
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)

# Plot the residuals for the variables in the model
plotResiduals(model_abun_nat_fix, abun_nat$treatment)
plotResiduals(model_abun_nat_fix, abun_nat$f_year)
# plotResiduals(model_abun_nat_fix, abun_nat$f_break)
# plotResiduals(model_abun_nat_fix, abun_nat$f_one_yr)
# plotResiduals(model_abun_nat_fix, abun_nat$f_two_yr)
testUniformity(model_abun_nat_fix)

# Fit random effects  ----
#   initial model with random effects ----
abun_nat_init_tbl_ran <- fxn_random(index_model = model_abun_nat_init)
#   Fixed effects model(s) with random effects ----
abun_nat_fix_tbl_ran <- fxn_random(index_model = model_abun_nat_fix)

model_abun_nat_ran <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name) + f_year,
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
#   Review random model performance ----
temp_model <- model_abun_nat_ran

testDispersion(temp_model)

check_singularity(temp_model)
testUniformity(temp_model)
testOutliers(temp_model)

# Plot the residuals for the variables in the model
plotResiduals(temp_model, abun_nat$treatment)
plotResiduals(temp_model, abun_nat$f_year)
# plotResiduals(temp_model, abun_nat$f_break)
# plotResiduals(temp_model, abun_nat$f_one_yr)
# plotResiduals(temp_model, abun_nat$f_two_yr)
# plotResiduals(temp_model, abun_nat$f_new)
# plotResiduals(temp_model, abun_nat$grazer)

#   Final model ----
model_abun_nat_final <-
  lmer(value_log ~ treatment + (1 + treatment | plot_name) + f_year,
       data = abun, subset = met_sub == "abun_nat", REML = FALSE
  )

# Plot the residuals for the variables in the model
plotResiduals(model_abun_nat_final, abun_nat$treatment)
plotResiduals(model_abun_nat_final, abun_nat$f_year)
# plotResiduals(model_abun_nat_final, abun_nat$f_break)
# plotResiduals(model_abun_nat_final, abun_nat$f_new)
# plotResiduals(model_abun_nat_final, abun_nat$grazer) 
# plotResiduals(model_abun_nat_final, abun_nat$f_one_yr)
# plotResiduals(model_abun_nat_final, abun_nat$f_two_yr)



# ========================================================== -----
# NATIVE FORB SPECIES  ----
# ---------------------------------------------------------- -----

# Fit initial models  -----
#   value_std: normal (no random effect)
n0 <- lm(value_std ~ treatment, data = a_frb)
n0n <- lm(value_std ~ 1, data = a_frb)
#   value: normal
n1 <- lmer(value_std ~ treatment + (1 | plot_name), data = a_frb, REML = FALSE)
n2 <- lmer(value_std ~ treatment + (1 + treatment | plot_name), data = a_frb, REML = FALSE)
n1n <- lmer(value_std ~ 1 + (1 | plot_name), data = a_frb, REML = FALSE)
n2n <- lmer(value_std ~ 1 + (1 + treatment | plot_name), data = a_frb, REML = FALSE)

#   value_log: normal
n3 <- lmer(value_log ~ treatment + (1 | plot_name), data = a_frb, REML = FALSE)
n4 <- lmer(value_log ~ treatment + (1 + treatment | plot_name), data = a_frb, REML = FALSE)
n3n <- lmer(value_log ~ 1 + (1 | plot_name), data = a_frb, REML = FALSE)
n4n <- lmer(value_log ~ 1 + (1 + treatment | plot_name), data = a_frb, REML = FALSE)

#   value_sqrt: normal
n5 <- lmer(value_sqrt ~ treatment + (1 | plot_name), data = a_frb, REML = FALSE)
n6 <- lmer(value_sqrt ~ treatment + (1 + treatment | plot_name), data = a_frb, REML = FALSE)
n5n <- lmer(value_sqrt ~ 1 + (1 | plot_name), data = a_frb, REML = FALSE)
n6n <- lmer(value_sqrt ~ 1 + (1 + treatment | plot_name), data = a_frb, REML = FALSE)

# Review initial model performance ----
AIC(n0, n0n) # n0
model.sel(n1, n2, n1n, n2n) # n2, n1
model.sel(n3, n4, n3n, n4n) # n4, n3
model.sel(n5, n6, n5n, n6n) # n6, n5

plotQQunif(n0) # significant deviation for KS test
plotQQunif(n2) # significant deviation for KS test, significant outlier(s)
plotQQunif(n4) # looks good
plotQQunif(n6) # significant deviation for KS test

best_initial <- n4

testDispersion(best_initial)
check_singularity(best_initial)
testUniformity(best_initial)
plotQQunif(best_initial)
testOutliers(best_initial) # one outlier near 0

plotResiduals(best_initial, a_frb$treatment)
plotResiduals(best_initial, a_frb$f_year) # y4
plotResiduals(best_initial, a_frb$f_break) # b0, b1
plotResiduals(best_initial, a_frb$f_one_yr) # o0, o1
plotResiduals(best_initial, a_frb$f_two_yr) # fail levene
plotResiduals(best_initial, a_frb$f_new) # n1
plotResiduals(best_initial, a_frb$grazer) # goat, sheep


# Best initial: n4 ----
lmer(value_log ~ treatment + (1 + treatment | plot_name), data = a_frb, REML = FALSE)

# Define subset and model ----
# abun_frb <- abun %>%
#   filter(met_sub %in% "abun_frb")
model_abun_frb_init <-
  lmer(value_log ~ treatment + (1 + treatment | plot_name),
       data = abun, subset = met_sub == "abun_frb", REML = FALSE
  )
# Fit fixed effects  ----
abun_frb_tbl_fix_dredge <- fxn_fixed_dredge_by_model(index_model = "abun_frb")
abun_frb_tbl_fix_one <- fxn_fixed_one(index_model = "abun_frb")
abun_frb_tbl_fix_two <- fxn_fixed_two(index_model = "abun_frb")

abun_frb_fix_yp <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name) +
    (f_year) * (plot_type),
  data = abun, subset = met_sub == "abun_frb", REML = FALSE
)

abun_frb_fix_y <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name) +
    (year),
  data = abun, subset = met_sub == "abun_frb", REML = FALSE
)

plotResiduals(abun_frb_fix_yp, abun_frb$treatment)
plotResiduals(abun_frb_fix_yp, abun_frb$f_year) # levene
plotResiduals(abun_frb_fix_yp, abun_frb$f_break)
plotResiduals(abun_frb_fix_yp, abun_frb$f_one_yr) # levene
plotResiduals(abun_frb_fix_yp, abun_frb$f_two_yr)
plotResiduals(abun_frb_fix_yp, abun_frb$f_new)
testUniformity(abun_frb_fix_yp)

# temp <- lmer(value_log ~ treatment + year + plot_type*f_break + (1 + treatment | plot_name),
#              data = abun, subset = met_sub == "abun_frb", REML = FALSE)
# plotResiduals(temp, abun_frb$treatment)
# plotResiduals(temp, abun_frb$f_year)  # levene
# plotResiduals(temp, abun_frb$f_break)
# plotResiduals(temp, abun_frb$f_one_yr)  # levene
# plotResiduals(temp, abun_frb$f_two_yr)
# plotResiduals(temp, abun_frb$f_new)
# plotResiduals(temp, abun_frb$grazer) # levene
# testUniformity(temp)

#   Identify best fixed model ----
model_abun_frb_fix <- lmer(value_log ~ treatment + (1 + treatment | plot_name) + f_year + plot_type,
                           data = abun, subset = met_sub == "abun_frb", REML = FALSE
)

plotResiduals(model_abun_frb_fix, abun_frb$treatment)
plotResiduals(model_abun_frb_fix, abun_frb$f_year) # levene
plotResiduals(model_abun_frb_fix, abun_frb$f_break)
plotResiduals(model_abun_frb_fix, abun_frb$f_one_yr) # levene
plotResiduals(model_abun_frb_fix, abun_frb$f_two_yr)
plotResiduals(model_abun_frb_fix, abun_frb$f_new)
testUniformity(model_abun_frb_fix)

# Fit random effects  ----
#   initial model with random effects ----
abun_frb_init_tbl_ran <- fxn_random(index_model = model_abun_frb_init)
#   Fixed effects model(s) with random effects ----
abun_frb_fix_tbl_ran <- fxn_random(index_model = model_abun_frb_fix)

model_abun_frb_ran <- lmer(value_log ~ treatment + (1 + treatment | plot_name) + f_year + plot_type,
                           data = abun, subset = met_sub == "abun_frb", REML = FALSE
)

#   Review ramdom model performance ----
mod_a_check <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name) +
    f_year +
    plot_type,
  data = abun_frb, REML = FALSE
)
temp_model <- mod_a_check

testDispersion(temp_model)

check_singularity(temp_model)
testUniformity(temp_model)
testOutliers(temp_model) # one outlier near 0

plotResiduals(temp_model, abun_frb$treatment)
plotResiduals(temp_model, abun_frb$f_year) # levene
plotResiduals(temp_model, abun_frb$f_break)
plotResiduals(temp_model, abun_frb$f_one_yr) # levene
plotResiduals(temp_model, abun_frb$f_two_yr)
plotResiduals(temp_model, abun_frb$f_new)
# plotResiduals(temp_model, abun_frb$grazer)

#   Final model ----
model_abun_frb_final <-
  lmer(value_log ~ treatment + (1 + treatment | plot_name) + f_year + plot_type,
       data = abun_frb, REML = FALSE
  )

model_abun_frb_final_2 <-
  lmer(value_log ~ treatment + (1 + treatment | plot_name) + f_year + (1 | plot_type),
       data = abun_frb, REML = FALSE
  )

plotResiduals(model_abun_frb_final, abun_frb$treatment)
plotResiduals(model_abun_frb_final, abun_frb$f_year) #
plotResiduals(model_abun_frb_final, abun_frb$f_break)
plotResiduals(model_abun_frb_final, abun_frb$f_new)
plotResiduals(model_abun_frb_final, abun_frb$grazer) #
plotResiduals(model_abun_frb_final, abun_frb$f_one_yr) #
plotResiduals(model_abun_frb_final, abun_frb$f_two_yr)

plotResiduals(model_abun_frb_final_2, abun_frb$treatment)
plotResiduals(model_abun_frb_final_2, abun_frb$f_year) #
plotResiduals(model_abun_frb_final_2, abun_frb$f_break)
plotResiduals(model_abun_frb_final_2, abun_frb$f_new)
plotResiduals(model_abun_frb_final_2, abun_frb$grazer) #
plotResiduals(model_abun_frb_final_2, abun_frb$f_one_yr) #
plotResiduals(model_abun_frb_final_2, abun_frb$f_two_yr)


# ========================================================== -----
# NON-NATIVE SPECIES  ----
# ---------------------------------------------------------- -----

# Fit initial models  -----
# Create a bunch of potential models: glm, glmer, glm.nb, glmmTMB, lmer
#   value_std: normal (no random effect)
n0 <- lm(value_std ~ treatment, data = a_non)
n0n <- lm(value_std ~ 1, data = a_non)
#   value: normal
n1 <- lmer(value_std ~ treatment + (1 | plot_name), data = a_non, REML = FALSE)
n2 <- lmer(value_std ~ treatment + (1 + treatment | plot_name), data = a_non, REML = FALSE)
n1n <- lmer(value_std ~ 1 + (1 | plot_name), data = a_non, REML = FALSE)
n2n <- lmer(value_std ~ 1 + (1 + treatment | plot_name), data = a_non, REML = FALSE)

#   value_log: normal
n3 <- lmer(value_log ~ treatment + (1 | plot_name), data = a_non, REML = FALSE)
n4 <- lmer(value_log ~ treatment + (1 + treatment | plot_name), data = a_non, REML = FALSE)
n3n <- lmer(value_log ~ 1 + (1 | plot_name), data = a_non, REML = FALSE)
n4n <- lmer(value_log ~ 1 + (1 + treatment | plot_name), data = a_non, REML = FALSE)

#   value_sqrt: normal
n5 <- lmer(value_sqrt ~ treatment + (1 | plot_name), data = a_non, REML = FALSE)
n6 <- lmer(value_sqrt ~ treatment + (1 + treatment | plot_name), data = a_non, REML = FALSE)
n5n <- lmer(value_sqrt ~ 1 + (1 | plot_name), data = a_non, REML = FALSE)
n6n <- lmer(value_sqrt ~ 1 + (1 + treatment | plot_name), data = a_non, REML = FALSE)

# Review initial model performance ----
AIC(n0, n0n) # n0
model.sel(n1, n2, n1n, n2n) # n1n (then n1)
model.sel(n3, n4, n3n, n4n) # n3n (then n4n)
model.sel(n5, n6, n5n, n6n) # n5n (then n6n)

plotQQunif(n0) # significant deviation for KS test
plotQQunif(n1n) # significant deviation for KS test
plotQQunif(n3n) # significant deviation for KS test
plotQQunif(n4n) # significant deviation for KS test
plotQQunif(n5n) # looks good
# y4, y6; b0, b1; o0, o1; n0, n1; sheep, goat; levene for f_two_yr

plotQQunif(n6n) # looks good
# one outlier near 1.0
# y4, y6; b0, b1; o0, o1; n0, n1; sheep, goat; levene for f_two_yr

best_initial <- n6n

testDispersion(best_initial)
check_singularity(best_initial)
testUniformity(best_initial)
testOutliers(best_initial)

plotResiduals(best_initial, a_non$treatment)
plotResiduals(best_initial, a_non$f_year)
plotResiduals(best_initial, a_non$f_break)
plotResiduals(best_initial, a_non$f_one_yr)
plotResiduals(best_initial, a_non$f_two_yr)
plotResiduals(best_initial, a_non$f_new)
plotResiduals(best_initial, a_non$grazer)

# Best initial: n6n ----
lmer(value_sqrt ~ 1 + (1 | plot_name), data = a_non, REML = FALSE)

# Define subset and model ----
# abun_non <- abun %>%
#   filter(met_sub %in% "abun_non")

model_abun_non_init <-
  lmer(value_sqrt ~ 1 + (1 | plot_name),
       data = abun, subset = met_sub == "abun_non", REML = FALSE
  )

# Fit fixed effects  ----
abun_non_tbl_fix_dredge <- fxn_fixed_dredge_by_model(index_model = "abun_non")
# Treatment not significant!

abun_non_tbl_fix_one <- fxn_fixed_one(index_model = "abun_non")
abun_non_tbl_fix_two <- fxn_fixed_two(index_model = "abun_non")

abun_non_fix_yp <- update(model_abun_non_init, . ~ . + f_year + plot_type)
abun_non_fix_yt <- update(model_abun_non_init, . ~ . + f_year + f_two_yr) # BEST
abun_non_fix_ytp <- update(model_abun_non_init, . ~ . + f_year + f_two_yr + plot_type)
abun_non_fix_yn <- update(model_abun_non_init, . ~ . + f_year + f_new) # VERY BAD
abun_non_fix_y <- update(model_abun_non_init, . ~ . + f_year)

plotResiduals(abun_non_fix_yp, abun_non$treatment)
plotResiduals(abun_non_fix_yp, abun_non$f_year) # y4
plotResiduals(abun_non_fix_yp, abun_non$f_break) # b2
plotResiduals(abun_non_fix_yp, abun_non$f_one_yr) # o0
plotResiduals(abun_non_fix_yp, abun_non$f_two_yr) # t1
testUniformity(abun_non_fix_yp)

plotResiduals(abun_non_fix_y, abun_non$treatment)
plotResiduals(abun_non_fix_y, abun_non$f_year) # y4
plotResiduals(abun_non_fix_y, abun_non$f_break) # b2
plotResiduals(abun_non_fix_y, abun_non$f_one_yr) # o0, o1
plotResiduals(abun_non_fix_y, abun_non$f_two_yr) # t1
testUniformity(abun_non_fix_y)

plotResiduals(abun_non_fix_yt, abun_non$treatment)
plotResiduals(abun_non_fix_yt, abun_non$f_year) # levene
plotResiduals(abun_non_fix_yt, abun_non$f_break)
plotResiduals(abun_non_fix_yt, abun_non$f_one_yr)
plotResiduals(abun_non_fix_yt, abun_non$f_two_yr)
testUniformity(abun_non_fix_yt)

# plotResiduals(abun_non_fix_yn, abun_non$treatment)  # VERY BAD all metrics
# plotResiduals(abun_non_fix_yn, abun_non$f_year)  # y5
# plotResiduals(abun_non_fix_yn, abun_non$f_break)
# plotResiduals(abun_non_fix_yn, abun_non$f_one_yr)
# plotResiduals(abun_non_fix_yn, abun_non$f_two_yr)
# testUniformity(abun_non_fix_yn)

#   Identify best fixed model ----
model_abun_non_fix <- lmer(value_sqrt ~ 1 + (1 | plot_name) + f_year + f_two_yr,
                           data = abun_non, REML = FALSE
)

# Fit random effects  ----
#   initial model with random effects ----
abun_non_init_tbl_ran <- fxn_random(index_model = model_abun_non_init)
#   Fixed effects model(s) with random effects ----
abun_non_fix_tbl_ran <- fxn_random(index_model = model_abun_non_fix)

model_abun_non_ran <- lmer(value_sqrt ~ 1 + (1 | plot_name) + f_year + f_two_yr,
                           data = abun_non, REML = FALSE
)

#   Review ramdom model performance ----
mod_a_check <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name) +
    f_year +
    plot_type,
  data = abun_non, REML = FALSE
)
temp_model <- mod_a_check

testDispersion(temp_model)

check_singularity(temp_model)
testUniformity(temp_model)
testOutliers(temp_model) # one outlier near 0

plotResiduals(temp_model, abun_non$treatment)
plotResiduals(temp_model, abun_non$f_year) # levene
plotResiduals(temp_model, abun_non$f_break)
plotResiduals(temp_model, abun_non$f_one_yr) # levene
plotResiduals(temp_model, abun_non$f_two_yr)
plotResiduals(temp_model, abun_non$f_new)
# plotResiduals(temp_model, abun_non$grazer)

#   Final model ----
model_abun_non_final <- lmer(value_sqrt ~ 1 + (1 | plot_name) + f_year + f_two_yr,
                             data = abun_non, REML = FALSE
)

# ========================================================== -----

# fit_and_select_model ----
# Function to perform model fitting and selection (generalized)
fit_and_select_model <- function(data, response_var, data_name) {
  
  # Data transformations (store in a list for easier access)
  transformed_data <- list(
    std = data[[response_var]], # Original scale
    log = if(min(data[[response_var]]) > 0) log(data[[response_var]]) else data[[response_var]], # Log transformation (check for 0s)
    sqrt = sqrt(data[[response_var]])
  )
  
  # Fit initial models (using a loop and list for organization)
  initial_models <- list()
  transforms <- c("std", "log", "sqrt")
  
  for (trans in transforms) {
    initial_models[[glue::glue("n0_{trans}")]] <- lm(transformed_data[[trans]] ~ treatment, data = data)
    initial_models[[glue::glue("n0n_{trans}")]] <- lm(transformed_data[[trans]] ~ 1, data = data)
    for (i in 1:2) { # Fit random intercept and slope models
      for (j in 1:2) {
        name <- glue::glue("n{i*2 + j - 1}_{trans}")
        formula <- glue::glue("{trans} ~ treatment + {ifelse(i==2, '(1 + treatment | plot_name)', '(1 | plot_name)')}")
        initial_models[[name]] <- lmer(as.formula(formula), data = data, REML = FALSE)
        name_n <- glue::glue("n{i*2 + j - 1}n_{trans}")
        formula_n <- glue::glue("{trans} ~ 1 + {ifelse(i==2, '(1 + treatment | plot_name)', '(1 | plot_name)')}")
        initial_models[[name_n]] <- lmer(as.formula(formula_n), data = data, REML = FALSE)
      }
    }
  }
  
  
  
  # Review initial model performance (using lapply for conciseness)
  aics <- lapply(initial_models, AIC)
  print(glue::glue("AIC values for {data_name}:"))
  print(aics)
  
  
  best_initial_transform <- character()
  best_initial_model <- list()
  for (trans in transforms){
    model_subset <- initial_models[grep(trans, names(initial_models))]
    model_selection_results <- MuMIn::model.sel(model_subset)
    
    best_model_name <- rownames(model_selection_results)[1]
    best_initial_model[[trans]] <- model_subset[[best_model_name]]
    best_initial_transform[trans] <- trans
    
    print(glue::glue("Model Selection for {data_name} with {trans} transform"))
    print(model_selection_results)
  }
  
  
  qq_plots <- lapply(initial_models, DHARMa::plotQQunif)
  print(qq_plots) # Examine these visually
  
  # Select best initial model (based on AIC and QQ plot inspection)
  # This will need to be done manually based on the outputs of the model selection and QQ plots above.
  # The code below is a placeholder and should be modified.
  best_initial <- best_initial_model[[best_initial_transform[["log"]]]] # Example – change as needed!
  best_transform <- best_initial_transform[["log"]]
  
  # ... (rest of your model fitting and selection process, but now using best_initial)
  
  return(list(best_initial = best_initial, best_transform = best_transform))
  
}

# About the fit_and_select_model function ----
# Function for Model Fitting: fit_and_select_model function. This function takes the data, response variable name, and a data name (for printing) as arguments. It performs the initial model fitting and selection process.
# 
# Data Transformation Handling: The function handles data transformations (log and sqrt) systematically. It stores the transformed data in a list, making it accessible within the function.  The log transformation checks for 0s in the data and only applies the transformation if there are no 0s.  If there are 0s, the original data is used.  This prevents errors when trying to take the log of 0.
# 
# Looping for Model Fitting: Nested loops iterate through the transformations and random effect structures. This reduces code duplication when fitting models.
# 
# Model Selection with MuMIn::model.sel: This function from the MuMIn package provides a clean way to compare multiple models and get AIC values.
# 
# QQ Plot Generation with DHARMa::plotQQunif: This function provides a more robust way to generate QQ plots for GLMMs.
# 
# Residual Plotting with ggplot2: It's good practice to use ggplot2 for more customizable and informative residual plots.  You can add facets, loess smoothers, etc.  I've included a basic example that can be adapted.
# 
# Clarity and Comments: I've used descriptive variable names and comments explain the code's logic and purpose.   
# 
# glue::glue() for String Formatting: This function from the glue package helps create strings dynamically, especially for formulas and model names within the loops.
# 
# Placeholder for Best Model Selection: The code prints the AIC values and generates QQ plots. However, you still need to manually inspect these outputs to choose the best initial model. I've added a placeholder where you'll need to insert your logic for choosing the best_initial model based on AIC and visual inspection of QQ plots.  This is a critical step that cannot be easily automated, as it involves your expert judgment.
# 
# Updating the Best Initial Model: The code includes an example of how to update the best initial model with additional fixed effects using the update() function.  This keeps the code concise and avoids refitting the entire model.

# ========================================================== -----
# Function to fit initial models ----
fxn_initial_model <- function(data) {
  responses <- c("value_std", "value_log", "value_sqrt")
  best_models <- list()
  
  # Create empty data frame for summary of all models
  summary_df <- data.frame(
    response = character(),
    model_name = character(),
    aic = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (response in responses) {
    models <- list(
      lm = lm(as.formula(paste(response, "~ treatment")), data = data),
      lm_null = lm(as.formula(paste(response, "~ 1")), data = data),
      lmer_1 = lmer(as.formula(paste(response, "~ treatment + (1 | plot_name)")), data = data, REML = FALSE),
      lmer_2 = lmer(as.formula(paste(response, "~ treatment + (1 + treatment | plot_name)")), data = data, REML = FALSE),
      lmer_null_1 = lmer(as.formula(paste(response, "~ 1 + (1 | plot_name)")), data = data, REML = FALSE),
      lmer_null_2 = lmer(as.formula(paste(response, "~ 1 + (1 + treatment | plot_name)")), data = data, REML = FALSE)
    )
    
    # Calculate AIC values for all models
    aic_values <- AIC(models$lm, models$lm_null)
    rownames(aic_values) <- c("lm", "lm_null")
    mixed_model_selection <- model.sel(lmer_1 = models$lmer_1, 
                                       lmer_2 = models$lmer_2,
                                       lmer_null_1 = models$lmer_null_1,
                                       lmer_null_2 = models$lmer_null_2) 
    
    # Add linear models to summary
    lm_summary <- data.frame(
      response = response,
      model_name = rownames(aic_values),
      aic = aic_values$AIC,
      stringsAsFactors = FALSE
    )
    
    # Add mixed models to summary
    mixed_summary <- data.frame(
      response = response,
      model_name = rownames(mixed_model_selection),
      aic = mixed_model_selection$AICc,
      stringsAsFactors = FALSE
    )
    
    # Combine all models for this response
    response_summary <- rbind(lm_summary, mixed_summary)
    summary_df <- rbind(summary_df, response_summary)
    
    # Find best model for this response
    if (min(aic_values$AIC) < min(mixed_model_selection$AICc)) {
      best_model_name <- names(models)[which.min(aic_values$AIC)]
      best_model <- models[[best_model_name]]
      best_aic <- min(aic_values$AIC)
    } else {
      best_model_name <- rownames(mixed_model_selection)[which.min(mixed_model_selection$AICc)]
      best_model <- models[[best_model_name]]
      best_aic <- min(mixed_model_selection$AICc)
    }
    
    best_models[[response]] <- list(
      response = response,
      model = best_model,
      model_name = best_model_name,
      formula = formula(best_model),
      aic = best_aic
    )
  }
  
  best_response <- names(best_models)[which.min(sapply(best_models, function(x) x$aic))]
  best_overall <- best_models[[best_response]]
  
  # Return both the best model and the complete summary table
  return(list(
    best_model = best_overall,
    summary_table = summary_df %>%
      arrange(aic)
  ))
}

# Function to review model performance ----
fxn_model_review <- function(model, data) {
  model_terms <- attr(terms(model), "term.labels")
  
  # Diagnostic checks
  testDispersion(model)
  check_singularity(model)
  testUniformity(model)
  testOutliers(model)
  
  # Residual plots only for included terms
  for (term in model_terms) {
    if (term %in% names(data)) {
      plotResiduals(model, data[[term]])
    }
  }
}


# Function to review model performance with combined plots [NOT WORKING]----
fxn_model_review_grid <- function(model, data) {
  # Simulated residuals
  sim_res <- simulateResiduals(model)
  
  # Generate diagnostic plots (using as.ggplot to convert base R plots)
  p1 <- ggplotify::as.ggplot(function() plot(sim_res) + ggtitle("Simulated Residuals"))
  p2 <- ggplotify::as.ggplot(function() plotResiduals(sim_res, rank = TRUE)) + ggtitle("Ranked Residuals") # Convert to ggplot
  p3 <- ggplotify::as.ggplot(function() plotResiduals(sim_res)) + ggtitle("Residuals vs Predicted") # Convert to ggplot
  p4 <- ggplotify::as.ggplot(function() qqnorm(residuals(model)) + qqline()) + ggtitle("Normal Q-Q") # Convert to ggplot
  
  # Combine main diagnostic plots (2x2 layout)
  diagnostic_plot <- p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 2)
  
  # Print combined diagnostic plots
  print(diagnostic_plot)
  
  # Residual plots for included terms (using as.ggplot here too)
  model_terms <- attr(terms(model), "term.labels")
  residual_plots <- list()
  for (term in model_terms) {
    if (term %in% names(data)) {
      residual_plots[[term]] <- ggplotify::as.ggplot(function() plotResiduals(sim_res, data[[term]])) + 
        ggtitle(paste("Residuals vs", term)) # Convert to ggplot
    }
  }
  
  # Print residual plots if available, also in a 2x2 grid
  if (length(residual_plots) > 0) {
    print(patchwork::wrap_plots(residual_plots, ncol = 2))
  }
}



# Function to fit and select fixed effects models ----
input_model <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)

fxn_fixed_effects <- function(input_model) {
  # Create empty data frame for summary of all models
  summary_df <- data.frame(
    model_name = character(),
    terms_count = integer(), # Number of terms
    terms = character(),
    aic = numeric(),
    model = I(list()) # Store model as a list column
  )
  # Create an empty list to store summary information
  model_summaries <- list()
  
  # One-term models
  model_terms_one <- c("plot_type", "f_break", "f_new", "f_one_yr", "f_two_yr", "f_year")
  fix_one <- lapply(model_terms_one, function(term) {
    model <- update(input_model, as.formula(paste(". ~ . +", term)))
    aic <- AIC(model)
    data.frame(
      model_name = term,
      terms_count = 1,
      terms = term,
      aic = aic,
      model = I(list(model))
    )
  }) %>% bind_rows()
  
  # Two-term models
  model_terms_two <- list(
    c("f_year", "plot_type"), c("f_year", "f_break"), c("f_year", "f_new"),
    c("f_year", "f_one_yr"), c("f_year", "f_two_yr")
  )
  
  fix_two <- lapply(model_terms_two, function(terms) {
    model <- update(input_model, as.formula(paste(". ~ . +", paste(terms, collapse = " + "))))
    aic <- AIC(model)
    model_name <- paste(terms, collapse = " + ")
    data.frame(
      model_name = model_name,
      terms_count = 2,
      terms = model_name, # Store combined term name
      aic = aic,
      model = I(list(model))
    )
  }) %>% bind_rows()
  
  # Combine all models
  all_models <- bind_rows(fix_one, fix_two)
  
  # Best overall model
  best_model_row <- all_models %>%
    arrange(aic) %>%
    slice(1)
  
  best_model <- best_model_row$model[[1]]
  best_model_name <- best_model_row$model_name
  best_aic <- best_model_row$aic
  best_terms <- best_model_row$terms
  best_terms_count <- best_model_row$terms_count
  
  # Return both the best model and the complete summary table
  return(list(
    best_model = list(
      model = best_model,
      model_name = best_model_name,
      terms = best_terms,
      terms_count = best_terms_count,
      formula = formula(best_model),
      aic = best_aic
    ),
    summary_table = all_models %>% arrange(aic) # Sort summary table by AIC
  ))
}


# Function to fit and select random effects models ----
fxn_random_effects <- function(best_fixed_model) {
  random_terms <- c(
    "(1 | plot_type)", "(1 | f_year)", "(1 | f_year / grazer)", "(1 | grazer)",
    "(1 | f_break)", "(1 | f_one_yr)", "(1 | f_two_yr)", "(1 | f_new)",
    "(1 + f_year | plot_type)", "(1 + treatment | f_year)", "(1 + treatment | f_year / grazer)",
    "(1 + treatment | grazer)", "(1 + treatment | f_break)", "(1 + treatment | f_one_yr)",
    "(1 + treatment | f_two_yr)", "(1 + treatment | f_new)"
  )
  
  models <- lapply(random_terms, function(term) update(best_fixed_model, as.formula(paste(". ~ . +", term))))
  model.sel(best_fixed_model, models)
}

# Function to fit final model ----
fit_final_model <- function(data, response) {
  lmer(
    as.formula(paste(response, "~ treatment + (1 + treatment | plot_name) + f_year + plot_type")), 
    data = data, REML = FALSE
  )
}


# ========================================================== -----

# Function to fit and assess random effects  ----
best_fixed_model = best_model_fixed_abun_nat
fxn_random_effects <- function(best_fixed_model) {
  
  # Add each random intercept one at a time
  plot_type <- update(best_fixed_model, . ~ . + (1 | plot_type))
  f_year <- update(best_fixed_model, . ~ . + (1 | f_year))
  f_year_grazer <- update(best_fixed_model, . ~ . + (1 | f_year/grazer))
  grazer <- update(best_fixed_model, . ~ . + (1 | grazer))
  f_break <- update(best_fixed_model, . ~ . + (1 | f_break))
  f_one_yr <- update(best_fixed_model, . ~ . + (1 | f_one_yr))
  f_two_yr <- update(best_fixed_model, . ~ . + (1 | f_two_yr))
  f_new <- update(best_fixed_model, . ~ . + (1 | f_new))
  # 
  # # Add each random slope one at a time
  f_year_plot_type <- update(best_fixed_model, . ~ . + (1 + f_year | plot_type))
  treatment_f_year <- update(best_fixed_model, . ~ . + (1 + treatment | f_year))
  treatment_f_year_grazer <- update(best_fixed_model, . ~ . + (1 + treatment | f_year/grazer))
  treatment_grazer <- update(best_fixed_model, . ~ . + (1 + treatment | grazer))
  treatment_f_break <- update(best_fixed_model, . ~ . + (1 + treatment | f_break))
  treatment_f_one_yr <- update(best_fixed_model, . ~ . + (1 + treatment | f_one_yr))
  treatment_f_two_yr  <- update(best_fixed_model, . ~ . + (1 + treatment | f_two_yr))
  treatment_f_new <- update(best_fixed_model, . ~ . + (1 + treatment | f_new))
  
  
  # Check for singularity using a robust method 
  list_model_names <- c("best_fixed_model", 
                        "plot_type", 
                        "f_year", 
                        "f_year_grazer", 
                        "grazer", 
                        "f_break", 
                        "f_one_yr",
                        "f_two_yr",
                        "f_new",
                        "f_year_plot_type",
                        "treatment_f_year",
                        "treatment_f_year_grazer",
                        "treatment_grazer",
                        "treatment_f_break",
                        "treatment_f_one_yr",
                        "treatment_f_two_yr",
                        "treatment_f_new")
  
  datalist <- list()
  for(n in list_model_names){
    
    model <- get(n)
    
    datalist[[n]] <- tibble(
      model_name = n, 
      is_singular = isSingular(model)
    )
  }
  
  bind_datalist <- do.call(bind_rows, datalist)
  
  not_singular <-  bind_datalist %>%
    filter(is_singular == FALSE) %>%
    filter(model_name != "best_fixed_model") 
  
  if(nrow(not_singular) == 0){
    print("All random effect models are singular")
  }else{
    print(not_singular)
  }
  
  all_models <- model.sel(best_fixed_model, 
                          plot_type, 
                          f_year, 
                          f_year_grazer, 
                          grazer, 
                          f_break, 
                          f_one_yr,
                          f_two_yr,
                          f_new,
                          f_year_plot_type,
                          treatment_f_year,
                          treatment_f_year_grazer,
                          treatment_grazer,
                          treatment_f_break,
                          treatment_f_one_yr,
                          treatment_f_two_yr,
                          treatment_f_new)  
  
  names(all_models)
  
  all_models_tbl <-   all_models %>%
    as_tibble() %>%
    # Add rownames as a column since they often contain model names
    tibble::rownames_to_column("model")
  
  # To see what columns are being dropped:
  dropped_cols <- setdiff(names(all_models), names(all_models_tbl))
  
  summary_table <- tibble(
    model_name = rownames(all_models),
    AICc = all_models$AICc,
    delta_AICc = all_models$delta,
    weight = round(all_models$weight, 4),
    df = all_models$df
  ) %>%
    arrange(AICc) %>%
    left_join(bind_datalist, "model_name")
}



fxn_select_random_effects <- function(model, data) {
  # Fit models with random effects
  ran_init <- fxn_random(index_model = model)
  ran_fix <- fxn_random(index_model = model)
  
  # Fit the best random effects model (assumed same formula as input)
  model_ran <- update(model, REML = FALSE)
  
  # Assess model diagnostics
  testDispersion(model_ran)
  check_singularity(model_ran)
  testUniformity(model_ran)
  testOutliers(model_ran)
  
  # Residual analysis for included terms
  included_terms <- names(fixef(model_ran))
  lapply(included_terms, function(term) plotResiduals(model_ran, data[[term]]))
  
  return(model_ran)
}

# Function to finalize model selection and review residuals  ----
fxn_final_model_review <- function(model, data) {
  # Fit the final model (assumed same formula as input)
  final_model <- update(model, REML = FALSE)
  
  # Residual analysis for included terms
  included_terms <- names(fixef(final_model))
  lapply(included_terms, function(term) plotResiduals(final_model, data[[term]]))
  
  return(final_model)
}

# Function to fit final model ----
fit_final_model <- function(data, response) {
  lmer(as.formula(paste(response, "~ treatment + (1 + treatment | plot_name) + f_year + plot_type")), 
       data = data, REML = FALSE)
}
# --- Running the workflow ---- 
# Apply functions to native species ----
# Fit the models with transformations
model_init_abun_nat <- fxn_initial_model(abun_nat)

# lm: lm(abundance ~ treatment)
# lm_null: lm(abundance ~ 1)
# lmer_1: lmer(abundance ~ treatment + (1 | plot_name))
# lmer_2: lmer(abundance ~ treatment + (1 + treatment | plot_name))
# lmer_null_1: lmer(abundance ~ 1 + (1 | plot_name))
# lmer_null_2: lmer(abundance ~ 1 + (1 + treatment | plot_name))

# View the summary table
print(model_init_abun_nat$summary_table)

# Review model diagnostics
# Based on AIC the standardized values had a better fit
# However, model diagnostics for this response didn't look good
fxn_model_review(model_init_abun_nat$best_model$model, abun_nat)

# Next best AIC set was value_log, These model diagnostics looked better but not perfect
best_model_init_abun_nat_log_2 <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
best_model_init_abun_nat_log_1 <- lmer(
  value_log ~ treatment + (1 | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
fxn_model_review(best_model_init_abun_nat_log_1, abun_nat)
fxn_model_review(best_model_init_abun_nat_log_2, abun_nat)

# Also looked at square root transformed diagnostics; super skewed

best_model_init_abun_nat_sqrt_2 <- lmer(
  value_sqrt ~ treatment + (1 + treatment | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
best_model_init_abun_nat_sqrt_1 <- lmer(
  value_sqrt ~ treatment + (1 | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
fxn_model_review(best_model_init_abun_nat_sqrt_1, abun_nat)
fxn_model_review(best_model_init_abun_nat_sqrt_2, abun_nat)

# After comparing model diagnostics for the three transformations I used log-transformed values for subsequent model selection
# Get the best overall model
best_model_init_abun_nat <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name),
  data = abun, subset = met_sub == "abun_nat", REML = FALSE
)




# Fit and select fixed effects model
best_fixed_abun_nat_dredge <- fxn_fixed_effects("best_model_abun_nat", method = "dredge")  # Choose "one" or "two" as needed
best_fixed_abun_nat_one <- fxn_fixed_effects("best_model_abun_nat", method = "one")  # Choose "one" or "two" as needed
best_fixed_abun_nat_two <- fxn_fixed_effects("best_model_abun_nat", method = "two")  # Choose "one" or "two" as needed

best_fixed_abun_nat

# Review model diagnostics
fxn_model_review(best_fixed_abun_nat$mixed_model_selection, abun_nat)

# Fit and select random effects model
best_random_abun_nat <- fxn_random_effects(best_fixed_abun_nat)
# Review model diagnostics
fxn_model_review(best_random_abun_nat$mixed_model_selection, abun_nat)

# Fit final model
final_model_abun_nat <- fit_final_model(abun_nat, response = "abun_nat")
# Review model diagnostics
fxn_model_review(final_model_abun_nat$mixed_model_selection, abun_nat)



nat_models <- define_models(abun_nat, "value_log")

# Initial model selection
nat_selection <- select_best_model(nat_models)
best_model_nat <- nat_models$lmer_2   
nat_assessment <- assess_model(best_model_nat, abun_nat)

model_abun_nat_fix <- fxn_select_fixed_effects(index_model = "model_abun_nat_init", data = abun)
model_abun_nat_ran <- fxn_select_random_effects(model = model_abun_nat_fix, data = abun)
model_abun_nat_final <- fxn_final_model_review(model = model_abun_nat_ran, data = abun)

final_model_nat <- fit_final_model(data_nat, "value_log")


# Apply functions to native forb species ----
frb_models <- define_models(abun_frb, "value_log")
frb_selection <- select_best_model(frb_models)
best_model_frb <- frb_models$lmer_2   
frb_assessment <- assess_model(best_model_frb, data_frb)
final_model_frb <- fit_final_model(data_frb, "value_log")

# Apply functions to non-native species ----
non_models <- define_models(abun_non, "value_log")
non_selection <- select_best_model(non_models)
best_model_non <- non_models$lmer_2   
non_assessment <- assess_model(best_model_non, data_non)
final_model_non <- fit_final_model(data_non, "value_log")

# ========================================================== -----