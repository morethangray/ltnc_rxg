# ========================================================== -----
# abundance model selection ----
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
