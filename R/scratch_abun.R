# ========================================================== -----
# abundance model selection ----
# ---------------------------------------------------------- -----
# Functions ----
#   fxn_fixed_dredge_by_model ----
fxn_fixed_dredge_by_model <- function(index_model) {
  input_model <- get(index_model)
  
  f0 <- update(input_model, . ~ . +
                 plot_type +
                 f_year +
                 # grazer +
                 f_break +
                 f_new +
                 f_one_yr +
                 f_two_yr)
  
  dd <- dredge(f0)
  
  aic_fixed <-
    subset(dd, delta < 6) %>%
    as_tibble(rownames = "model_name") %>%
    mutate(
      rank = 1:n(),
      input_model = index_model
    ) %>%
    clean_names() %>%
    dplyr::relocate(input_model,
                    rank,
                    treatment,
                    f_one_yr,
                    f_two_yr,
                    f_break,
                    f_new,
                    plot_type,
                    f_year,
                    # grazer,
                    delta,
                    weight,
                    aicc = ai_cc,
                    df,
                    model_name
    )
}

#   fxn_fixed_one ----
# index_subset <- "n_val_i"
fxn_fixed_one <- function(index_subset) {
  input_model <- get(paste0(index_subset, "_red"))
  
  f01 <- update(input_model, . ~ . + plot_type)
  # f02 <- update(input_model, . ~ . + grazer)
  f03 <- update(input_model, . ~ . + f_break)
  f04 <- update(input_model, . ~ . + f_new)
  f05 <- update(input_model, . ~ . + f_one_yr)
  f06 <- update(input_model, . ~ . + f_two_yr)
  f07 <- update(input_model, . ~ . + f_year)
  
  model.sel(
    f01,
    # f02,
    f03,
    f04,
    f05, f06, f07
  ) %>%
    as_tibble(rownames = "model_name") %>%
    clean_names() %>%
    mutate(
      input_model = index_subset,
      rank = 1:n()
    ) %>%
    dplyr::relocate(input_model,
                    rank,
                    treatment,
                    f_one_yr,
                    f_two_yr,
                    f_break,
                    f_new,
                    plot_type,
                    f_year,
                    # grazer,
                    delta,
                    weight,
                    aicc = ai_cc,
                    df,
                    model_name
    )
}
#   fxn_fixed_two ----
fxn_fixed_two <- function(index_subset) {
  input_model <- get(paste0(index_subset, "_red"))
  
  f08 <- update(input_model, . ~ . + f_year + plot_type)
  # f09 <- update(input_model, . ~ . + f_year + grazer)
  f10 <- update(input_model, . ~ . + f_year + f_break)
  f11 <- update(input_model, . ~ . + f_year + f_new)
  f12 <- update(input_model, . ~ . + f_year + f_one_yr)
  f13 <- update(input_model, . ~ . + f_year + f_two_yr)
  
  model.sel(
    f08,
    # f09,
    f10, f11, f12, f13
  ) %>%
    as_tibble(rownames = "model_name") %>%
    clean_names() %>%
    mutate(
      input_model = index_subset,
      rank = 1:n()
    ) %>%
    dplyr::relocate(input_model,
                    rank,
                    treatment,
                    f_one_yr,
                    f_two_yr,
                    f_break,
                    f_new,
                    plot_type,
                    f_year,
                    # grazer,
                    delta,
                    weight,
                    aicc = ai_cc,
                    df,
                    model_name
    )
}


#   fxn_random ----
fxn_random <- function(index_model) {
  best_fixed <- index_model
  
  # Add each random intercept one at a time
  m11 <- update(best_fixed, . ~ . + (1 | plot_type))
  m12 <- update(best_fixed, . ~ . + (1 | f_year))
  m13 <- update(best_fixed, . ~ . + (1 | f_year / grazer))
  m14 <- update(best_fixed, . ~ . + (1 | grazer))
  m15 <- update(best_fixed, . ~ . + (1 | f_break))
  m16 <- update(best_fixed, . ~ . + (1 | f_one_yr))
  m17 <- update(best_fixed, . ~ . + (1 | f_two_yr))
  m18 <- update(best_fixed, . ~ . + (1 | f_new))
  #
  # # Add each random slope one at a time
  m21 <- update(best_fixed, . ~ . + (1 + f_year | plot_type))
  m22 <- update(best_fixed, . ~ . + (1 + treatment | f_year))
  m23 <- update(best_fixed, . ~ . + (1 + treatment | f_year / grazer))
  m24 <- update(best_fixed, . ~ . + (1 + treatment | grazer))
  m25 <- update(best_fixed, . ~ . + (1 + treatment | f_break))
  m26 <- update(best_fixed, . ~ . + (1 + treatment | f_one_yr))
  m27 <- update(best_fixed, . ~ . + (1 + treatment | f_two_yr))
  m28 <- update(best_fixed, . ~ . + (1 + treatment | f_new))
  
  model.sel(
    best_fixed,
    m11,
    m12, m13, m14,
    m15, m16, m17, m18,
    m21, m22, m23, m24,
    m25, m26, m27, m28
  )
}
# ---------------------------------------------------------- -----
# Fit reduced models  -----
# Create a bunch of potential models: glm, glmer, glm.nb, glmmTMB, lmer

#   value_std: normal (no random effect)
n0 <- lm(value_std ~ treatment, data = abun_nat)
n0n <- lm(value_std ~ 1, data = abun_nat)
#   value: normal
n1 <- lmer(value_std ~ treatment + (1 | plot_name), data = abun_nat, REML = FALSE)
n2 <- lmer(value_std ~ treatment + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)
n1n <- lmer(value_std ~ 1 + (1 | plot_name), data = abun_nat, REML = FALSE)
n2n <- lmer(value_std ~ 1 + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)

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

# Review reduced model performance ----
AIC(n0, n0n) # n0
model.sel(n1, n2, n1n, n2n) # n2, n1
model.sel(n3, n4, n3n, n4n) # n4, n3
model.sel(n5, n6, n5n, n6n) # n6, n5

plotQQunif(n0) # significant deviation for KS test
plotQQunif(n2) # significant deviation for KS test
plotQQunif(n4) # looks good
plotQQunif(n6) # significant deviation for KS test

best_reduced <- n4

testDispersion(best_reduced)
check_singularity(best_reduced)
testUniformity(best_reduced)
plotQQunif(best_reduced)
testOutliers(best_reduced)

plotResiduals(best_reduced)
plotResiduals(best_reduced, abun_nat$treatment)
plotResiduals(best_reduced, abun_nat$f_year) # y4
plotResiduals(best_reduced, abun_nat$f_break)
plotResiduals(best_reduced, abun_nat$f_one_yr)
plotResiduals(best_reduced, abun_nat$f_two_yr)
plotResiduals(best_reduced, abun_nat$f_new)
plotResiduals(best_reduced, abun_nat$grazer) # grazer

# Best reduced: n4 ----
lmer(value_log ~ treatment + (1 + treatment | plot_name), data = abun_nat, REML = FALSE)

# Define subset and model ----
abun_nat <- abun %>%
  filter(met_sub %in% "abun_nat")

model_abun_nat_red <-
  lmer(value_log ~ treatment + (1 + treatment | plot_name),
       data = abun, subset = met_sub == "abun_nat", REML = FALSE
  )

# Fit fixed effects  ----
abun_nat_tbl_fix_dredge <- fxn_fixed_dredge_by_model(index_subset = "abun_nat")
abun_nat_tbl_fix_one <- fxn_fixed_one(index_subset = "abun_nat")
abun_nat_tbl_fix_two <- fxn_fixed_two(index_subset = "abun_nat")


abun_nat_graz <- lmer(value_log ~ treatment + (1 + treatment | plot_name) + grazer,
                      data = abun_nat, REML = FALSE
)

plotResiduals(abun_nat_graz, abun_nat$treatment)
plotResiduals(abun_nat_graz, abun_nat$f_year)
plotResiduals(abun_nat_graz, abun_nat$f_break)
plotResiduals(abun_nat_graz, abun_nat$f_one_yr)
plotResiduals(abun_nat_graz, abun_nat$f_two_yr)
# plotResiduals(abun_nat_graz, abun_nat$grazer)
testUniformity(abun_nat_graz)

abun_nat_yr <- lmer(
  value_log ~ treatment + (1 + treatment | plot_name) +
    f_year +
    plot_type,
  data = abun_nat, REML = FALSE
)

plotResiduals(abun_nat_yr, abun_nat$treatment)
plotResiduals(abun_nat_yr, abun_nat$f_year)
plotResiduals(abun_nat_yr, abun_nat$f_break)
plotResiduals(abun_nat_yr, abun_nat$f_one_yr)
plotResiduals(abun_nat_yr, abun_nat$f_two_yr)
testUniformity(abun_nat_yr)

#   Identify best fixed model ----
model_abun_nat_fix <- lmer(value_log ~ treatment + (1 + treatment | plot_name) + f_year + plot_type,
                           data = abun, subset = met_sub == "abun_nat", REML = FALSE
)

plotResiduals(model_abun_nat_fix, abun_nat$treatment)
plotResiduals(model_abun_nat_fix, abun_nat$f_year)
plotResiduals(model_abun_nat_fix, abun_nat$f_break)
plotResiduals(model_abun_nat_fix, abun_nat$f_one_yr)
plotResiduals(model_abun_nat_fix, abun_nat$f_two_yr)
testUniformity(model_abun_nat_fix)

# Fit random effects  ----
#   Reduced model with random effects ----
abun_nat_red_tbl_ran <- fxn_random(index_model = model_abun_nat_red)
#   Fixed effects model(s) with random effects ----
abun_nat_fix_tbl_ran <- fxn_random(index_model = model_abun_nat_fix)

model_abun_nat_ran <- lmer(value_log ~ treatment + (1 + treatment | plot_name) + f_year + plot_type,
                           data = abun, subset = met_sub == "abun_nat", REML = FALSE
)
#   Review random model performance ----
temp_model <- model_abun_nat_ran

testDispersion(temp_model)

check_singularity(temp_model)
testUniformity(temp_model)
testOutliers(temp_model)

plotResiduals(temp_model, abun_nat$treatment)
plotResiduals(temp_model, abun_nat$f_year)
plotResiduals(temp_model, abun_nat$f_break)
plotResiduals(temp_model, abun_nat$f_one_yr)
plotResiduals(temp_model, abun_nat$f_two_yr)
plotResiduals(temp_model, abun_nat$f_new)
# plotResiduals(temp_model, abun_nat$grazer)

#   Final model ----
model_abun_nat_final <-
  lmer(value_log ~ treatment + (1 + treatment | plot_name) + f_year + plot_type,
       data = abun, subset = met_sub == "abun_nat", REML = FALSE
  )

plotResiduals(model_abun_nat_final, abun_nat$treatment)
plotResiduals(model_abun_nat_final, abun_nat$f_year)
plotResiduals(model_abun_nat_final, abun_nat$f_break)
plotResiduals(model_abun_nat_final, abun_nat$f_new)
plotResiduals(model_abun_nat_final, abun_nat$grazer) #
plotResiduals(model_abun_nat_final, abun_nat$f_one_yr)
plotResiduals(model_abun_nat_final, abun_nat$f_two_yr)


# ========================================================== -----