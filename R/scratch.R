# # Appendix {#sec-Appendix}
# 
# **Diagnostic plots for excluded fixed effect models fit to native species abundance**
#   
#   The following are the diagnostic plots for models with lower AIC to illustrate their issues with uniformity.
# 
# Model diagnostics for the model with the lowest AIC: f_one_yr
# 
# ```{r}
# # fxn_lmer_model_review(models_fixed_abun_nat$best_model$model, abun_nat)
# model_fixed_abun_nat_1 <- update(best_model_init_abun_nat, . ~ . + f_one_yr)
# 
# fxn_lmer_model_review(model_fixed_abun_nat_1, abun_nat)
# 
# ```
# 
# <br>
#   
#   Model diagnostics for the model with the second lowest AIC: f_break
# 
# ```{r}
# model_fixed_abun_nat_2 <- update(best_model_init_abun_nat, . ~ . + f_break)
# 
# fxn_lmer_model_review(model_fixed_abun_nat_2, abun_nat)
# 
# ```
# 
# <br>
#   
#   Model diagnostics for the model with the third lowest AIC: f_new
# 
# ```{r}
# model_fixed_abun_nat_3 <- update(best_model_init_abun_nat, . ~ . + f_new)
# 
# fxn_lmer_model_review(model_fixed_abun_nat_3, abun_nat)
# 
# 
# ```
# 
# <br>
#   
#   Model diagnostics for the model with the fourth lowest AIC: plot_type
# 
# ```{r}
# model_fixed_abun_nat_4 <- update(best_model_init_abun_nat, . ~ . + plot_type)
# 
# fxn_lmer_model_review(model_fixed_abun_nat_4, abun_nat)
# 
# 
# ```
# 
