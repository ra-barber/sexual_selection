

################################################################################
###### SQRT plots #####


mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sqrt(family_ss))) + 
  categorical_settings + ylim(c(0, 2.7)) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_se, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = sqrt(sexual_mean - (1.96*sexual_se)), 
                    ymax = sqrt(sexual_mean + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_se, inherit.aes = FALSE,
             aes(x=migration_binary, y = sqrt(sexual_mean)), size = 3.6) + theme(axis.text.y = element_blank())

tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sqrt(family_ss))) + 
  categorical_settings +ylim(c(0, 2.7)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sqrt(sexual_mean - (1.96*sexual_se)), 
                    ymax = sqrt(sexual_mean + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sqrt(sexual_mean)), size = 3.6) + theme(axis.text.y = element_blank())

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=sqrt(family_ss))) + 
  categorical_settings +ylim(c(0, 2.7)) +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_se, inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = sqrt(sexual_mean - (1.96*sexual_se)), 
                    ymax = sqrt(sexual_mean + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_se, inherit.aes = FALSE,
             aes(x=territoriality_binary, y = sqrt(sexual_mean)), size = 3.6) + theme(axis.text.y = element_blank())

# Make some nice plots.
temp_plot <- family_seasonality_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + ylim(c(0,2.7)) +
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())

side_plots_sqrt <- ggarrange(tro_plot + ylab("Sexual selection"), 
                             mig_plot,
                             terr_plot + ylab("Sexual selection"), 
                             temp_plot,
                             nrow = 2, ncol = 2, align = "h",
                             widths = c(1.2,1))
side_plots_sqrt
side_plots


temp_plot <- family_seasonality_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") +
  #  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
  #             aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_smooth(data = temp_plot_predictions, aes(x = temp_log, y = sqrt(estimate__)), method = "lm", linetype = "dashed", linewidth = 1, se = FALSE, colour = "black") + 
  #geom_smooth(data = temp_plot_predictions, aes(x = temp_log, y = sqrt(lower__)), method = "lm", linetype = "dotted", linewidth = 1, se = FALSE, colour = "black") + 
  #geom_smooth(data = temp_plot_predictions, aes(x = temp_log, y = sqrt(upper__)), method = "lm", linetype = "dotted", linewidth = 1, se = FALSE, colour = "black") + 
  theme(axis.text.y = element_blank()) + ylim(c(0,2.7)) + coord_cartesian(xlim = c(3.3,7.2), clip = "on")
temp_plot_predictions$lower__[temp_plot_predictions$lower__ < 0] <- 0







# 
# ###############################################################################
#                 ##### Extract model estimates #####
# 
# terr_estimates <- convert_ordinal_back(terr_uncen_all_data[[1]])
# tro_estimates <- convert_ordinal_back(tro_uncen_all_data[[1]])
# mig_estimates <- convert_ordinal_back(mig_uncen_all_data[[1]])
# 












fit <- fruit_example_model


post <- as_draws_df(fit) %>% 
  select(
    starts_with("b_")
  )

posterior_interval(fit)
?posterior_interval
postEst <- post %>% 
  mutate(
    EstControl = 1*pnorm(`b_Intercept[1]`) + 
      2*(pnorm(`b_Intercept[2]`) - pnorm(`b_Intercept[1]`)) +
      3*(pnorm(`b_Intercept[3]`) - pnorm(`b_Intercept[2]`)) +
      4*(pnorm(`b_Intercept[4]`) - pnorm(`b_Intercept[3]`)) +
      5*(1 - pnorm(`b_Intercept[4]`)),
    EstTreatment = 1*pnorm(`b_Intercept[1]` - b_trop_non_tropnon_trop) + 
      2*(pnorm(`b_Intercept[2]` - b_trop_non_tropnon_trop) - pnorm(`b_Intercept[1]` - b_trop_non_tropnon_trop)) +
      3*(pnorm(`b_Intercept[3]` - b_trop_non_tropnon_trop) - pnorm(`b_Intercept[2]` - b_trop_non_tropnon_trop)) +
      4*(pnorm(`b_Intercept[4]` - b_trop_non_tropnon_trop) - pnorm(`b_Intercept[3]` - b_trop_non_tropnon_trop)) +
      5*(1 - pnorm(`b_Intercept[4]` - b_trop_non_tropnon_trop))
  )


# Could try this on the standard error bit of conditional output.
options(scipen = 999)

tro_intervals <- posterior_interval(tro_example_model, variable = "b_", regex = TRUE)

test_intervals <- tro_intervals %>% as.data.frame() %>%  mutate(low_ci = pnorm(`2.5%`), high_ci = pnorm(`97.5%`))

mean(test_intervals$low_ci[1:4]) # 0.3

mean(test_intervals$high_ci[1:4]) # 0.58
model <- mig_uncen_all_data[[1]]
convert_ordinal_back <- function(model){
  summ <- summary(model)
  model_coeffs <- summ$fixed
  
  # Convert to normal scale using pnorm.
  model_coeffs$est_norm <- pnorm(model_coeffs$Estimate)
  model_coeffs$err_norm <- pnorm(model_coeffs$Est.Error)
  model_coeffs$low_norm <- pnorm(model_coeffs$`l-95% CI`)
  model_coeffs$upp_norm <- pnorm(model_coeffs$`u-95% CI`)
  
  inter_est <- mean(model_coeffs$est_norm[1:4])
  inter_err <- mean(model_coeffs$err_norm[1:4])
  inter_low <- mean(model_coeffs$low_norm[1:4])
  inter_upp <- mean(model_coeffs$upp_norm[1:4])
  
  effect_est <- model_coeffs$est_norm[5]
  effect_err <- model_coeffs$err_norm[5]
  effect_low <- model_coeffs$low_norm[5]
  effect_upp <- model_coeffs$upp_norm[5]
  
  error_bar_data <- data.frame(predictor = c("Intercept", "Effect"), 
                               estimate = c(inter_est, effect_est+inter_est),
                               lower = c(inter_low, effect_low+inter_est),
                               upper = c(inter_upp, effect_upp+inter_est))
  return(error_bar_data)
}

# Extract model coeffcients.
summ <- summary(tro_example_model)
model_coeffs <- summ$fixed

# Convert to normal scale using pnorm.
model_coeffs$est_norm <- pnorm(model_coeffs$Estimate)
model_coeffs$err_norm <- pnorm(model_coeffs$Est.Error)
model_coeffs$low_norm <- pnorm(model_coeffs$`l-95% CI`)
model_coeffs$upp_norm <- pnorm(model_coeffs$`u-95% CI`)

inter_est <- mean(model_coeffs$est_norm[1:4])
inter_err <- mean(model_coeffs$est_err[1:4])
inter_low <- mean(model_coeffs$low_norm[1:4])
inter_upp <- mean(model_coeffs$upp_norm[1:4])

effect_est <- model_coeffs$est_norm[5]
effect_err <- model_coeffs$est_err[5]
effect_low <- model_coeffs$low_norm[5]
effect_upp <- model_coeffs$upp_norm[5]

error_bar_data <- data.frame(trophic_level_binary = c("Secondary", "Primary"), 
           estimate = c(inter_est, effect_est),
           lower = c(inter_low, effect_low),
           upper = c(inter_upp, effect_upp))

trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = error_bar_data, inherit.aes = FALSE,
                aes(x=predictor, ymin = lower, 
                    ymax = upper),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = error_bar_data, inherit.aes = FALSE,
             aes(x=predictor, y = estimate), size = 3.6)



model_coeffs <- summ$fixed

# Convert to normal scale using pnorm.
model_coeffs$est_norm <- pnorm(model_coeffs$Estimate)
model_coeffs$err_norm <- pnorm(model_coeffs$Est.Error)
model_coeffs$low_norm <- pnorm(model_coeffs$`l-95% CI`)
model_coeffs$upp_norm <- pnorm(model_coeffs$`u-95% CI`)

inter_est <- mean(model_coeffs$est_norm[1:4])
inter_err <- mean(model_coeffs$est_err[1:4])
inter_low <- mean(model_coeffs$low_norm[1:4])
inter_upp <- mean(model_coeffs$upp_norm[1:4])

effect_est <- model_coeffs$est_norm[5]
effect_err <- model_coeffs$est_err[5]
effect_low <- model_coeffs$low_norm[5]
effect_upp <- model_coeffs$upp_norm[5]

error_bar_data <- data.frame(predictor = c("Secondary", "Primary"), 
                             estimate = c(inter_est, effect_est),
                             lower = c(inter_low, effect_low),
                             upper = c(inter_upp, effect_upp))

pnorm(0.6)
model_coeffs %>% 
  mutate(
    EstControl = 1*pnorm(`b_Intercept[1]`) + 
      2*(pnorm(`b_Intercept[2]`) - pnorm(`b_Intercept[1]`)) +
      3*(pnorm(`b_Intercept[3]`) - pnorm(`b_Intercept[2]`)) +
      4*(pnorm(`b_Intercept[4]`) - pnorm(`b_Intercept[3]`)) +
      5*(1 - pnorm(`b_Intercept[4]`)),
    EstTreatment = 1*pnorm(`b_Intercept[1]` - b_trop_non_tropnon_trop) + 
      2*(pnorm(`b_Intercept[2]` - b_trop_non_tropnon_trop) - pnorm(`b_Intercept[1]` - b_trop_non_tropnon_trop)) +
      3*(pnorm(`b_Intercept[3]` - b_trop_non_tropnon_trop) - pnorm(`b_Intercept[2]` - b_trop_non_tropnon_trop)) +
      4*(pnorm(`b_Intercept[4]` - b_trop_non_tropnon_trop) - pnorm(`b_Intercept[3]` - b_trop_non_tropnon_trop)) +
      5*(1 - pnorm(`b_Intercept[4]` - b_trop_non_tropnon_trop))
  )

# Calculate intercept estimate.
model_coeffs$est_norm[1] + 2*(model_coeffs$est_norm[2] - model_coeffs$est_norm[1]) + 3*(model_coeffs$est_norm[3] - model_coeffs$est_norm[2]) + 4*(model_coeffs$est_norm[4] - model_coeffs$est_norm[3]) + 5*(1 - model_coeffs$est_norm[4])  

# Calculate intercept effect size.
model_coeffs$err_norm[1] + 2*(model_coeffs$err_norm[2] - model_coeffs$err_norm[1]) + 3*(model_coeffs$err_norm[3] - model_coeffs$err_norm[2]) + 4*(model_coeffs$err_norm[4] - model_coeffs$err_norm[3]) + 5*(1 - model_coeffs$err_norm[4])  

model_coeffs$low_norm[1] + 2*(model_coeffs$low_norm[2] - model_coeffs$low_norm[1]) + 3*(model_coeffs$low_norm[3] - model_coeffs$low_norm[2]) + 4*(model_coeffs$low_norm[4] - model_coeffs$low_norm[3]) + 5*(1 - model_coeffs$low_norm[4])  
model_coeffs$upp_norm[1] + 2*(model_coeffs$upp_norm[2] - model_coeffs$upp_norm[1]) + 3*(model_coeffs$upp_norm[3] - model_coeffs$upp_norm[2]) + 4*(model_coeffs$upp_norm[4] - model_coeffs$upp_norm[3]) + 5*(1 - model_coeffs$upp_norm[4])  


model_coeffs$est_norm[1] + 2*(model_coeffs$est_norm[2] - model_coeffs$est_norm[1]) + 3*(model_coeffs$est_norm[3] - model_coeffs$est_norm[2]) + 4*(model_coeffs$est_norm[4] - model_coeffs$est_norm[3]) + 5*(1 - model_coeffs$est_norm[4])  

tro_est <- 1*pnorm(model_coeffs$Estimate[1] - model_coeffs$Estimate[5]) + 
  2*(pnorm(model_coeffs$Estimate[2] - model_coeffs$Estimate[5]) - pnorm(model_coeffs$Estimate[1] - model_coeffs$Estimate[5])) +
  3*(pnorm(model_coeffs$Estimate[3] - model_coeffs$Estimate[5]) - pnorm(model_coeffs$Estimate[2] - model_coeffs$Estimate[5])) +
  4*(pnorm(model_coeffs$Estimate[4] - model_coeffs$Estimate[5]) - pnorm(model_coeffs$Estimate[3] - model_coeffs$Estimate[5])) +
  5*(1 - pnorm(model_coeffs$Estimate[4] - model_coeffs$Estimate[5]))


tro_err <- 1*pnorm(model_coeffs$Est.Error[1] - model_coeffs$Est.Error[5]) + 
  2*(pnorm(model_coeffs$Est.Error[2] - model_coeffs$Est.Error[5]) - pnorm(model_coeffs$Est.Error[1] - model_coeffs$Est.Error[5])) +
  3*(pnorm(model_coeffs$Est.Error[3] - model_coeffs$Est.Error[5]) - pnorm(model_coeffs$Est.Error[2] - model_coeffs$Est.Error[5])) +
  4*(pnorm(model_coeffs$Est.Error[4] - model_coeffs$Est.Error[5]) - pnorm(model_coeffs$Est.Error[3] - model_coeffs$Est.Error[5])) +
  5*(1 - pnorm(model_coeffs$Est.Error[4] - model_coeffs$Est.Error[5]))


tro_est-1
(tro_err-1)*1.96
1.2*1.96

pnorm(-10.48)

# Raw values from posterior interval.
# Lower CI
-1.029926e+01   
-5.075882e+00   
 -1.884825e+00  
 5.978222e+00 
 
# Upper CI
 -0.13555793
 3.52903200
 7.76350775
 26.61281750
 
 mean(c(pnorm(-1.029926e+01), pnorm(-5.075882e+00), pnorm(-1.884825e+00), pnorm(5.978222e+00 )))
 
 # Lower CI
 lo_int_1 <- pnorm(-1.029926e+01)   
 lo_int_2 <- pnorm(-5.075882e+00)   
 lo_int_3 <- pnorm(-1.884825e+00)  
 lo_int_4 <- pnorm(5.978222e+00) 
 
 1*lo_int_1 + 2*(lo_int_2 - lo_int_1) + 3*(lo_int_3 - lo_int_2) + 4*(lo_int_4 - lo_int_3) + 4*(1 - lo_int_4)
 
pnorm(0.2)

cat_test <- cat[[1]]

library(tidyr)

cat_test %>% pivot_wider(names_from = effect2__, values_from = se__)

cat_trop <- cat_test %>% filter(trop_non_trop == "trop")


trop_se_prob_1 <- pnorm(cat_trop[cat_trop$cats__==1,"se__"])
trop_se_prob_2 <- pnorm(cat_trop[cat_trop$cats__==2,"se__"])
trop_se_prob_3 <- pnorm(cat_trop[cat_trop$cats__==3,"se__"])
trop_se_prob_4 <- pnorm(cat_trop[cat_trop$cats__==4,"se__"])
trop_se_prob_5 <- pnorm(cat_trop[cat_trop$cats__==5,"se__"])

1*trop_se_prob_1 + 2*(trop_se_prob_2 - trop_se_prob_1) + 3*(trop_se_prob_3 - trop_se_prob_2) + 4*(trop_se_prob_4 - trop_se_prob_3) + 4*(trop_se_prob_5 - trop_se_prob_4) + 5*(1 - trop_se_prob_5)


1*pnorm(cat_trop[cat_trop$cats__==1,"se__"])



mean(postEst$EstControl)
mean(postEst$EstTreatment)
standard_error(postEst$EstControl)

?conditional_effects
qt(0.975, 8000)*(sd(postEst$EstControl)/sqrt(1000))
fruit_example_model

cat_2[[1]]

pnorm(0.03)
conditional_effects(fruit_example_model, conditons =conditions)[[1]]

2.72+0.639*1.96
2.72-0.639*1.96



conditions <- data.frame(trop_non_trop = c(0, 1))














# Try boxplots.

fit <- tro_example_model


tro_draws <- as_draws_df(fit) %>% 
  select(
    starts_with("b_")
  )


tro_est_test <- tro_draws %>% 
  mutate(
    EstControl = 1*pnorm(`b_Intercept[1]`) + 
      2*(pnorm(`b_Intercept[2]`) - pnorm(`b_Intercept[1]`)) +
      3*(pnorm(`b_Intercept[3]`) - pnorm(`b_Intercept[2]`)) +
      4*(pnorm(`b_Intercept[4]`) - pnorm(`b_Intercept[3]`)) +
      5*(1 - pnorm(`b_Intercept[4]`)),
    EstTreatment = 1*pnorm(`b_Intercept[1]` - b_trophic_level_binaryPrimary) + 
      2*(pnorm(`b_Intercept[2]` - b_trophic_level_binaryPrimary) - pnorm(`b_Intercept[1]` - b_trophic_level_binaryPrimary)) +
      3*(pnorm(`b_Intercept[3]` - b_trophic_level_binaryPrimary) - pnorm(`b_Intercept[2]` - b_trophic_level_binaryPrimary)) +
      4*(pnorm(`b_Intercept[4]` - b_trophic_level_binaryPrimary) - pnorm(`b_Intercept[3]` - b_trophic_level_binaryPrimary)) +
      5*(1 - pnorm(`b_Intercept[4]` - b_trophic_level_binaryPrimary))
  )

draws_longer <- tro_est_test %>% select(EstControl, EstTreatment) %>%  pivot_longer(cols = c(EstControl, EstTreatment), names_to = "level", values_to = "draws")

draws_longer %>% ggplot(aes(x = level, y = draws-1)) + geom_violin(fill = NA)













