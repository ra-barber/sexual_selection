###############################################################################
                        ##### Script Title #####
###############################################################################

# This script does something new. Pretty sick eh.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(tictoc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(brms)
library(stringr)

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                  #### Read in models ####

# Read in the new files.
seasonality_files <- list.files(path = c("Results/Models/Seasonality_sample/"), full.names = T, include.dirs = FALSE, recursive = FALSE)
migration_files <- list.files(path = c("Results/Models/Migration_sample/"), full.names = T, include.dirs = FALSE, recursive = FALSE)

migration_files <- migration_files[-2]
migration_files <- migration_files[-5]

# Read in models.
all_seasonal_models <- combine_brms("/all", seasonality_files)
high_seasonal_models <- combine_brms("/high", seasonality_files)

all_migration_models <- combine_brms("/all", migration_files)
high_migration_models <- combine_brms("/high", migration_files)


all_mig_1 <- readRDS(migration_files[[1]])
all_mig_10 <- readRDS(migration_files[[2]])
all_mig_2 <- readRDS(migration_files[[3]])

combine_models(all_mig_1, all_mig_2, check_data = FALSE)
?combine_models

high_models <- combine_brms("/high")
gc()

saveRDS(object = all_models, file =  "Results/Models/Combined_models/Multivariate/multivariate_all_models.rds")
saveRDS(object = high_models, file =  "Results/Models/Combined_models/Multivariate/multivariate_high_models.rds")

###############################################################################
                       #### Model R2 ####

library(bayestestR)
library(bayesplot)


brms_forest(all_seasonal_models)
brms_forest(all_migration_models)







# Bayesian McKelvey-Zavoina R2 ------------------------------------------------
# Bayesian version of McKelvey and Zavoina's pseudo-R2 for binary and ordinal
# brms models (McKelvey and Zavoina, 1975). See also Gelman et al. (2018). This
# pseudo-R2 closely approximates the R2 that would have been obtained if a
# linear model had have been run on observations of the continuous latent
# variable underlying the discrete responses (Veall and Zimmermann, 1992; Hagle
# and Mitchell, 1992; Veall and Zimmermann, 1994).
Bayes_R2_MZ <- function(fit, ...) {
  y_pred <- fitted(fit, scale = "linear", summary = FALSE, ...)
  var_fit <- apply(y_pred, 1, var)
  if (fit$formula$family$family == "cumulative" ||
      fit$formula$family$family == "bernoulli") {
    if (fit$formula$family$link == "probit" || 
        fit$formula$family$link == "probit_approx") {
      var_res <- 1
    }
    else if (fit$formula$family$link == "logit") {
      var_res <- pi^2 / 3 
    }
  } 
  else {
    sum_fit <- summary(fit)
    sig_res <- sum_fit$spec_pars["sigma", "Estimate"]
    var_res <- sig_res^2
  } 
  R2_MZ <- var_fit / (var_fit + var_res)
  print(
    data.frame(
      Estimate = mean(R2_MZ), 
      Est.Error = sd(R2_MZ), 
      "l-95% CI" = quantile(R2_MZ, 0.025),
      "u-95% CI" = quantile(R2_MZ, 0.975),
      row.names = "Bayes_R2_MZ", 
      check.names = FALSE), 
    digits = 2)
}









all_r_squared <- bayes_R2(all_models)
all_r_squared
high_r_squared <- bayes_R2(high_models)
high_r_squared


###############################################################################
                     #### Plot conditional effects ####
# 
# 
# pdf("Results/Models/Combined_models/cond_plots_all_nocat.pdf")
# plot(conditional_effects(all_models, categorical = FALSE), ask = FALSE)
# dev.off()
# 
# pdf("Results/Models/Combined_models/cond_plots_high_nocat.pdf")
# plot(conditional_effects(high_models, categorical = FALSE), ask = FALSE)
# dev.off()
# 
# pdf("Results/Models/Combined_models/cond_plots_all_cat.pdf")
# plot(conditional_effects(all_models, categorical = TRUE), ask = FALSE)
# dev.off()
# 
# pdf("Results/Models/Combined_models/cond_plots_high_cat.pdf")
# plot(conditional_effects(high_models, categorical = TRUE), ask = FALSE)
# dev.off()

###############################################################################
                           #### Section 3 ####


###############################################################################
                           #### Section 4 ####

###############################################################################
                           #### Section 5 ####

###############################################################################
                           #### Section 6 ####


###############################################################################
                           #### Section 7 ####


###############################################################################
                           #### Section 8 ####

# Look Rob, you've had your fun with the sectioning. 
# They'll be no more sectioning today.


###############################################################################
                             #### END ####
###############################################################################


###############################################################################
              #### All the stuff I'm afraid to delete ####


