###############################################################################
                ##### Combine multivariate models #####
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
interaction_files <- list.files(path = c("Results/Models/Interaction_longer/"), full.names = T, include.dirs = FALSE, recursive = FALSE)
migration_files <- list.files(path = c("Results/Models/Migration/"), full.names = T, include.dirs = FALSE, recursive = FALSE)
seasonality_files <- list.files(path = c("Results/Models/Seasonality/"), full.names = T, include.dirs = FALSE, recursive = FALSE)

# Create all models.
all_interaction_models <- combine_brms(model_pattern = "/all", file_pathways = interaction_files)
all_migration_models <- combine_brms(model_pattern = "/all", file_pathways = migration_files)
all_seasonality_models <- combine_brms(model_pattern = "/all", file_pathways = seasonality_files)

# Create high models.
high_interaction_models <- combine_brms(model_pattern = "/high", file_pathways = interaction_files)
high_migration_models <- combine_brms(model_pattern = "/high", file_pathways = migration_files)
high_seasonality_models <- combine_brms(model_pattern = "/high", file_pathways = seasonality_files)

gc()

saveRDS(object = all_interaction_models, file =  "Results/Models/Combined_models/Multivariate/interaction_all_models.rds")
saveRDS(object = all_migration_models, file =  "Results/Models/Combined_models/Multivariate/migration_all_models.rds")
saveRDS(object = all_seasonality_models, file =  "Results/Models/Combined_models/Multivariate/seasonality_all_models.rds")

saveRDS(object = high_interaction_models, file =  "Results/Models/Combined_models/Multivariate/interaction_high_models.rds")
saveRDS(object = high_migration_models, file =  "Results/Models/Combined_models/Multivariate/migration_high_models.rds")
saveRDS(object = high_seasonality_models, file =  "Results/Models/Combined_models/Multivariate/seasonality_high_models.rds")



###############################################################################
                       #### Model R2 ####



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









Bayes_R2_MZ(all_interaction_models)
Bayes_R2_MZ(all_seasonality_models)
Bayes_R2_MZ(all_migration_models)

Bayes_R2_MZ(high_interaction_models)
Bayes_R2_MZ(high_seasonality_models)
Bayes_R2_MZ(high_migration_models)


###############################################################################
                             #### END ####
###############################################################################



