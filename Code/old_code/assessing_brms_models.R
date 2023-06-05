###############################################################################
                   # Assess the fit of brms models  #
###############################################################################


# Packages to load.
library(magrittr)
library(tictoc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(brms)

# Clear the workspace.
rm(list=ls())

# Functions.
source("Code/functions.R")


###############################################################################
                  #### Read in the models #####

# Read in the model.
standard_all <- readRDS("Z:/home/sexual_selection/Results/Models/37_all.rds")
standard_high <- readRDS("Z:/home/sexual_selection/Results/Models/37_high.rds")





sum <- summary(standard_high)

sum$fixed
summary(standard_high)

library(bayesplot)

# 
# # Plot the model.
plot(all_brms_model)
plot(standard_high)
# 
# # Check the posterior distribution.
pp_check(standard_high)
pp_check(centered_all)
# # Look at the pairs.
# #pairs(brms_model)
# 
library(bayesplot)
# 
# # Do mcmc_areas.
mcmc_areas(standard_high, regex_pars = "^b_[a-z]", prob = 0.95)
mcmc_areas(centered_all, regex_pars = "^b_[a-z]", prob = 0.95)
mcmc_areas(no_env_all, regex_pars = "^b_[a-z]", prob = 0.95)

hyp <- "sd_tree_tip__Intercept^2 / (sd_tree_tip__Intercept^2 + disc^2) = 0"

(hyp_2 <- hypothesis(standard_all, hyp, class = NULL))
(hyp_2 <- hypothesis(no_env_all, hyp, class = NULL))


?brm
