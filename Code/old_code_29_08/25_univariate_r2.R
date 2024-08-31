###############################################################################
               # Assess the r2 of univariate models  #
###############################################################################


# R squared from models are pasted at the bottom of the script. 

# Packages to load.
library(magrittr)
library(tictoc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(brms)
library(bayesplot)
library(ggdist)
library(ggbeeswarm)
library(stringr)
library(bayestestR)
library(tidybayes)

# Clear the workspace.
rm(list=ls())

# Functions.
source("Code/functions.R")


# Read in the models.
first_half <- "Z:/home/sexual_selection/Results/Models/Consensus/Univariate/ordinal_"

# Read in models using raw data.
tro_model <- readRDS(paste0(first_half, "tro_all.rds"))
mig_model <- readRDS(paste0(first_half, "mig_all.rds"))
terr_model <- readRDS(paste0(first_half, "terr_all.rds"))
temp_model <- readRDS(paste0(first_half, "temp_z_all.rds"))

high_tro_model <- readRDS(paste0(first_half, "tro_high.rds"))
high_mig_model <- readRDS(paste0(first_half, "mig_high.rds"))
high_terr_model <- readRDS(paste0(first_half, "terr_high.rds"))
high_temp_model <- readRDS(paste0(first_half, "temp_z_high.rds"))

# Read in the multivariate models.
multi_model <- readRDS("Z:/home/sexual_selection/Results/Models/Consensus/Multivariate/ordinal_all.rds")
high_multi_model <- readRDS("Z:/home/sexual_selection/Results/Models/Consensus/Multivariate/ordinal_high.rds")





library(tictoc)
tic()
tro_marg_r2 <- marginal_R2_MZ(tro_model)
toc()

tic()
mig_marg_r2 <- marginal_R2_MZ(mig_model)
toc()

terr_marg_r2 <- marginal_R2_MZ(terr_model)
temp_marg_r2 <- marginal_R2_MZ(temp_model)

tic()
high_tro_marg_r2 <- marginal_R2_MZ(high_tro_model)
toc()

tic()
high_mig_marg_r2 <- marginal_R2_MZ(high_mig_model)
toc()

high_terr_marg_r2 <- marginal_R2_MZ(high_terr_model)
high_temp_marg_r2 <- marginal_R2_MZ(high_temp_model)

# Multivariate marginal.
multi_marg_r2 <- marginal_R2_MZ(multi_model)
high_multi_marg_r2 <- marginal_R2_MZ(high_multi_model)



tic()
tro_full_r2 <- Bayes_R2_MZ(tro_model)
toc()

tic()
mig_full_r2 <- Bayes_R2_MZ(mig_model)
toc()

terr_full_r2 <- Bayes_R2_MZ(terr_model)
temp_full_r2 <- Bayes_R2_MZ(temp_model)

tic()
high_tro_full_r2 <- Bayes_R2_MZ(high_tro_model)
toc()

tic()
high_mig_full_r2 <- Bayes_R2_MZ(high_mig_model)
toc()

high_terr_full_r2 <- Bayes_R2_MZ(high_terr_model)
high_temp_full_r2 <- Bayes_R2_MZ(high_temp_model)

################################################################################
                       ###### Reported stats #####

## Marginal stats.

tro_marg_r2 <- marginal_R2_MZ(tro_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ    0.022     0.015   0.0011    0.059

> mig_marg_r2 <- marginal_R2_MZ(mig_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.12     0.019    0.084     0.16

> terr_marg_r2 <- marginal_R2_MZ(terr_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.34     0.046     0.25     0.42

> temp_marg_r2 <- marginal_R2_MZ(temp_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.31     0.033     0.25     0.38


> high_tro_marg_r2 <- marginal_R2_MZ(high_tro_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ    0.029     0.022  0.00057    0.086

> high_mig_marg_r2 <- marginal_R2_MZ(high_mig_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.11     0.022    0.071     0.16

> high_terr_marg_r2 <- marginal_R2_MZ(high_terr_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.46     0.053     0.36     0.57

> high_temp_marg_r2 <- marginal_R2_MZ(high_temp_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.38     0.044      0.3     0.47


# Full R squared.

> tro_full_r2 <- Bayes_R2_MZ(tro_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.93    0.0089     0.91     0.94

> mig_full_r2 <- Bayes_R2_MZ(mig_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.92    0.0086     0.91     0.94

> terr_full_r2 <- Bayes_R2_MZ(terr_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.92     0.009      0.9     0.94

> temp_full_r2 <- Bayes_R2_MZ(temp_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.93    0.0081     0.91     0.94



> high_tro_full_r2 <- Bayes_R2_MZ(high_tro_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.95    0.0081     0.93     0.96

> high_mig_full_r2 <- Bayes_R2_MZ(high_mig_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.94    0.0082     0.92     0.96

> high_terr_full_r2 <- Bayes_R2_MZ(high_terr_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.94    0.0087     0.92     0.96

> high_temp_full_r2 <- Bayes_R2_MZ(high_temp_model)
            Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.95    0.0073     0.93     0.96


# Multivariate models.
> multi_marg_r2 <- marginal_R2_MZ(multi_model)
Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.45     0.032     0.38     0.51
> high_multi_marg_r2 <- marginal_R2_MZ(high_multi_model)
Estimate Est.Error l-95% CI u-95% CI
Bayes_R2_MZ     0.54     0.036     0.47     0.61



