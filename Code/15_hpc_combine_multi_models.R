###############################################################################
                ##### Combine multivariate models #####
###############################################################################

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
model_files <- list.files(path = c("Results/Models/Multivariate/"), 
                          full.names = T, include.dirs = FALSE, 
                          recursive = FALSE)

# Function to extract and combine brms models.
combine_brms <- function(model_pattern){
  pattern_matched <- model_files %>% str_subset(pattern = model_pattern)
  if (length(pattern_matched) == 1){
    brms_models <- readRDS(pattern_matched)
  } else {
    
    for (x in 1:length(pattern_matched)){
      if (x ==1){
        brms_model1 <- readRDS(pattern_matched[x])
        print(summary(brms_model1))
      } else if (x == 2){
        brms_model2 <- readRDS(pattern_matched[x])
        print(summary(brms_model2))
        brms_models <- combine_models(brms_model1, brms_model2, check_data = FALSE)
        rm(brms_model1, brms_model2)
      } else {
        brms_model <- readRDS(pattern_matched[x])
        print(summary(brms_model))
        brms_models <- combine_models(brms_models, brms_model, check_data = FALSE)
        rm(brms_model)
      }
      print(x)
    }
  }
  return(brms_models)
}

# Read in models.
all_models <- combine_brms("/all")
high_models <- combine_brms("/high")
gc()

saveRDS(object = all_models, 
        file =  "Results/Models/Combined_models/Multivariate/all_models.rds")
saveRDS(object = high_models, 
        file =  "Results/Models/Combined_models/Multivariate/high_models.rds")


all_r_squared <- Bayes_R2_MZ(all_models)
all_r_squared
high_r_squared <- Bayes_R2_MZ(high_models)
high_r_squared


###############################################################################
                             #### END ####
###############################################################################