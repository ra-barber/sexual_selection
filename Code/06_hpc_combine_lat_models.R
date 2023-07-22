###############################################################################
                ##### Combine brms models for latitude #####
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


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Read in the new files.
model_files <- list.files(path = c("Results/Models/Latitude/"), 
                          full.names = T, include.dirs = FALSE, recursive = FALSE)
# model_files <- lat_files

#model_files <- lat_files
model_names <- model_files %>% str_extract("[a-z]*_[a-z]*_[a-z]*")
model_list <- unique(model_names)

# Do it manually to change frugnect name (don't add extra underscore next time).
model_list <- c(
  #"all_centered_all", "all_centered_high",
  "all_uncentered_all", "all_uncentered_high",
  #"certainty_centered_all", "certainty_centered_high",   
  "certainty_uncentered_all", "certainty_uncentered_high",
  #"frug_nect_centered_all", "frug_nect_centered_high",
  #"frug_nect_uncentered_all", "frug_nect_uncentered_high",
  #"frugivore_centered_all", "frugivore_centered_high",
  "frugivore_uncentered_all", "frugivore_uncentered_high",
  #"invertivore_centered_all", "invertivore_centered_high",
  "invertivore_uncentered_all", "invertivore_uncentered_high",
  #"nectarivore_centered_all", "nectarivore_centered_high",
  #"nectarivore_uncentered_all", "nectarivore_uncentered_high",
  #"primary_centered_all", "primary_centered_high",
  "primary_uncentered_all", "primary_uncentered_high",
  #"secondary_centered_all", "secondary_centered_high",
  "secondary_uncentered_all", "secondary_uncentered_high")


# Extract model name.
current_model_name <- model_list[[array_number]]

# Read in model.
current_model <- combine_brms(current_model_name, file_pathways = model_files) 
gc()

# Save the combined model.
file_pathway <- paste0("Results/Models/Combined_models/Latitude/", current_model_name, "_models.rds")
saveRDS(object = current_model, file =  file_pathway)

# model_pattern = current_model_name
# # Function to extract and combine brms models.
# combine_brms <- function(model_pattern, file_pathways = lat_files){
#   pattern_matched <- file_pathways %>% str_subset(pattern = model_pattern)
#   if (length(pattern_matched) == 1){
#     brms_models <- readRDS(pattern_matched)
#   } else {
#     
#     for (x in 1:length(pattern_matched)){
#       if (x ==1){
#         brms_model1 <- readRDS(pattern_matched[x])
#         print(summary(brms_model1))
#       } else if (x == 2){
#         brms_model2 <- readRDS(pattern_matched[x])
#         print(summary(brms_model2))
#         brms_models <- combine_models(brms_model1, brms_model2, check_data = FALSE)
#         rm(brms_model1, brms_model2)
#       } else {
#         brms_model <- readRDS(pattern_matched[x])
#         print(summary(brms_model))
#         brms_models <- combine_models(brms_models, brms_model, check_data = FALSE)
#         rm(brms_model)
#       }
#       print(x)
#     }
#   }
#   return(brms_models)
# }
# 



###############################################################################
             ##### Extract r squared and model predictions #####


# Extra model predictions for plotting models.
current_model_predictions <- conditional_effects(current_model)[[1]]

# Extract r-squared as most time consuming bit.
current_model_r_squared <- Bayes_R2_MZ(current_model)

# All details.
current_model_deets <- list(current_model, current_model_predictions, current_model_r_squared)

# Export model details for future use.
file_pathway <- paste0("Results/Models/Combined_models/Latitude/", current_model_name, "_data.rds")
saveRDS(current_model_deets, file = file_pathway)



###############################################################################
                             #### END ####
###############################################################################


