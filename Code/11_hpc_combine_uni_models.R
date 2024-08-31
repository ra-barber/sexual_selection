###############################################################################
             ##### Combine univarate models from the HPC #####
###############################################################################

# Combine univariate models as an array job into a single folder on the HPC:
# Combined_models/Univariate


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
model_files <- list.files(path = c("Results/Models/Univariate/"), 
                          full.names = T, include.dirs = FALSE, recursive = FALSE)

#model_files <- lat_files
model_names <- model_files %>% str_extract("[a-z]*_[a-z]*_[a-z]*")
model_list <- unique(model_names)

# Extract model name.
current_model_name <- model_list[[array_number]]

# Read in model.
current_model <- combine_brms(current_model_name, file_pathways = model_files) 
gc()

# Save the combined model.
file_pathway <- paste0("Results/Models/Combined_models/Univariate/", current_model_name, "_models.rds")
saveRDS(object = current_model, file =  file_pathway)


###############################################################################
             ##### Extract r squared and model predictions #####


# Extra model predictions for plotting models.
current_model_predictions <- conditional_effects(current_model)[[1]]

# Extract r-squared as most time consuming bit.
current_model_r_squared <- Bayes_R2_MZ(current_model)

# All details.
current_model_deets <- list(current_model, current_model_predictions, current_model_r_squared)

# Export model details for future use.
file_pathway <- paste0("Results/Models/Combined_models/Univariate/", current_model_name, "_data.rds")
saveRDS(current_model_deets, file = file_pathway)


###############################################################################
                             #### END ####
###############################################################################