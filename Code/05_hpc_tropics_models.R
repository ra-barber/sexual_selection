###############################################################################
                 # Binary tropics models on the cluster  #
###############################################################################

# This script runs Bayesian ordinal models on the high performance cluster. It 
# uses 32 cores in parallel, but runs each model sequentially. It uses 100 chains
# per model, to be comparable with phylogenetic models used later in the pipeline
# for evolutionary regression.


# Packages to load.
library(magrittr)
#library(tictoc)
library(caper)
library(dplyr)
library(effectsize)
#library(car)
#library(ggplot2)   # Check you can remove these packages.
#library(ggpubr) 
library(phytools)
library(brms)
library(graph4lg)
library(janitor)
library(readxl)

# Set the seed.
set.seed(1993)

################################################################################
                    #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Models that need both high cert and full data.
model_type <- c("allbirds","primary", "fruit", "secondary", "invert", 
                "mig", "no_mig", "terr", "no_terr")

# Centered or uncentered.
center <- c("uncentered", "centered")

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(center, data_type, model_type)

# Extra array specific info.
center <- all_combos[array_number, 1] %>% as.character()
data_type <- all_combos[array_number, 2] %>% as.character()
model_type <- all_combos[array_number, 3] %>% as.character()



###############################################################################
                       #### Read in the data #####

# Functions.
source("Code/functions.R")

# Read in the life history traits.
model_data <- read_excel("Data/sexual_selection_dataset_04_09.xlsx", sheet = 2, na = "NA") %>% 
  clean_names()
model_data$tree_tip <- gsub(" ", "_", model_data$scientific_name_bird_tree)

###############################################################################
              #### Prepare predictor variables ######

# Prepare response variables.
model_data$sexual_selection <- model_data$sexual_selection + 1

# Add the tropical non-tropical models.
model_data$trop_non_trop <- NA
model_data$trop_non_trop[abs(model_data$latitude) < 23.43624] <- "trop"
model_data$trop_non_trop[is.na(model_data$trop_non_trop)] <- "non_trop"



###############################################################################
           ### Run models on all birds (not shown in publication) ##


# all_model <- trop_brms_model(data_set = model_data)
# saveRDS(all_model, "Results/Models/Nonphy_models/Tropics/all_model.rds")
# 
# high_model <- trop_brms_model(data_set = model_data %>% filter(data_certainty > 2))
# saveRDS(high_model, "Results/Models/Nonphy_models/Tropics/high_model.rds")
# 



# Filter for high certainty if it's for the conservative analysis.
if (data_type == "high"){
  model_data %<>% filter(data_certainty > 2)
}

# Filter for dietary data.
primary_data <- model_data %>% filter(trophic_level_binary == "Primary")
secondary_data <- model_data %>% filter(trophic_level_binary == "Secondary")
fruit_data <- model_data %>% filter(trophic_niche == "Frugivore")
invert_data <- model_data %>% filter(trophic_niche == "Invertivore")

# Filter for eco roles.
mig_data <- model_data %>% filter(migration_binary == "Strong")
non_mig_data <- model_data %>% filter(migration_binary == "Weak")
terr_data <- model_data %>% filter(territoriality_binary == "Territorial")
non_terr_data <- model_data %>% filter(territoriality_binary == "Non-territorial")

# Pull out the specific data partition.
datasets <- list(model_data, primary_data, fruit_data, secondary_data, invert_data,
                 mig_data, non_mig_data, terr_data, non_terr_data)
names(datasets) <-  c("allbirds", "primary", "fruit", "secondary", "invert",
                  "mig", "no_mig", "terr", "no_terr")
model_data <- datasets[[model_type]]

# Center predictor after filtering.
model_data %<>% mutate(
  trop_non_trop = relevel(as.factor(trop_non_trop), ref = "trop"),
  trop_non_trop_c = center_categorical(trop_non_trop))

# Get the predictor if it's centered or un-centered. 
# We don't used the centered results in the publication, as i think centering categorical preditors is 
# still too novel that it adds an extra layer of confusion for those unfamiliar. 
if (center =="uncentered"){
  model_predictor <- "trop_non_trop"
} else {
  model_predictor <- "trop_non_trop_c"
}



# ###############################################################################
#                 #### Run main publication models ######
# 
# # Tropical brms function.
# trop_brms_model <- function(data_set = model_data, response = "sexual_selection",  
#                             predictor = "trop_non_trop", family = "cumulative"){
#   
#   # Center categorical predictors.
#   data_set %<>%  mutate(
#     trop_non_trop = relevel(as.factor(trop_non_trop), ref = "trop"),
#     trop_non_trop_c = center_categorical(trop_non_trop))
#   
#   # Create model formula.
#   model_formula <- formula(paste0(response, " ~ ", predictor))
#   
#   # Add the model family.
#   if (family == "cumulative"){
#     brms_formula <- brmsformula(model_formula, family = cumulative())
#   } else {
#     brms_formula <- brmsformula(model_formula, family = bernoulli())
#   }
#   
#   # Priors
#   linear_priors <- c(prior(normal(0,1), class="Intercept"),
#                      prior(normal(0,1), class="b"))
#   # Run brms models.
#   brm(
#     brms_formula,
#     data = data_set,
#     prior = linear_priors,
#     iter = 10000,
#     warmup = 5000,
#     chains = 100,
#     thin = 20,
#     cores = 32,
#     init = 0,
#     #file = model_pathway,
#     normalize = FALSE,
#     backend = "cmdstanr",
#     threads = threading(16)
#   )
# }

# # Run the models 
# tropical_model <- trop_brms_model(predictor = model_predictor)

# Create model formula.
model_formula <- formula(paste0(model_predictor, " ~ ", predictor))

# Add the model family.
brms_formula <- brmsformula(model_formula, family = cumulative())

# Priors.
linear_priors <- c(prior(normal(0,1), class="Intercept"),
                   prior(normal(0,1), class="b"))

# Run brms models.
tropical_model <- brm(
  brms_formula,
  data = model_data,
  prior = linear_priors,
  iter = 10000,
  warmup = 5000,
  chains = 100,
  thin = 20,
  cores = 32,
  init = 0,
  normalize = FALSE,
  backend = "cmdstanr",
  threads = threading(16))

# Model pathway.
first_half <- "Results/Models/Nonphy_models/Tropics/"
model_pathway <- paste0(first_half, model_type, "_", center, "_", data_type, "_model.rds")

# Export the model.
saveRDS(tropical_model, model_pathway)



#################################################################################
                           ##### End #######







# 
# # Run the niche models.
# primary_model <- trop_brms_model(data_set = primary_data)
# secondary_model <- trop_brms_model(data_set = secondary_data)
# fruit_model <- trop_brms_model(data_set = fruit_data)
# invert_model <- trop_brms_model(data_set = invert_data)
# 
# mig_model <- trop_brms_model(data_set = mig_data)
# non_mig_model <- trop_brms_model(data_set = non_mig_data)
# terr_model <- trop_brms_model(data_set = terr_data)
# non_terr_model <- trop_brms_model(data_set = non_terr_data)
# 
# # Export the models.
# saveRDS(primary_model, "Results/Models/Nonphy_models/Tropics/primary_model.rds")
# saveRDS(secondary_model, "Results/Models/Nonphy_models/Tropics/secondary_model.rds")
# saveRDS(fruit_model, "Results/Models/Nonphy_models/Tropics/fruit_model.rds")
# saveRDS(invert_model, "Results/Models/Nonphy_models/Tropics/invert_model.rds")
# 
# saveRDS(mig_model, "Results/Models/Nonphy_models/Tropics/mig_model.rds")
# saveRDS(non_mig_model, "Results/Models/Nonphy_models/Tropics/non_mig_model.rds")
# saveRDS(terr_model, "Results/Models/Nonphy_models/Tropics/terr_model.rds")
# saveRDS(non_terr_model, "Results/Models/Nonphy_models/Tropics/non_terr_model.rds")
# 
# # Run the models with standardised predictors.
# primary_model <- trop_brms_model(data_set = primary_data, predictor = "trop_non_trop_c")
# secondary_model <- trop_brms_model(data_set = secondary_data, predictor = "trop_non_trop_c")
# fruit_model <- trop_brms_model(data_set = fruit_data, predictor = "trop_non_trop_c")
# invert_model <- trop_brms_model(data_set = invert_data, predictor = "trop_non_trop_c")
# 
# mig_model <- trop_brms_model(data_set = mig_data, predictor = "trop_non_trop_c")
# non_mig_model <- trop_brms_model(data_set = non_mig_data, predictor = "trop_non_trop_c")
# terr_model <- trop_brms_model(data_set = terr_data, predictor = "trop_non_trop_c")
# non_terr_model <- trop_brms_model(data_set = non_terr_data, predictor = "trop_non_trop_c")
# 
# 
# # Export the models.
# saveRDS(primary_model, "Results/Models/Nonphy_models/Tropics/centered_primary_model.rds")
# saveRDS(secondary_model, "Results/Models/Nonphy_models/Tropics/centered_secondary_model.rds")
# saveRDS(fruit_model, "Results/Models/Nonphy_models/Tropics/centered_fruit_model.rds")
# saveRDS(invert_model, "Results/Models/Nonphy_models/Tropics/centered_invert_model.rds")
# 
# saveRDS(mig_model, "Results/Models/Nonphy_models/Tropics/centered_mig_model.rds")
# saveRDS(non_mig_model, "Results/Models/Nonphy_models/Tropics/centered_non_mig_model.rds")
# saveRDS(terr_model, "Results/Models/Nonphy_models/Tropics/centered_terr_model.rds")
# saveRDS(non_terr_model, "Results/Models/Nonphy_models/Tropics/centered_non_terr_model.rds")
# 
# 
# ################################################################################
#                    ##### High certainty models #########
# 
# # Filter for high cert.
# model_data %<>% filter(data_certainty > 2)
# 
# # Filter for dietary data.
# primary_data <- model_data %>% filter(trophic_level_binary == "Primary")
# secondary_data <- model_data %>% filter(trophic_level_binary == "Secondary")
# fruit_data <- model_data %>% filter(trophic_niche == "Frugivore")
# invert_data <- model_data %>% filter(trophic_niche == "Invertivore")
# 
# # Filter for eco roles.
# mig_data <- model_data %>% filter(migration_binary == "Strong")
# non_mig_data <- model_data %>% filter(migration_binary == "Weak")
# terr_data <- model_data %>% filter(territoriality_binary == "Territory")
# non_terr_data <- model_data %>% filter(territoriality_binary == "No territory")
# 
# 
# # Run the niche models.
# primary_model <- trop_brms_model(data_set = primary_data)
# secondary_model <- trop_brms_model(data_set = secondary_data)
# fruit_model <- trop_brms_model(data_set = fruit_data)
# invert_model <- trop_brms_model(data_set = invert_data)
# 
# mig_model <- trop_brms_model(data_set = mig_data)
# non_mig_model <- trop_brms_model(data_set = non_mig_data)
# terr_model <- trop_brms_model(data_set = terr_data)
# non_terr_model <- trop_brms_model(data_set = non_terr_data)
# 
# # Export the models.
# saveRDS(primary_model, "Results/Models/Nonphy_models/Tropics/high_primary_model.rds")
# saveRDS(secondary_model, "Results/Models/Nonphy_models/Tropics/high_secondary_model.rds")
# saveRDS(fruit_model, "Results/Models/Nonphy_models/Tropics/high_fruit_model.rds")
# saveRDS(invert_model, "Results/Models/Nonphy_models/Tropics/high_invert_model.rds")
# 
# saveRDS(mig_model, "Results/Models/Nonphy_models/Tropics/high_mig_model.rds")
# saveRDS(non_mig_model, "Results/Models/Nonphy_models/Tropics/high_non_mig_model.rds")
# saveRDS(terr_model, "Results/Models/Nonphy_models/Tropics/high_terr_model.rds")
# saveRDS(non_terr_model, "Results/Models/Nonphy_models/Tropics/high_non_terr_model.rds")
# 
# # Run the models with standardised predictors.
# primary_model <- trop_brms_model(data_set = primary_data, predictor = "trop_non_trop_c")
# secondary_model <- trop_brms_model(data_set = secondary_data, predictor = "trop_non_trop_c")
# fruit_model <- trop_brms_model(data_set = fruit_data, predictor = "trop_non_trop_c")
# invert_model <- trop_brms_model(data_set = invert_data, predictor = "trop_non_trop_c")
# 
# mig_model <- trop_brms_model(data_set = mig_data, predictor = "trop_non_trop_c")
# non_mig_model <- trop_brms_model(data_set = non_mig_data, predictor = "trop_non_trop_c")
# terr_model <- trop_brms_model(data_set = terr_data, predictor = "trop_non_trop_c")
# non_terr_model <- trop_brms_model(data_set = non_terr_data, predictor = "trop_non_trop_c")
# 
# 
# # Export the models.
# saveRDS(primary_model, "Results/Models/Nonphy_models/Tropics/high_centered_primary_model.rds")
# saveRDS(secondary_model, "Results/Models/Nonphy_models/Tropics/high_centered_secondary_model.rds")
# saveRDS(fruit_model, "Results/Models/Nonphy_models/Tropics/high_centered_fruit_model.rds")
# saveRDS(invert_model, "Results/Models/Nonphy_models/Tropics/high_centered_invert_model.rds")
# 
# saveRDS(mig_model, "Results/Models/Nonphy_models/Tropics/high_centered_mig_model.rds")
# saveRDS(non_mig_model, "Results/Models/Nonphy_models/Tropics/high_centered_non_mig_model.rds")
# saveRDS(terr_model, "Results/Models/Nonphy_models/Tropics/high_centered_terr_model.rds")
# saveRDS(non_terr_model, "Results/Models/Nonphy_models/Tropics/high_centered_non_terr_model.rds")
# 






