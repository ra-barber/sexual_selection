###############################################################################
                     # Latitudinal brms models on cluster  #
###############################################################################

# Packages to load.
library(magrittr)
library(tictoc)
library(caper)
library(dplyr)
library(effectsize)
library(car)
library(ggplot2)
library(ggpubr)
library(phytools)
library(brms)
library(graph4lg)
library(readxl)
library(janitor)

# Set the seed.
set.seed(1993)


################################################################################
                 #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Models that need both high cert and full data.
model_type <- c("allbirds", "certainty", "primary", "fruit", "secondary", "invert", 
                "mig", "no_mig", "terr", "no_terr")

# Centered or uncentered.
center <- c("uncentered")

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(center, data_type, model_type)

# Models that don't need 2nd high cert sensitivity model.
model_type <- c("highcert", "prim_allterr", "prim_yearterr", "sec_allterr", "sec_yearterr")
data_type <- c("all")

# Add extra models.
extra_combos <- expand.grid(center, data_type, model_type)
all_combos <- rbind(all_combos, extra_combos)

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
model_data$abs_lat <- abs(model_data$latitude)

# Filter for high cert.
if (data_type == "high"){
  model_data %<>% filter(data_certainty > 2)
}


###############################################################################
              #### Prepare predictor variables / data ######

# Prepare response variables.
model_data$sexual_selection <- model_data$sexual_selection + 1

# Trophic niche model data.
primary_data <- model_data %>% filter(trophic_level_binary == "Primary")
secondary_data <- model_data %>% filter(trophic_level_binary == "Secondary")
fruit_data <- model_data %>% filter(trophic_niche == "Frugivore")
invert_data <- model_data %>% filter(trophic_niche == "Invertivore")

# Filter for eco roles.
mig_data <- model_data %>% filter(migration_binary == "Strong")
no_mig_data <- model_data %>% filter(migration_binary == "Weak")
terr_data <- model_data %>% filter(territoriality_binary == "Territory")
no_terr_data <- model_data %>% filter(territoriality_binary == "No territory")

# Highest certainty data.
hi_cert_data <- model_data %>% filter(data_certainty > 3)

# Create a list of datasets for easy reference.
all_datasets <- list(model_data, primary_data, fruit_data, secondary_data, invert_data,
                     mig_data, no_mig_data, terr_data, no_terr_data, hi_cert_data)
names(all_datasets) <- c("allbirds", "primary", "fruit", "secondary", "invert", 
                         "mig", "no_mig", "terr", "no_terr", "highcert")

###############################################################################
                #### Define brms model function. ######

lat_brms_model <- function(data_set = model_data, response = "sexual_selection",  
                           predictor = "abs_lat", family = "cumulative"){
  
  # Scale continuous predictors to two SD.
  data_set %<>% mutate(
    centroid_z = standardize(sqrt(abs_lat), two_sd = TRUE))
  
  # Create model formula.
  model_formula <- formula(paste0(response, " ~ ", predictor))
  
  # Add the model family.
  if (family == "cumulative"){
    brms_formula <- brmsformula(model_formula, family = cumulative())
  } else {
    brms_formula <- brmsformula(model_formula, family = bernoulli())
  }
  
  # Priors
  linear_priors <- c(prior(normal(0,1), class="Intercept"),
                     prior(normal(0,1), class="b"))
  # Run brms models.
  brm(
    brms_formula,
    data = data_set,
    prior = linear_priors,
    iter = 10000,
    warmup = 5000,
    chains = 100,
    thin = 20,
    cores = 32,
    init = 0,
    #file = model_pathway,
    normalize = FALSE,
    backend = "cmdstanr"
  )
}

###############################################################################
              #### Run main publication models ######


# Run the straightforward sexual selection models first.
if (model_type %in% names(all_datasets)){
  latitude_model <- lat_brms_model(data_set = all_datasets[[model_type]])
}

# Run the global data certainty analysis, with certainty as the response variable.
if (model_type == "certainty"){
  latitude_model <- lat_brms_model(response = "data_certainty")
}

# Make dummy territoriality metrics.
model_data$terr_dummy <- "0"
model_data$year_terr_dummy <- "0"
model_data %>% mutate(
  terr_dummy = replace(terr_dummy, territoriality_binary == "Non-Territorial", "1"),
  year_terr_dummy = replace(year_terr_dummy, territoriality == "Strong", "1"))
model_data$terr_dummy %<>% as.numeric()
model_data$year_terr_dummy %<>% as.numeric()

# Run territoriality models.
if (model_type == "prim_allterr"){
  latitude_model <- lat_brms_model(response = "terr_dummy", data_set = primary_data, family = "bernoulli")
}
if (model_type == "prim_yearterr"){
  latitude_model <- lat_brms_model(response = "year_terr_dummy", data_set = primary_data, family = "bernoulli")
}
if (model_type == "sec_allterr"){
  latitude_model <- lat_brms_model(response = "terr_dummy", data_set = secondary_data, family = "bernoulli")
}
if (model_type == "sec_yearterr"){
  latitude_model <- lat_brms_model(response = "year_terr_dummy", data_set = secondary_data, family = "bernoulli")
}

# Model pathway.
first_half <- "Results/Models/Nonphy_models/Latitude/"
model_pathway <- paste0(first_half, model_type, "_", data_type, "_model.rds")

# Export the model.
saveRDS(latitude_model, model_pathway)
