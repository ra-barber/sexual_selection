###############################################################################
                 # Latitudinal brms models on cluster  #
###############################################################################

# Script to run brms latitudinal models on the cluster.

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
model_type <- c("allbirds",  "primary", "fruit", "secondary", "invert", 
                "mig", "nomig", "terr", "noterr")

# Response type.
response_type <- c("ordinal")

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(response_type, data_type, model_type)

# Add the extra sensitivity analysis on highest certainty. 
model_type <- c("highcert")
data_type <- c("all")
cert_combos <- expand.grid(response_type, data_type, model_type)

# And the specific analysis on certainty. 
cert_combos <- rbind(cert_combos, data.frame(Var1 = "ordinal", 
                                             Var2 = "all", Var3 = "certainty"))

# Models that only use logistsic, unrelated to sexual selection certainty.
model_type <- c("primterr", "primyearterr", "secterr", "secyearterr")
response_type <- c("logistic")
extra_combos <- expand.grid(response_type, data_type, model_type)

# Combine all possible combinations.
all_combos <- rbind(all_combos, cert_combos, extra_combos)

# Extra array specific info.
response_type <- all_combos[array_number, 1] %>% as.character()
data_type <- all_combos[array_number, 2] %>% as.character()
model_type <- all_combos[array_number, 3] %>% as.character()


###############################################################################
                       #### Read in the data #####

# Functions.
source("Code/functions.R")

# Read in the life history traits.
model_data <- read_ss_data()

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
nomig_data <- model_data %>% filter(migration_binary == "Weak")
terr_data <- model_data %>% filter(territoriality_binary == "Territorial")
noterr_data <- model_data %>% filter(territoriality_binary == "Non-territorial")

# Highest certainty data.
hi_cert_data <- model_data %>% filter(data_certainty > 3)

# Create a list of datasets for easy reference.
all_datasets <- list(model_data, primary_data, fruit_data, secondary_data, invert_data,
                     mig_data, nomig_data, terr_data, noterr_data, hi_cert_data)
names(all_datasets) <- c("allbirds", "primary", "fruit", "secondary", "invert", 
                         "mig", "nomig", "terr", "noterr", "highcert")

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
    warmup = 8000,
    chains = 4,
    cores = 8,
    file = model_pathway,
    normalize = FALSE,
    control = list(adapt_delta = 0.95),
    threads = threading(2),
    backend = "cmdstanr"
  )
}

###############################################################################
              #### Run latitudinal models ######


# Model pathway.
first_half <- "Results/Models/Latitudinal/Latitude/"
model_pathway <- paste0(first_half, response_type, "_", 
                        model_type, "_", data_type, ".rds")


# Run the straightforward sexual selection models first.
if (model_type %in% names(all_datasets)){
  model_data <- all_datasets[[model_type]]
    lat_brms_model(model_data)
}

# Run the global data certainty analysis, with certainty as the response variable.
if (model_type == "certainty"){
  lat_brms_model(response = "data_certainty")
}

# Run territoriality models.
if (model_type == "primterr"){
  lat_brms_model(response = "terr_dummy", 
                 data_set = primary_data, family = "bernoulli")
}
if (model_type == "primyearterr"){
  lat_brms_model(response = "year_terr_dummy", 
                 data_set = primary_data, family = "bernoulli")
}
if (model_type == "secterr"){
  lat_brms_model(response = "terr_dummy", 
                 data_set = secondary_data, family = "bernoulli")
}
if (model_type == "secyearterr"){
  lat_brms_model(response = "year_terr_dummy", 
                 data_set = secondary_data, family = "bernoulli")
}

