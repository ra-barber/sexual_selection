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

# Set the seed.
set.seed(1993)


###############################################################################
                       #### Read in the data #####

# Functions.
source("Code/functions.R")

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)
model_data$abs_lat <- abs(model_data$complete_latitude)

###############################################################################
              #### Prepare predictor variables ######

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE))

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1

###############################################################################
                #### Define brms model function. ######

lat_brms_model <- function(data_set = model_data, response = "sexual_score",  
                           predictor = "abs_lat", family = "cumulative"){
  
  
  # Scale continuous predictors to two SD.
  data_set %<>% mutate(
    centroid_z = standardize(centroid_sqrt, two_sd = TRUE))
  
  # Create model formula.
  model_formula <- formula(paste0(response, " ~ ", predictor))
  
  # Add the model family.
  if (family == "cumulative"){
    brms_formula <- brmsformula(model_formula, family = cumulative(threshold = "equidistant"))
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
    iter = 50000,
    warmup = 25000,
    chains = 2,
    thin = 20,
    cores = 32,
    init = 0,
    #file = model_pathway,
    normalize = FALSE,
    backend = "cmdstanr",
    threads = threading(16)
  )
}

###############################################################################
              #### Run main publication models ######


# Filter for primary and secondary data.
primary_data <- model_data %>% filter(trophic_binary == "Primary")
secondary_data <- model_data %>% filter(trophic_binary == "Secondary")
fruit_data <- model_data %>% filter(trophic_niche == "Frugivore")
invert_data <- model_data %>% filter(trophic_niche == "Invertivore")
gran_data <- model_data %>% filter(trophic_niche == "Granivore")

# Run the main paper models.
allbirds_model <- lat_brms_model(response = "sexual_score")
cert_model <- lat_brms_model(response = "cert_reverse")
primary_model <- lat_brms_model(data_set = primary_data)
secondary_model <- lat_brms_model(data_set = secondary_data)
fruit_model <- lat_brms_model(data_set = fruit_data)
invert_model <- lat_brms_model(data_set = invert_data)
gran_model <- lat_brms_model(data_set = gran_data)

# Export the models.
saveRDS(allbirds_model, "Results/Models/Nonphy_models/Latitude/allbirds_model.rds")
saveRDS(cert_model, "Results/Models/Nonphy_models/Latitude/cert_model.rds")
saveRDS(primary_model, "Results/Models/Nonphy_models/Latitude/primary_model.rds")
saveRDS(secondary_model, "Results/Models/Nonphy_models/Latitude/secondary_model.rds")
saveRDS(fruit_model, "Results/Models/Nonphy_models/Latitude/fruit_model.rds")
saveRDS(invert_model, "Results/Models/Nonphy_models/Latitude/invert_model.rds")
saveRDS(gran_model, "Results/Models/Nonphy_models/Latitude/gran_model.rds")


# Run the models with standardised predictors.
allbirds_model <- lat_brms_model(response = "sexual_score", predictor = "centroid_z")
cert_model <- lat_brms_model(response = "cert_reverse", predictor = "centroid_z")
primary_model <- lat_brms_model(data_set = primary_data, predictor = "centroid_z")
secondary_model <- lat_brms_model(data_set = secondary_data, predictor = "centroid_z")
fruit_model <- lat_brms_model(data_set = fruit_data, predictor = "centroid_z")
invert_model <- lat_brms_model(data_set = invert_data, predictor = "centroid_z")
gran_model <- lat_brms_model(data_set = gran_data, predictor = "centroid_z")

# Export the models.
saveRDS(allbirds_model, "Results/Models/Nonphy_models/Latitude/centered_allbirds_model.rds")
saveRDS(cert_model, "Results/Models/Nonphy_models/Latitude/centered_cert_model.rds")
saveRDS(primary_model, "Results/Models/Nonphy_models/Latitude/centered_primary_model.rds")
saveRDS(secondary_model, "Results/Models/Nonphy_models/Latitude/centered_secondary_model.rds")
saveRDS(fruit_model, "Results/Models/Nonphy_models/Latitude/centered_fruit_model.rds")
saveRDS(invert_model, "Results/Models/Nonphy_models/Latitude/centered_invert_model.rds")
saveRDS(gran_model, "Results/Models/Nonphy_models/Latitude/centered_gran_model.rds")


###############################################################################
            #### Run supplementary materials models ######


# Filter for primary and secondary data.
med_cert_data <- model_data %>% filter(cert_reverse > 2)
hi_cert_data <- model_data %>% filter(cert_reverse > 3)

# Run the supplementary models.
med_cert_model <- lat_brms_model(response = "sexual_score", data_set = med_cert_data)
hi_cert_model <- lat_brms_model(response = "sexual_score", data_set = hi_cert_data)
terr_primary_model <- lat_brms_model(response = "terr_dummy", data_set = primary_data, family = "bernoulli")
yearterr_primary_model <- lat_brms_model(response = "year_terr_dummy", data_set = primary_data, family = "bernoulli")
terr_secondary_model <- lat_brms_model(response = "terr_dummy", data_set = secondary_data, family = "bernoulli")
yearterr_secondary_model <- lat_brms_model(response = "year_terr_dummy", data_set = secondary_data, family = "bernoulli")


# Export the models.
saveRDS(med_cert_model, "Results/Models/Nonphy_models/Latitude/med_cert_model.rds")
saveRDS(hi_cert_model, "Results/Models/Nonphy_models/Latitude/hi_cert_model.rds")
saveRDS(terr_primary_model, "Results/Models/Nonphy_models/Latitude/terr_primary_model.rds")
saveRDS(yearterr_primary_model, "Results/Models/Nonphy_models/Latitude/yearterr_primary_model.rds")
saveRDS(terr_secondary_model, "Results/Models/Nonphy_models/Latitude/terr_secondary_model.rds")
saveRDS(yearterr_secondary_model, "Results/Models/Nonphy_models/Latitude/yearterr_secondary_model.rds")

# Run the supplementary models.
centered_med_cert_model <- lat_brms_model(response = "sexual_score", data_set = med_cert_data, predictor = "centroid_z")
centered_hi_cert_model <- lat_brms_model(response = "sexual_score", data_set = hi_cert_data, predictor = "centroid_z")
centered_terr_primary_model <- lat_brms_model(response = "terr_dummy", data_set = primary_data, family = "bernoulli", predictor = "centroid_z")
centered_yearterr_primary_model <- lat_brms_model(response = "year_terr_dummy", data_set = primary_data, family = "bernoulli", predictor = "centroid_z")
centered_terr_secondary_model <- lat_brms_model(response = "terr_dummy", data_set = secondary_data, family = "bernoulli", predictor = "centroid_z")
centered_yearterr_secondary_model <- lat_brms_model(response = "year_terr_dummy", data_set = secondary_data, family = "bernoulli", predictor = "centroid_z")

# Export the models.
saveRDS(centered_med_cert_model, "Results/Models/Nonphy_models/Latitude/centered_med_cert_model.rds")
saveRDS(centered_hi_cert_model, "Results/Models/Nonphy_models/Latitude/centered_hi_cert_model.rds")
saveRDS(centered_terr_primary_model, "Results/Models/Nonphy_models/Latitude/centered_terr_primary_model.rds")
saveRDS(centered_yearterr_primary_model, "Results/Models/Nonphy_models/Latitude/centered_yearterr_primary_model.rds")
saveRDS(centered_terr_secondary_model, "Results/Models/Nonphy_models/Latitude/centered_terr_secondary_model.rds")
saveRDS(centered_yearterr_secondary_model, "Results/Models/Nonphy_models/Latitude/centered_yearterr_secondary_model.rds")



###############################################################################
         #### Run extended data models for territory + migration ####


# Filter for eco roles.
mig_data <- model_data %>% filter(migration_binary == "Strong")
no_mig_data <- model_data %>% filter(migration_binary == "Weak")
terr_data <- model_data %>% filter(territory_binary == "Territory")
no_terr_data <- model_data %>% filter(territory_binary == "No territory")

# Run the models.
mig_model <- lat_brms_model(data_set = mig_data)
no_mig_model <- lat_brms_model(data_set = no_mig_data)
terr_model <- lat_brms_model(data_set = terr_data)
no_terr_model <- lat_brms_model(data_set = no_terr_data)

# Export the models.
saveRDS(mig_model, "Results/Models/Nonphy_models/Latitude/mig_model.rds")
saveRDS(no_mig_model, "Results/Models/Nonphy_models/Latitude/no_mig_model.rds")
saveRDS(terr_model, "Results/Models/Nonphy_models/Latitude/terr_model.rds")
saveRDS(no_terr_model, "Results/Models/Nonphy_models/Latitude/no_terr_model.rds")

# Run the models with standardised predictors.

mig_model <- lat_brms_model(data_set = mig_data, predictor = "centroid_z")
no_mig_model <- lat_brms_model(data_set = no_mig_data, predictor = "centroid_z")
terr_model <- lat_brms_model(data_set = terr_data, predictor = "centroid_z")
no_terr_model <- lat_brms_model(data_set = no_terr_data, predictor = "centroid_z")

# Export the models.
saveRDS(mig_model, "Results/Models/Nonphy_models/Latitude/centered_mig_model.rds")
saveRDS(no_mig_model, "Results/Models/Nonphy_models/Latitude/centered_no_mig_model.rds")
saveRDS(terr_model, "Results/Models/Nonphy_models/Latitude/centered_terr_model.rds")
saveRDS(no_terr_model, "Results/Models/Nonphy_models/Latitude/centered_no_terr_model.rds")



