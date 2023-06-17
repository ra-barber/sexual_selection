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


###############################################################################
              #### Prepare predictor variables ######

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1

# Add the tropical non-tropical models.
model_data$trop_non_trop <- NA
model_data$trop_non_trop[abs(model_data$complete_latitude) < 23.43624] <- "trop"
model_data$trop_non_trop[is.na(model_data$trop_non_trop)] <- "non_trop"

# Tropical brms function.
trop_brms_model <- function(data_set = model_data, response = "sexual_score",  
                            predictor = "trop_non_trop", family = "cumulative"){
  
  # Center categorical predictors.
  data_set %<>%  mutate(
    trop_non_trop = relevel(as.factor(trop_non_trop), ref = "trop"),
    trop_non_trop_c = center_categorical(trop_non_trop))
  
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


# Filter for dietary data.
primary_data <- model_data %>% filter(trophic_binary == "Primary")
secondary_data <- model_data %>% filter(trophic_binary == "Secondary")
fruit_data <- model_data %>% filter(trophic_niche == "Frugivore")
invert_data <- model_data %>% filter(trophic_niche == "Invertivore")

# Filter for eco roles.
mig_data <- model_data %>% filter(migration_binary == "Strong")
non_mig_data <- model_data %>% filter(migration_binary == "Weak")
terr_data <- model_data %>% filter(territory_binary == "Territory")
non_terr_data <- model_data %>% filter(territory_binary == "No territory")


# Run the niche models.
primary_model <- trop_brms_model(data_set = primary_data)
secondary_model <- trop_brms_model(data_set = secondary_data)
fruit_model <- trop_brms_model(data_set = fruit_data)
invert_model <- trop_brms_model(data_set = invert_data)

mig_model <- trop_brms_model(data_set = mig_data)
non_mig_model <- trop_brms_model(data_set = non_mig_data)
terr_model <- trop_brms_model(data_set = terr_data)
non_terr_model <- trop_brms_model(data_set = non_terr_data)

# Export the models.
saveRDS(primary_model, "Results/Models/Nonphy_models/Tropics/primary_model.rds")
saveRDS(secondary_model, "Results/Models/Nonphy_models/Tropics/secondary_model.rds")
saveRDS(fruit_model, "Results/Models/Nonphy_models/Tropics/fruit_model.rds")
saveRDS(invert_model, "Results/Models/Nonphy_models/Tropics/invert_model.rds")

saveRDS(mig_model, "Results/Models/Nonphy_models/Tropics/mig_model.rds")
saveRDS(non_mig_model, "Results/Models/Nonphy_models/Tropics/non_mig_model.rds")
saveRDS(terr_model, "Results/Models/Nonphy_models/Tropics/terr_model.rds")
saveRDS(non_terr_model, "Results/Models/Nonphy_models/Tropics/non_terr_model.rds")

# Run the models with standardised predictors.
primary_model <- trop_brms_model(data_set = primary_data, predictor = "trop_non_trop_c")
secondary_model <- trop_brms_model(data_set = secondary_data, predictor = "trop_non_trop_c")
fruit_model <- trop_brms_model(data_set = fruit_data, predictor = "trop_non_trop_c")
invert_model <- trop_brms_model(data_set = invert_data, predictor = "trop_non_trop_c")

mig_model <- trop_brms_model(data_set = mig_data, predictor = "trop_non_trop_c")
non_mig_model <- trop_brms_model(data_set = non_mig_data, predictor = "trop_non_trop_c")
terr_model <- trop_brms_model(data_set = terr_data, predictor = "trop_non_trop_c")
non_terr_model <- trop_brms_model(data_set = non_terr_data, predictor = "trop_non_trop_c")


# Export the models.
saveRDS(primary_model, "Results/Models/Nonphy_models/Tropics/centered_primary_model.rds")
saveRDS(secondary_model, "Results/Models/Nonphy_models/Tropics/centered_secondary_model.rds")
saveRDS(fruit_model, "Results/Models/Nonphy_models/Tropics/centered_fruit_model.rds")
saveRDS(invert_model, "Results/Models/Nonphy_models/Tropics/centered_invert_model.rds")

saveRDS(mig_model, "Results/Models/Nonphy_models/Tropics/centered_mig_model.rds")
saveRDS(non_mig_model, "Results/Models/Nonphy_models/Tropics/centered_non_mig_model.rds")
saveRDS(terr_model, "Results/Models/Nonphy_models/Tropics/centered_terr_model.rds")
saveRDS(non_terr_model, "Results/Models/Nonphy_models/Tropics/centered_non_terr_model.rds")