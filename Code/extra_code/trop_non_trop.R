###############################################################################
                        ##### Script Title #####
###############################################################################

# This script does something new. Pretty sick eh.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(skimr)
library(tictoc)
library(stringr)
library(caper)
library(dplyr)
library(janitor)
library(ggpubr)

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                             #### Data ####

# Read in some data.
data <- read.csv("")

# Read in a tree.
tree <- read.tree("")

# Read in a raster.
raster <- raster("")


###############################################################################
                           #### Section 1 ####

model_data$trop_non_trop <- NA

model_data$trop_non_trop[abs(model_data$complete_latitude) < 23.43624] <- "trop"
model_data$trop_non_trop[is.na(model_data$trop_non_trop)] <- "non_trop"
model_data %>% filter(trophic_niche == "Frugivore") %>% group_by(trop_non_trop) %>% summarise(mean_ss = mean(sexual_score))

# model_data$trop_non_trop[abs(model_data$complete_latitude) < 20] <- "trop"
# model_data$trop_non_trop[is.na(model_data$trop_non_trop)] <- "non_trop"
# 
# 
# model_data$trop_non_trop <- NA
# model_data$trop_non_trop[abs(model_data$complete_latitude) < 15] <- "trop"
# model_data$trop_non_trop[is.na(model_data$trop_non_trop)] <- "non_trop"
# model_data %>% filter(trophic_niche == "Frugivore") %>% count(trop_non_trop)
# model_data %>% filter(trophic_niche == "Frugivore") %>% group_by(trop_non_trop) %>% summarise(mean_ss = mean(sexual_score))

model_data %<>% filter(trophic_niche == "Frugivore")

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

trop_brms_model <- function(data_set = model_data, response = "sexual_score",  
                           predictor = "trop_non_trop", family = "cumulative"){
  
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
    iter = 10000,
    warmup = 2500,
    chains = 2,
    thin = 20,
    cores = 8,
    init = 0,
    #file = model_pathway,
    normalize = FALSE,
    backend = "cmdstanr",
    threads = threading(4)
  )
}

frug_model <- trop_brms_model()


###############################################################################
                           #### Section 2 ####


###############################################################################
                           #### Section 3 ####


###############################################################################
                           #### Section 4 ####

###############################################################################
                           #### Section 5 ####

###############################################################################
                           #### Section 6 ####


###############################################################################
                           #### Section 7 ####


###############################################################################
                           #### Section 8 ####

# Look Rob, you've had your fun with the sectioning. 
# They'll be no more sectioning today.


###############################################################################
                             #### END ####
###############################################################################


###############################################################################
              #### All the stuff I'm afraid to delete ####


