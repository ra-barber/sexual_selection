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

# Clear the workspace.
rm(list=ls())


################################################################################
                  #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Create types for hpc jobs.
tree_number <- 1:50

# Model (In order of size and therefore speed)
model_type <- c("frugivore", "primary", "invertivore", "secondary", "all", "certainty")

# Centered or uncentered.
center <- c("uncentered")

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(tree_number, center, data_type, model_type)

# Tree type.
tree_number <- all_combos[array_number, 1] %>% as.numeric()

# Center type.
center <- all_combos[array_number, 2] %>% as.character()

# Data type.
data_type <- all_combos[array_number, 3] %>% as.character()

# Model type.
model_type <- all_combos[array_number, 4] %>% as.character()

# Put working directory as correct one.
#setwd("..")

###############################################################################
                       #### Read in the data #####


# Functions.
source("Code/functions.R")

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[tree_number]]

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)

# Remove missing data for certain models.
#model_data %<>% filter(!is.na(fed_sqrt))

# Filter for high cert.
if (data_type == "high"){
  model_data %<>% filter(sexual_certainty < 3)
}

# Filter for trophic level.
if (model_type == "primary"){
  model_data %<>% filter(trophic_binary == "Primary")
}
if (model_type == "secondary"){
  model_data %<>% filter(trophic_binary == "Secondary")
}
if (model_type == "frugivore"){
  model_data %<>% filter(trophic_niche == "Frugivore")
}
if (model_type == "invertivore"){
  model_data %<>% filter(trophic_niche == "Invertivore")
}
if (model_type == "nectarivore"){
  model_data %<>% filter(trophic_niche == "Invertivore")
}
if (model_type == "frug_nect"){
  model_data %<>% filter(trophic_niche == "Invertivore")
}

# Drop tips on the tree.
model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Make a covariance matrix, and order data the same.
model_covar <- ape::vcv.phylo(model_tree)

# Reorder the matrix so it's random, to maximise parallel processing speed.
mat_order <- sample(1:nrow(model_covar), size = nrow(model_covar), replace = FALSE)
model_covar <- reorder_mat(model_covar, rownames(model_covar)[mat_order])
row.names(model_data) <- model_data$tree_tip
model_data <- model_data[row.names(model_covar),]


###############################################################################
              #### Prepare predictor variables ######

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE),
  abs_lat = abs(complete_latitude)
)

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1


###############################################################################
                #### Set model formula ######

# Centered model formula.
if (center == "centered"){
  model_formula <- "sexual_score ~ centroid_z + (1|gr(tree_tip, cov=A))"
} else {
  model_formula <- "sexual_score ~ abs_lat + (1|gr(tree_tip, cov=A))"
}

# For sexual certainty.
if (model_type == "certainty"){
  if (center == "centered"){
    model_formula <- "cert_reverse ~ centroid_z + (1|gr(tree_tip, cov=A))"
  } else {
    model_formula <- "cert_reverse ~ abs_lat + (1|gr(tree_tip, cov=A))"
  }
}


# brms formula.
brms_formula <- brmsformula(model_formula, family = cumulative())

# Simple models.
model_pathway <- paste0("Results/Models/Latitude/", model_type, "_", center, "_", data_type, "_", tree_number, ".rds") 

#library(standist) ~ for visualising priors.
# # Add un-informative priors.
normal_priors <- c(prior(normal(0,1), class="Intercept"),
                   prior(normal(0,1), class="b"),
                   prior(gamma(2,1), "sd")) # Gamma 2,1 probs seems to work well given all previous models end up with values around 1

# Report time before starting.
Sys.time()

# Run brms models.
brms_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  prior = normal_priors,
  iter = 2500,
  warmup = 500,
  chains = 2,
  thin = 20,
  cores = 32,
  init = 0,
  file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  threads = threading(16)
)

