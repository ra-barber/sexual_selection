###############################################################################
                  # Logistic Latitudinal brms models on cluster  #
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
core_number <- as.numeric(Sys.getenv("NUM_CORES"))
core_number

# Create types for hpc jobs.
tree_number <- 1

# Model (In order of size and therefore speed)
model_type <- c("fruit", "primary","mig", "no_terr", "invert", "secondary", 
                "no_mig", "terr", "allbirds", "certainty")

# Centered or uncentered.
center <- c("uncentered", "centered")
center <- c("uncentered")

# Set the data types.
data_type <- c("all", "high")

# Response type.
response_type <- c("ordinal", "logistic")

# Expand the grid.
all_combos <- expand.grid(tree_number, center, data_type, response_type, model_type)

# Tree type.
tree_number <- all_combos[array_number, 1] %>% as.numeric()

# Center type.
center <- all_combos[array_number, 2] %>% as.character()

# Data type.
data_type <- all_combos[array_number, 3] %>% as.character()

# Response type.
response_type <- all_combos[array_number, 4] %>% as.character()

# Model type.
model_type <- all_combos[array_number, 5] %>% as.character()


###############################################################################
                     #### Read in the data #####


# Functions.
source("Code/functions.R")

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_consensus_tree.tre")

# Read in the life history traits.
model_data <- read_ss_data()

# Filter for high cert.
if (data_type == "high"){
  model_data %<>% filter(data_certainty > 2)
}

# Trophic niche model data.
primary_data <- model_data %>% filter(trophic_level_binary == "Primary")
secondary_data <- model_data %>% filter(trophic_level_binary == "Secondary")
fruit_data <- model_data %>% filter(trophic_niche == "Frugivore")
invert_data <- model_data %>% filter(trophic_niche == "Invertivore")

# Filter for eco roles.
mig_data <- model_data %>% filter(migration_binary == "Strong")
no_mig_data <- model_data %>% filter(migration_binary == "Weak")
terr_data <- model_data %>% filter(territoriality_binary == "Territorial")
no_terr_data <- model_data %>% filter(territoriality_binary == "Non-territorial")

# Create a list of datasets for easy reference.
all_datasets <- list(model_data, primary_data, fruit_data, secondary_data, invert_data,
                     mig_data, no_mig_data, terr_data, no_terr_data)
names(all_datasets) <- c("allbirds", "primary", "fruit", "secondary", "invert", 
                         "mig", "no_mig", "terr", "no_terr")

# Pull out the correct dataset for the model, and overwrite model data.
if (model_type != "certainty"){
  model_data <-  all_datasets[[model_type]]
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
  centroid_z = standardize(abs_lat, two_sd = TRUE))

# Prepare response variables.
model_data$sexual_selection <- model_data$sexual_selection + 1


###############################################################################
                   #### Set model formula ######

# Predictor variable.
if (center == "centered"){
  predictor <- "centroid_z"
} else {
  predictor <- "abs_lat"
}

# Response variable.
if (response_type == "ordinal"){
  response <- "sexual_selection"
} else {
  response <- "sexual_binary"
}
if (model_type == "certainty"){
  response <- "data_certainty"
}

# # Centered model formula.
# if (center == "centered"){
#   model_formula <- "sexual_selection ~ centroid_z + (1|gr(tree_tip, cov=A))"
# } else {
#   model_formula <- "sexual_selection ~ abs_lat + (1|gr(tree_tip, cov=A))"
# }
# 
# # For sexual certainty.
# if (model_type == "certainty"){
#   if (center == "centered"){
#     model_formula <- "data_certainty ~ centroid_z + (1|gr(tree_tip, cov=A))"
#   } else {
#     model_formula <- "data_certainty ~ abs_lat + (1|gr(tree_tip, cov=A))"
#   }
# }

# Paste it all together.
model_formula <- paste0(response, " ~ ", predictor, " + (1|gr(tree_tip, cov=A))")

# Response variable.
if (response_type == "ordinal"){
  brms_formula <- brmsformula(model_formula, family = cumulative())
} else {
  brms_formula <- brmsformula(model_formula, family = bernoulli())
}
if (model_type == "certainty"){
  brms_formula <- brmsformula(model_formula, family = cumulative())
}

# # brms formula.
# brms_formula <- brmsformula(model_formula, family = cumulative())
# brms_formula <- brmsformula(model_formula, family = bernoulli())


# Simple models.
#model_pathway <- paste0("Results/Models/Latitude/", model_type, "_", center, "_", data_type, "_", tree_number, ".rds") 

model_pathway <- paste0("Results/Models/Consensus/Latitude/", response_type, "_", model_type, "_", center, "_", data_type, "_", tree_number, ".rds") 


# #library(standist) ~ for visualising priors.
# # # Add un-informative priors.
# normal_priors <- c(prior(normal(0,1), class="Intercept"),
#                    prior(normal(0,1), class="b"),
#                    prior(gamma(2,1), "sd")) # Gamma 2,1 probs seems to work well given all previous models end up with values around 1
# 
# normal_priors <- c(prior(normal(0,10), class="Intercept"),
#                    prior(normal(0,10), class="b")) 
# 

# Report time before starting.
Sys.time()

# Run brms models.
brms_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  #prior = normal_priors,
  iter = 10000,
  warmup = 8000,
  chains = 4,
  #thin = 2,
  cores = core_number,
  #init = 0,
  file = model_pathway,
  control = list(adapt_delta = 0.95),
  normalize = FALSE,
  backend = "cmdstanr",
  threads = threading(core_number/4)
)

