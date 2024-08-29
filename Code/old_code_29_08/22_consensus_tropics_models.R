###############################################################################
                   #  Tropical brms models on cluster  #
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

# Model (In order of size and therefore speed).
model_type <- c("fruit", "primary", "mig", "invert", "secondary", "terr",
                "no_mig",  "no_terr")

# Set the data certainty.
data_type <- c("all", "high")

# Response type.
response_type <- c("ordinal", "logistic")

# Expand the grid.
all_combos <- expand.grid(data_type, response_type, model_type)

# Extract array job variables.
data_type <- all_combos[array_number, 1] %>% as.character()
response_type <- all_combos[array_number, 2] %>% as.character()
model_type <- all_combos[array_number, 3] %>% as.character()


###############################################################################
                       #### Read in the data #####


# Functions.
source("Code/functions.R")

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_consensus_tree.tre")

# Read in the life history traits.
model_data <- read_ss_data()

# Filter for data certainty.
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

# Pull out the correct dataset for the model, and overwrite model data.
all_datasets <- list(model_data, primary_data, fruit_data, secondary_data, invert_data,
                     mig_data, no_mig_data, terr_data, no_terr_data)
names(all_datasets) <- c("allbirds", "primary", "fruit", "secondary", "invert", 
                         "mig", "no_mig", "terr", "no_terr")
model_data <-  all_datasets[[model_type]]

# Create covariance matrix.
model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))
model_covar <- ape::vcv.phylo(model_tree)

# Reorder the matrix so it's random, to maximise parallel processing speed.
mat_order <- sample(1:nrow(model_covar), size = nrow(model_covar), replace = FALSE)
model_covar <- reorder_mat(model_covar, rownames(model_covar)[mat_order])
row.names(model_data) <- model_data$tree_tip
model_data <- model_data[row.names(model_covar),]


###############################################################################
                #### Prepare predictor variables ######

# Prepare response variables.
model_data$sexual_selection <- model_data$sexual_selection + 1

# Add the tropical non-tropical models.
model_data$trop_non_trop <- NA
model_data$trop_non_trop[abs(model_data$latitude) < 23.43624] <- "trop"
model_data$trop_non_trop[is.na(model_data$trop_non_trop)] <- "non_trop"


###############################################################################
                    #### Run brms models ######


# Tropical brms function.
trop_brms_model <- function(data_set = model_data, response = "sexual_selection",  
                            predictor = "trop_non_trop", family = "cumulative"){
  
  # Center categorical predictors.
  data_set %<>%  mutate(
    trop_non_trop = relevel(as.factor(trop_non_trop), ref = "trop"))
  
  # Create model formula.
  model_formula <- formula(paste0(response, " ~ ", predictor, " + (1|gr(tree_tip, cov=A))"))
  
  # Add the model family.
  if (family == "cumulative"){
    brms_formula <- brmsformula(model_formula, family = cumulative())
  } else {
    brms_formula <- brmsformula(model_formula, family = bernoulli())
  }
  
  # Run brms models.
  brm(
    brms_formula,
    data = data_set,
    data2 = list(A=model_covar),
    #prior = linear_priors,
    iter = 10000,
    warmup = 8000,
    chains = 4,
    cores = 32,
    #init = 0,
    file = model_pathway,
    control = list(adapt_delta = 0.95),
    normalize = FALSE,
    backend = "cmdstanr",
    threads = threading(8)
  )
}

# Create the model pathway to save.
model_pathway <- paste0("Results/Models/Consensus/Tropics/", response_type, "_", model_type, "_", data_type, ".rds") 

# Run the niche models.
if (response_type == "ordinal"){
  trop_brms_model()
} else {
  trop_brms_model(response = "sexual_binary", family = "bernoulli")
}
