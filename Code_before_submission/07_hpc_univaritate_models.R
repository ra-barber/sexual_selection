###############################################################################
               # Running univariate phylogenetic models  #
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
# library(posterior)
# library(bayesplot)


# Clear the workspace.
rm(list=ls())


################################################################################
                  #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Create types for hpc jobs.
tree_number <- 1:50

# Set the data types.
data_type <- c("all", "high")

# Create the response variable.
response <- c("temp", "mig", "tro", "terr")

# Centered or uncentered.
center <- c("centered", "uncentered")

# Expand the grid.
all_combos <- expand.grid(response, center, tree_number, data_type)

# Reponse.
response <- all_combos[array_number, 1] %>% as.character()

# Centered.
center <- all_combos[array_number, 2] %>% as.character()

# Data type.
data_type <- all_combos[array_number, 3] %>% as.character()

# Tree type.
tree_number <- all_combos[array_number, 4] %>% as.numeric()

# Time for array 1 - 400
# 72 hours

# Time for array 401 - 800
# 48 Hours




###############################################################################
                       #### Read in the data #####


# Functions.
source("Code/functions.R")

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[tree_number]]

# Read in the life history traits.
model_data <- read_ss_data()

# Filter for high cert.
if (data_type == "high"){
  model_data %<>% filter(data_certainty > 2)
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

# Set as factor, then re-level for appropriate reference group.
model_data %<>% mutate(
  territoriality_binary = relevel(as.factor(territoriality_binary), ref = "Non-territorial"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_level_binary = relevel(as.factor(trophic_level_binary), ref = "Secondary"))

# Center categorical predictors.
model_data %<>% mutate(
  terr_bi_c = center_categorical(territoriality_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_level_binary))

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  temp_log = log(seasonality),
  temp_seasonality_z = standardize(temp_log, two_sd = TRUE))

# Prepare response variables.
model_data$sexual_selection <- model_data$sexual_selection + 1


###############################################################################
                #### Set model formula ######

# If response variable is centered.
if (center == "centered"){
  centered_responses <- c("temp_seasonality_z", "migration_bi_c", "trophic_level_c", "terr_bi_c")
  names(centered_responses) <- c("temp", "mig", "tro", "terr")
  model_response <- centered_responses[[response]]
} else {
  centered_responses <- c("temp_log", "migration_binary", "trophic_level_binary", "territoriality_binary")
  names(centered_responses) <- c("temp", "mig", "tro", "terr")
  model_response <- centered_responses[[response]]
}

# Model formula.
model_formula <- paste0("sexual_selection ~ ", model_response, " + (1|gr(tree_tip, cov=A))")
 
# brms formula.
brms_formula <- brmsformula(model_formula, family = cumulative())

# Simple models.
model_pathway <- paste0("Results/Models/Univariate/", response, "_", center, "_", data_type, "_", tree_number, ".rds") 

#library(standist) ~ for visualising priors.
# # Add un-informative priors.
normal_priors <- c(prior(normal(0,1), class="Intercept"),
                   prior(normal(0,1), class="b"),
                   prior(gamma(2,1), "sd")) # Gamma 2,1 probs seems to work well given all previous models end up with values around 1

  # Run brms models.
  brms_model <- brm(
    brms_formula,
    data = model_data,
    data2 = list(A=model_covar),
    prior = normal_priors,
    iter = 10000,
    warmup = 5000,
    chains = 2,
    thin = 20,
    cores = 16,
    init = 0,
    file = model_pathway,
    normalize = FALSE,
    backend = "cmdstanr",
    threads = threading(8))
