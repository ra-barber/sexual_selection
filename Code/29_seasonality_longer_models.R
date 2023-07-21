###############################################################################
                     # Simple brms models on cluster  #
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
tree_number <- 1:10

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(tree_number, data_type)

# Tree type.
tree_number <- all_combos[array_number, 1] %>% as.numeric()

# Data type.
data_type <- all_combos[array_number, 2] %>% as.character()


###############################################################################
                       #### Read in the data #####


# Functions.
source("Code/functions.R")

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[tree_number]]

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)

# Filter for high cert.
if (data_type == "high"){
  model_data %<>% filter(sexual_certainty < 3)
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
  territory_binary = relevel(as.factor(territory_binary), ref = "No territory"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_binary = relevel(as.factor(trophic_binary), ref = "Secondary"),
  devo_binary = relevel(as.factor(devo_mode_wang), ref = "altricial"),
  devo_mode = relevel(as.factor(devo_mode), ref = "altrical"),
  wang_edited = relevel(as.factor(wang_edited), ref = "altricial"),
  diet_binary = relevel(as.factor(diet_binary), ref = "Non-frug-nect")
)

# Center categorical predictors.
model_data %<>% mutate(
  terr_bi_c = center_categorical(territory_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary),
  devo_bi_c = center_categorical(devo_binary),
  devo_ed_c = center_categorical(wang_edited),
  devo_mode_c = center_categorical(devo_mode),
  diet_c = center_categorical(diet_binary)
)

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  body_mass_log_z = standardize(body_mass_log, two_sd = TRUE),
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE),
  temp_seasonality_z = standardize(temp_log, two_sd = TRUE),
  npp_z = standardize(npp_sqrt, two_sd = TRUE),
  gen_z = standardize(gen_log, two_sd = TRUE),
  dens_z = standardize(dens_sqrt, two_sd = TRUE),
  fed_z = standardize(fed_sqrt, two_sd = TRUE),
  chick_z = standardize(chick_pc1, two_sd = TRUE),
  chick_sqrt_z = standardize(chick_sqrt, two_sd = TRUE)
)

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1


###############################################################################
                #### Set model formula ######


model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + trophic_level_c + 
terr_bi_c*temp_seasonality_z + migration_bi_c*temp_seasonality_z + trophic_level_c*temp_seasonality_z + (1|gr(tree_tip, cov=A))"


# brms formula.
brms_formula <- brmsformula(model_formula, family = cumulative(), decomp = "QR")

# Simple models.
model_pathway <- paste0("Results/Models/Seasonality/", data_type, "_", tree_number, ".rds") 

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
    cores = 32,
    init = 0,
    file = model_pathway,
    normalize = FALSE,
    backend = "cmdstanr",
    #control = list(adapt_delta = 0.6),
    threads = threading(16),
  )
