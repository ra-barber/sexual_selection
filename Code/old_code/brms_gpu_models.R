###############################################################################
# Phylo models for quick exploration  #
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
library(phyr)

# Clear the workspace.
rm(list=ls())

# Functions.
source("../Code/functions.R")


################################################################################
#### Set up array iteration ####


# Get the array number from the job script.
#array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number <- 4

# Create responses for hpc jobs.
responses <- c("sexual_score", "sexual_binary")

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(responses, data_type)

# Response.
response <- all_combos[array_number, 1] %>% as.character()

# Data type.
data_type <- all_combos[array_number, 2] %>% as.character()

# Put working directory as correct one.
setwd("..")



###############################################################################
#### Read in the data #####


# Read in the tree.
model_tree <- read.tree("Data/Trees/model_trees.tre")[[1]]

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$bird_tree_name)

if (data_type == "high"){
  model_data %<>% filter(sexual_certainty < 3)
}

# Drop tips on the tree.
model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Make a covariance matrix.
model_covar <- ape::vcv.phylo(model_tree)

# Make sure the covariance matrix and data are in the same order.
row.names(model_data) <- model_data$tree_tip
model_data <- model_data[row.names(model_covar),]

###############################################################################
#### Prepare predictor variables ######

# colnames(model_data)

# Set as factor, then re-level for appropriate reference group.
model_data %<>% mutate(
  territory_binary = relevel(as.factor(territory_binary), ref = "No territory"),
  habitat_binary = relevel(as.factor(habitat_binary), ref = "Open"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_binary = relevel(as.factor(trophic_binary), ref = "Secondary")
)

# # Center categorical predictors.
# model_data %<>% mutate(
#   terr_year_c = center_categorical(territory_year_round),
#   terr_bi_c = center_categorical(territory_binary),
#   habitat_bi_c = center_categorical(habitat_binary),
#   migration_bi_c = center_categorical(migration_binary),
#   trophic_level_c = center_categorical(trophic_binary)
# )
# 

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  body_mass_log_z = standardize(body_mass_log, two_sd = TRUE),
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE),
  range_z = standardize(range_log, two_sd = TRUE),
  temp_stab_z = standardize(log(bio4), two_sd = TRUE),
  rain_stab_z = standardize(bio15, two_sd = TRUE),
  npp_z = standardize(npp, two_sd = TRUE)
)



# Prepare response variables.
model_data$sexual_binary %<>% as.factor() %>% relevel(ref = "Weak")
model_data$sexual_binary %<>% as.numeric()
model_data$sexual_binary <- model_data$sexual_binary -1
model_data$sexual_score <- model_data$sexual_score + 1



library(brms)


# Brms formula options.

if(response == "sexual_binary"){
  
  
  
  brms_formula <- brmsformula(sexual_binary ~  territory_binary + trophic_binary + centroid_z + habitat_binary + 
                                migration_binary + 
                                body_mass_log_z + trophic_binary*centroid_z + temp_stab_z + 
                                rain_stab_z + npp_z + 
                                (1|gr(tree_tip, cov=A)), 
                              family = bernoulli(),
                              decomp = "QR")
} else {
  brms_formula <- brmsformula(sexual_score ~  territory_binary + trophic_binary + centroid_z + habitat_binary + 
                                migration_binary + 
                                body_mass_log_z + trophic_binary*centroid_z + temp_stab_z + 
                                rain_stab_z + npp_z + 
                                (1|gr(tree_tip, cov=A)), 
                              family = cumulative(),
                              decomp = "QR")
}

model_pathway <- paste0("Results/Models/brms/", response, "_", data_type, ".rds")


# Add un-informative priors.
normal_priors <- c(prior(normal(0,5),class="Intercept"),
                   prior(normal(0,5),class="b"),
                   prior(student_t(3,0,20), "sd"))

tic()
brms_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  prior = normal_priors,
  iter = 2500,
  warmup = 1000,
  chains = 4,
  cores = 32,
  init = 0,
  #control = list(adapt_delta=0.99, max_treedepth = 12), # no control for now to see how it goes.
  #file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  opencl = opencl(c(0, 0))
  #threads = threading(8),
)
toc()


summary(brms_model)

# # Plot the conditional effects.
# #plot(conditional_effects(brms_model), points=TRUE)
# 
# # Plot the model.
# plot(brms_model)
# 
# # Check the posterior distribution.
# pp_check(brms_model)
# 
# # Look at the pairs.
# #pairs(brms_model)
# 
# library(bayesplot)
# 
# # Do mcmc_areas.
# mcmc_areas(brms_model, regex_pars = "^b_[a-z]", prob = 0.89)

