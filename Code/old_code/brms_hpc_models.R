###############################################################################
              # Testing culumative brms models on cluster  #
###############################################################################
#dyn.load("C:/Windows/System32/OpenCL.dll", FALSE, TRUE)

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

# Clear the workspace.
rm(list=ls())
#setwd("Manuscript/")
# Functions.
source("../Code/functions.R")


################################################################################
                      #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number <- 1

# Create types for hpc jobs.
model_type <- c("standard", "centered", "no_env", "centered_no_env")

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(model_type, data_type)

# Model type.
model_type <- all_combos[array_number, 1] %>% as.character()

# Data type.
data_type <- all_combos[array_number, 2] %>% as.character()

# Put working directory as correct one.
setwd("..")
#setwd("Sexual_selection/")


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

# Center categorical predictors.
model_data %<>% mutate(
  terr_bi_c = center_categorical(territory_binary),
  habitat_bi_c = center_categorical(habitat_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary)
)

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
model_data$sexual_score <- model_data$sexual_score + 1

# #test sample size.
x <- 500
test_data <- model_data[sample(1:nrow(model_data), x),]

# Drop tips on the tree.
test_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, test_data$tree_tip))

# Make a covariance matrix.
test_covar <- ape::vcv.phylo(test_tree)

# Make sure the covariance matrix and data are in the same order.
row.names(test_data) <- test_data$tree_tip
test_data <- test_data[row.names(test_covar),]



# Create the formula based on the options.
if(model_type == "standard"){
  model_formula <- formula(sexual_score ~  territory_binary + trophic_binary + centroid_z + habitat_binary + 
                            migration_binary + 
                            body_mass_log_z + trophic_binary*centroid_z + temp_stab_z + 
                            rain_stab_z + npp_z + 
                            (1|gr(tree_tip, cov=A)))
}

if(model_type == "centered"){
  model_formula <- formula(sexual_score ~  terr_bi_c + trophic_level_c + centroid_z + habitat_bi_c + 
                            migration_bi_c + 
                            body_mass_log_z + trophic_level_c*centroid_z + temp_stab_z + 
                            rain_stab_z + npp_z + 
                            (1|gr(tree_tip, cov=A)))
}
if(model_type == "no_env"){
  model_formula <- formula(sexual_score ~  territory_binary + trophic_binary + centroid_z + 
                            habitat_binary + migration_binary + 
                            body_mass_log_z + trophic_binary*centroid_z +  
                            (1|gr(tree_tip, cov=A)))
}
if(model_type == "centered_no_env"){
  model_formula <- formula(sexual_score ~  terr_bi_c + trophic_level_c + centroid_z + 
                             habitat_bi_c + migration_bi_c + 
                             body_mass_log_z + trophic_level_c*centroid_z +  
                             (1|gr(tree_tip, cov=A)))
}


# Create the brms formula.
brms_formula <- brmsformula(model_formula, family = cumulative(), decomp = "QR")

# Create a unique model pathway for each model.
model_pathway <- paste0("Results/Models/", model_type, "_", data_type, ".rds")

# Add un-informative priors.
normal_priors <- c(prior(normal(0,5), class="Intercept"),
                   prior(normal(0,5), class="b"),
                   prior(student_t(3,0,20), "sd"))

tic()
brms_model2 <- brm(
  brms_formula,
  data = test_data,
  data2 = list(A=test_covar),
  prior = normal_priors,
  iter = 1000,
  warmup = 500,
  chains = 2,
  thin = 10,
  cores = 2,
  init = 0,
  #file = model_pathway, 
  normalize = FALSE,
  backend = "cmdstanr"
  #threads = threading(4),
  #algorithm = "fixed_param"
)
toc()

# Summary
summary(brms_model)
tic()
brms_model2 <- brm(
  body_mass_log_z ~ centroid_z,
  data = test_data,
  iter = 1000,
  warmup = 500,
  chains = 2,
  thin = 10,
  cores = 2,
  backend = "rstan",
  algorithm = "meanfield"
)
# 
toc()
library(rstan)
rstan::stan_version()
install.packages("StanHeaders")
install.packages("rstan")
library(rstantools)

library(cmdstanr)
cmdstanr::check_cmdstan_toolchain()
cmdstanr::check_cmdstan_toolchain(fix = TRUE) 
remotes::install_git("https://github.com/hsbadr/rstan", subdir = "StanHeaders", ref= "develop")
