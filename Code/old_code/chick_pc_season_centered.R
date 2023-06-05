###############################################################################
                     # Testing brms models on cluster  #
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

#setwd("Manuscript/")
#setwd("Sexual_selection/Manuscript/")
# Functions.
source("../Code/functions.R")


################################################################################
                  #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Create types for hpc jobs.
tree_number <- 1:3

# Set the data types.
#data_type <- c("all", "high")
data_type <- c("all")

# Expand the grid.
all_combos <- expand.grid(tree_number, data_type)

# Model type.
tree_number <- all_combos[array_number, 1] %>% as.numeric()

# Data type.
data_type <- all_combos[array_number, 2] %>% as.character()

# Put working directory as correct one.
setwd("..")


###############################################################################
                       #### Read in the data #####


# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[tree_number]]

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)

# Remove missing data for certain models.
#model_data %<>% filter(!is.na(devo_mode))
#model_data %<>% filter(chick_pc_ascores == "RealData")

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
  wang_edited = relevel(as.factor(wang_edited), ref = "altricial")
)

# Center categorical predictors.
model_data %<>% mutate(
  terr_bi_c = center_categorical(territory_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary),
  devo_bi_c = center_categorical(devo_binary),
  devo_ed_c = center_categorical(wang_edited),
  devo_mode_c = center_categorical(devo_mode)
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
  chick_z = standardize(chick_pc1, two_sd = TRUE)
)

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1

# Do the model formula.
# model_formula <- "sexual_score ~ gen_z + chick_z + migration_binary + territory_binary +
#                   trophic_binary*temp_seasonality_z +
#                   trophic_binary*npp_z + trophic_binary*territory_binary +
#                   chick_z*temp_seasonality_z +
#                   chick_z*npp_z + chick_z*territory_binary +
#                   (1|gr(tree_tip, cov=A))"

# model_formula <- "sexual_score ~ gen_z + chick_z + migration_binary + territory_binary +
#                   trophic_binary*temp_seasonality_z + trophic_binary*territory_binary +
#                   chick_z*temp_seasonality_z + chick_z*territory_binary +
#                   (1|gr(tree_tip, cov=A))"

model_formula <- "sexual_score ~ gen_z + chick_z + migration_bi_c + terr_bi_c +
                  trophic_level_c*temp_seasonality_z + trophic_level_c*terr_bi_c +
                  (1|gr(tree_tip, cov=A))"


# model_formula <- "sexual_score ~ gen_z + chick_z + migration_binary + territory_binary +
#                   trophic_binary*temp_seasonality_z +
#                   trophic_binary*npp_z + trophic_binary*territory_binary +
#                   (1|gr(tree_tip, cov=A))"

# model_formula <- "sexual_score ~ gen_z + chick_z + migration_bi_c + terr_bi_c +
#                   trophic_level_c*temp_seasonality_z +
#                   trophic_level_c*npp_z + trophic_level_c*terr_bi_c +
#                   (1|gr(tree_tip, cov=A))"


# brms formula.
brms_formula <- brmsformula(model_formula, family = cumulative(), decomp = "QR")

# Create a unique model pathway for each model.

# Changing ref group.
#model_pathway <- paste0("Results/Models/chick_pc", data_type, "_", tree_number, ".rds")                         ~~~~~Done~~~~~~
#model_pathway <- paste0("Results/Models/chick_pc_centered", data_type, "_", tree_number, ".rds")                 ~~~~~Done~~~~~~
#model_pathway <- paste0("Results/Models/chick_pc_inter", data_type, "_", tree_number, ".rds")                ~~~~~Done~~~~~~
#model_pathway <- paste0("Results/Models/chick_pc_season", data_type, "_", tree_number, ".rds")             ~~~~~Already Done as model 49~~~~~~
#model_pathway <- paste0("Results/Models/chick_pc_season_inter_centered", data_type, "_", tree_number, ".rds")      ~~~~~Done~~~~~~
model_pathway <- paste0("Results/Models/chick_pc_season_centered", data_type, "_", tree_number, ".rds") 
#model_pathway <- paste0("Results/Models/chick_pc_everything", data_type, "_", tree_number, ".rds") 
#model_pathway <- paste0("Results/Models/chick_pc_everything_centered", data_type, "_", tree_number, ".rds") 


# # Add un-informative priors.
normal_priors <- c(prior(normal(0,5), class="Intercept"),
                   prior(normal(0,5), class="b"),
                   prior(student_t(3,0,20), "sd"))

# Run brms models.
brms_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  prior = normal_priors,
  iter = 10000,
  warmup = 5000,
  chains = 2,
  thin = 10,
  cores = 32,
  init = 0,
  file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  threads = threading(16),
)




# #test sample size.
# x <- 2000
# test_data <- model_data[sample(1:nrow(model_data), x),]
# 
# # Drop tips on the tree.
# test_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, test_data$tree_tip))
# 
# # Make a covariance matrix.
# test_covar <- ape::vcv.phylo(test_tree)
# 
# # Make sure the covariance matrix and data are in the same order.
# row.names(test_data) <- test_data$tree_tip
# test_data <- test_data[row.names(test_covar),]
# 
# tic()
# brms_model <- brm(
#   brms_formula,
#   data = test_data,
#   data2 = list(A=test_covar),
#   prior = normal_priors,
#   iter = 1000,
#   warmup = 500,
#   chains = 2,
#   #thin = 0,
#   cores = 8,
#   init = 0,
#   #file = model_pathway,
#   normalize = FALSE,
#   backend = "cmdstanr",
#   threads = threading(4),
# )
# toc()
# 
# summary(brms_model)
# pp_check(brms_model)
# 
# library(bayesplot)
# # Plot model results with 95% confidence intervals.
# mcmc_areas(brms_model, regex_pars = "^b_[a-z]", prob = 0.95)




