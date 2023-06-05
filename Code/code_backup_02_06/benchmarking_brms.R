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

# Functions.
source("Code/functions.R")


################################################################################
                  #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Create types for hpc jobs.
tree_number <- 1

# Set the data types.
#data_type <- c("all", "high")
data_type <- c("all")

# Model type.
model_type <- 2

# We don't want to use cratio/sratio, because it assumes it can only move in 
# one direction. I.e. can't become less monogamous again. Likewise acat only does
# the probability between adjacent groups. So cumulative is the best shout.

# Family.    
family <- c("cumulative")

# Grain sizes.
grain_size <- c(100, 200, 400, 800, 1200, 1600)

# Threads.
#thread_number <- c(4, 8, 16)
thread_number <- c(4)

# Expand the grid.
all_combos <- expand.grid(tree_number, data_type, model_type, family, grain_size, thread_number)

# Model type.
tree_number <- all_combos[array_number, 1] %>% as.numeric()

# Data type.
data_type <- all_combos[array_number, 2] %>% as.character()

# Model type.
model_type <- all_combos[array_number, 3] %>% as.numeric()

# Family
family <- all_combos[array_number, 4] %>% as.character()

# Grain size.
grain_size <- all_combos[array_number, 5] %>% as.numeric() 

# Thread number.
thread_number <- all_combos[array_number, 6] %>% as.numeric()

# Put working directory as correct one.
setwd("..")
#setwd("Sexual_selection/")


###############################################################################
                       #### Read in the data #####


# Read in the tree.
model_tree <- read.tree("Data/Trees/model_trees.tre")[[tree_number]]

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$bird_tree_name)

if (data_type == "high"){
  model_data %<>% filter(sexual_certainty < 3)
}

# Drop tips on the tree.
model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Make a covariance matrix, and order data the same.
model_covar <- ape::vcv.phylo(model_tree)

# Reorder the matrix so it's random, to maximise parallel processing speed.
set.seed(1993)
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
  trophic_binary = relevel(as.factor(trophic_binary), ref = "Secondary")
)

# Center categorical predictors.
model_data %<>% mutate(
  terr_bi_c = center_categorical(territory_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary)
)

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  body_mass_log_z = standardize(body_mass_log, two_sd = TRUE),
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE),
  temp_seasonality_z = standardize(temp_log, two_sd = TRUE),
  npp_z = standardize(npp_sqrt, two_sd = TRUE),
  gen_z = standardize(gen_log, two_sd = TRUE),
)

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1

model_formula <- "sexual_score ~ gen_z + migration_binary + territory_binary + 
                  trophic_binary*temp_seasonality_z + 
                  trophic_binary*npp_z + trophic_binary*territory_binary + 
                  (1|gr(tree_tip, cov=A))"

# Family.
brms_formula <- brmsformula(model_formula, family = cumulative(), decomp = "QR")

# Create a unique model pathway for each model.
model_pathway <- paste0("Results/Models/Test/", family, "_", thread_number, "_", grain_size, "_", tree_number, ".rds")


tic()
brms_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  #prior = normal_priors,
  iter = 200,
  warmup = 50,
  chains = 2,
  #thin = 10,
  cores = thread_number*2,
  init = 0,
  #file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  threads = threading(thread_number, grainsize = grain_size, static = TRUE),
  seed = 1999
)
toc()

# Create and export a dataframe of times for threads.
model_threads <- data.frame(array_number = array_number, 
                            family = family,
                            tree_number = tree_number,
                            grain_size = grain_size,
                            thread_number = thread_number,
                            elapsed = mean(rstan::get_elapsed_time(brms_model$fit)))

# Append the table.
write.table(model_threads, "benchmarking.csv", sep = ",", append = TRUE, 
            row.names = FALSE, col.names=!file.exists("benchmarking.csv"))