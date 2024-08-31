###############################################################################
                   ##### HPC sex metric models #####
###############################################################################

# This script runs the analysis comparing sexual selection scores against 
# alternative measures of sexual selection.


# Packages to load.
library(magrittr)
library(caper)
library(dplyr)
library(effectsize)
library(car)
library(ggplot2)
library(ggpubr)
library(phytools)
library(brms)
library(graph4lg)
library(janitor)
library(readxl)

# Clear the workspace.
rm(list=ls())

# Read in the functions. 
source("Code/functions.R")


################################################################################
                  #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Create types for hpc jobs.
tree_number <- 1:50

# Model (In order of size and therefore speed)
model_type <- c("testes", "bateman", "oss", "sensoss")

# Centered or uncentered.
center <- c("centered", "uncentered")

# Expand the grid.
all_combos <- expand.grid(tree_number, center, model_type)

# Extract array variables.
tree_number <- all_combos[array_number, 1] %>% as.numeric()
center <- all_combos[array_number, 2] %>% as.character()
model_type <- all_combos[array_number, 3] %>% as.character()


## Array job spec info ##

# Array num: 1 - 100 = 16 cores for 2 hours
# Array num: 101 - 300 = 16 cores for 30 minutes
# Array number 301 - 600 8 cores for 30 minutes

# Set cores.
if (array_number %in% 1:300){
  core_number <- 16
  thread_number <- 8
}
if (array_number %in% 301:600){
  core_number <- 8
  thread_number <- 4
}

###############################################################################
                        #### Read in data ####

# Read in the sex metric data.
pathway <- "Data/supplementary_dataset_1.xlsx"
metric_data <- read_excel(pathway, sheet = 4, na = "NA") %>% clean_names()
birdtree_data <-  read_excel(pathway, sheet = 2, na = "NA") %>% clean_names()
model_data <- left_join(metric_data, birdtree_data)
model_data$tree_tip <- gsub(" ", "_", model_data$scientific_name_bird_tree)

# Select testes mass.
testes_data <- model_data %>% 
  dplyr::select(scientific_name_bird_tree, tree_tip, sexual_selection, data_certainty, 
                sex_role_reversal, residual_testes_mass) %>% na.omit()
testes_data %<>% filter(data_certainty > 2)

# Select Bateman gradient metrics.
bateman_data <- model_data %>% 
  dplyr::select(scientific_name_bird_tree, tree_tip, sexual_selection, data_certainty, 
                bateman_gradient, bateman_individuals) %>% na.omit()

# Select OSS metrics.
oss_data <- model_data %>% 
  dplyr::select(scientific_name_bird_tree, tree_tip, sexual_selection, data_certainty, 
                sex_role_reversal, opportunity_for_sexual_selection, 
                oss_individuals, oss_estimated)  %>% na.omit()

# Sensitivity analysis.
oss_data %<>% filter(opportunity_for_sexual_selection != 0)

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[tree_number]]

# Drop tips on the tree.
testes_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, testes_data$tree_tip))
bateman_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, bateman_data$tree_tip))
oss_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, oss_data$tree_tip))

# Make a covariance matrix, and order data the same.
testes_covar <- ape::vcv.phylo(testes_tree)
bateman_covar <- ape::vcv.phylo(bateman_tree)
oss_covar <- ape::vcv.phylo(oss_tree)


###############################################################################
                      #### Predictors ####


# Scale continuous predictors to two SD.
testes_data %<>% mutate(
  testes_z = standardize(residual_testes_mass, two_sd = TRUE))

bateman_data %<>% mutate(
  bateman_z = standardize(bateman_gradient, two_sd = TRUE))

oss_data %<>% mutate(
  oss_z = standardize(opportunity_for_sexual_selection, two_sd = TRUE),
  oss_log = log(opportunity_for_sexual_selection + 1),
  oss_log_z = standardize(oss_log, two_sd = TRUE))

# Add 1 to sexual selection score so it can be used with ordinal models.
testes_data$sexual_selection <- testes_data$sexual_selection + 1
bateman_data$sexual_selection <- bateman_data$sexual_selection + 1
oss_data$sexual_selection <- oss_data$sexual_selection + 1


###############################################################################
                     #### Run brms models #####


# Formula constants.
resp <- "sexual_selection ~ "
random <- " + (1|gr(tree_tip, cov=A))"

# Make the formula.
if (center == "uncentered"){
teste_formula <- brmsformula(paste0(resp, "residual_testes_mass", random), family = cumulative())
bateman_male_formula <- brmsformula(paste0(resp, "bateman_gradient", random), family = cumulative())
oss_formula <- brmsformula(paste0(resp, "oss_log", random), family = cumulative())

} else {
teste_formula <- brmsformula(paste0(resp, "testes_z", random), family = cumulative())
bateman_male_formula <- brmsformula(paste0(resp, "bateman_z", random), family = cumulative())
oss_formula <- brmsformula(paste0(resp, "oss_log_z", random), family = cumulative())
}
  
# Model pathway
model_pathway <- paste0("Results/Models/Sex_metrics/", model_type, "_", center, 
                        "_", data_type, "_", tree_number, ".rds") 

# Run models.
if (model_type == "testes"){
  teste_model <- brm(
    teste_formula, data = testes_data, data2 = list(A=testes_covar),
    iter = 10000, warmup = 5000, chains = 2, thin = 20, cores = core_number, init = 0, 
    control = list(adapt_delta = 0.99), file = model_pathway, threads = thread_number,
    normalize = FALSE, backend = "cmdstanr")
}
if (model_type == "bateman"){
bateman_male_model <- brm(
  bateman_male_formula, data = bateman_data, data2 = list(A=bateman_covar),
  iter = 10000, warmup = 5000, chains = 2, thin = 20, cores = core_number, init = 0, 
  control = list(adapt_delta = 0.99), file = model_pathway, threads = thread_number,
  normalize = FALSE, backend = "cmdstanr")
}
if (model_type %in% c("oss", "sensoss")){
oss_model <- brm(
  oss_formula, data = oss_data, data2 = list(A=oss_covar),
  iter = 10000, warmup = 5000, chains = 2, thin = 20, cores = core_number, init = 0, 
  control = list(adapt_delta = 0.99), file = model_pathway, threads = thread_number,
  normalize = FALSE, backend = "cmdstanr")
}



###############################################################################
                              #### END ####
###############################################################################