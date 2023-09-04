###############################################################################
                   ##### HPC sex metric models #####
###############################################################################

# This script runs the analysis comparing sexual selection scores against 
# alternative measures of sexual selection.

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
model_type <- c("testes", "bateman", "oss")

# Centered or uncentered.
center <- c("centered", "uncentered")

# Set the data types.
data_type <- c("phylo", "linear")

# Expand the grid.
all_combos <- expand.grid(tree_number, center, model_type, data_type)

# Tree type.
tree_number <- all_combos[array_number, 1] %>% as.numeric()

# Center type.
center <- all_combos[array_number, 2] %>% as.character()

# Model type.
model_type <- all_combos[array_number, 3] %>% as.character()

# Data type.
data_type <- all_combos[array_number, 4] %>% as.character()


## Array job spec info ##

# Array num: 1 - 100 = 16 cores for 2 hours
# Array num: 101 - 300 = 16 cores for 30 minutes
# Array number 301 - 600 8 cores for 30 minutes

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
metric_data <- read_excel("Data/sexual_selection_dataset_23_08.xlsx", sheet = 4, na = "NA") %>% clean_names()
birdtree_data <-  read_excel("Data/sexual_selection_dataset_23_08.xlsx", sheet = 2, na = "NA") %>% clean_names()
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
  testes_z = standardize(residual_testes_mass, two_sd = TRUE)
)

bateman_data %<>% mutate(
  bateman_z = standardize(bateman_gradient, two_sd = TRUE)
)

oss_data %<>% mutate(
  oss_z = standardize(opportunity_for_sexual_selection, two_sd = TRUE),
  oss_log = log(opportunity_for_sexual_selection + 1),
  oss_log_z = standardize(oss_log, two_sd = TRUE)
)

# Add 1 to sexual selection score so it can be used with ordinal models.
testes_data$sexual_selection <- testes_data$sexual_selection + 1
bateman_data$sexual_selection <- bateman_data$sexual_selection + 1
oss_data$sexual_selection <- oss_data$sexual_selection + 1



###############################################################################
                   #### Run brms models #####

# Make the formula.

# For phylo models.
if (data_type == "phylo"){
  
  # For un-centered models.
if (center == "uncentered"){
teste_formula <- brmsformula("sexual_selection ~ residual_testes_mass + (1|gr(tree_tip, cov=A))", family = cumulative())
bateman_male_formula <- brmsformula("sexual_selection ~ bateman_gradient + (1|gr(tree_tip, cov=A))", family = cumulative())
oss_formula <- brmsformula("sexual_selection ~ oss_log + (1|gr(tree_tip, cov=A))", family = cumulative())

  # For centered models.
} else {
teste_formula <- brmsformula("sexual_selection ~ testes_z + (1|gr(tree_tip, cov=A))", family = cumulative())
bateman_male_formula <- brmsformula("sexual_selection ~ bateman_z + (1|gr(tree_tip, cov=A))", family = cumulative())
oss_formula <- brmsformula("sexual_selection ~ oss_log_z + (1|gr(tree_tip, cov=A))", family = cumulative())
}
  
  
# For non-phylo models.
} else {
  if (center == "uncentered"){
    teste_formula <- brmsformula("sexual_selection ~ residual_testes_mass", family = cumulative())
    bateman_male_formula <- brmsformula("sexual_selection ~ bateman_gradient", family = cumulative())
    oss_formula <- brmsformula("sexual_selection ~ oss_log", family = cumulative())
  } else {
    teste_formula <- brmsformula("sexual_selection ~ testes_z", family = cumulative())
    bateman_male_formula <- brmsformula("sexual_selection ~ bateman_z", family = cumulative())
    oss_formula <- brmsformula("sexual_selection ~ oss_log_z", family = cumulative())
  }
}

# Model pathway
model_pathway <- paste0("Results/Models/Sex_metrics/", model_type, "_", center, "_", data_type, "_", tree_number, ".rds") 

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
if (model_type == "oss"){
oss_model <- brm(
  oss_formula, data = oss_data, data2 = list(A=oss_covar),
  iter = 10000, warmup = 5000, chains = 2, thin = 20, cores = core_number, init = 0, 
  control = list(adapt_delta = 0.99), file = model_pathway, threads = thread_number,
  normalize = FALSE, backend = "cmdstanr")
}






# 
# 
# ###############################################################################
#                      #### Read in the data #####
# 
# # Read in some data.
# model_data <- read.csv("Data/both_testes_datasets_01_08.csv")
# 
# # Filter for high cert.
# model_data %<>% filter(data_certainty > 2)
# 
# model_data <- model_data[1:250,]
# 
# # Read in the tree.
# model_tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]
# 
# # Drop tips on the tree.
# model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))
# 
# # Make a covariance matrix, and order data the same.
# model_covar <- ape::vcv.phylo(model_tree)
# 
# # Reorder the matrix so it's random, to maximise parallel processing speed.
# mat_order <- sample(1:nrow(model_covar), size = nrow(model_covar), replace = FALSE)
# model_covar <- reorder_mat(model_covar, rownames(model_covar)[mat_order])
# row.names(model_data) <- model_data$tree_tip
# model_data <- model_data[row.names(model_covar),]
# 
# 
# ###############################################################################
#               #### Prepare predictor variables ######
# 
# # Scale continuous predictors to two SD.
# model_data %<>% mutate(
#   testes_z = standardize(residual_testes_mass, two_sd = TRUE))
# 
# # Prepare response variables.
# model_data$sexual_selection <- model_data$sexual_selection + 1
# 
# 
# ###############################################################################
#                     #### Set model formula ######
# 
# # model formula.
# full_formula <- "sexual_selection ~ testes_z + (1|gr(tree_tip, cov=A))"
# phylo_formula <- "sexual_selection ~ 1 + (1|gr(tree_tip, cov=A))"
# pred_formula <- "sexual_selection ~ testes_z"
# 
# # brms formula.
# full_formula %<>%  brmsformula(family = cumulative())
# phylo_formula %<>%  brmsformula(family = cumulative())
# pred_formula %<>%  brmsformula(family = cumulative())
# 
# # Simple models.
# model_pathway <- paste0("Results/Models/Testes/teste_model_1.rds") 
# 
# # Report time before starting.
# Sys.time()
# 
# # Run brms models.
# full_model <- brm(
#   full_formula, data = model_data, data2 = list(A=model_covar),
#   iter = 2000, warmup = 500, chains = 2, thin = 20,
#   cores = 32, init = 0, #file = model_pathway, 
#   normalize = FALSE,
#   backend = "cmdstanr", threads = threading(16))
# 
# # Run brms models.
# phy_model <- brm(
#   phylo_formula, data = model_data, data2 = list(A=model_covar),
#   iter = 1000, warmup = 500, chains = 2, thin = 20,
#   cores = 32, init = 0, #file = model_pathway, 
#   normalize = FALSE,
#   backend = "cmdstanr", threads = threading(16))
# 
# library(performance)
# 
# marginal_R2_MZ(full_model)
# 
# 
# 
# r2_bayes(full_model)
# ?r2_bayes
# R2.lik(full_model, phy_model)
# 
# bayes_R2(full_model)
# bayes_R2(full_model, resp = "testes_z")
# summary(full_model)
# bayes_R2(phy_model)
# 
# r2_loo(full_model)
# 
# y_pred <- fitted(full_model, scale = "linear", summary = FALSE, re_formula = NA)
# var_fit <- apply(y_pred, 1, var)
# 
# var_res <- pi^2 / 3 
# 
# R2_MZ <- var_fit / (var_fit + var_res)
# 
# mean(R2_MZ)
# 
#  # 0.89
#  # 0.15 %
# library(rr2)
# 
# 
# 
# Est.Error = sd(R2_MZ)
# 
#     "l-95% CI" = quantile(R2_MZ, 0.025),
#     "u-95% CI" = quantile(R2_MZ, 0.975),
#     row.names = "Bayes_R2_MZ", 
#     check.names = FALSE), 
#   digits = 2)
# 
# 
# # Bayesian McKelvey-Zavoina R2 ------------------------------------------------
# # Bayesian version of McKelvey and Zavoina's pseudo-R2 for binary and ordinal
# # brms models (McKelvey and Zavoina, 1975). See also Gelman et al. (2018). This
# # pseudo-R2 closely approximates the R2 that would have been obtained if a
# # linear model had have been run on observations of the continuous latent
# # variable underlying the discrete responses (Veall and Zimmermann, 1992; Hagle
# # and Mitchell, 1992; Veall and Zimmermann, 1994).
# Bayes_R2_MZ <- function(fit, ...) {
#   y_pred <- fitted(fit, scale = "linear", summary = FALSE, ...)
#   var_fit <- apply(y_pred, 1, var)
#   
#   if (fit$formula$family$family == "cumulative" ||
#       fit$formula$family$family == "bernoulli") {
#     if (fit$formula$family$link == "probit" || 
#         fit$formula$family$link == "probit_approx") {
#       var_res <- 1
#     }
#     else if (fit$formula$family$link == "logit") {
#       var_res <- pi^2 / 3 
#     }
#   } 
#   else {
#     sum_fit <- summary(fit)
#     sig_res <- sum_fit$spec_pars["sigma", "Estimate"]
#     var_res <- sig_res^2
#   } 
#   R2_MZ <- var_fit / (var_fit + var_res)
#   print(
#     data.frame(
#       Estimate = mean(R2_MZ), 
#       Est.Error = sd(R2_MZ), 
#       "l-95% CI" = quantile(R2_MZ, 0.025),
#       "u-95% CI" = quantile(R2_MZ, 0.975),
#       row.names = "Bayes_R2_MZ", 
#       check.names = FALSE), 
#     digits = 2)
# }
# 
# 
# 
# 
# 
# ###############################################################################
#                            #### Section 2 ####
# 
# ###############################################################################
# ###### R squared data ########
# 
# 
# # Load r-squared package.
# library(rr2)
# 
# # Lambda zero object for measuring effect of phylogeny.
# no_phylo <- data.frame(lambda = 0)
# 
# # Object to save results in. # Need to add trait and tree in as well.
# r_squared_data <- data.frame(total_R = NA, phy_R = NA, pred_R = NA)
# 
# # Get formula.
# formula <- trait_model$formula
# 
# # Get lambda.
# lambda_opt <- data.frame(lambda = trait_model$optpar)
# 
# # Get data.
# model_data <- cbind(trait_model$y, trait_model$X[,2:6])
# colnames(model_data)[1] <- model_trait
# #colnames(model_data)[2:3] <- c("sex_bin", "territory_bi") # doesn't matter if center categorical.
# model_data %<>% as.data.frame()
# 
# # Do a model with lambda set to zero.
# no_lambda_model <- phylolm(formula, data=model_data, phy = model_tree, model = "lambda", starting.value = no_phylo, 
#                            lower.bound = no_phylo$lambda, upper.bound = no_phylo$lambda)
# 
# # Do a phylo model with intercept only.
# no_predictor_model <- phylolm(paste(model_trait,"~ 1"), data=model_data, phy = model_tree, model = "lambda", starting.value = lambda_opt, 
#                               lower.bound = lambda_opt$lambda, upper.bound = lambda_opt$lambda)
# 
# # Calculate full R squared. Includes predictors, phylogeny and intercept/mean.
# r_squared_data[1,1] <- R2.lik(trait_model)
# 
# # Get phylo model R squared.
# r_squared_data[1,2] <- R2.lik(trait_model, no_lambda_model)
# 
# # Get predictor R squared.
# r_squared_data[1,3] <- R2.lik(trait_model, no_predictor_model)
# 
# # Add in array data.
# r_squared_data$trait <- trait
# r_squared_data$tree_number <- tree_number
# r_squared_data$data_type <- data_type
# 
# # Inspect results.
# r_squared_data
# 
# # Make a pathway name.
# file_pathway <- paste0("../Results/Phylolm/", data_type, "_", trait, "_", tree_number, "_absolute_dimorphism_r_squared_data.csv")
# 
# # Export the data.
# write.csv(r_squared_data, file_pathway, row.names = FALSE)
# 
# 
# 
# # # Export the data.
# # write.table(r_squared_data, "../Results/absolute_dimorphism_r_squared_data.csv", 
# #             sep = ",", append = TRUE, row.names = FALSE, 
# #             col.names=!file.exists("../Results/absolute_dimorphism_r_squared_data.csv"))
# 
# 
# ################################################################################
# #### End ####
# 
# 
# 
# ###############################################################################
#                            #### Section 3 ####
# 
# 
# ###############################################################################
#                            #### Section 4 ####
# 
# ###############################################################################
#                            #### Section 5 ####
# 
# ###############################################################################
#                            #### Section 6 ####
# 
# 
# ###############################################################################
#                            #### Section 7 ####
# 
# 
# ###############################################################################
#                            #### Section 8 ####
# 
# # Look Rob, you've had your fun with the sectioning. 
# # They'll be no more sectioning today.
# 
# 
# ###############################################################################
#                              #### END ####
# ###############################################################################
# 
# 
# ###############################################################################
#               #### All the stuff I'm afraid to delete ####
# 
# 
