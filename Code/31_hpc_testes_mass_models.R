###############################################################################
                   ##### HPC teste mass models #####
###############################################################################

# This script does something new. Pretty sick eh.


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

# Read in the functions. 
source("Code/functions.R")

###############################################################################
                     #### Read in the data #####

# Read in some data.
model_data <- read.csv("Data/both_testes_datasets_01_08.csv")

# Filter for high cert.
model_data %<>% filter(data_certainty > 2)

model_data <- model_data[1:250,]

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]

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
  testes_z = standardize(residual_testes_mass, two_sd = TRUE))

# Prepare response variables.
model_data$sexual_selection <- model_data$sexual_selection + 1


###############################################################################
                    #### Set model formula ######

# model formula.
full_formula <- "sexual_selection ~ testes_z + (1|gr(tree_tip, cov=A))"
phylo_formula <- "sexual_selection ~ 1 + (1|gr(tree_tip, cov=A))"
pred_formula <- "sexual_selection ~ testes_z"

# brms formula.
full_formula %<>%  brmsformula(family = cumulative())
phylo_formula %<>%  brmsformula(family = cumulative())
pred_formula %<>%  brmsformula(family = cumulative())

# Simple models.
model_pathway <- paste0("Results/Models/Testes/teste_model_1.rds") 

# Report time before starting.
Sys.time()

# Run brms models.
full_model <- brm(
  full_formula, data = model_data, data2 = list(A=model_covar),
  iter = 2000, warmup = 500, chains = 2, thin = 20,
  cores = 32, init = 0, #file = model_pathway, 
  normalize = FALSE,
  backend = "cmdstanr", threads = threading(16))

# Run brms models.
phy_model <- brm(
  phylo_formula, data = model_data, data2 = list(A=model_covar),
  iter = 1000, warmup = 500, chains = 2, thin = 20,
  cores = 32, init = 0, #file = model_pathway, 
  normalize = FALSE,
  backend = "cmdstanr", threads = threading(16))

library(performance)

marginal_R2_MZ(full_model)



r2_bayes(full_model)
?r2_bayes
R2.lik(full_model, phy_model)

bayes_R2(full_model)
bayes_R2(full_model, resp = "testes_z")
summary(full_model)
bayes_R2(phy_model)

r2_loo(full_model)

y_pred <- fitted(full_model, scale = "linear", summary = FALSE, re_formula = NA)
var_fit <- apply(y_pred, 1, var)

var_res <- pi^2 / 3 

R2_MZ <- var_fit / (var_fit + var_res)

mean(R2_MZ)

 # 0.89
 # 0.15 %
library(rr2)



Est.Error = sd(R2_MZ)

    "l-95% CI" = quantile(R2_MZ, 0.025),
    "u-95% CI" = quantile(R2_MZ, 0.975),
    row.names = "Bayes_R2_MZ", 
    check.names = FALSE), 
  digits = 2)


# Bayesian McKelvey-Zavoina R2 ------------------------------------------------
# Bayesian version of McKelvey and Zavoina's pseudo-R2 for binary and ordinal
# brms models (McKelvey and Zavoina, 1975). See also Gelman et al. (2018). This
# pseudo-R2 closely approximates the R2 that would have been obtained if a
# linear model had have been run on observations of the continuous latent
# variable underlying the discrete responses (Veall and Zimmermann, 1992; Hagle
# and Mitchell, 1992; Veall and Zimmermann, 1994).
Bayes_R2_MZ <- function(fit, ...) {
  y_pred <- fitted(fit, scale = "linear", summary = FALSE, ...)
  var_fit <- apply(y_pred, 1, var)
  
  if (fit$formula$family$family == "cumulative" ||
      fit$formula$family$family == "bernoulli") {
    if (fit$formula$family$link == "probit" || 
        fit$formula$family$link == "probit_approx") {
      var_res <- 1
    }
    else if (fit$formula$family$link == "logit") {
      var_res <- pi^2 / 3 
    }
  } 
  else {
    sum_fit <- summary(fit)
    sig_res <- sum_fit$spec_pars["sigma", "Estimate"]
    var_res <- sig_res^2
  } 
  R2_MZ <- var_fit / (var_fit + var_res)
  print(
    data.frame(
      Estimate = mean(R2_MZ), 
      Est.Error = sd(R2_MZ), 
      "l-95% CI" = quantile(R2_MZ, 0.025),
      "u-95% CI" = quantile(R2_MZ, 0.975),
      row.names = "Bayes_R2_MZ", 
      check.names = FALSE), 
    digits = 2)
}





###############################################################################
                           #### Section 2 ####

###############################################################################
###### R squared data ########


# Load r-squared package.
library(rr2)

# Lambda zero object for measuring effect of phylogeny.
no_phylo <- data.frame(lambda = 0)

# Object to save results in. # Need to add trait and tree in as well.
r_squared_data <- data.frame(total_R = NA, phy_R = NA, pred_R = NA)

# Get formula.
formula <- trait_model$formula

# Get lambda.
lambda_opt <- data.frame(lambda = trait_model$optpar)

# Get data.
model_data <- cbind(trait_model$y, trait_model$X[,2:6])
colnames(model_data)[1] <- model_trait
#colnames(model_data)[2:3] <- c("sex_bin", "territory_bi") # doesn't matter if center categorical.
model_data %<>% as.data.frame()

# Do a model with lambda set to zero.
no_lambda_model <- phylolm(formula, data=model_data, phy = model_tree, model = "lambda", starting.value = no_phylo, 
                           lower.bound = no_phylo$lambda, upper.bound = no_phylo$lambda)

# Do a phylo model with intercept only.
no_predictor_model <- phylolm(paste(model_trait,"~ 1"), data=model_data, phy = model_tree, model = "lambda", starting.value = lambda_opt, 
                              lower.bound = lambda_opt$lambda, upper.bound = lambda_opt$lambda)

# Calculate full R squared. Includes predictors, phylogeny and intercept/mean.
r_squared_data[1,1] <- R2.lik(trait_model)

# Get phylo model R squared.
r_squared_data[1,2] <- R2.lik(trait_model, no_lambda_model)

# Get predictor R squared.
r_squared_data[1,3] <- R2.lik(trait_model, no_predictor_model)

# Add in array data.
r_squared_data$trait <- trait
r_squared_data$tree_number <- tree_number
r_squared_data$data_type <- data_type

# Inspect results.
r_squared_data

# Make a pathway name.
file_pathway <- paste0("../Results/Phylolm/", data_type, "_", trait, "_", tree_number, "_absolute_dimorphism_r_squared_data.csv")

# Export the data.
write.csv(r_squared_data, file_pathway, row.names = FALSE)



# # Export the data.
# write.table(r_squared_data, "../Results/absolute_dimorphism_r_squared_data.csv", 
#             sep = ",", append = TRUE, row.names = FALSE, 
#             col.names=!file.exists("../Results/absolute_dimorphism_r_squared_data.csv"))


################################################################################
#### End ####



###############################################################################
                           #### Section 3 ####


###############################################################################
                           #### Section 4 ####

###############################################################################
                           #### Section 5 ####

###############################################################################
                           #### Section 6 ####


###############################################################################
                           #### Section 7 ####


###############################################################################
                           #### Section 8 ####

# Look Rob, you've had your fun with the sectioning. 
# They'll be no more sectioning today.


###############################################################################
                             #### END ####
###############################################################################


###############################################################################
              #### All the stuff I'm afraid to delete ####


