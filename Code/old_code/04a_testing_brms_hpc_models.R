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

setwd("Manuscript/")
#setwd("Sexual_selection/Manuscript/")
# Functions.
source("../Code/functions.R")


################################################################################
                  #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number <- 1

# Create types for hpc jobs.
tree_number <- 1:10

# Set the data types.
data_type <- c("all", "high")

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

if (data_type == "high"){
  model_data %<>% filter(sexual_certainty < 3)
}

# For testing devo mode.
#model_data %<>% tidyr::drop_na(devo_mode)

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
  devo_binary = relevel(as.factor(devo_mode_wang), ref = "altricial")
)

# Center categorical predictors.
model_data %<>% mutate(
  terr_bi_c = center_categorical(territory_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary),
  devo_wang_c = center_categorical(devo_binary)
)

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  body_mass_log_z = standardize(body_mass_log, two_sd = TRUE),
  centroid_z = standardize(centroid_sqrt, two_sd = TRUE),
  temp_seasonality_z = standardize(temp_log, two_sd = TRUE),
  npp_z = standardize(npp_sqrt, two_sd = TRUE),
  gen_z = standardize(gen_log, two_sd = TRUE),
  dens_z = standardize(dens_sqrt, two_sd = TRUE),
  fed_z = standardize(fed_sqrt, two_sd = TRUE)
)

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1

# Do the model formula. (removed density from models)
model_formula <- "sexual_score ~ (body_mass_log_z + migration_binary + territory_binary +
                  trophic_binary*temp_seasonality_z + 
                  trophic_binary*npp_z + trophic_binary*territory_binary)*devo_binary + 
                  (1|gr(tree_tip, cov=A))"

# brms formula.
brms_formula <- brmsformula(model_formula, family = cumulative(), decomp = "QR")

# Create a unique model pathway for each model.
model_pathway <- paste0("Results/Models/dens", data_type, "_", tree_number, ".rds")

# # Add un-informative priors.
normal_priors <- c(prior(normal(0,5), class="Intercept"),
                   prior(normal(0,5), class="b"),
                   prior(student_t(3,0,20), "sd"))

# Run brms models.
tic()
brms_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  prior = normal_priors,
  iter = 1000,
  warmup = 500,
  chains = 2,
  #thin = 10,
  cores = 4,
  init = 0,
  #file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  threads = threading(2),
)
toc()

summary(brms_model)
 
model_data$devo_mode %>% as.factor() %>% as.numeric() %>% cor(y = model_data$gen_log)

cor(model_data$fed_sqrt, y = model_data$gen_log)

library(GGally)

test <- model_data %>% dplyr::select(body_mass_log_z, gen_z, devo_binary) 
test$devo_binary %<>% as.factor() %>% as.numeric()

test %>% ggcorr(label = TRUE, method =  c("everything","spearman"))




#test sample size.
x <- 500
test_data <- model_data[sample(1:nrow(model_data), x),]

# Drop tips on the tree.
test_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, test_data$tree_tip))

# Make a covariance matrix.
test_covar <- ape::vcv.phylo(test_tree)

# Make sure the covariance matrix and data are in the same order.
row.names(test_data) <- test_data$tree_tip
test_data <- test_data[row.names(test_covar),]

tic()
brms_model <- brm(
  brms_formula,
  data = test_data,
  data2 = list(A=test_covar),
  prior = normal_priors,
  iter = 1000,
  warmup = 500,
  chains = 2,
  #thin = 0,
  cores = 8,
  init = 0,
  #file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  threads = threading(4),
)
toc()

summary(brms_model)
pp_check(brms_model)

library(bayesplot)
# Plot model results with 95% confidence intervals.
mcmc_areas(brms_model, regex_pars = "^b_[a-z]", prob = 0.95)


################################################################################
             #### testing models using phylom ######

library(phylolm)

# Using wang and kimball developmental mode.
model_formula <- "sexual_score ~ devo_wang_c + gen_z + migration_binary + territory_binary +
                  trophic_binary*temp_seasonality_z + 
                  trophic_binary*npp_z + trophic_binary"
test_model <- phylolm(formula = model_formula, data = model_data, phy = model_tree)
summary(test_model)

model_formula <- "sexual_score ~ wang_edited + migration_binary + territory_binary +
                  trophic_binary*temp_seasonality_z + 
                  trophic_binary*npp_z + trophic_binary"
test_model <- phylolm(formula = model_formula, data = model_data, phy = model_tree)
summary(test_model)


# Using greisser version.
model_formula <- "sexual_score ~ devo_mode + gen_z + migration_binary + territory_binary +
                  trophic_binary*temp_seasonality_z + 
                  trophic_binary*npp_z + trophic_binary"
test_model <- phylolm(formula = model_formula, data = model_data, phy = model_tree)
summary(test_model)

# Using cooney version.
model_formula <- "sexual_score ~ cooney_dev + gen_z + migration_binary + territory_binary +
                  trophic_binary*temp_seasonality_z + 
                  trophic_binary*npp_z + trophic_binary"
test_model <- phylolm(formula = model_formula, data = model_data, phy = model_tree)
summary(test_model)

# Using time fed.
model_formula <- "sexual_score ~ fed_sqrt + gen_z + migration_binary + territory_binary +
                  trophic_binary*temp_seasonality_z + 
                  trophic_binary*npp_z + trophic_binary"
test_model <- phylolm(formula = model_formula, data = model_data, phy = model_tree)
summary(test_model)


model_formula <- "sexual_score ~ chick_pc1 + migration_binary + territory_binary +
                  trophic_binary*temp_seasonality_z + 
                  trophic_binary*npp_z + trophic_binary"
test_model <- phylolm(formula = model_formula, data = model_data, phy = model_tree)
summary(test_model)



model_formula <- "sexual_score ~ wang_edited + migration_binary + wang_edited*territory_binary +
                  trophic_binary*temp_seasonality_z +
                  trophic_binary*npp_z + trophic_binary*territory_binary + wang_edited*temp_seasonality_z + wang_edited*npp_z"
test_model <- phylolm(formula = model_formula, data = model_data, phy = model_tree)
summary(test_model)


################################################################################
                    #### testing models using sensiphy #####


# Extract relevant devo data.
devo_data <- model_data %>% dplyr::select(tree_tip, order, family, sexual_score, devo_mode_wang, devo_mode, devo_mode_full, cooney_dev, fed_sqrt)


#test sample size.
x <- 500
test_data <- model_data[sample(1:nrow(model_data), x),]
test_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, test_data$tree_tip))
test_covar <- ape::vcv.phylo(test_tree)
row.names(test_data) <- test_data$tree_tip
test_data <- test_data[row.names(test_covar),]


base_model <- phylolm(sexual_score ~ devo_wang_c, data = test_data, phy = test_tree)
summary(base_model)


library(sensiPhy)

# Using clades doesn't work apparently.
influ_clade <- clade_phylm(sexual_score ~ devo_wang_c, phy = test_tree, 
                     data = test_data, track=TRUE, clade.col = test_data$family)
summary(influ_clade)
sensi_plot(influ, "Emberizidae")


# Try looking at most influential species.
influ_2 <- influ_phylm(sexual_score ~ devo_wang_c, phy = test_tree, 
              data = test_data, track=TRUE, cutoff = 1)
summary(influ_2)
sensi_plot(influ_2, 2)

# Extract intercept and estimate species and check.
influ_estimate <- summary(influ_2)$Estimate
colnames(influ_estimate)[1] <- "tree_tip"

influ_intercept <- summary(influ_2)$Intercept
colnames(influ_intercept)[1] <- "tree_tip"

left_join(influ_estimate, devo_data) %>% filter(DIFestimate > 0)
left_join(influ_intercept, devo_data) %>% filter(DIFintercept < 0)


## Try changing these developmental modes as they're tricky.

model_data %>% filter(family == "Ardeidae")

changed_data <- model_data
changed_data$devo_mode_wang[changed_data$family == "Alcidae"] <- "precocial"
changed_data$devo_mode_wang[changed_data$family == "Dromadidae"] <- "altricial"
changed_data$devo_mode_wang[changed_data$family == "Ardeidae"] <- "altricial"

model_data %<>% filter(chick_pc_ascores == "RealData")

all_model <- phylolm(sexual_score ~ chick_pc1, data = model_data, phy = model_tree)
summary(all_model)
change_model <- phylolm(sexual_score ~ devo_mode_wang, data = changed_data, phy = model_tree)
summary(change_model)

## Try with a glm. ##

model_data$sex_binary <- model_data$sexual_score
model_data$sex_binary[model_data$sex_binary < 4] <- 0
model_data$sex_binary[model_data$sex_binary > 0] <- 1

base_glm <- phyloglm(sex_binary ~ devo_wang_c, data = test_data, phy = test_tree)
summary(base_glm)

influ_3 <- influ_phyglm(sex_binary ~ devo_wang_c, phy = test_tree, 
                        data = test_data, track=TRUE)
summary(influ_3)
influ_glm_estimate <- summary(influ_3)$Intercept
colnames(influ_glm_estimate)[1] <- "tree_tip"
left_join(influ_glm_estimate, devo_data)




## Try looking at time fed. ##

influ_4 <- influ_phylm(sexual_score ~ fed_sqrt, phy = test_tree, 
                       data = test_data, track=TRUE, cutoff = 1)

summary(influ_4)

not_rhea <- !test_data$tree_tip %in% "Rhea_americana"

test_data<- test_data[not_rhea,]

base_model <- phylolm(sexual_score ~ fed_sqrt, data = test_data, phy = test_tree)
summary(base_model)

influ_4estimate <- summary(influ_4)$Estimate
colnames(influ_4estimate)[1] <- "tree_tip"

influ_4intercept <- summary(influ_4)$Intercept
colnames(influ_4intercept)[1] <- "tree_tip"

left_join(influ_4estimate, devo_data) %>% filter(DIFestimate > 0)
left_join(influ_4intercept, devo_data) %>% filter(DIFintercept < 0)


