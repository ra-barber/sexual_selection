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

setwd("Manuscript/")

# Functions.
source("../Code/functions.R")


################################################################################
                  #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Create types for hpc jobs.
model_type <- 54:69

# Set the data types.
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(model_type, data_type)

# Model type.
model_type <- all_combos[array_number, 1] %>% as.numeric()

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

# Make a covariance matrix, and order data the same.
model_covar <- ape::vcv.phylo(model_tree)
row.names(model_data) <- model_data$tree_tip
model_data <- model_data[row.names(model_covar),]

#table(model_data$sexual_score)

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
  temp_seasonality_z = standardize(temp_log, two_sd = TRUE),
  rain_seasonality_z = standardize(bio15, two_sd = TRUE),
  npp_z = standardize(npp_sqrt, two_sd = TRUE),
  gen_z = standardize(gen_log, two_sd = TRUE),
  survival_z = standardize(adult_survival, two_sd = TRUE),
  breeding_z = standardize(breed_log, two_sd = TRUE),
  longevity_z = standardize(long_log, two_sd = TRUE),
  hatch_z = standardize(hatchling_pc1, two_sd = TRUE),
  chick_z = standardize(chick_pc1, two_sd = TRUE)
)

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1


# library(GGally)
# model_data %>% select(body_mass_log_z, survival_z, longevity_z, hatch_z, gen_z) %>% ggcorr(label = TRUE)




# 
# test <- lm(sexual_score ~ territory_binary*trophic_binary + migration_binary + trophic_binary*temp_seasonality_z +
#                   trophic_binary + trophic_binary*centroid_z, data = model_data)
# 
# 
# car::vif(test)

# # 
# library(phylolm)
# test <- phylolm(sexual_score ~ territory_binary*trophic_binary + migration_binary +
#                   trophic_binary + trophic_binary*centroid_z, data = model_data, phy = model_tree)

# summary(test)
# 
# test <- phylolm(sexual_score ~ territory_binary*trophic_binary + migration_binary +
#                   trophic_binary + trophic_binary*temp_seasonality_z + hatch_z, data = model_data, phy = model_tree)
# 
# summary(test)


# Predictors used in all models.
basic_preds <- "territory_binary + migration_binary + trophic_binary "

# Model one. 
extra_preds_1 <- "+ trophic_binary*centroid_z"

# Model two.
extra_preds_2 <- "+ body_mass_log_z + trophic_binary*centroid_z"

# Model threes.
extra_preds_3 <- "+ trophic_binary*temp_seasonality_z"
extra_preds_4 <- "+ trophic_binary*rain_seasonality_z"
extra_preds_5 <- "+ trophic_binary*npp_z"

# Model fours.
extra_preds_6 <- "+ gen_z"
extra_preds_7 <- "+ survival_z"
extra_preds_8 <- "+ longevity_z"

# Model fives.
extra_preds_9 <- "+ body_mass_log_z + trophic_binary*temp_seasonality_z"
extra_preds_10 <- "+ body_mass_log_z + trophic_binary*rain_seasonality_z"
extra_preds_11 <- "+ body_mass_log_z + trophic_binary*npp_z"

# Model sixes.
extra_preds_12 <- "+ trophic_binary*centroid_z + gen_z"
extra_preds_13 <- "+ trophic_binary*centroid_z + survival_z"
extra_preds_14 <- "+ trophic_binary*centroid_z + longevity_z"

# Without trophic interaction.
extra_preds_15 <- "+ temp_seasonality_z"
extra_preds_16 <- "+ rain_seasonality_z"
extra_preds_17 <- "+ npp_z"

# And with body mass too.
extra_preds_18 <- "+ body_mass_log_z + temp_seasonality_z"
extra_preds_19 <- "+ body_mass_log_z + rain_seasonality_z"
extra_preds_20 <- "+ body_mass_log_z + npp_z"

# No interaction.
extra_preds_21 <- "+ centroid_z"
extra_preds_22  <- "+ body_mass_log_z"
extra_preds_23  <- "+ centroid_z + body_mass_log_z"

# Full models.
extra_preds_24  <- "+ trophic_binary*temp_seasonality_z + trophic_binary*npp_z + body_mass_log_z"
extra_preds_25  <- "+ trophic_binary*temp_seasonality_z + trophic_binary*npp_z + survival_z + longevity_z"

# Last one.
extra_preds_26  <- "+ trophic_binary*centroid_z + survival_z*longevity_z"

# Extra ones now.
extra_preds_27 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*territory_binary + hatch_z"
extra_preds_28 <- "+ trophic_binary*temp_seasonality_z + hatch_z"

# Do ones with chick hatching instead.
extra_preds_29 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*territory_binary + chick_z"
extra_preds_30 <- "+ trophic_binary*temp_seasonality_z + chick_z"

# See how generation length affects hatching and chick success.
extra_preds_31 <- "+ trophic_binary*temp_seasonality_z + chick_z + gen_z"
extra_preds_32 <- "+ trophic_binary*temp_seasonality_z + hatch_z + gen_z"

# Ones with latitude and seasonality together.
extra_preds_33 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*territory_binary + trophic_binary*centroid_z"
extra_preds_34 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*territory_binary + hatch_z"

extra_preds_35 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*territory_binary + hatch_z"
extra_preds_36 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*territory_binary + hatch_z"

# The three good models.
extra_preds_37 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*centroid_z + gen_z"
extra_preds_38 <- "+ trophic_binary*temp_seasonality_z + gen_z"
extra_preds_39 <- "+ trophic_binary*centroid_z + gen_z"

# Lods of extra models.
extra_preds_40 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*centroid_z + gen_z + hatch_z"
extra_preds_41 <- "+ trophic_binary*temp_seasonality_z + gen_z + hatch_z"
extra_preds_42 <- "+ trophic_binary*centroid_z + gen_z + hatch_z"

extra_preds_43 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*centroid_z + gen_z + chick_z"
extra_preds_44 <- "+ trophic_binary*temp_seasonality_z + gen_z + chick_z"
extra_preds_45 <- "+ trophic_binary*centroid_z + gen_z + chick_z"

extra_preds_46 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*centroid_z + trophic_binary*territory_binary + gen_z + hatch_z"
extra_preds_47 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*territory_binary + gen_z + hatch_z"
extra_preds_48 <- "+ trophic_binary*centroid_z + trophic_binary*territory_binary + gen_z + hatch_z"

extra_preds_49 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*centroid_z + trophic_binary*territory_binary + gen_z + chick_z"
extra_preds_50 <- "+ trophic_binary*temp_seasonality_z + trophic_binary*territory_binary + gen_z + chick_z"
extra_preds_51 <- "+ trophic_binary*centroid_z + trophic_binary*territory_binary + gen_z + chick_z"

extra_preds_52 <- "+ temp_seasonality_z + trophic_binary*centroid_z + gen_z + hatch_z"
extra_preds_53 <- "+ temp_seasonality_z + gen_z + hatch_z"

# Last models.
extra_preds_54 <- "+ gen_z + centroid_z"
#extra_preds_54 <- "+ gen_z + trophic_binary*centroid_z"
extra_preds_55 <- "+ gen_z + npp_z + temp_seasonality_z"
extra_preds_56 <- "+ gen_z + trophic_binary*npp_z + trophic_binary*temp_seasonality_z"

#extra_preds_54 <- "+ centroid_z"
#extra_preds_54 <- "+ trophic_binary*centroid_z"
extra_preds_57 <- "+ npp_z + temp_seasonality_z"
extra_preds_58 <- "+ trophic_binary*npp_z + trophic_binary*temp_seasonality_z"

extra_preds_59 <- "+ chick_z + centroid_z"
extra_preds_60 <- "+ chick_z + trophic_binary*centroid_z"
extra_preds_61 <- "+ chick_z + npp_z + temp_seasonality_z"
extra_preds_62 <- "+ chick_z + trophic_binary*npp_z + trophic_binary*temp_seasonality_z"

extra_preds_63 <- "+ chick_z + gen_z + centroid_z"
#extra_preds_62 <- "+ chick_z + gen_z + trophic_binary*centroid_z"
extra_preds_64 <- "+ chick_z + gen_z + npp_z + temp_seasonality_z"
extra_preds_65 <- "+ chick_z + gen_z + trophic_binary*npp_z + trophic_binary*temp_seasonality_z"

extra_preds_66 <- "+ chick_z + trophic_binary*npp_z"
extra_preds_67 <- "+ gen_z + trophic_binary*npp_z"
extra_preds_68 <- "+ chick_z + gen_z + trophic_binary*npp_z"
extra_preds_69 <- "+ gen_z + npp_z"


# Assemble list of models.
extra_pred_list <- list(extra_preds_1, 
                        extra_preds_2,
                        extra_preds_3,
                        extra_preds_4,
                        extra_preds_5,
                        extra_preds_6,
                        extra_preds_7,
                        extra_preds_8,
                        extra_preds_9,
                        extra_preds_10,
                        extra_preds_11,
                        extra_preds_12,
                        extra_preds_13,
                        extra_preds_14,
                        extra_preds_15,
                        extra_preds_16,
                        extra_preds_17,
                        extra_preds_18,
                        extra_preds_19,
                        extra_preds_20,
                        extra_preds_21,
                        extra_preds_22,
                        extra_preds_23,
                        extra_preds_24,
                        extra_preds_25,
                        extra_preds_26,
                        extra_preds_27,
                        extra_preds_28,
                        
                        extra_preds_29,
                        extra_preds_30,
                        extra_preds_31,
                        extra_preds_32,
                        extra_preds_33,
                        extra_preds_34,
                        extra_preds_35,
                        extra_preds_36,
                        extra_preds_37,
                        extra_preds_38,
                        extra_preds_39,
                        extra_preds_40,
                        
                        extra_preds_41,
                        extra_preds_42,
                        extra_preds_43,
                        extra_preds_44,
                        extra_preds_45,
                        extra_preds_46,
                        extra_preds_47,
                        extra_preds_48,
                        extra_preds_49,
                        
                        extra_preds_50,
                        extra_preds_51,
                        extra_preds_52,
                        extra_preds_53,
                        
                        extra_preds_54,
                        extra_preds_55,
                        extra_preds_56,
                        extra_preds_57,
                        extra_preds_58,
                        extra_preds_59,
                        
                        extra_preds_60,
                        extra_preds_61,
                        extra_preds_62,
                        extra_preds_63,
                        extra_preds_64,
                        extra_preds_65,
                        extra_preds_66,
                        extra_preds_67,
                        extra_preds_68,
                        extra_preds_69
                        )

# Assemble formula.
model_formula <- paste0("sexual_score ~ ", basic_preds, extra_pred_list[model_type], " + (1|gr(tree_tip, cov=A))")

# Create the brms formula.
brms_formula <- brmsformula(model_formula, family = cumulative(), decomp = "QR")

# Create a unique model pathway for each model.
model_pathway <- paste0("Results/Models/", model_type, "_", data_type, ".rds")

# Add un-informative priors.
normal_priors <- c(prior(normal(0,5), class="Intercept"),
                   prior(normal(0,5), class="b"),
                   prior(student_t(3,0,20), "sd"))

tic()
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
toc()

# Create pathway for plots.
plot_pathway <- paste0("Results/Plots/", model_type, "_", data_type, ".tiff")

library(bayesplot)

# Save the bayes plot.
tiff(plot_pathway, res = 150, width = 2000, height = 2000)

# Plot model results with 95% confidence intervals.
mcmc_areas(brms_model, regex_pars = "^b_[a-z]", prob = 0.95)

# Turn device off.
dev.off()


# 
# 
# #test sample size.
# x <- 3000
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
