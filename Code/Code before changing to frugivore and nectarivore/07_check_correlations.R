###############################################################################
                        ##### Check correlations #####
###############################################################################

# This script does something new. Pretty sick eh.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(skimr)
library(tictoc)
library(stringr)
library(caper)
library(dplyr)
library(janitor)
library(ggpubr)
library(effectsize)

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                             #### Data ####

# Read in some data.
model_data <- read.csv("Data/sexual_traits.csv")

# Read in a tree.
tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]


###############################################################################
                    #### Transform data as usual ####

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
  chick_z = standardize(chick_pc1, two_sd = TRUE),
  chick_sqrt_z = standardize(chick_sqrt, two_sd = TRUE)
)

###############################################################################
                    #### Check raw correlations ####

continuous_variables <- model_data %>% dplyr::select(temp_seasonality_z, npp_z, chick_sqrt_z)
library(GGally)

ggcorr(continuous_variables, label = TRUE, method =c("everything","pearson"), label_round = 3)

lat_variables <- model_data %>% dplyr::select(centroid_z, temp_seasonality_z, 
                                              terr_bi_c, migration_bi_c, trophic_level_c)

 
ggcorr(lat_variables, label = TRUE, method =c("everything","pearson"), label_round = 3)
ggcorr(lat_variables, label = TRUE, method =c("everything","spearman"), label_round = 3)



###############################################################################
                       #### Check threads ####

# Uncentered.
model_formula <- "sexual_score ~ territory_binary + migration_binary + 
    trophic_binary + chick_z + trophic_binary*temp_seasonality_z + 
    trophic_binary*territory_binary + trophic_binary*npp_z + 
    chick_z*temp_seasonality_z + chick_z*territory_binary + chick_z*npp_z"

# Centered.
model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + chick_sqrt_z + trophic_level_c*temp_seasonality_z + 
    trophic_level_c*terr_bi_c + trophic_level_c*npp_z + 
    chick_sqrt_z*temp_seasonality_z + chick_sqrt_z*terr_bi_c + chick_sqrt_z*npp_z"

# No NPP
model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + chick_z + trophic_level_c*temp_seasonality_z + 
    trophic_level_c*terr_bi_c + 
    chick_z*temp_seasonality_z + chick_z*terr_bi_c"

# No NPP centerered.
model_formula <- "sexual_score ~ territory_binary + migration_binary + 
    trophic_binary + chick_z + trophic_binary*temp_seasonality_z + 
    trophic_binary*territory_binary + 
    chick_z*temp_seasonality_z + chick_z*territory_binary"

# With chick sqrt.
model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + trophic_level_c*temp_seasonality_z  + 
    chick_sqrt_z*temp_seasonality_z"

# Simple.
model_formula <- "sexual_score ~ terr_bi_c + migration_bi_c + 
    trophic_level_c + chick_z + trophic_level_c*temp_seasonality_z"



library(car)
test_model <- lm(model_formula, data = model_data)

vif(test_model)
vif(test_model, type = "predictor")
sqrt(3)

summary(test_model)

model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)
row.names(model_data) <- model_data$tree_tip

phylolm(model_formula, model_data, tree) %>% summary()


###############################################################################
                           #### Section 4 ####

library(phylolm)


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


