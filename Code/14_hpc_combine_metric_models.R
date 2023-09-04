###############################################################################
                ##### Combine multivariate models #####
###############################################################################

# This script does something new. Pretty sick eh.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(tictoc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(brms)
library(stringr)

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                  #### Read in models ####

# Read in the new files.
model_files <- list.files(path = c("Results/Models/Sex_metrics/"), full.names = T, include.dirs = FALSE, recursive = FALSE)

# Bateman Models.
bateman_centered_linear  <- combine_brms("/bateman_centered_linear", model_files)
bateman_centered_phylo  <- combine_brms("/bateman_centered_phylo", model_files)

bateman_uncentered_linear  <- combine_brms("/bateman_uncentered_linear", model_files)
bateman_uncentered_phylo  <- combine_brms("/bateman_uncentered_phylo", model_files)

# OSS models.
oss_centered_linear  <- combine_brms("/oss_centered_linear", model_files)
oss_centered_phylo  <- combine_brms("/oss_centered_phylo", model_files)

oss_uncentered_linear  <- combine_brms("/oss_uncentered_linear", model_files)
oss_uncentered_phylo  <- combine_brms("/oss_uncentered_phylo", model_files)

# Testes models. 
testes_centered_linear  <- combine_brms("/testes_centered_linear", model_files)
testes_centered_phylo <- combine_brms("/testes_centered_phylo", model_files)

testes_uncentered_linear <- combine_brms("/testes_uncentered_linear", model_files)
testes_uncentered_phylo <- combine_brms("/testes_uncentered_phylo", model_files)

gc()


# Export the models.
saveRDS(object = bateman_centered_linear, file =  "Results/Models/Combined_models/Sex_metrics/bateman_centered_linear_models.rds")
saveRDS(object = bateman_centered_phylo, file =  "Results/Models/Combined_models/Sex_metrics/bateman_centered_phylo_models.rds")

saveRDS(object = bateman_uncentered_linear, file =  "Results/Models/Combined_models/Sex_metrics/bateman_uncentered_linear_models.rds")
saveRDS(object = bateman_uncentered_phylo, file =  "Results/Models/Combined_models/Sex_metrics/bateman_uncentered_phylo_models.rds")

# OSS models.
saveRDS(object = oss_centered_linear, file =  "Results/Models/Combined_models/Sex_metrics/oss_centered_linear_models.rds")
saveRDS(object = oss_centered_phylo, file =  "Results/Models/Combined_models/Sex_metrics/oss_centered_phylo_models.rds")

saveRDS(object = oss_uncentered_linear, file =  "Results/Models/Combined_models/Sex_metrics/oss_uncentered_linear_models.rds")
saveRDS(object = oss_uncentered_phylo, file =  "Results/Models/Combined_models/Sex_metrics/oss_uncentered_phylo_models.rds")

# Testes models. 
saveRDS(object = testes_centered_linear, file =  "Results/Models/Combined_models/Sex_metrics/testes_centered_linear_models.rds")
saveRDS(object = testes_centered_phylo, file =  "Results/Models/Combined_models/Sex_metrics/testes_centered_phylo_models.rds")

saveRDS(object = testes_uncentered_linear, file =  "Results/Models/Combined_models/Sex_metrics/testes_uncentered_linear_models.rds")
saveRDS(object = testes_uncentered_phylo, file =  "Results/Models/Combined_models/Sex_metrics/testes_uncentered_phylo_models.rds")


###############################################################################
                             #### END ####
###############################################################################


