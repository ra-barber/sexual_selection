###############################################################################
                   ##### Combine metric models #####
###############################################################################

# This script combines models run on the cluster.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
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
model_files <- list.files(path = c("Results/Models/Sex_metrics/"), 
                          full.names = T, include.dirs = FALSE, recursive = FALSE)

# Bateman Models.
bateman_centered_phylo  <- combine_brms("/bateman_centered_phylo", model_files)
bateman_uncentered_phylo  <- combine_brms("/bateman_uncentered_phylo", model_files)

# OSS models.
oss_centered_phylo  <- combine_brms("/oss_centered_phylo", model_files)
oss_uncentered_phylo  <- combine_brms("/oss_uncentered_phylo", model_files)

# OSS models.
oss_centered_phylo  <- combine_brms("/oss_centered_phylo", model_files)
oss_uncentered_phylo  <- combine_brms("/oss_uncentered_phylo", model_files)
sensoss_centered_phylo  <- combine_brms("/oss_centered_phylo", model_files)

# Testes models. 
testes_centered_phylo <- combine_brms("/testes_centered_phylo", model_files)
testes_uncentered_phylo <- combine_brms("/testes_uncentered_phylo", model_files)

gc()

# Export models.
first_half <- "Z:/home/sexual_selection/Results/Models/Metrics/"
export_model <- function(model, name){
  pathway <- paste0(first_half, name, "_phylo_models.rds")
  saveRDS(object = model, file = pathway)
}


# Export the models.
export_model(bateman_centered_phylo, "bateman_centered")
export_model(bateman_uncentered_phylo, "bateman_uncentered")

# OSS models.
export_model(oss_centered_phylo, "oss_centered")
export_model(oss_uncentered_phylo, "oss_uncentered")
export_model(sensoss_centered_phylo, "sensoss_centered")

# Testes models. 
export_model(testes_centered_phylo, "testes_centered")
export_model(testes_uncentered_phylo, "testes_uncentered")



###############################################################################
                             #### END ####
###############################################################################