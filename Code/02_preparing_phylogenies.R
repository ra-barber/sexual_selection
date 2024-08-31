###############################################################################
                       # Make some random trees #
###############################################################################


# Clean the environment.
rm(list=ls())

# Packages to load.
library(caper)
library(phytools)



################################################################################
                ##### Read the prum trees #####


# Set seed and read in trees.
set.seed(1993)
prum_trees <- read.tree("Data/Trees/Prum_merge_hackett_stage2_1000trees.tre")

# Get a set of random numbers to subset the trees.
prum_numbers <- sample(1:length(prum_trees), 50)

# Extract the multiphylo trees.
prum_model_trees <- prum_trees[c(prum_numbers)]

# Export.
write.tree(prum_model_trees, "Data/Trees/prum_trees.tre")

# Create a consensus tree.
prum_consensus_tree <- consensus.edges(prum_model_trees)

# Export the trees. 
write.tree(prum_consensus_tree, "Data/Trees/prum_consensus_tree.tre")


################################################################################
                              #### End ####
################################################################################


# tree_pathway <- "../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett"
# 
# # Read in all available Hackett back trees.
# bird_tree1 <- read.tree(paste0(tree_pathway, "1.tre"))
# bird_tree2 <- read.tree(paste0(tree_pathway, "2.tre"))
# bird_tree3 <- read.tree(paste0(tree_pathway, "3.tre"))
# bird_tree4 <- read.tree(paste0(tree_pathway, "4.tre"))
# bird_tree5 <- read.tree(paste0(tree_pathway, "5.tre"))
# bird_tree6 <- read.tree(paste0(tree_pathway, "6.tre"))
# bird_tree7 <- read.tree(paste0(tree_pathway, "7.tre"))
# bird_tree8 <- read.tree(paste0(tree_pathway, "8.tre"))
# bird_tree9 <- read.tree(paste0(tree_pathway, "9.tre"))
# bird_tree10 <- read.tree(paste0(tree_pathway, "10.tre"))
# 
# # Create a new multi-phylo object.
# all_trees <- c(bird_tree1, bird_tree2, bird_tree3, bird_tree4, bird_tree5, 
#                bird_tree6, bird_tree7, bird_tree8, bird_tree9, bird_tree10)
# 
# # Set the seed.
# set.seed(1993)

# # Get a set of random numbers to subset the trees.
# model_numbers <- sample(1:length(all_trees), 50)
# 
# # Two extra set of trees. One for plotting and one for diagnostics.
# diagnostic_numbers <- sample(1:length(all_trees), 50)
# plot_numbers <- sample(1:length(all_trees), 50)
# 
# # Extract the multiphylo trees.
# model_trees <- all_trees[c(model_numbers)]
# diagnostic_trees <- all_trees[c(diagnostic_numbers)]
# plot_trees <- all_trees[c(plot_numbers)]

# # Export the trees. 
# write.tree(model_trees, "Data/Trees/model_trees.tre")
# write.tree(diagnostic_trees, "Data/Trees/diagnostic_trees.tre")
# write.tree(plot_trees, "Data/Trees/plot_trees.tre")
# 
# 
# # Read in 50 random trees and then make a consensus.
# model_trees <- read.tree("Data/Trees/model_trees.tre")
# model_consensus_tree <- consensus.edges(model_trees, method="least.squares")
# model_consensus_tree <- consensus.edges(model_trees)
# 
# 
# # Consensus tree.
# model_consensus_tree <- consensus.edges(model_trees, p = 0.5)