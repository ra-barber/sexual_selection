###############################################################################
                      # Generate random phylogenies #
###############################################################################

# This script generates random phylogenies used in all downstream analysis.


# Clean the environment.
rm(list=ls())

# Packages to load.
library(caper)
library(phytools)


################################################################################
                       #### Generate phylogeny ####


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