###############################################################################
                     ##### Make a consensus tree #####
###############################################################################

# This script does something new. Pretty sick eh.


# Clean the environment.
rm(list=ls())

# Load packages.
library(caper)
library(phytools)

###############################################################################
                             #### Data ####

# Read in 50 random trees and then make a consensus.
prum_trees <- read.tree("Data/Trees/prum_trees.tre")

prum_consensus_tree <- consensus.edges(prum_trees)

# Export the trees. 
write.tree(prum_consensus_tree, "Data/Trees/prum_consensus_tree.tre")


###############################################################################
                             #### END ####
###############################################################################