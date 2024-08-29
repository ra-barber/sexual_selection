###############################################################################
                       # Make some random trees #
###############################################################################

# Clean the environment.
rm(list=ls())

# Packages to load.
library(caper)

# Read in all available Hackett back trees.
bird_tree1 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett1.tre")
bird_tree2 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett2.tre")
bird_tree3 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett3.tre")
bird_tree4 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett4.tre")
bird_tree5 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett5.tre")
bird_tree6 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett6.tre")
bird_tree7 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett7.tre")
bird_tree8 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett8.tre")
bird_tree9 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett9.tre")
bird_tree10 <- read.tree("../../Phd_data/Trees/Jetz_trees/Hackett/BirdzillaHackett10.tre")

# Create a new multi-phylo object.
all_trees <- c(bird_tree1, bird_tree2, bird_tree3, bird_tree4, bird_tree5, 
               bird_tree6, bird_tree7, bird_tree8, bird_tree9, bird_tree10)

# Set the seed.
set.seed(1993)

# Get a set of random numbers to subset the trees.
model_numbers <- sample(1:length(all_trees), 50)

# Two extra set of trees. One for plotting and one for diagnostics.
diagnostic_numbers <- sample(1:length(all_trees), 50)
plot_numbers <- sample(1:length(all_trees), 50)

# Extract the multiphylo trees.
model_trees <- all_trees[c(model_numbers)]
diagnostic_trees <- all_trees[c(diagnostic_numbers)]
plot_trees <- all_trees[c(plot_numbers)]

# Export the trees. 
write.tree(model_trees, "Data/Trees/model_trees.tre")
write.tree(diagnostic_trees, "Data/Trees/diagnostic_trees.tre")
write.tree(plot_trees, "Data/Trees/plot_trees.tre")


# Read in 50 random trees and then make a consensus.
model_trees <- read.tree("Data/Trees/model_trees.tre")
library(phytools)
model_consensus_tree <- consensus.edges(model_trees, method="least.squares")
model_consensus_tree <- consensus.edges(model_trees)


# Consensus tree.
model_consensus_tree <- consensus.edges(model_trees, p = 0.5)


################################################################################
                ##### Read the prum trees #####


set.seed(1993)
prum_trees <- read.tree("Data/Trees/Prum_merge_hackett_stage2_1000trees.tre")

# Get a set of random numbers to subset the trees.
prum_numbers <- sample(1:length(prum_trees), 50)

# Extract the multiphylo trees.
prum_model_trees <- prum_trees[c(prum_numbers)]

# Export.
write.tree(prum_model_trees, "Data/Trees/prum_trees.tre")



################################################################################
                            ##### Ebird tree #####

set.seed(1993)
ebird_trees <- read.tree("Data/Trees/ebird_24may.tre")

plot(ebird_trees)


################################################################################
                              #### End ####
################################################################################