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


################################################################################
                              #### End ####
################################################################################