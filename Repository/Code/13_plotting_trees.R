###############################################################################
                             # Plot trees #
###############################################################################

# A script to create a publication ready phylogenetic tree.

# Packages to load.
library(magrittr)
library(caper)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(phytools)
library(stringr)
library(ggtree)
library(ggtreeExtra)

# Clear the workspace.
rm(list=ls())

# Read in the functions. 
source("Code/functions.R")


################################################################################
                            #### Data ####


# Read in modelling data.
model_data <- read_ss_data()
row.names(model_data) <- model_data$tree_tip

# Select the traits we want to model.
full_data <- model_data %>% dplyr::select(
  # Taxonomy.
  scientific_name_bird_tree, tree_tip, order_bird_tree, family_bird_tree,
  # Sexual selection.                                   
  sexual_selection, data_certainty)

# Check the data.
skimr::skim(full_data)
row.names(full_data) <- full_data$tree_tip

# Read in the tree. Pick the same one for now.
plot_tree <- read.tree("Data/Trees/prum_consensus_tree.tre")
plot_tree

# Drop tips to create trees.
plot_tree <- drop.tip(plot_tree, setdiff(plot_tree$tip.label, full_data$tree_tip))

# Add node ID to dataframes.
full_data$node <-  nodeid(plot_tree, full_data$tree_tip)

# Read in clade function and assign.
source("Code/clade_function.R")
full_data %<>% assign_prum_clades()

# Create padded tip labels so they're all the same length.
full_data$tip_label <- str_pad(full_data$family_bird_tree, 
                               max(nchar(full_data$family_bird_tree)), 
                               side="right", pad=" ")


################################################################################
                      #### Summarise data ####


# Create a function to group by different traits.
full_group_clade <- function(data, clade, tree = plot_tree){
  # Group by and summarise.
  data %<>% group_by(!!! syms(clade)) %>% 
    summarise(tree_tip = first(tree_tip),
              higher_clade = first(higher_clade),
              mean_score = mean(sexual_selection))
  # Drop tips to create trees.
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$tree_tip))
  # Add node ID to dataframes.
  data$node <-  nodeid(tree, data$tree_tip)
  return(data)
}

# Family.
full_family_data <- full_group_clade(full_data, "family_bird_tree")

# function to join tree and data.
data_tree <- function(tree = plot_tree, data){
  # Drop tips to create trees.
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$tree_tip))
  # Join trees with data.
  full_join(tree, data, by = "node")
}

# Create the trees with data.
full_family_tree <- data_tree(plot_tree, full_family_data)


################################################################################
                       #### Plot settings ####


# Reorder colours for plotting.
prum_clade_colours <- c("#7494EA", # Palaeognathae
                        "#8D0801", # Galloanserae
                        "#77AD78", # Strisores
                        "#5941A9", # Columbaves
                        "#988F2A", # Gruiformes
                        "#136F63", # Aequorlitornithes
                        "#AA4465", # Opisthocomiformes
                        "#C8B8DB", # Accipitriformes
                        "#05299E", # Strigiformes
                        "#EB6424", # Coraciimorphae
                        "#E0CA3C") # Australaves

# Add names to colours.
names(prum_clade_colours) <- c("Palaeognathae", "Galloanserae", "Strisores",
                          "Columbaves", "Gruiformes", "Aequorlitornithes", 
                          "Opisthocomiformes", "Accipitriformes", "Strigiformes",
                          "Coraciimorphae", "Australaves")

# Change legend order.
full_family_tree@data$higher_clade %<>% 
  factor(levels = c("Palaeognathae", "Opisthocomiformes", "Galloanserae", 
                    "Accipitriformes", "Strisores", "Strigiformes",
                    "Columbaves", "Coraciimorphae", "Gruiformes", 
                    "Australaves", "Aequorlitornithes"))

# Make a function to add images to the tips of a tree.
tree_image <- function(node_num, pathway, image_size = i_size, nudge_right = 0, nudge_up = 1){
  geom_tiplab(aes(subset=node == {{node_num}}, image=pathway),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size = image_size, 
              nudge_x = nudge_right, nudge_y = nudge_up)
}


################################################################################ 
                   #### Plot the tree #####


# Saved variables for adding images to the tree.
outer <- -5
inner <- -20
i_size <- 0.02 

# Plot a circular tree with highlighted clades and navy bars.
family_plot <- ggtree(full_family_tree, layout="fan", 
                      size = 0.5, colour = "grey", alpha = 0.5) +
  # Add the sexual selection bars.
  geom_fruit(geom=geom_bar, mapping=aes(y=node, x=mean_score), pwidth=0.25,
             orientation="y", offset = 0.01, stat="identity", 
             fill="navy", colour="navy", width=0.4) + 
  # Add the highlighting for clades.
  geom_fruit(geom=geom_bar, mapping=aes(y=node, x=1, fill = higher_clade),
             pwidth=0.35, colour = NA, orientation="y", offset = -0.605,
             stat="identity", width=1, alpha = 0.5) + 
  # Colours.
  scale_fill_manual(values = prum_clade_colours, name = NULL, na.value = NA) +
  scale_colour_manual(values = prum_clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical", ncol = 2)) + 
  # Add black tip labels for family, offset inside the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family_bird_tree), fontface = 2, hjust =1) + 
  
  # Makes the tree the right size/zoom.
  xlim(0,92) +
  theme(text = element_text(size = 12), 
        legend.position = c(0.555, 0.42),
        legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 11.5), 
        legend.title = element_text(size = 12), 
        legend.background = element_rect(colour = NA, fill = NA, linetype='solid'),
        legend.spacing.y = unit(0.2, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(l = -100, b = -185, t = -175, r = -100))

# Add the clade images.
figure_tree <- family_plot +
  tree_image(193, "Data/phylopic/ostrich.png", i_size*1.15, inner+1, 1) +
  tree_image(186, "Data/phylopic/duck.png", i_size*1.15, inner+1, 4) +
  tree_image(173, "Data/phylopic/hummingbird.png", i_size*1.05, inner+2.5, 1) +
  tree_image(164, "Data/phylopic/Grus canadensis 2.png", i_size*1.1, inner+3, 0.5) +
  tree_image(131, "Data/phylopic/jacana.png", i_size*1.1, inner+3.5, -2) + 
  tree_image(129, "Data/phylopic/buttonquail.png", i_size*1.1, inner+4, 4) +
  tree_image(156, "Data/phylopic/albatross.png", i_size*1.3, inner+3, -1) +
  tree_image(143, "Data/phylopic/nycticorax.png", i_size, inner+3, -1) +
  tree_image(121, "Data/phylopic/accipiter.png", i_size*1.6, inner+2, -1) + 
  tree_image(115, "Data/phylopic/hoopoe.png", i_size, inner+2, -1) + 
  tree_image(103, "Data/phylopic/toucan.png", i_size*1.3, inner+1, 5) +
  tree_image(98, "Data/phylopic/amazona.png", i_size*1, inner+1, 3) +
  tree_image(87, "Data/phylopic/huet.png", i_size*1.1, inner+1, 3) +
  tree_image(83, "Data/phylopic/menura.png", i_size*2, inner-1, 11) +
  tree_image(50, "Data/phylopic/corvus.png", i_size*1.5, inner+3, 3) + 
  tree_image(47, "Data/phylopic/paradisaeidae.png", i_size*1.4, inner+5, 12) + 
  tree_image(57, "Data/phylopic/euryceros.png", i_size*1.05, inner+4, 7) + 
  tree_image(39, "Data/phylopic/swallow.png", i_size*1.8, inner-1, 7) + 
  tree_image(24, "Data/phylopic/sturnus.png", i_size*1.2, inner+2, 4) + 
  tree_image(7, "Data/phylopic/bullfinch.png", i_size*1.25, inner-1, 0)

# Save the tree.
ggsave("Plots/Trees/figure_prum_tree.png", dpi = 900, width = 10, height = 10)

# Read in the side plots.
side_plots_preds <- readRDS("Plots/Trees/side_plots.rds")

# Group together side plots and tree.
both_plots <- ggarrange(figure_tree, side_plots_preds, widths = c(3,2), nrow = 1,
                        labels = c("a", ""), 
                        font.label = list(size = 30))

# Export.
ggsave("Figures/Fig5_large.tiff", dpi = 600, width = 17, height = 10, compression = "lzw")


################################################################################
                           ###### End ########
################################################################################