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
library(wesanderson)
library(shades)
library(janitor)

# Clear the workspace.
rm(list=ls())

getwd()
read.csv("")

################################################################################
                            #### Data ####


# Read in modelling data.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- model_data$birdtree_name %>% str_replace(" ", "_")
row.names(model_data) <- model_data$tree_tip

# Read in all the data, including the species we don't model.
#bird_traits <- read.csv("Data/birdtree_ecotraits_24_01_2023.csv") %>% clean_names()
bird_traits <- read.csv("Data/birdtree_ecotraits_09_05_2023.csv") %>% clean_names()

# Select the traits we want to model.
full_data <- bird_traits %>% dplyr::select(
  # Taxonomy.
  birdtree_name, order, family,
  # Sexual selection.                                   
  sexual_score, sexual_certainty,
  # Social Selection.
  territory)

# Check the data.
skimr::skim(full_data)

full_data$tree_tip <- full_data$birdtree_name %>% str_replace(replacement = "_", pattern = " ")
row.names(full_data) <- full_data$tree_tip

# Read in the tree. Pick the same one for now.
plot_tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]
plot_tree

# Drop tips to create trees.
plot_tree <- drop.tip(plot_tree, setdiff(plot_tree$tip.label, full_data$tree_tip))

# Add node ID to dataframes.
model_data$node <-  nodeid(plot_tree, model_data$tree_tip)
full_data$node <-  nodeid(plot_tree, full_data$tree_tip)

# Read in clade function and assign.
source("Code/clade_function.R")
model_data %<>% assign_clades()
full_data %<>% assign_clades()

full_data %<>% assign_prum_clades()

# Create padded tip labels so they're all the same length.
full_data$tip_label <- str_pad(full_data$family, max(nchar(full_data$family)), side="right", pad=" ")


################################################################################
                      #### Summarise data ####

# Create a function to group by different traits.
group_clade <- function(data, clade, tree = plot_tree){
  # Group by and summarise.
  data %<>% group_by(!!! syms(clade)) %>% 
    summarise(tree_tip = first(tree_tip),
              higher_clade = first(higher_clade),
              mean_score = mean(sexual_score))
  # Drop tips to create trees.
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$tree_tip))
  # Add node ID to dataframes.
  data$node <-  nodeid(tree, data$tree_tip)
  return(data)
}

# Create a function to group by different traits.
full_group_clade <- function(data, clade, tree = plot_tree){
  # Group by and summarise.
  data %<>% group_by(!!! syms(clade)) %>% 
    summarise(tree_tip = first(tree_tip),
              higher_clade = first(higher_clade),
              mean_score = mean(sexual_score))
  # Drop tips to create trees.
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$tree_tip))
  # Add node ID to dataframes.
  data$node <-  nodeid(tree, data$tree_tip)
  return(data)
}

# Family.
family_data <- group_clade(model_data, "family")
full_family_data <- full_group_clade(full_data, "family")


################################################################################
                #### Create trees with data ####

# function to join tree and data.
data_tree <- function(tree = plot_tree, data){
  # Drop tips to create trees.
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$tree_tip))
  # Join trees with data.
  full_join(tree, data, by = "node")
}

# Create the trees with data.
species_tree <- data_tree(plot_tree, model_data)
family_tree <- data_tree(plot_tree, family_data)
full_family_tree <- data_tree(plot_tree, full_family_data)



################################################################################
                       #### Make Palettes ####

# Zissou palette. 
wes_pal <- wes_palette("Zissou1", 100, type = "continuous")

# Reorder colours for plotting.
clade_colours <- c(#NA,   #"#FFFFFF" , # Blank space for legend title.
                   "#7494EA", # Palaeognathae
                   "#8D0801", # Galloanseres
                   "#77AD78", # Strisores
                   "#5941A9", # Columbimorphae
                   "#988F2A",# Otidimorphae
                   "#AA4465", # Gruiformes
                   "#05299E", # Mirandornithes
                   "#C8B8DB", # Charadriifromes
                   "#DBA159", # Eurpygimorphae
                   "#136F63", # Aequornithes
                   "#A89B8C", # Opisthocomiformes
                   "#EB6424", # Afroaves
                   "#E0CA3C") # Australaves

# Add the names so it plots in order.
names(clade_colours) <- c(#"", 
                          "Palaeognathae", "Galloanseres", "Strisores", 
                          "Columbimorphae", "Otidimorphae", "Gruiformes",
                           "Mirandornithes", "Charadriiformes", "Eurpygimorphae",
                           "Aequornithes", "Opisthocomiformes",
                          "Afroaves", "Australaves")

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

names(prum_clade_colours) <- c("Palaeognathae", "Galloanserae", "Strisores",
                          "Columbaves", "Gruiformes", "Aequorlitornithes", 
                          "Opisthocomiformes", "Accipitriformes", "Strigiformes",
                          "Coraciimorphae", "Australaves")




# # Extra colours:
# #"#136F63", # 1. Dark green
# "#E0CA3C", # 2. Citrine (dark matt yellow)
# #"#390040", # 3. Dark purple
# "#690375", # 3. Patriarch (Lighter Purple)
# "#44CCFF", # 4. Light Blue
# "#3D405B", # 4. Dark matt blue (Independence)
# #"#05299E", #4. International Klein Blue (Dark Blue)
# "#7494EA", # 5. Cornflower (matt blue)
# "#C8B8DB", # 6. Thistle col (Light purple)
# "#77AD78", # 7. Forest green crayola (light-ish matt green)
# "#DBA159", # 8. Earth yellow (kinda mustardy yellow colour)
# "#AA4465", # 9. Irresistible (Matt pink)
# "#A89B8C", # 10. Grullo (Brownish grey)
# "#B8B8D1", # 10. Lavender Grey
# "#8F5C38", # 11. Coyote Brown
# #"#8D0801" # 11/ Dark Red
#"#EB6424" # Brighter orange colour.


# Make a function to add images to the tips of a tree.
tree_image <- function(node_num, pathway, image_size = i_size, nudge_right = 0, nudge_up = 1){
  geom_tiplab(aes(subset=node == {{node_num}}, image=pathway),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size = image_size, 
              nudge_x = nudge_right, nudge_y = nudge_up)
}


###############################################################################
                    #### Create the family tree ####


# Plot a circular tree with highlighted clades and navy bars.
family_plot <- ggtree(full_family_tree, layout="fan", #open.angle = 5, 
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
  scale_fill_manual(values = prum_clade_colours, name = "Clade", na.value = NA,
                   limits = names(prum_clade_colours)) +
  scale_colour_manual(values = prum_clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 7, title.position = "top")) + 
  # Add black tip labels for family, offset inside the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1) + 
  
  # Makes the tree the right size/zoom.
  xlim(0,92) +
  theme(text = element_text(size = 12, face = "bold"), 
        #legend.title=element_text(vjust = -5.5, hjust = 0.11), 
        legend.position = c(0.55,0.42), 
        legend.direction = "horizontal", 
        
        #legend.title.align = 2,
        legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(colour = NA, fill = NA, linetype='solid'),
        legend.spacing.y = unit(0.2, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.margin = margin(-9, 6, 6, 6),
        plot.margin = margin(l = -100, b = -185, t = -175, r = -100))
  


# Saved variables for adding images to the tree.
outer <- -5
inner <- -20
i_size <- 0.03 

# Plot with tree with bird images.
figure_tree <- family_plot +
  ## Add the clade images.
  tree_image(193, "Data/phylopic/ostrich.png", i_size*1.15, inner+1, 1) +
  tree_image(186, "Data/phylopic/duck.png", i_size*1.15, inner+1, 4) +
  tree_image(173, "Data/phylopic/hummingbird.png", i_size*1.05, inner+2.5, 0) +
  tree_image(164, "Data/phylopic/Grus canadensis 2.png", i_size*1.1, inner+3, 1) +
  tree_image(131, "Data/phylopic/jacana.png", i_size*1.1, inner+3.5, -2) + #public
  tree_image(129, "Data/phylopic/buttonquail.png", i_size*1.1, inner+4, 4) +
  tree_image(156, "Data/phylopic/albatross.png", i_size*1.3, inner+3, 0) +
  tree_image(143, "Data/phylopic/nycticorax.png", i_size, inner+3, 0) +
  tree_image(121, "Data/phylopic/accipiter.png", i_size*1.6, inner+1, 0) + #public
  tree_image(115, "Data/phylopic/hoopoe.png", i_size, inner+2, -1) + #public
  tree_image(103, "Data/phylopic/toucan.png", i_size*1.3, inner+1, 5) +
  tree_image(98, "Data/phylopic/amazona.png", i_size*1, inner+1, 3) +
  tree_image(87, "Data/phylopic/huet.png", i_size*1.1, inner+1, 3) +
  tree_image(83, "Data/phylopic/menura.png", i_size*2, inner-1, 11) +
  tree_image(50, "Data/phylopic/corvus.png", i_size*1.5, inner+3, 3) + #public
  tree_image(47, "Data/phylopic/paradisaeidae.png", i_size*1.4, inner+5, 12) + #public
  tree_image(57, "Data/phylopic/euryceros.png", i_size*1.05, inner+4, 7) + #public
  tree_image(39, "Data/phylopic/swallow.png", i_size*1.8, inner-1, 7) + #public
  tree_image(24, "Data/phylopic/sturnus.png", i_size*1.2, inner+2, 4) + #public
  tree_image(7, "Data/phylopic/bullfinch.png", i_size*1.25, inner-1, 0) #public
  
# Save the tree.
ggsave("Plots/Trees/figure_prum_tree.png", dpi = 900, width = 10, height = 10)


###############################################################################
       #### Put it all together with script to make side plots ####


# Group together side plots and tree.
both_plots <- ggarrange(figure_tree, side_plots, widths = c(3,2), nrow = 1,
                        labels = c("a", ""), 
                        font.label = list(size = 28))

# Export.
ggsave("Plots/Trees/tree_and_plots.png", dpi = 900, width = 17, height = 10)
ggsave("Plots/Trees/tree_and_plots.pdf", width = 17, height = 10)


################################################################################ 
                   #### Make a non-bold version #####

# Plot a circular tree with highlighted clades and navy bars.
family_plot <- ggtree(full_family_tree, layout="fan", #open.angle = 5, 
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
  scale_fill_manual(values = prum_clade_colours, name = "Clade", na.value = NA,
                    limits = names(prum_clade_colours)) +
  scale_colour_manual(values = prum_clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 7, title.position = "top")) + 
  # Add black tip labels for family, offset inside the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1) + 
  
  # Makes the tree the right size/zoom.
  xlim(0,92) +
  theme(text = element_text(size = 12), 
        #legend.title=element_text(vjust = -5.5, hjust = 0.11), 
        legend.position = c(0.54,0.415),     # original
        #legend.position = c(0.57,0.43), 
        legend.direction = "horizontal", 
        
        #legend.title.align = 2,
        legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 11.5), 
        legend.title = element_text(size = 12), 
        legend.background = element_rect(colour = NA, fill = NA, linetype='solid'),
        legend.spacing.y = unit(0.2, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.margin = margin(-9, 6, 6, 6),
        plot.margin = margin(l = -100, b = -185, t = -175, r = -100))
figure_tree <- family_plot +
  ## Add the clade images.
  tree_image(193, "Data/phylopic/ostrich.png", i_size*1.15, inner+1, 1) +
  tree_image(186, "Data/phylopic/duck.png", i_size*1.15, inner+1, 4) +
  tree_image(173, "Data/phylopic/hummingbird.png", i_size*1.05, inner+2.5, 0) +
  tree_image(164, "Data/phylopic/Grus canadensis 2.png", i_size*1.1, inner+3, 1) +
  tree_image(131, "Data/phylopic/jacana.png", i_size*1.1, inner+3.5, -2) + #public
  tree_image(129, "Data/phylopic/buttonquail.png", i_size*1.1, inner+4, 4) +
  tree_image(156, "Data/phylopic/albatross.png", i_size*1.3, inner+3, 0) +
  tree_image(143, "Data/phylopic/nycticorax.png", i_size, inner+3, 0) +
  tree_image(121, "Data/phylopic/accipiter.png", i_size*1.6, inner+1, 0) + #public
  tree_image(115, "Data/phylopic/hoopoe.png", i_size, inner+2, -1) + #public
  tree_image(103, "Data/phylopic/toucan.png", i_size*1.3, inner+1, 5) +
  tree_image(98, "Data/phylopic/amazona.png", i_size*1, inner+1, 3) +
  tree_image(87, "Data/phylopic/huet.png", i_size*1.1, inner+1, 3) +
  tree_image(83, "Data/phylopic/menura.png", i_size*2, inner-1, 11) +
  tree_image(50, "Data/phylopic/corvus.png", i_size*1.5, inner+3, 3) + #public
  tree_image(47, "Data/phylopic/paradisaeidae.png", i_size*1.4, inner+5, 12) + #public
  tree_image(57, "Data/phylopic/euryceros.png", i_size*1.05, inner+4, 7) + #public
  tree_image(39, "Data/phylopic/swallow.png", i_size*1.8, inner-1, 7) + #public
  tree_image(24, "Data/phylopic/sturnus.png", i_size*1.2, inner+2, 4) + #public
  tree_image(7, "Data/phylopic/bullfinch.png", i_size*1.25, inner-1, 0) #public


# Group together side plots and tree.
both_plots <- ggarrange(figure_tree, side_plots, widths = c(3,2), nrow = 1,
                        labels = c("a", ""), 
                        font.label = list(size = 30))

# Export.
ggsave("Plots/Trees/nobold_tree_and_plots.png", dpi = 900, width = 17, height = 10)
ggsave("Plots/Trees/nobold_tree_and_plots.pdf", width = 17, height = 10)





###############################################################################
                  #### Add clade labels ####

# Plot a circular tree with highlighted clades and navy bars.
family_plot <- ggtree(full_family_tree, layout="fan", #open.angle = 180, 
                     size = 0.5, colour = "grey", alpha = 0.5) +
  # Add the sexual selection bars.
  geom_fruit(geom=geom_bar, mapping=aes(y=node, x=mean_score), pwidth=0.25,
             orientation="y", offset = 0.1, stat="identity", 
             fill="navy", colour="navy", width=0.4) + 
  # Add the highlighting for clades.
  geom_fruit(geom=geom_bar, mapping=aes(y=node, x=1, fill = higher_clade),
             pwidth=0.05, colour = NA, orientation="y", offset = -0.325,
             stat="identity", width=1, alpha = 0.5) + 
  # Colours.
  scale_fill_manual(values = prum_clade_colours, name = "Clade", na.value = NA,
                    limits = names(prum_clade_colours)) +
  scale_colour_manual(values = prum_clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 7, title.position = "top")) + 
  # Add black tip labels for family, offset inside the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1) + 
  
  # Makes the tree the right size/zoom.
  xlim(0,99) +
  theme(text = element_text(size = 12, face = "bold"), 
        #legend.title=element_text(vjust = -5.5, hjust = 0.11), 
        legend.position = "none", 
        legend.direction = "horizontal", 
        #legend.title.align = 2,
        legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(colour = NA, fill = NA, linetype='solid'),
        legend.spacing.y = unit(0.2, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(l = -100, b = -100, t = -90, r = -100))
family_plot

# Make a function to add images to the tips of a tree.
tree_label <- function(node_num, clade_label, label_angle = 0){
  geom_cladelabel(node=node_num, label=clade_label,
                  angle = label_angle, hjust=-0.1, align = TRUE,
                  geom = "text", offset = 1.5, barsize = NA)
}


node_positions <- c(383, # Palaeognathae
                    376, # Galloanseres
                    369, # Strisores
                    365, # Columbimorphae
                    367, # Otidimorphae
                    360, # Gruiformes
                    341, # Mirandornithes
                    325, # Charadriiformes
                    358, # Eurpygimorphae
                    343, # Aequornithes
                    123, # Opisthocomiformes
                    303, # Afroaves
                    204) # Australaves

node_positions <- c(383, # Palaeognathae
                    376, # Galloanserae
                    369, # Strisores
                    364, # Columbaves
                    360, # Gruiformes
                    323, # Aequorlitornithes
                    123, # Opisthocomiformes
                    322, # Accipitriformes
                    321, # Strigiformes
                    303, # Coraciimorphae
                    204) # Australaves









clade_names <- c("Palaeognathae", "Galloanseres", "Strisores", 
                 "Columbimorphae", "Otidimorphae", "Gruiformes",
                 "Mirandornithes", "Charadriiformes", "Eurpygimorphae",
                 "Aequornithes", "Opisthocomiformes",
                 "Afroaves", "Australaves")

clade_names <- c("1", "Galloanseres", "Strisores", 
                 "4", "5", "Gruiformes",
                 "6", "Charadriiformes", "7",
                 "Aequornithes", "8",
                 "Afroaves", "Australaves")
clade_names <- c("Palaeognathae", "Galloanserae", "Strisores",
  "Columbaves", "Gruiformes", "Aequorlitornithes", 
  "Opisthocomiformes", "Accipitriformes", "Strigiformes",
  "Coraciimorphae", "Australaves")

clade_names <- c("1", "Galloanserae", "Strisores",
                 "2", "3", "Aequorlitornithes", 
                 "4", "5", "6",
                 "Coraciimorphae", "Australaves")


angles <- c(90, # Palaeognathae
            110, # Galloanserae
            125, # Strisores
            90, # Columbaves
            90, # Gruiformes
            180, # Aequorlitornithes
            90, # Opisthocomiformes
            90, # Accipitriformes
            90, # Strigiformes
            255, # Coraciimorphae
            180) # Australaves

hjusts <- c(0.5, # Palaeognathae
            0.5, # Galloanserae
            0.5, # Strisores
            0.5, # Columbaves
            0.5, # Gruiformes
            0.35, # Aequorlitornithes
            0.5, # Opisthocomiformes
            0.5, # Accipitriformes
            0.5, # Strigiformes
            0.5, # Coraciimorphae
            0.5) # Australaves


clade_label_info <- data.frame(node = node_positions, 
                               clade = clade_names,
                               angles = angles,
                               horizontal = FALSE,
                               hjust = hjusts,
                               vjust = -0.1)

family_plot + geom_cladelab(data=clade_label_info,
                            mapping=aes(
                              node=node,
                              label=clade,
                              #horizontal=horizontal,
                              hjust=hjust,
                              vjust=vjust,
                              #angle = angles
                              ),
                            align = TRUE,
                            angle = "auto",
                            horizontal = FALSE,
                            fontsize = 3, offset.text = 0.7,
                            fontface = 2)
  

  
# family_plot +
#   geom_cladelabel(node=204, label="Australaves", angle = 0, hjust=-0.1, 
#                   geom = "text", offset = 1.5, barsize = NA) +
#   geom_cladelabel(node=303, label="Afroaves",  
#                   angle = 58, hjust=-0.1, offset = 1.5, barsize = NA)
?geom_cladelabel2


family_plot +
  geom_cladelabel(node=204, label="Australaves", 
                  angle = 0, hjust=-0.1, color = clade_colours["Australaves"], 
                  geom = "text", offset = 1, barsize = 5) +
  geom_cladelab(node=303, label="Afroaves", 
                  angle = 0, hjust=-0.1, color = clade_colours["Afroaves"], 
                  geom = "text", offset.text = 2, barsize = 5)
  
family_plot +
  geom_cladelabel(node=204, label="Australaves", angle = 0, hjust=-0.1, 
                  geom = "text", offset = 1, barsize = NA) +
  geom_cladelabel(node=303, label="Afroaves",  
                  angle = 55, hjust=-0.1, offset = 1, barsize = NA, extend = TRUE)





ggtree(full_family_tree, layout="fan", #open.angle = 5, 
        size = 0.5, colour = "grey", alpha = 0.5, aes(group = higher_clade))  +
  geom_cladelabel(node=204, label="Australaves", 
                  angle = 0, hjust=-0.1, color = clade_colours["Australaves"], 
                  geom = "text", offset = 1, barsize = 5, extend = TRUE) +
  geom_cladelabel(node=303, label="Afroaves",  
                angle = 55, hjust=-0.1, color = clade_colours["Afroaves"], 
                geom = "text", offset = 1, offset.text = 1, barsize = 5, extend = TRUE)









?geom_cladelabel

family_plot + geom_nodepoint(aes(label = node, col = higher_clade), size = 0.5)
?geom_label2
family_data
MCRA

nodes <- family_data %>% select(higher_clade, node) %>% 
  filter(higher_clade == "Australaves") %>% pull(node)

nodes_2 <- family_data %>% select(higher_clade, node) %>% 
  filter(higher_clade == "Afroaves") %>% pull(node)

MRCA(full_family_tree, nodes)
MRCA(full_family_tree, nodes_2)

names(clade_colours) <- c(#"", 
  "Palaeognathae", "Galloanseres", "Strisores", 
  "Columbimorphae", "Otidimorphae", "Gruiformes",
  "Mirandornithes", "Charadriiformes", "Eurpygimorphae",
  "Aequornithes", "Opisthocomiformes",
  "Afroaves", "Australaves")
###############################################################################
                   #### Arc style tree ####


# Plot a circular tree with highlighted clades and navy bars.
family_plot <-ggtree(full_family_tree, layout="fan", open.angle = 180, 
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
  scale_fill_manual(values = clade_colours, name = "Clade", na.value = NA,
                    limits = names(clade_colours)) +
  scale_colour_manual(values = clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 7, title.position = "top")) + 
  # Add black tip labels for family, offset inside the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1) + 
  
  # Makes the tree the right size/zoom.
  xlim(0,92) +
  theme(text = element_text(size = 12, face = "bold"), 
        #legend.title=element_text(vjust = -5.5, hjust = 0.11), 
        legend.position = c(0.4,0.55), 
        legend.direction = "horizontal", 
        
        #legend.title.align = 2,
        legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 10, face = "bold"), 
        legend.background = element_rect(colour = NA, fill = NA, linetype='solid'),
        legend.spacing.y = unit(0.2, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
family_plot


ggtree(full_family_tree, layout="fan", open.angle = 10, 
       size = 0.5, colour = "grey", alpha = 0.5) +
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1) +
  geom_cladelabel(node=77, label="Clade Label 1", angle = 120, hjust=-0.1) 
?MCRA
full_family_data
?geom_cladelabel

library(magick)
f <- tempfile(fileext=".pdf")
ggsave(filename = f, plot = family_plot, width = 17, height = 10)
x <- image_read(f, density=300)
x <- image_trim(x)
image_ggplot(x)



# Saved variables for adding images to the tree.
outer <- -5
inner <- -20
i_size <- 0.02 


# Plot with tree with bird images.
figure_tree <- family_plot +
  ## Add the clade images.
  tree_image(193, "Data/phylopic/ostrich.png", i_size*1.15, inner+1, 1) +
  tree_image(186, "Data/phylopic/duck.png", i_size*1.15, inner+1, 4) +
  tree_image(173, "Data/phylopic/hummingbird.png", i_size*1.05, inner+2.5, 0) +
  tree_image(164, "Data/phylopic/Grus canadensis 2.png", i_size*1.1, inner+3, 1) +
  tree_image(131, "Data/phylopic/jacana.png", i_size*1.1, inner+3.5, -2) + #public
  tree_image(129, "Data/phylopic/buttonquail.png", i_size*1.1, inner+4, 4) +
  tree_image(156, "Data/phylopic/albatross.png", i_size*1.3, inner+3, 0) +
  tree_image(143, "Data/phylopic/nycticorax.png", i_size, inner+3, 0) +
  tree_image(121, "Data/phylopic/accipiter.png", i_size*1.6, inner+1, 0) + #public
  tree_image(115, "Data/phylopic/hoopoe.png", i_size, inner+2, -1) + #public
  tree_image(103, "Data/phylopic/toucan.png", i_size*1.3, inner+1, 5) +
  tree_image(98, "Data/phylopic/amazona.png", i_size*1, inner+1, 3) +
  tree_image(87, "Data/phylopic/huet.png", i_size*1.1, inner+1, 3) +
  tree_image(83, "Data/phylopic/menura.png", i_size*2, inner-1, 11) +
  tree_image(50, "Data/phylopic/corvus.png", i_size*1.5, inner+3, 3) + #public
  tree_image(47, "Data/phylopic/paradisaeidae.png", i_size*1.4, inner+5, 12) + #public
  tree_image(57, "Data/phylopic/euryceros.png", i_size*1.05, inner+4, 7) + #public
  tree_image(39, "Data/phylopic/swallow.png", i_size*1.8, inner-1, 7) + #public
  tree_image(24, "Data/phylopic/sturnus.png", i_size*1.2, inner+2, 4) + #public
  tree_image(7, "Data/phylopic/bullfinch.png", i_size*1.25, inner-1, 0) #public
figure_tree
