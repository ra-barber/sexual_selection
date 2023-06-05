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



################################################################################
                            #### Data ####


# Read in modelling data.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- model_data$birdtree_name %>% str_replace(" ", "_")
row.names(model_data) <- model_data$tree_tip

# Read in all the data, including the species we don't model.
bird_traits <- read.csv("Data/birdtree_ecotraits_24_01_2023.csv") %>% clean_names()

# Select the traits we want to model.
full_data <- bird_traits %>% dplyr::select(
  # Taxonomy.
  birdtree_name, order, family,
  # Sexual selection.                                   
  sexual_score, sexual_certainty,
  # Social Selection.
  territory, mass, complete_latitude)

# Check the data.
skimr::skim(full_data)

full_data$tree_tip <- full_data$birdtree_name %>% str_replace(replacement = "_", pattern = " ")
row.names(full_data) <- full_data$tree_tip

# Read in the tree. Pick the same one for now.
plot_tree <- read.tree("Data/Trees/prum_trees.tre")
#prum_tree <- read.tree("Data/Trees/new_Prum_merge_hackett_stage2_1000trees.tre")
#plot_tree_2 <- prum_tree[1:50]
#plot_tree_2[[1]]$tip.label
# 
# pdf("Plots/Trees/2nd_test.pdf", width = 500, height = 500)
# plot(consensus_tree)
# dev.off()
# for (x in 1:50){
#   plot_tree[[x]]$tip.label <- rev(plot_tree[[x]]$tip.label)
# }

# prum_tree <- read.tree("Data/Trees/prum_trees.tre")
# plot_tree <- prum_tree
# Make a consensus tree for plotting.
consensus_tree <- consensus(plot_tree_2, check.labels = TRUE)
consensus_tree <- plot_tree[[1]]
#consensus_tree <- consensus(plot_tree, check.labels = TRUE, p = 0.5)
#plot_tree$tip.label <- rev(plot_tree$tip.label)

plot(consensus_tree)

# Drop tips to create trees.
full_tree <- drop.tip(consensus_tree, setdiff(consensus_tree$tip.label, full_data$tree_tip))
consensus_tree <- drop.tip(consensus_tree, setdiff(consensus_tree$tip.label, model_data$tree_tip))

# Add node ID to dataframes.
model_data$node <-  nodeid(consensus_tree, model_data$tree_tip)
full_data$node <-  nodeid(full_tree, full_data$tree_tip)

# Make sexual selection certainty just for 1s.
model_data$cert_binary <- model_data$sexual_certainty
model_data$cert_binary[model_data$cert_binary > 1] <- 0
full_data$cert_binary <- full_data$sexual_certainty
full_data$cert_binary[full_data$cert_binary > 1] <- 0

# Create a genus column.
model_data$genus <- str_split(model_data$birdtree_name, pattern = " ", simplify = TRUE)[,1]

# Read in clade function and assign.
source("Code/clade_function.R")
model_data %<>% assign_clades()
full_data %<>% assign_clades()

full_data$tip_label <- str_pad(full_data$family, max(nchar(full_data$family)), side="right", pad=" ")


################################################################################
                      #### Summarise data ####

# Tree for plotting for quick editing.
plot_tree <- consensus_tree

# Create a function to group by different traits.
group_clade <- function(data, clade, tree = plot_tree){
  # Group by and summarise.
  data %<>% group_by(!!! syms(clade)) %>% 
    summarise(tree_tip = first(tree_tip),
              higher_clade = first(higher_clade),
              mean_score = mean(sexual_score),
              mean_certainty = mean(sexual_certainty),
              cert_binary = mean(cert_binary),
              trophic_binary  = mean(as.numeric(as.factor(trophic_binary)))-1,
              migration_binary  = mean(as.numeric(as.factor(migration_binary)))-1,
              gen_length = mean(gen_length))
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
              #tip_label = first(tip_label),
              higher_clade = first(higher_clade),
              mean_score = mean(sexual_score),
              mean_certainty = mean(sexual_certainty),
              cert_binary = mean(cert_binary))
  # Drop tips to create trees.
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$tree_tip))
  # Add node ID to dataframes.
  data$node <-  nodeid(tree, data$tree_tip)
  return(data)
}

# Get genus averages.
genus_data <- group_clade(model_data, c("family", "genus"))

# Family.
family_data <- group_clade(model_data, "family")
cert_family_data <- group_clade(filter(model_data, sexual_certainty < 2), "family")
full_family_data <- full_group_clade(full_data, "family", full_tree)
prum_family_data <- full_group_clade(full_data, "family", full_tree)

# Order.
order_data <- group_clade(model_data, "order")

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
genus_tree <- data_tree(plot_tree, genus_data)
family_tree <- data_tree(plot_tree, family_data)
cert_family_tree <- data_tree(plot_tree, cert_family_data)
order_tree <- data_tree(plot_tree, order_data)

full_family_tree <- data_tree(full_tree, full_family_data)
prum_family_tree <- data_tree(full_tree, prum_family_data)


################################################################################
                       #### Make Palettes ####


# Colour blind palette.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2,4,6,7)]

# Zissou palette. 
wes_pal <- wes_palette("Zissou1", 100, type = "continuous")

# Colour palette for 12 clades.
clade_colours <- c(
           "#136F63", # 1. Dark green
           "#E0CA3C", # 2. Citrine (dark matt yellow)
           "#390040", # 3. Dark purple
           "#05299E", # 4. International Klein Blue (Dark Blue)
           "#AA4465", # 5. Irresistible (Matt pink)
           "#C8B8DB", # 6. Thistle col (Light purple)
           "#77AD78", # 7. Forest green crayola (light-ish matt green)
           "#DBA159", # 8. Earth yellow (kinda mustardy yellow colour)
           "#A89B8C", # 9. Grullo (Brownish grey)
           "#7494EA", # 10. Cornflower (matt blue)
           "#8D0801", # 11/ Dark Red
           "#EB6424") # 12.Brighter orange colour.

# Reorder colours for plotting.
clade_colours <- c("#7494EA","#8D0801", "#77AD78","#A89B8C",
                            "#05299E",   "#AA4465", "#390040",
                            "#C8B8DB",   "#136F63",   "#DBA159",
                             "#EB6424", "#E0CA3C")

# Add the names so it plots in order.
names(clade_colours) <- c("Palaeognathae", "Galloanseres", "Strisores", 
                          "Eurpygimorphae", "Mirandornithes", "Columbimorphae",
                          "Otidimorphae", "Gruimorphae", "Aequornithes", "Opisthocomiformes",
                          "Afroaves", "Australaves")

# Extra colours:
# "#136F63", # 1. Dark green
# "#E0CA3C", # 2. Citrine (dark matt yellow)
# "#390040", # 3. Dark purple
# "#690375", # 3. Patriarch (Lighter Purple)
# "#44CCFF", # 4. Light Blue
# "#3D405B", # 4. Dark matt blue (Independence)
# "#05299E", #4. International Klein Blue (Dark Blue)
# "#7494EA", # 5. Cornflower (matt blue)
# "#C8B8DB", # 6. Thistle col (Light purple)
# "#77AD78", # 7. Forest green crayola (light-ish matt green)
# "#DBA159", # 8. Earth yellow (kinda mustardy yellow colour)
# "#AA4465", # 9. Irresistible (Matt pink)
# "#A89B8C", # 10. Grullo (Brownish grey)
# "#B8B8D1", # 10. Lavender Grey
# "#8F5C38", # 11. Coyote Brown
# "#8D0801" # 11/ Dark Red
# "#EB6424" # Brighter orange colour.


# Make a function to add images to the tips of a tree.
tree_image <- function(node_num, pathway, image_size = i_size, nudge_right = 0, nudge_up = 1){
  geom_tiplab(aes(subset=node == {{node_num}}, image=pathway),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size = image_size, 
              nudge_x = nudge_right, nudge_y = nudge_up)
}


###############################################################################
                    #### Family trees ####

# Make one with bars,
ggtree(family_tree, layout="fan", open.angle = 10,   size = 1) + 
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab2(size = 2.5, colour = "black", offset =0.1,
              aes(label = family), fontface = 2) +
  
  theme(legend.title=element_blank(), 
        legend.position = c(0.65,0.47), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(1.7, "cm"), legend.key.height = unit(0.16, "cm"), 
        legend.text = element_text(size = 10), legend.margin = NULL, 
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = 5, t = 5, r = -40)) +

geom_fruit(geom=geom_bar,
           mapping=aes(y=node, x=mean_score),
           pwidth=0.1,
           orientation="y", offset = 0.5,
           stat="identity", fill="navy", colour="navy", 
           width=0.4)
ggsave("Plots/Trees/family_sex_tree_bars.png", dpi = 900, width = 10, height = 10)



# Make one with bars,
ggtree(family_tree, layout="fan", open.angle = 10, size = 0.5, colour = "grey") + 
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  
  theme(legend.title=element_blank(), 
        legend.position = c(0.65,0.1), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 10), legend.margin = NULL, 
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = 5, t = 5, r = -40)) +
  
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score, fill = higher_clade, colour = higher_clade),
             pwidth=0.25,
             orientation="y", offset = 0.01,
             stat="identity", #fill="navy", colour="navy", 
             width=0.4) + scale_fill_manual(values = clade_colours) +
  scale_colour_manual(values = clade_colours) +
 geom_highlight(data = order_data, mapping = aes(node=node, fill = higher_clade),
               type = "gradient", align = "left", gradient.length.out = 2)
  # geom_fruit(geom=geom_bar,
  #            mapping=aes(y=node, x=1, fill = higher_clade, colour = higher_clade),
  #            pwidth=0.35,colour = NA,
  #            orientation="y", offset = -0.605,
  #            stat="identity", #fill="navy", colour="navy", 
  #            width=1, alpha = 0.5)

ggsave("Plots/Trees/family_clade_highlighted.png", dpi = 900, width = 10, height = 10)

# Plot a circular tree with highlighted clades and navy bars.
family_plot <-ggtree(full_family_tree, layout="fan", #open.angle = 5, 
         size = 0.5, colour = "grey") 

# + 
  # Makes the tree the right size/zoom.
  #xlim(0,37)

family_plot
ggsave("Plots/Trees/test_tree.png", dpi = 1200, width = 25, height = 25, limitsize = FALSE, )
ggsave("Plots/Trees/test_tree.jpg", width = 80, height = 80, limitsize = FALSE)



# Rotate tips to allign clades.
family_plot <- flip(family_plot, 333, 328)  # Rotate hummingbirds.
family_plot <- flip(family_plot, 306, 320)  # Rotate hummingbirds.
family_plot <- flip(family_plot, 306, 190)  # Rotate hummingbirds.

# Saved variables for adding images to the tree.
outer <- -5
inner <- -8
i_size <- 0.03 



# Plot a circular tree with highlighted clades and navy bars.
family_plot +
  # Theme.
  theme(legend.title=element_blank(), 
        legend.position = c(0.56,0.435), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 9), legend.margin = NULL, 
        legend.background = element_blank(),
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = -75, t = -35, r = -40)) +
  
  # Add the sexual selection bars.
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.25,
             orientation="y", offset = 0.01,
             stat="identity", fill="navy", colour="navy", 
             width=0.4) + 
  scale_fill_manual(values = clade_colours) +
  scale_colour_manual(values = clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 6)) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=1, fill = higher_clade),
             pwidth=0.35,colour = NA,
             orientation="y", offset = -0.605,
             stat="identity",
             width=1, alpha = 0.5) + geom_text(aes(label=node), size = 1) +
  ## Add the clade images.
    tree_image(2, "Data/phylopic/ostrich.png", i_size, outer, 2) +
    tree_image(14, "Data/phylopic/duck.png", i_size, inner, 2) +
    tree_image(20, "Data/phylopic/hummingbird.png", i_size, outer, 0) +
    tree_image(16, "Data/phylopic/kagu.png", i_size, inner, 1) +
    tree_image(165, "Data/phylopic/flamingo.png", i_size*0.9, outer, 2) +
    tree_image(192, "Data/phylopic/pigeon.png", i_size, inner, 0) +
    tree_image(167, "Data/phylopic/bustard.png", i_size, outer, -2) +  # Have to give credit for this one on phylopic.
    tree_image(169, "Data/phylopic/rail.png", i_size, inner, 0) +
  tree_image(25, "Data/phylopic/buttonquail.png", i_size, inner, 0) +
  
    tree_image(178, "Data/phylopic/albatross.png", i_size*1.1, inner, 2) +
    tree_image(182, "Data/phylopic/nycticorax.png", i_size, inner, 2) +
  
    #tree_image(183, "Data/phylopic/pelecan.png", i_size, outer, -2) +
    #tree_image(189, "Data/phylopic/cormorant.png", i_size, outer, -2) +
    tree_image(193, "Data/phylopic/hoazin.png", i_size*1.1, inner+4, 1.5) + # Have to give credit for this one on phylopic.
    tree_image(36, "Data/phylopic/charadrius.png", i_size*1.1, outer, 0) +
     tree_image(56, "Data/phylopic/toucan.png", i_size*1.3, inner, 0) +
  
    tree_image(62, "Data/phylopic/accipiter.png", i_size*1.6, inner, 0) +
  tree_image(65, "Data/phylopic/amazona.png", i_size*1, inner, 0) +
  tree_image(80, "Data/phylopic/menura.png", i_size*2, inner, 0) +
  tree_image(112, "Data/phylopic/corvus.png", i_size*2, inner, 0) +
  
    tree_image(148, "Data/phylopic/bullfinch.png", i_size*1.2, outer, 0)

ggsave("Plots/Trees/family_clade_highlighted_2.png", dpi = 900, width = 10, height = 10)




# Plot with tiplabs.
family_plot +
  # Theme.
  theme(legend.title=element_blank(), 
        legend.position = c(0.56,0.435), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 9), legend.margin = NULL, 
        legend.background = element_blank(),
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = -75, t = -35, r = -40)) +
  
  # Add the sexual selection bars.
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.25,
             orientation="y", offset = 0.01,
             stat="identity", fill="navy", colour="navy", 
             width=0.4) + 
  scale_fill_manual(values = clade_colours) +
  scale_colour_manual(values = clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 6)) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=1, fill = higher_clade),
             pwidth=0.35,colour = NA,
             orientation="y", offset = -0.605,
             stat="identity",
             width=1, alpha = 0.5) +   geom_text(aes(label=node), size = 1) +
  ## Add the clade images.00
  tree_image(2, "Data/phylopic/ostrich.png", i_size, inner, 2) +
  tree_image(14, "Data/phylopic/duck.png", i_size, inner, 2) +
  tree_image(20, "Data/phylopic/hummingbird.png", i_size, inner, 0) +
 # tree_image(16, "Data/phylopic/kagu.png", i_size, inner, 1) +
  #tree_image(165, "Data/phylopic/flamingo.png", i_size*0.9, outer, 2) +
  tree_image(192, "Data/phylopic/pigeon.png", i_size, inner, 0) +
  #tree_image(167, "Data/phylopic/bustard.png", i_size, inner, -2) +  # Have to give credit for this one on phylopic.
  #tree_image(169, "Data/phylopic/rail.png", i_size, inner, 0) +
  tree_image(25, "Data/phylopic/buttonquail.png", i_size, inner, 0) +
  
  tree_image(178, "Data/phylopic/albatross.png", i_size*1.1, inner, 2) +
  tree_image(182, "Data/phylopic/nycticorax.png", i_size, inner, 2) +
  
  #tree_image(183, "Data/phylopic/pelecan.png", i_size, outer, -2) +
  #tree_image(189, "Data/phylopic/cormorant.png", i_size, outer, -2) +
 # tree_image(193, "Data/phylopic/hoazin.png", i_size*1.1, inner+4, 1.5) + # Have to give credit for this one on phylopic.
  tree_image(36, "Data/phylopic/charadrius.png", i_size*1.3, inner, -2) +
  tree_image(56, "Data/phylopic/toucan.png", i_size*1.3, inner, 0) +
  
  tree_image(62, "Data/phylopic/accipiter.png", i_size*1.6, inner, 0) + #public
  tree_image(65, "Data/phylopic/amazona.png", i_size*1, inner, 0) +
  tree_image(80, "Data/phylopic/menura.png", i_size*2, inner-2, 2) +
  tree_image(97, "Data/phylopic/mohoua.png", i_size*1.5, inner, 2) + #public
  tree_image(112, "Data/phylopic/corvus.png", i_size*1.5, inner, 2) + #public
  tree_image(125, "Data/phylopic/swallow.png", i_size*1.8, inner, -4) + #public
  tree_image(132, "Data/phylopic/sturnus.png", i_size*1.2, inner, -2) + #public
  
  #tree_image(144, "Data/phylopic/sugarbird.png", i_size*1.5, inner, 2) + #public
  
  tree_image(148, "Data/phylopic/bullfinch.png", i_size*1.2, inner, 0) + #public
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1)

ggsave("Plots/Trees/family_clade_highlighted_3.png", dpi = 900, width = 10, height = 10)





# Plot with tiplabs.
family_plot +
  # Theme.
  theme(legend.title=element_blank(), 
        legend.position = c(0.5,0.5), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 9), 
        legend.background = element_rect(colour = 'grey', fill = 'white', linetype='solid'),
        legend.spacing.y = unit(0.2, "cm"),
        legend.margin = margin(6, 6, 6, 6),
        #legend.box.margin = margin(b = 200, t = 20),
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = -75, t = -35, r = -40)) +
  
  # Add the sexual selection bars.
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.25,
             orientation="y", offset = 0.01,
             stat="identity", fill="navy", colour="navy", 
             width=0.4) + 
  scale_fill_manual(values = clade_colours) +
  scale_colour_manual(values = clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 6)) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=1, fill = higher_clade),
             pwidth=0.35,colour = NA,
             orientation="y", offset = -0.605,
             stat="identity",
             width=1, alpha = 0.5) +  # geom_text(aes(label=node), size = 1) +
  ## Add the clade images.00
  tree_image(2, "Data/phylopic/ostrich.png", i_size, inner, 2) +
  tree_image(14, "Data/phylopic/duck.png", i_size, inner, 2) +
  tree_image(20, "Data/phylopic/hummingbird.png", i_size, inner, 0) +
  # tree_image(16, "Data/phylopic/kagu.png", i_size, inner, 1) +
  #tree_image(165, "Data/phylopic/flamingo.png", i_size*0.9, outer, 2) +
  tree_image(192, "Data/phylopic/pigeon.png", i_size, inner, 0) +
  #tree_image(167, "Data/phylopic/bustard.png", i_size, inner, -2) +  # Have to give credit for this one on phylopic.
  #tree_image(169, "Data/phylopic/rail.png", i_size, inner, 0) +
  tree_image(25, "Data/phylopic/buttonquail.png", i_size, inner, 0) +
  
  tree_image(178, "Data/phylopic/albatross.png", i_size*1.1, inner, 2) +
  tree_image(182, "Data/phylopic/nycticorax.png", i_size, inner, 2) +
  
  #tree_image(183, "Data/phylopic/pelecan.png", i_size, outer, -2) +
  #tree_image(189, "Data/phylopic/cormorant.png", i_size, outer, -2) +
  # tree_image(193, "Data/phylopic/hoazin.png", i_size*1.1, inner+4, 1.5) + # Have to give credit for this one on phylopic.
  tree_image(36, "Data/phylopic/charadrius.png", i_size*1.3, inner, -2) +
  tree_image(56, "Data/phylopic/toucan.png", i_size*1.3, inner, 0) +
  
  tree_image(62, "Data/phylopic/accipiter.png", i_size*1.6, inner, 0) + #public
  tree_image(65, "Data/phylopic/amazona.png", i_size*1, inner, 0) +
  tree_image(80, "Data/phylopic/menura.png", i_size*2, inner-2, 2) +
  tree_image(97, "Data/phylopic/mohoua.png", i_size*1.5, inner, 2) + #public
  tree_image(112, "Data/phylopic/corvus.png", i_size*1.5, inner, 2) + #public
  tree_image(125, "Data/phylopic/swallow.png", i_size*1.8, inner, -4) + #public
  tree_image(132, "Data/phylopic/sturnus.png", i_size*1.2, inner, -2) + #public
  
  #tree_image(144, "Data/phylopic/sugarbird.png", i_size*1.5, inner, 2) + #public
  
  tree_image(148, "Data/phylopic/bullfinch.png", i_size*1.2, inner, 0) + #public
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1)



ggsave("Plots/Trees/family_clade_highlighted_4.png", dpi = 900, width = 10, height = 10)


# Plot with tiplabs.
family_plot +
  # Theme.
  theme(legend.title=element_blank(), 
        legend.position = c(0.5,0.425), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 9), 
        legend.background = element_rect(colour = 'grey', fill = 'white', linetype='solid'),
        legend.spacing.y = unit(0.2, "cm"),
        legend.margin = margin(6, 6, 6, 6),
        #legend.box.margin = margin(b = 200, t = 20),
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = -75, t = -35, r = -40)) +
  
  # Add the sexual selection bars.
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.25,
             orientation="y", offset = 0.01,
             stat="identity", fill="navy", colour="navy", 
             width=0.4) + 
  scale_fill_manual(values = clade_colours) +
  scale_colour_manual(values = clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 6)) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=1, fill = higher_clade),
             pwidth=0.35,colour = NA,
             orientation="y", offset = -0.605,
             stat="identity",
             width=1, alpha = 0.5) +

  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1) + geom_text(aes(label=node), size = 1)

ggsave("Plots/Trees/prum_highlighted.png", dpi = 900, width = 10, height = 10)


# Saved variables for adding images to the tree.
outer <- -5
inner <- -20
i_size <- 0.03 



# Plot with tiplabs.
family_plot +
  # Theme.
  theme(legend.title=element_blank(), 
        legend.position = c(0.5,0.5), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size = 9), 
        legend.background = element_rect(colour = 'grey', fill = 'white', linetype='solid'),
        legend.spacing.y = unit(0.2, "cm"),
        legend.margin = margin(6, 6, 6, 6),
        #legend.box.margin = margin(b = 200, t = 20),
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = -75, t = -35, r = -40)) +
  
  # Add the sexual selection bars.
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.25,
             orientation="y", offset = 0.01,
             stat="identity", fill="navy", colour="navy", 
             width=0.4) + 
  scale_fill_manual(values = clade_colours, limits = names(clade_colours)) +
  scale_colour_manual(values = clade_colours) +
  guides(fill = guide_legend(byrow = TRUE, nrow = 6)) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=1, fill = higher_clade),
             pwidth=0.35,colour = NA,
             orientation="y", offset = -0.605,
             stat="identity",
             width=1, alpha = 0.5) +  # geom_text(aes(label=node), size = 1) +
  ## Add the clade images.00
  tree_image(193, "Data/phylopic/ostrich.png", i_size, inner, 1) +
  
  tree_image(186, "Data/phylopic/duck.png", i_size, inner, 2) +
  tree_image(173, "Data/phylopic/hummingbird.png", i_size, inner, 0) +
  # tree_image(16, "Data/phylopic/kagu.png", i_size, inner, 1) +
  #tree_image(165, "Data/phylopic/flamingo.png", i_size*0.9, outer, 2) +
  tree_image(167, "Data/phylopic/pigeon.png", i_size, inner, 3) +
  #tree_image(167, "Data/phylopic/bustard.png", i_size, inner, -2) +  # Have to give credit for this one on phylopic.
  #tree_image(169, "Data/phylopic/rail.png", i_size, inner, 0) +
  tree_image(129, "Data/phylopic/buttonquail.png", i_size, inner, 2) +
  
  tree_image(156, "Data/phylopic/albatross.png", i_size*1.1, inner, 0) +
  tree_image(143, "Data/phylopic/nycticorax.png", i_size, inner, 0) +
  
  #tree_image(183, "Data/phylopic/pelecan.png", i_size, outer, -2) +
  #tree_image(189, "Data/phylopic/cormorant.png", i_size, outer, -2) +
  # tree_image(193, "Data/phylopic/hoazin.png", i_size*1.1, inner+4, 1.5) + # Have to give credit for this one on phylopic.
  tree_image(135, "Data/phylopic/charadrius.png", i_size*1.3, inner, 4) +
  tree_image(103, "Data/phylopic/toucan.png", i_size*1.3, inner, 0) +
  
  tree_image(121, "Data/phylopic/accipiter.png", i_size*1.6, inner, 0) + #public
  tree_image(98, "Data/phylopic/amazona.png", i_size*1, inner, 0) +
  tree_image(83, "Data/phylopic/menura.png", i_size*2, inner-2, 2) +
  tree_image(65, "Data/phylopic/mohoua.png", i_size*1.5, inner, 2) + #public
  tree_image(50, "Data/phylopic/corvus.png", i_size*1.5, inner, 2) + #public
  tree_image(39, "Data/phylopic/swallow.png", i_size*1.8, inner, -4) + #public
  tree_image(24, "Data/phylopic/sturnus.png", i_size*1.2, inner, -2) + #public
  
  #tree_image(144, "Data/phylopic/sugarbird.png", i_size*1.5, inner, 2) + #public
  
  tree_image(7, "Data/phylopic/bullfinch.png", i_size*1.2, inner, 0) + #public
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset = -0.05,
              aes(label = family), fontface = 2, hjust =1)


ggsave("Plots/Trees/prum_highlighted_2.png", dpi = 900, width = 10, height = 10)




prum_family_data %>% dplyr::select(family, higher_clade, node) %>% as.data.frame()

full_family_data %>% filter(node == 92)

###############################################################################
                    #### Try an order tree  ####

# Create basic plot.
order_plot <- ggtree(order_tree, layout="fan", open.angle = 10,   size = 0.2) + 
  geom_tiplab(size = 2.5, colour = "black", offset = 1, aes(label = order), fontface = 2)

# Create heat map.
plot_data <- as.data.frame(order_data[,3])
rownames(plot_data) <- order_data$tree_tip
gheatmap(order_plot, plot_data, width=.3, offset = 5) + scale_fill_viridis_c()

# Create one with clades highlighted.
ggtree(order_tree, layout="rectangular", size = 0.5, ladderize = FALSE) + 
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset =0.1,
              aes(label = order), fontface = 2) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=sqrt(mean_score), fill = higher_clade, colour = higher_clade),
             pwidth=0.3,
             orientation="y", offset = 0.5,
             stat="identity", 
             width=0.4) +
  scale_fill_manual(values = clade_colours) +
  scale_colour_manual(values = clade_colours) +
  geom_highlight(data = order_data, mapping = aes(node=node, fill = higher_clade),
                 type = "gradient", align = "left", gradient.length.out = 2) +
  theme(legend.key.width = unit(0.4, "cm"), legend.key.height = unit(1, "cm"), 
        legend.text = element_text(size = 10), legend.position = "none",
        text = element_text(face = "bold")) + geom_text(aes(label=node), hjust=-.3)
ggsave("Plots/Trees/order_rectangle.png", dpi = 900, width = 6, height = 10)


###############################################################################
           #### Order tree with rearranged clades ####


# Create the basic tree first.
node_plot <- ggtree(order_tree, layout="rectangular", size = 0.5, 
                    ladderize = FALSE, alpha = 0.8, colour = "darkgrey")

# Reorder the plots if p = 0.5 for consensus.
# node_plot <- flip(node_plot, 45, 57)
# node_plot <- flip(node_plot, 48, 23)
# node_plot <- flip(node_plot, 71, 8)
# node_plot <- flip(node_plot, 60, 36)
# node_plot <- flip(node_plot, 69, 67)
# node_plot <- flip(node_plot, 40, 39)
# node_plot <- flip(node_plot, 66, 29)
# node_plot <- flip(node_plot, 69, 67)
# node_plot <- flip(node_plot, 66, 59)
#  node_plot <- flip(node_plot, 58, 67)
#  node_plot <- flip(node_plot, 69, 58)
#  node_plot <- flip(node_plot, 69, 67)
#  node_plot <- flip(node_plot, 40, 39)
#  node_plot <- flip(node_plot, 69, 67)
#  node_plot <- flip(node_plot, 40, 39)
 
# Reorder for p = 1
node_plot <- flip(node_plot, 45, 54) 
node_plot <- flip(node_plot, 23, 47) # Sorts out landbird clades.
node_plot <- flip(node_plot, 8, 65) # swaps round hummingbirds
node_plot <- flip(node_plot, 62, 55) # move hoazin out the way
node_plot <- flip(node_plot, 40, 39)
node_plot <- flip(node_plot, 40, 55) # link up pigeons
node_plot <- flip(node_plot, 36, 56)
node_plot <- flip(node_plot, 36, 61)

# Saved variables for adding images to the tree.
right <- -3.3
left <- -4.5
i_size <- 0.05 

# Create a tree plot without legend first.
ordered_plot <- node_plot + 
  
  # Add the coloured bars.
  geom_fruit(geom=geom_bar, 
             mapping = aes(y=node, x=mean_score, fill = higher_clade, colour = higher_clade),
             pwidth=0.3, orientation="y", offset = 0.01, stat="identity", width=0.4) +
  
  # Add the colour scale.
  scale_fill_manual(values = clade_colours, breaks = rev(names(clade_colours)), name = "Clade") +
  scale_colour_manual(values = clade_colours, breaks = rev(names(clade_colours)), name = "Clade") +
  guides(fill = guide_legend(byrow = TRUE), colour = guide_legend(byrow = TRUE)) +
  
  # Highlight the clades.
  geom_highlight(data = order_data, mapping = aes(node=node, fill = higher_clade),
                 type = "gradient", align = "left", gradient.length.out = 2) +
  
  # Remove legend from first plot.
  theme(legend.position = "none") +
    
  ## Add the clade images.
  tree_image(2, "Data/phylopic/ostrich.png", i_size, left, 1) +
  tree_image(7, "Data/phylopic/duck.png", i_size, right, -0.5) +
  tree_image(9, "Data/phylopic/hummingbird.png", i_size, left, 0.5) +
  tree_image(8, "Data/phylopic/kagu.png", i_size, right, 0.5) +
  tree_image(25, "Data/phylopic/flamingo.png", i_size*0.9, left, 0.5) +
  tree_image(38, "Data/phylopic/pigeon.png", i_size, right, 0) +
  tree_image(27, "Data/phylopic/bustard.png", i_size, left, 0) +  # Have to give credit for this one on phylopic.
  tree_image(29, "Data/phylopic/rail.png", i_size, right, 0) +
  tree_image(34, "Data/phylopic/pelecan.png", i_size, left, -1.5) +
  tree_image(39, "Data/phylopic/hoazin.png", i_size*1.1, right, 0) + # Have to give credit for this one on phylopic.
  tree_image(11, "Data/phylopic/charadrius.png", i_size*1.1, left, 0) +
  tree_image(18, "Data/phylopic/accipiter.png", i_size*1.6, right, -2.5) +
  tree_image(22, "Data/phylopic/bullfinch.png", i_size*1.2, left, -1)
  
# Export the tree.
ggsave("Plots/Trees/3rd_rearranged_order_rectangle.png", dpi = 900, width = 6, height = 6)

# Add legend information.
legend_plot <- ordered_plot + theme(text = element_text(face = "bold"),
                      legend.key.width = unit(0.4, "cm"), 
                      legend.key.height = unit(0.1, "cm"),
                      #legend.key=element_rect(colour="black"),
                      legend.text = element_text(size = 10, face = "plain"), 
                      legend.position = c(0.12, 0.78),
                      legend.background = element_rect(fill='transparent'),
                      legend.spacing.y = unit(0.2, 'cm')) 

# Export the tree with legend.
ggsave("Plots/Trees/legend_rearranged_order_rectangle.png", dpi = 900, width = 6, height = 6)


# Code to create a legend with multiple rows/columns.
legend_settings <- list(scale_colour_manual(values = family_colours, 
                                            breaks = legend_order,
                                            guide = guide_legend(title.position = "top", 
                                                                 nrow = 2,
                                                                 override.aes = list(size = 65, shape = 16))),
                        theme(text = element_text(face = "bold"),
                              legend.text = element_text(size=rel(0.9)),
                              legend.title = element_text(size=rel(1.1)),
                              legend.title.align = 0.5,
                              legend.position = c(0.55, 0.3), 
                              legend.justification = "center", 
                              legend.direction = "horizontal",
                              legend.margin = margin(0, 600, 500, 0),
                              #legend.text.align = 0.5,
                              legend.key.height = unit(5, 'cm'),
                              #legend.key.size = unit(1, "cm")
                        ))




###############################################################################
                 #### Adding trait data to the plot ####

library(ggnewscale)


# Trying to add diet.
# Maybe have two squares that's green and red for primary vs secondary.
# Then have the alpha value for each one be related to their percentage in the clade.
node_plot + 
  
  geom_fruit(geom=geom_bar,
             mapping = aes(y=node, x = 0.1, alpha=trophic_binary), fill = "forestgreen",
             pwidth=0.05, orientation="y", offset = 0.025,
             stat="identity", width=0.5) +
geom_fruit(geom=geom_bar,
           mapping = aes(y=node, x = 0.1, alpha=migration_binary), fill = "darkred",
           pwidth=0.05, orientation="y", offset = 0.025,
           stat="identity", width=0.5) +
  scale_alpha_continuous(range = c(0.1,1)) +
  new_scale("alpha") +
  
  
  geom_fruit(geom=geom_bar,
             mapping = aes(y=node, x = 0.1, alpha=log(gen_length)), fill = "darkblue",
             pwidth=0.05, orientation="y", offset = 0.025,
             stat="identity", width=0.5) +

  scale_alpha_continuous(range = c(0.1,1)) +

  
  # Add the colour scheme.
  # scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  # scale_fill_gradientn(colours = opacity(wes_pal, 5))
  
  
  
  geom_fruit(geom=geom_bar,
                       mapping = aes(y=node, x=sqrt(mean_score), 
                                     fill = clade, colour = clade),
                       pwidth=0.3, orientation="y", offset = 0.01,
                       stat="identity", width=0.4) +
  scale_fill_manual(values = clade_colours) +
  scale_colour_manual(values = clade_colours) +
  geom_highlight(data = order_data, mapping = aes(node=node, fill = clade),
                 type = "gradient", align = "left", gradient.length.out = 2) +
  theme(legend.key.width = unit(0.4, "cm"), legend.key.height = unit(1, "cm"), 
        legend.text = element_text(size = 10), legend.position = "none",
        text = element_text(face = "bold")) +
  
  
  geom_tiplab(aes(subset=node == 2, image="Data/phylopic/ostrich.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -5, nudge_y = 1,
  ) + 
  geom_tiplab(aes(subset=node == 22, image="Data/phylopic/bullfinch.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -5, nudge_y = -1,
  ) +
  geom_tiplab(aes(subset=node == 7, image="Data/phylopic/duck.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -4, nudge_y = -0.5,
  ) +
  geom_tiplab(aes(subset=node == 9, image="Data/phylopic/hummingbird.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -5, nudge_y = 0.5,
  )  + 
  geom_tiplab(aes(subset=node == 8, image="Data/phylopic/kagu.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -4, nudge_y = 0.5,
  )  +
  geom_tiplab(aes(subset=node == 25, image="Data/phylopic/flamingo.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -5, nudge_y = 0.5
  )  + geom_tiplab(aes(subset=node == 38, image="Data/phylopic/pigeon.png"),  
                   inherit.aes = FALSE, geom="image", 
                   linetype = "solid",  size= 0.05, nudge_x = -4, nudge_y = 0
  ) +
  geom_tiplab(aes(subset=node == 18, image="Data/phylopic/accipiter.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.075, nudge_x = -4, nudge_y = -2.5
  ) + 
  geom_tiplab(aes(subset=node == 11, image="Data/phylopic/charadrius.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -5, nudge_y = 0
  ) +
  # Have to give credit for this one on phylopic.
  geom_tiplab(aes(subset=node == 27, image="Data/phylopic/bustard.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -5, nudge_y = 0
  ) +
  geom_tiplab(aes(subset=node == 34, image="Data/phylopic/pelecan.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -5, nudge_y = -1.5
  ) +
  geom_tiplab(aes(subset=node == 29, image="Data/phylopic/rail.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -4, nudge_y = 0
  ) +
  # Have to give credit for this one on phylopic.
  geom_tiplab(aes(subset=node == 39, image="Data/phylopic/hoazin.png"),  
              inherit.aes = FALSE, geom="image", 
              linetype = "solid",  size= 0.05, nudge_x = -4, nudge_y = 0
  )

ggsave("Plots/Trees/traits_order_rectangle.png", dpi = 900, width = 6, height = 6)


# Code to make a legend. 

test_plot <- ggtree(test_tree, layout="rectangular", size = 0.5, ladderize = FALSE) + 
  # Add black tip labels for family, offset slightly from the tree.
  
  geom_tiplab(size = 1.5, colour = "black", offset =2.5,
              aes(label = order), fontface = 2) +
  
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=sqrt(mean_score), fill = clade, colour = clade),
             pwidth=0.3,
             orientation="y", offset = 0.025,
             stat="identity", 
             width=0.4) +
  scale_fill_manual(values = clade_colours,
                    guide = guide_legend(title.position = "top", 
                                         nrow = 6,
                                         override.aes = list(size = 3, shape = 16))) +
  scale_colour_manual(values = clade_colours) +
  guides(color = FALSE) +
  geom_highlight(data = order_data, mapping = aes(node=node, fill = clade),
                 type = "gradient", align = "left", gradient.length.out = 2, show.legend = FALSE) +
  theme(legend.text = element_text(size = 10),
        text = element_text(face = "bold"), 
        legend.position = c(0.1, 0.725))









ggsave("Plots/Trees/order_rectangle.png", dpi = 900, width = 8, height = 10)

order_data$order

test_orders <- c("ACCIPITRIFORMES","ANSERIFORMES","APODIFORMES","APTERYGIFORMES","BUCEROTIFORMES")

test_tree <- order_tree@phylo


test_tree <- drop.tip(test_tree, setdiff(test_tree$tip.label, order_data$tree_tip[1:5]))

t1 <- order_data$tree_tip[order_data$order == "ACCIPITRIFORMES"]
t2 <- order_data$tree_tip[order_data$order == "ANSERIFORMES"]
t3 <- order_data$tree_tip[order_data$order == "APODIFORMES"]
t4 <- order_data$tree_tip[order_data$order == "APTERYGIFORMES"]
t5 <- order_data$tree_tip[order_data$order == "BUCEROTIFORMES"]

y <- ape::rotateConstr(test_tree, as.character(c(t1, t4, t5, t2, t3)))

p1=ggtree(test_tree, ladderize=F) + geom_tiplab(offset = -1)
p2=ggtree(y, ladderize=F) + geom_tiplab(offset = -1)
cowplot::plot_grid(p1, p2, ncol=2)

ape::rotateConstr(test_tree, c("Apteryx_australis",
                               "Aix_galericulata","Abeillia_abeillei",
                               "Aceros_cassidix","Accipiter_albogularis"))

y <- ape::rotateConstr(test_tree, c("Aix_galericulata", "Apteryx_australis",
                               "Abeillia_abeillei",
                               "Aceros_cassidix","Accipiter_albogularis"))

order_for_plotting<- c("STRUTHIONIFORMES", "RHEIFORMES", "APTERYGIFORMES", 
                       "CASUARIIFORMES", "TINAMIFORMES","GALLIFORMES", 
                       "ANSERIFORMES","EURYPYGIFORMES", "PHAETHONTIFORMES", "APODIFORMES", "CAPRIMULGIFORMES",
                        "COLUMBIFORMES",  "MESITORNITHIFORMES", "PTEROCLIDIFORMES",
                       "MUSOPHAGIFORMES", "OTIDIFORMES", "CUCULIFORMES","GRUIFORMES",
                       "GAVIIFORMES", "SPHENISCIFORMES", "PROCELLARIIFORMES", "CICONIIFORMES", 
                       "SULIFORMES", "PELECANIFORMES","OPISTHOCOMIFORMES", 
                       "CHARADRIIFORMES", "PHOENICOPTERIFORMES","PODICIPEDIFORMES",
                       "ACCIPITRIFORMES", "STRIGIFORMES", "COLIIFORMES", 
                       "LEPTOSOMIFORMES",
                       "TROGONIFORMES", "BUCEROTIFORMES", "CORACIIFORMES", 
                       "PICIFORMES","CARIAMIFORMES",  "FALCONIFORMES", 
                       "PSITTACIFORMES", "PASSERIFORMES")

tip_labs <- order_data$tree_tip
names(tip_labs) <- order_data$order

as.character(tip_labs[order_for_plotting])

tip_order <- order_data$tree_tip[order_data$order %in% order_for_plotting]

y <- ape::rotateConstr(test_tree, tip_order)

order_tree <- data_tree(y, order_data)

ggtree(order_tree) + geom_la








ggtree(order_tree, layout = "circular") + geom_tiplab(aes(label = order))
  
  
  
order_tree@phylo$tip.label




?rotateConstr
orders


order_tree@data

order_data

test_tree <- groupOTU(order_tree, grp, "Species")

legend_settings <- list(scale_colour_manual(values = family_colours, 
                                            breaks = legend_order,
                                            guide = guide_legend(title.position = "top", 
                                                                 nrow = 2,
                                                                 override.aes = list(size = 65, shape = 16))),
                        coord_cartesian(ylim = c(-7,4)),
                        scale_x_discrete(labels = c("Low", "High")),
                        theme_classic(base_size = 150),
                        theme(text = element_text(face = "bold"),
                              legend.text = element_text(size=rel(0.9)),
                              legend.title = element_text(size=rel(1.1)),
                              legend.title.align = 0.5,
                              legend.position = c(0.55, 0.3), 
                              legend.justification = "center", 
                              legend.direction = "horizontal",
                              legend.margin = margin(0, 600, 500, 0),
                              #legend.text.align = 0.5,
                              legend.key.height = unit(5, 'cm'),
                              #legend.key.size = unit(1, "cm")
                        ))







+ geom_hilight(data=d, aes(node=node, fill=type),
               type = "roundrect")


#create a dataframe of just dichromatism, changing row names to families to match tree tips
dichro <- plot_data
dichro<- as.data.frame(dichro)
dichro$Jetz_Name <- NULL
row.names(dichro) <- as.character(dichro$Family)
dichro$Family <- NULL
dichro$Dichro <- round(dichro$Dichro)

#make a heat map of dichromatism
gheat <- gheatmap(p=circle, data=dichro, offset = -10, width=0.1, colnames = FALSE) + 
  scale_fill_gradient(low = "navy", high = "gold", breaks=c(0,2,4,6,8,10)) +
  geom_tiplab(offset = 7, cex=5) + 
  theme(legend.title=element_blank(), legend.position = "bottom", 
        legend.key.width = unit(2, "cm"), legend.key.height = unit(0.08, "cm"), #
        legend.text = element_text(size = 12))
gheat




################################################################################
               #### Plot Jetz phylogeny ####

# Saved settings.
ggplot_settings <- list(
  
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)),
                        
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 2.5, colour = "black", offset = 1, 
               aes(label = family), fontface = 2)
  )


# Family tree with coloured tips.
score_plot <- ggtree(family_tree_data, layout="fan", 
                   open.angle = 10,   size = 1,# aes(colour=mean_score), 
                      ) + 
  # Makes the tree the right size/zoom.
  xlim(0,60) +
  # Add the tip points.
  #geom_tippoint(mapping=aes(colour=mean_score), size=2.5, show.legend=TRUE) +
  # Add the settings.
  ggplot_settings +
  
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.38,
             orientation="y", offset = 0.6,
             stat="identity", fill="navy", colour="navy", width=0.2) +

theme(legend.title=element_blank(), 
      legend.position = c(0.65,0.47), legend.direction = "horizontal", legend.title.align = 1,
      legend.key.width = unit(1.7, "cm"), legend.key.height = unit(0.16, "cm"), 
      legend.text = element_text(size = 10), legend.margin = NULL, 
      text = element_text(face = "bold"), plot.margin = margin(l = -80, b = -80, t = -80, r = -80))
  

score_plot

ggsave("Plots/Trees/family_sex_tree.tiff", dpi = 600, width = 10, height = 10)

# Family tree with coloured tips.
cert_plot <- ggtree(family_tree_data, layout="fan", 
                     open.angle = 10,   size = 1, aes(colour=binary_cert)) + 
  # Add the tip points.
  geom_tippoint(mapping=aes(colour=binary_cert), size=2.5, show.legend=TRUE) +
  # Add the settings.
  ggplot_settings 




# Family tree with coloured tips.
both_plot <- ggtree(family_tree_data, layout="fan", 
                     open.angle = 10,   size = 1, aes(colour=cert_binary)) + 
  # Add the tip points.
  #geom_tippoint(mapping=aes(colour=cert_binary), size=2.5, show.legend=TRUE) +
  # Add the settings.
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  scale_fill_gradientn(colours = opacity(wes_pal, 5)) +
# Add black tip labels for family, offset slightly from the tree.
geom_tiplab(size = 2.5, colour = "black", offset = 3.5 ,
            aes(label = family), fontface = 2) +
  
  theme(legend.title=element_blank(), 
        legend.position = c(0.65,0.48), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.08, "cm"), 
        legend.text = element_text(size = 9), legend.margin = NULL, 
        plot.margin=margin(-40,-80,-40,-80)) +
  
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score, fill = cert_binary, colour=cert_binary),
             pwidth=0.1,
             orientation="y", offset = 0.01,
             stat="identity",# fill="navy", colour="navy", 
             width=0.4)

both_plot

ggsave("Plots/Trees/family_sex_tree.tiff", dpi = 600, width = 10, height = 10)


both_plot <- ggtree(family_tree_data, layout="fan", 
                    open.angle = 10,   size = 1, aes(colour=cert_binary)) + 
  # Add the tip points.
  #geom_tippoint(mapping=aes(colour=cert_binary), size=2.5, show.legend=TRUE) +
  # Add the settings.
  ggplot_settings +
  
  theme(legend.title=element_blank(), 
        legend.position = c(0.65,0.48), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.08, "cm"), 
        legend.text = element_text(size = 9), legend.margin = NULL, plot.margin=margin(-40,-80,-40,-80)) +
  
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.1,
             orientation="y", offset = 0.5,
             stat="identity", fill="navy", colour="navy", width=0.4)

both_plot

################################################################################
                      #### Order trees ####



both_plot <- ggtree(order_tree, layout="fan", 
                    open.angle = 10,   size = 1, aes(colour=cert_binary)) + 
  # Add the tip points.
  #geom_tippoint(mapping=aes(colour=cert_binary), size=2.5, show.legend=TRUE) +
  # Add the settings.
 # ggplot_settings +
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 2.5, colour = "black", offset = 1,  align = TRUE,
              aes(label = order), fontface = 2) +
  
  theme(legend.title=element_blank(), 
        legend.position = c(0.65,0.48), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.08, "cm"), 
        legend.text = element_text(size = 9), legend.margin = NULL, plot.margin=margin(-40,-80,-40,-80)) +
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.38,
             orientation="y", offset = 0.6,
             stat="identity", fill="navy", colour="navy", width=0.2)

both_plot


both_plot <- ggtree(order_tree, size = 1, aes(colour=cert_binary)) + 
  # Add the tip points.
  #geom_tippoint(mapping=aes(colour=cert_binary), size=2.5, show.legend=TRUE) +
  # Add the settings.
  # ggplot_settings +
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 2.5, colour = "black", offset = 0.15,  align = FALSE,
              aes(label = order), fontface = 2) +
  
  theme(#legend.title=element_blank(), 
        text = element_text(face = "bold"),
        legend.position = c(0.25,0.95), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(1.8, "cm"), legend.key.height = unit(0.08, "cm"), 
        legend.text = element_text(size = 9), legend.margin = NULL) +
        #plot.margin=margin(-40,-40,-40,-40)) +
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5), 
                        guide = guide_legend(title.position = "top"),
                        name = "Certainty Percentage") +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.18,
             orientation="y", offset = 0.23,
             stat="identity", fill="navy", colour="navy", width=0.5)

both_plot


################################################################################
                      #### Species trees ####


ggtree(species_tree_data, layout="fan", 
       open.angle = 10,   size = 1, aes(colour=sexual_score)) + 
  # Add the tip points.
  #geom_tippoint(mapping=aes(colour=cert_binary), size=2.5, show.legend=TRUE) +
  # Add the settings.
  # ggplot_settings +
  scale_color_gradientn(colours = opacity(wes_pal, 5), 
                        guide = guide_legend(title.position = "top"),
                        name = "Certainty Percentage") +
  theme(legend.position = c(0.65,0.48), 
        legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(0.8, "cm"), 
        legend.key.height = unit(0.08, "cm"), 
        legend.text = element_text(size = 9), 
        legend.margin = NULL, plot.margin=margin(-40,-80,-40,-80)) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=sexual_score),
             pwidth=0.18,
             orientation="y", offset = 0.23,
             stat="identity", fill="navy", colour="navy", width=0.5)
  # Add the colour scheme.
