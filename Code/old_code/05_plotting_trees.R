###############################################################################
                   # Plot sexual selection trees #
###############################################################################


# Packages to load.
library(magrittr)
library(tictoc)
library(caper)
library(dplyr)
library(effectsize)
library(car)
library(ggplot2)
library(ggpubr)
library(phytools)
library(stringr)
library(ggtree)
library(ggtreeExtra)


# Clear the workspace.
rm(list=ls())





################################################################################
                            #### Data ####


# Read in order and family data from jetz data file.
phylo_data <- read.csv("Data/sexual_traits.csv")

# Insert underscore to match to tree and merge.
phylo_data$tree_tip <- phylo_data$bird_tree_name %>% str_replace(" ", "_")
row.names(phylo_data) <- phylo_data$tree_tip

# Read in the tree. Pick the same one for now.
plot_tree <- read.tree("Data/Trees/plot_trees.tre")

# Make a consensus tree for plotting.
consensus_tree <- consensus(plot_tree, check.labels = TRUE)
consensus_tree <- consensus(plot_tree, check.labels = TRUE, p = 0.5)

# Make sexual selection certainty just for 1s.
phylo_data$cert_binary <- phylo_data$sexual_certainty
phylo_data$cert_binary[phylo_data$cert_binary > 1] <- 0

# Create a genus column.
phylo_data$genus <- str_split(phylo_data$bird_tree_name, pattern = " ", simplify = TRUE)[,1]


################################################################################
                    #### Add clade info ####

# Orders for PALAEOGNATHAE (in brain paper)
Palaeognathae <- c("STRUTHIONIFORMES", "RHEIFORMES", "APTERYGIFORMES", 
                   "CASUARIIFORMES", "TINAMIFORMES")

# Orders for Galloanseres. (in brain paper)
Galloanseres <- c("GALLIFORMES", "ANSERIFORMES")

# Orders for aequornithes (in brain paper)
Aequornithes <- c("GAVIIFORMES", "SPHENISCIFORMES", "PROCELLARIIFORMES", "CICONIIFORMES", 
                  "SULIFORMES", "PELECANIFORMES", "OPISTHOCOMIFORMES")

# Telluraves. (Core land birds)
Telluraves <- c("ACCIPITRIFORMES", "STRIGIFORMES", "COLIIFORMES", "LEPTOSOMIFORMES",
                "TROGONIFORMES", "BUCEROTIFORMES", "CORACIIFORMES", "PICIFORMES")

# Australaves (Austrailian radiation)
Australaves <- c("CARIAMIFORMES",  "FALCONIFORMES", "PSITTACIFORMES", "PASSERIFORMES")

# Strisores. (in brain paper)
Strisores <- c("APODIFORMES", "CAPRIMULGIFORMES")

# Columbimorphae (in brain paper)
Columbimorphae <- c("COLUMBIFORMES",  "MESITORNITHIFORMES", "PTEROCLIDIFORMES")

# Mirandornithes
Mirandornithes <- c("PHOENICOPTERIFORMES","PODICIPEDIFORMES")

# Eurpygimorphae
Eurpygimorphae <- c("EURYPYGIFORMES", "PHAETHONTIFORMES")

# Otidimorphae (in brain paper)
Otidimorphae <- c("MUSOPHAGIFORMES", "OTIDIFORMES", "CUCULIFORMES")

# orphaned groups
Charadriiformes <- c("CHARADRIIFORMES")
Gruiformes <- c("GRUIFORMES")


order_tree_data@data$clade <- NA

order_tree_data@data$clade[order_tree_data@data$order %in% Palaeognathae] <- "Palaeognathae"
order_tree_data@data$clade[order_tree_data@data$order %in% Galloanseres] <- "Galloanseres"
order_tree_data@data$clade[order_tree_data@data$order %in% Aequornithes] <- "Aequornithes"
order_tree_data@data$clade[order_tree_data@data$order %in% Telluraves] <- "Telluraves"
order_tree_data@data$clade[order_tree_data@data$order %in% Australaves] <- "Australaves"
order_tree_data@data$clade[order_tree_data@data$order %in% Strisores] <- "Strisores"
order_tree_data@data$clade[order_tree_data@data$order %in% Columbimorphae] <- "Columbimorphae"
order_tree_data@data$clade[order_tree_data@data$order %in% Mirandornithes] <- "Mirandornithes"
order_tree_data@data$clade[order_tree_data@data$order %in% Eurpygimorphae] <- "Eurpygimorphae"
order_tree_data@data$clade[order_tree_data@data$order %in% Otidimorphae] <- "Otidimorphae"
order_tree_data@data$clade[order_tree_data@data$order %in% Charadriiformes] <- "Charadriiformes"
order_tree_data@data$clade[order_tree_data@data$order %in% Gruiformes] <- "Gruiformes"


################################################################################
                      #### Summarise data ####

# Create a function to group by different traits.
group_clade <- function(data, clade){
  data %>% group_by(!!! syms(clade)) %>% 
    summarise(tree_tip = first(tree_tip),
              mean_score = mean(sexual_score),
              mean_certainty = mean(sexual_certainty),
              cert_binary = mean(cert_binary),
              trophic_binary  = mean(as.numeric(as.factor(trophic_binary))))
}

# Get genus averages.
genus_data <- group_clade(phylo_data, c("family", "genus"))

# Family.
family_data <- group_clade(phylo_data, "family")
cert_family_data <- group_clade(filter(phylo_data, sexual_certainty < 2), "family")

# Order.
order_data <- group_clade(phylo_data, "family")

################################################################################
                #### Create trees with data ####

# function to join tree and data.
data_tree <- function(tree, data){
  # Drop tips to create trees.
  drop_tree <- drop.tip(tree, setdiff(tree$tip.label, data$tree_tip))
  # Add node ID to dataframes.
  data$node <-  nodeid(drop_tree, data$tree_tip)
  # Join trees with data.
  full_join(drop_tree, data, by = "node")
}

# Tree for plotting for quick editing.
plot_tree <- consensus_tree

# Create the trees with data.
species_tree <- data_tree(plot_tree, phylo_data)
genus_tree <- data_tree(plot_tree, genus_data)
family_tree <- data_tree(plot_tree, family_data)
cert_family_tree <- data_tree(plot_tree, cert_family_data)
order_tree <- data_tree(plot_tree, order_data)





################################################################################
                 #### Prepare for plots ####


library(wesanderson)
library(shades)

# Colour blind palette.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2,4,6,7)]

# Zissou palette. 
wes_pal <- wes_palette("Zissou1", 100, type = "continuous")



###############################################################################
   #### Make a standard family tree like the size dimorphism one ####

# Family tree with coloured tips.
score_plot <- ggtree(family_tree, layout="fan", 
                     open.angle = 10,   size = 1, aes(colour=mean_score)) + 
  
  # Makes the tree the right size/zoom.
  xlim(0,32) +
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +

  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 2.5, colour = "black", offset = 1, 
            aes(label = family), fontface = 2) +
  # Add the tip points.
  geom_tippoint(mapping=aes(colour=mean_score), size=2.5, show.legend=TRUE) +
  
  theme(legend.title=element_blank(), 
        legend.position = c(0.65,0.47), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(1.7, "cm"), legend.key.height = unit(0.16, "cm"), 
        legend.text = element_text(size = 10), legend.margin = NULL, 
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = 5, t = 5, r = -40))
  
score_plot 
ggsave("Plots/Trees/family_sex_tree.png", dpi = 900, width = 10, height = 10)


# Make one with bars,
ggtree(family_tree_data, layout="fan", 
                     open.angle = 10,   size = 1) + 
  
  # Makes the tree the right size/zoom.
  #xlim(0,80) +
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



ggtree(family_tree_data, layout="rectangular", size = 0.5) + 
  
  # Makes the tree the right size/zoom.
  #xlim(0,80) +
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset =0.1,
               aes(label = family), fontface = 2) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score),
             pwidth=0.3,
             orientation="y", offset = 0.7,
             stat="identity", fill="navy", colour="navy", 
             width=0.4)
ggsave("Plots/Trees/family_rectangle.png", dpi = 900, width = 2, height = 10)


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




###############################################################################
                  #### Try a genus level tree ####

# Family tree with coloured tips.
genus_plot <- ggtree(genus_tree_data, layout="fan", 
                     open.angle = 10,   size = 0.2, colour = "grey"
                     #aes(colour=mean_score)
                      ) + 

  # Makes the tree the right size/zoom.
  #xlim(0,650) +
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  scale_fill_gradientn(colours = opacity(wes_pal, 5)) +
  
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=mean_score, fill=mean_score, colour=mean_score),
             pwidth=0.05,
             orientation="y", offset = 0.01,
             stat="identity", width=0.01) +
  
  theme(legend.title=element_blank(), 
        legend.position = c(0.65,0.47), legend.direction = "horizontal", legend.title.align = 1,
        legend.key.width = unit(1.7, "cm"), legend.key.height = unit(0.16, "cm"), 
        legend.text = element_text(size = 10), legend.margin = NULL, 
        text = element_text(face = "bold"), 
        plot.margin = margin(l = -40, b = -40, t = -40, r = -40))

genus_plot
ggsave("Plots/Trees/genus_sex_tree.png", dpi = 900, width = 10, height = 10)


plot_data <- as.data.frame(genus_data[,4])
rownames(plot_data) <- genus_data$tree_tip
#plot_data[order_tree$tip.label,]
gheatmap(genus_plot, plot_data, width=.3, offset = 5)

###############################################################################
#### Try an order tree  ####

order_plot <- ggtree(order_tree_data, layout="fan", 
       open.angle = 10,   size = 0.2
       #aes(colour=mean_score)
) + geom_tiplab(size = 2.5, colour = "black", offset = 1, 
                aes(label = order), fontface = 2)
plot_data <- as.data.frame(order_data[,3])
rownames(plot_data) <- order_data$tree_tip
#plot_data[order_tree$tip.label,]
gheatmap(order_plot, plot_data, width=.3, offset = 5) + scale_fill_viridis_c()

hist(sqrt(order_data$mean_score), breaks = 200)

ggtree(order_tree_data, layout="rectangular", size = 0.5) + 
  
  # Makes the tree the right size/zoom.
  #xlim(0,80) +
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset =0.1,
                aes(label = order), fontface = 2) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=sqrt(mean_score)),
             pwidth=0.3,
             orientation="y", offset = 0.5,
             stat="identity", fill="navy", colour="navy", 
             width=0.4) +
  geom_hilight(data = order_data, mapping = aes(node=node),
               type = "roundrect")

ggsave("Plots/Trees/order_rectangle.png", dpi = 900, width = 4, height = 6)


# Orders for PALAEOGNATHAE (in brain paper)
Palaeognathae <- c("STRUTHIONIFORMES", "RHEIFORMES", "APTERYGIFORMES", 
                   "CASUARIIFORMES", "TINAMIFORMES")

# Orders for Galloanseres. (in brain paper)
Galloanseres <- c("GALLIFORMES", "ANSERIFORMES")


# Orders for aequornithes (in brain paper)
Aequornithes <- c("GAVIIFORMES", "SPHENISCIFORMES", "PROCELLARIIFORMES", "CICONIIFORMES", 
  "SULIFORMES", "PELECANIFORMES", "OPISTHOCOMIFORMES")

# Telluraves.
Telluraves <- c("ACCIPITRIFORMES", "STRIGIFORMES", "COLIIFORMES", "LEPTOSOMIFORMES",
  "TROGONIFORMES", "BUCEROTIFORMES", "CORACIIFORMES", "PICIFORMES")

# Australaves
Australaves <- c("CARIAMIFORMES",  "FALCONIFORMES", "PSITTACIFORMES", "PASSERIFORMES")

# Strisores. (in brain paper)
Strisores <- c("APODIFORMES", "CAPRIMULGIFORMES")

# Columbimorphae (in brain paper)
Columbimorphae <- c("COLUMBIFORMES",  "MESITORNITHIFORMES", "PTEROCLIDIFORMES")

# Mirandornithes
Mirandornithes <- c("PHOENICOPTERIFORMES","PODICIPEDIFORMES")

# Eurpygimorphae
Eurpygimorphae <- c("EURYPYGIFORMES", "PHAETHONTIFORMES")

# Otidimorphae (in brain paper)
Otidimorphae <- c("MUSOPHAGIFORMES", "OTIDIFORMES", "CUCULIFORMES")

# orphaned groups
Charadriiformes <- c("CHARADRIIFORMES")
Gruiformes <- c("GRUIFORMES")


order_tree_data@data$clade <- NA

order_tree_data@data$clade[order_tree_data@data$order %in% Palaeognathae] <- "Palaeognathae"
order_tree_data@data$clade[order_tree_data@data$order %in% Galloanseres] <- "Galloanseres"
order_tree_data@data$clade[order_tree_data@data$order %in% Aequornithes] <- "Aequornithes"
order_tree_data@data$clade[order_tree_data@data$order %in% Telluraves] <- "Telluraves"
order_tree_data@data$clade[order_tree_data@data$order %in% Australaves] <- "Australaves"
order_tree_data@data$clade[order_tree_data@data$order %in% Strisores] <- "Strisores"
order_tree_data@data$clade[order_tree_data@data$order %in% Columbimorphae] <- "Columbimorphae"
order_tree_data@data$clade[order_tree_data@data$order %in% Mirandornithes] <- "Mirandornithes"
order_tree_data@data$clade[order_tree_data@data$order %in% Eurpygimorphae] <- "Eurpygimorphae"
order_tree_data@data$clade[order_tree_data@data$order %in% Otidimorphae] <- "Otidimorphae"
order_tree_data@data$clade[order_tree_data@data$order %in% Charadriiformes] <- "Charadriiformes"
order_tree_data@data$clade[order_tree_data@data$order %in% Gruiformes] <- "Gruiformes"
 


ggtree(order_tree_data, layout="rectangular", size = 0.5, ladderize = FALSE) + 
  
  # Makes the tree the right size/zoom.
  #xlim(0,80) +
  # Add the colour scheme.
  scale_color_gradientn(colours = opacity(wes_pal, 5)) +
  
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset =0.1,
              aes(label = order), fontface = 2) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=sqrt(mean_score)),
             pwidth=0.3,
             orientation="y", offset = 0.5,
             stat="identity", fill="navy", colour="navy", 
             width=0.4) +
  scale_fill_manual(values = family_colours) +
  geom_highlight(data = order_data, mapping = aes(node=node, fill = clade),
               type = "gradient", align = "left", gradient.length.out = 2) +
  theme(legend.key.width = unit(0.4, "cm"), legend.key.height = unit(1, "cm"), 
        legend.text = element_text(size = 10), 
        text = element_text(face = "bold"))

ggsave("Plots/Trees/order_rectangle.png", dpi = 900, width = 6, height = 10)


# Colour palette for 11 families.
family_colours <- c(
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
           "#EB6424") # 11.Brighter orange colour.



# Colour palette for 11 families.
family_colours <- c(
  "#136F63", # 1. Dark green
           "#E0CA3C", # 2. Citrine (dark matt yellow)
           "#390040", # 3. Dark purple
           "#690375", # 3. Patriarch (Lighter Purple)
           "#44CCFF", # 4. Light Blue
           "#3D405B", # 4. Dark matt blue (Independence)
           "#05299E", #4. International Klein Blue (Dark Blue)
           "#7494EA", # 5. Cornflower (matt blue)
           "#C8B8DB", # 6. Thistle col (Light purple)
           "#77AD78", # 7. Forest green crayola (light-ish matt green)
           "#DBA159", # 8. Earth yellow (kinda mustardy yellow colour)
           "#AA4465", # 9. Irresistible (Matt pink)
           "#A89B8C", # 10. Grullo (Brownish grey)
           "#B8B8D1", # 10. Lavender Grey
           "#8F5C38", # 11. Coyote Brown
           "#8D0801" # 11/ Dark Red
           "#EB6424" # Brighter orange colour.
)



ggtree(order_tree_data, layout="rectangular", size = 0.5, ladderize = FALSE) + 
  # Add black tip labels for family, offset slightly from the tree.
  geom_tiplab(size = 1.5, colour = "black", offset =0.1,
              aes(label = order), fontface = 2) +
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=sqrt(mean_score), fill = clade, colour = clade),
             pwidth=0.3,
             orientation="y", offset = 0.5,
             stat="identity", 
             width=0.4) +
  scale_fill_manual(values = family_colours) +
  scale_colour_manual(values = family_colours) +
  geom_highlight(data = order_data, mapping = aes(node=node, fill = clade),
                 type = "gradient", align = "left", gradient.length.out = 2) +
  theme(legend.key.width = unit(0.4, "cm"), legend.key.height = unit(1, "cm"), 
        legend.text = element_text(size = 10), legend.position = "none",
        text = element_text(face = "bold"))


legend_order <- c("Cotingidae",  "Rhinocryptidae", "Dendrocolaptidae", "Scleruridae", 
                  "Furnariidae", "Thamnophilidae", "Grallariidae", "Tityridae", 
                  "Pipridae", "Tyrannidae", "Pipromorphidae")


ggtree(order_tree_data, layout="rectangular", size = 0.5, ladderize = FALSE) + 
  # Add black tip labels for family, offset slightly from the tree.
  geom_fruit(geom=geom_bar,
             mapping=aes(y=node, x=sqrt(mean_score), fill = clade, colour = clade),
             pwidth=0.3,
             orientation="y", offset = 0.025,
             stat="identity", 
             width=0.4) +
  scale_fill_manual(values = family_colours,
                    guide = guide_legend(title.position = "top", 
                                         nrow = 6,
                                         override.aes = list(size = 3, shape = 16))) +
  scale_colour_manual(values = family_colours) +
  guides(color = FALSE) +
  geom_highlight(data = order_data, mapping = aes(node=node, fill = clade),
                 type = "gradient", align = "left", gradient.length.out = 2, show.legend = FALSE) +
  theme(legend.text = element_text(size = 10),
        text = element_text(face = "bold"), 
        legend.position = c(0.1, 0.725))

ggsave("Plots/Trees/order_rectangle.png", dpi = 900, width = 8, height = 10)



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





# # Leftover neonaves
# c("APODIFORMES", "CAPRIMULGIFORMES", "CHARADRIIFORMES", "COLUMBIFORMES", 
#   "CUCULIFORMES", "EURYPYGIFORMES", "GRUIFORMES", "MESITORNITHIFORMES", 
#   "MUSOPHAGIFORMES", "OTIDIFORMES", "PHAETHONTIFORMES", "PHOENICOPTERIFORMES", 
#   "PODICIPEDIFORMES", "PTEROCLIDIFORMES")




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



both_plot <- ggtree(order_tree_data, layout="fan", 
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


both_plot <- ggtree(order_tree_data, size = 1, aes(colour=cert_binary)) + 
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
