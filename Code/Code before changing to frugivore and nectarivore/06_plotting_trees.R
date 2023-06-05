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
plot_tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]

# Make a consensus tree for plotting.
#consensus_tree <- consensus(plot_tree_2, check.labels = TRUE)
consensus_tree <- plot_tree

# Drop tips to create trees.
full_tree <- drop.tip(plot_tree, setdiff(plot_tree$tip.label, full_data$tree_tip))

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
clade_colours <- c(#NA,   #"#FFFFFF" , # Blank space for legend title.
                   "#7494EA", # Palaeognathae
                   "#8D0801", # Galloanseres
                   "#77AD78", # Strisores
                   "#5941A9", # Columbimorphae
                   "#988F2A",# Otidimorphae
                 #  "#390040", 
                  # "#472836" 
                  # "#502F4C"
                  # "#70587C",
                  # "#883955",
                
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
family_plot <-ggtree(full_family_tree, layout="fan", #open.angle = 5, 
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
  
  
  ### Code to check how the birds line up.
  
  
  # geom_fruit(geom=geom_bar, mapping=aes(y=node, x=1),
  #            pwidth=0.005, fill = "black", colour = "black", orientation="y", offset = -0.205,
  #            stat="identity", width=1, alpha = 0.5) +
  # 
  # geom_fruit(geom=geom_bar, mapping=aes(y=node, x=1),
  #            pwidth=0.0025, fill = "black", colour = "black", orientation="y", offset = -0.075,
  #            stat="identity", width=1, alpha = 0.5) +
  # 
  # geom_fruit(geom=geom_bar, mapping=aes(y=node, x=1),
  #            pwidth=0.0025, fill = "black", colour = "black", orientation="y", offset = 0.125,
  #            stat="identity", width=1, alpha = 0.5) +
  
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
  #tree_image(99, "Data/phylopic/peregrin.png", i_size*1.2, inner, 3) + #public
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


################################################################################
             ##### Run family level brms analysis #####

# Average traits for family.
var_data <- model_data %>% group_by(family) %>% 
  summarise(higher_clade = first(higher_clade),
            tree_tip = first(tree_tip),
            sexual_mean = mean(sexual_score),
            clade_sum = length(sexual_score),
            seasonal_mean = mean(temp_log),
            devo_mean = mean(chick_pc1),
            trophic_mean = mean(as.numeric(as.factor(trophic_binary)))-1,
            migration_mean  = mean(as.numeric(as.factor(migration_binary)))-1,
            territory_mean  = mean(as.numeric(as.factor(territory_binary)))-1)


# Make a covariance matrix, and order data the same.
row.names(var_data) <- var_data$tree_tip
model_covar <- ape::vcv.phylo(family_tree@phylo)
var_data <- var_data[row.names(model_covar),]



# Scale continuous predictors to two SD.
library(effectsize)
var_data %<>% mutate(seasonal_mean_z = standardize(seasonal_mean, two_sd = TRUE),
                     devo_mean_z = standardize(devo_mean, two_sd = TRUE))

# Set priors.
library(brms)
normal_priors <- c(prior(normal(0,5), class="Intercept"),
                   prior(normal(0,5), class="b"),
                   prior(gamma(2,1), "sd"))

# Run the brms models.
devo_lognormal <- brm(log(sexual_mean+1) ~ devo_mean_z + (1|gr(tree_tip, cov=A)), 
                      data = var_data, 
                      data2 = list(A=model_covar),
                      family = hurdle_lognormal(), 
                      prior = normal_priors,
                      iter = 5000, chains = 2, cores = 2, init = 0,
                      backend = "cmdstanr")
seas_lognormal <- brm(log(sexual_mean+1) ~ seasonal_mean_z + (1|gr(tree_tip, cov=A)), 
                      data = var_data, 
                      data2 = list(A=model_covar),
                      family = hurdle_lognormal(), 
                      prior = normal_priors,
                      iter = 5000, chains = 2, cores = 2, init = 0,
                      backend = "cmdstanr")

# Look at model summaries.
summary(devo_lognormal)
summary(seas_lognormal)

# Extract predictions for models.
seasonal_preds <- conditional_effects(seas_lognormal)[[1]] 
devo_preds <- conditional_effects(devo_lognormal)[[1]] 


###############################################################################
              #### Create continuous side plots ####

# Create settings for scatter plots.
point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                      size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                       labs(x = "", y = NULL,  colour = "Clade"),
                       scale_colour_manual(values = clade_colours),
                       scale_y_continuous(breaks = c(0, 1, 2)),
                       theme_classic(base_size = 20),
                       theme(text = element_text(face = "bold"),
                             legend.position = "none",
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.title = element_text(size = rel(0.8))))

# Make some nice plots.
seasonal_plot <- var_data %>% ggplot(aes(x=seasonal_mean_z, y=log(sexual_mean + 1))) + point_settings +
  xlab("Seasonality") + #theme(axis.text.y = element_blank()) +
  geom_ribbon(data = seasonal_preds, aes(ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2)  + 
  geom_line(data = seasonal_preds, aes(y = estimate__), linetype = "dashed", linewidth = 1) 
seasonal_plot

devo_plot <- var_data %>% ggplot(aes(x=devo_mean_z, y=log(sexual_mean + 1))) + point_settings +
  xlab("Chick maturity") +
  geom_ribbon(data = devo_preds, aes(ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2)  + 
  geom_line(data = devo_preds, aes(y = estimate__), linetype = "dashed", linewidth = 1) 

###############################################################################
             #### Create categorical side plots ####

# Group by family and trait for binary traits.
average_family <- function(varible){
  var_data <- model_data %>% dplyr::group_by(family, !!! syms(varible)) %>% 
    summarise(higher_clade = first(higher_clade),
              sexual_mean = mean(sexual_score),
              clade_size = length(sexual_score))
  var_data %>% group_by(family) %>% summarise(clade_sum = sum(clade_size)) %>% right_join(var_data)
}

# Get territory and migration means.
teritory_data <- average_family("territory_binary")

model_data %>% group_by(family, trophic_binary)  %>%   summarise(higher_clade = first(higher_clade),
                                                              sexual_mean = mean(sexual_score),
                                                              clade_size = length(sexual_score))

migration_data <- average_family("migration_binary")
trophic_data <- average_family("trophic_binary")

# Calculate overall summary statistics.
teritory_data$log_sex <- log(teritory_data$sexual_mean + 1)
migration_data$log_sex <- log(migration_data$sexual_mean + 1)
trophic_data$log_sex <- log(trophic_data$sexual_mean + 1)

terr_means <- teritory_data %>% group_by(territory_binary) %>% summarise(log_sex_mean = mean(log_sex),
                                                           log_sex_mean_se = sd(log_sex)/sqrt(length(log_sex)))
mig_means <- migration_data %>% group_by(migration_binary) %>% summarise(log_sex_mean = mean(log_sex),
                                                                         log_sex_mean_se = sd(log_sex)/sqrt(length(log_sex)))
tro_means <- trophic_data %>% group_by(trophic_binary) %>% summarise(log_sex_mean = mean(log_sex),
                                                                         log_sex_mean_se = sd(log_sex)/sqrt(length(log_sex)))



# terr_means <- summarySE(teritory_data, measurevar="log_sex", groupvars=c("territory_binary"))
# mig_means <- summarySE(migration_data, measurevar="log_sex", groupvars=c("migration_binary"))


# Settings for binary traits.
plot_settings <- list(geom_point(aes(group=family, colour = higher_clade, 
                                     size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
                                 position = position_jitter(0.2)),
  #geom_boxplot(col = "black",fill = NA, outlier.shape = NA, width = 0.1, lwd = 1),
  labs(x = "", y = NULL,  colour = "Clade"),
  scale_colour_manual(values = clade_colours),
  scale_y_continuous(breaks = c(0, 1, 2, 3)),
  theme_classic(base_size = 20),
  theme(text = element_text(face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title.x = element_text(size = rel(0.8)),
        axis.title.y = element_text(size = rel(0.8))))

# Territory plot.
terr_plot <- teritory_data %>% ggplot(aes(x=territory_binary, y=log(sexual_mean + 1))) + 
  plot_settings + scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality")  + theme(axis.text.y = element_blank()) +
  geom_errorbar(data = terr_means, inherit.aes = FALSE,
                aes(x=territory_binary, ymin = log_sex_mean - log_sex_mean_se, ymax = log_sex_mean + log_sex_mean_se),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_means, inherit.aes = FALSE,
             aes(x=territory_binary, y = log_sex_mean, size = 4))

# Migration plots.
#migration_data$migration_binary %>% factor(levels = c("Weak", "Strong"))
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=log(sexual_mean + 1))) + 
  plot_settings + scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") + theme(axis.text.y = element_blank()) +
  geom_errorbar(data = mig_means, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = log_sex_mean - log_sex_mean_se, ymax = log_sex_mean + log_sex_mean_se),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_means, inherit.aes = FALSE,
             aes(x=migration_binary, y = log_sex_mean), size = 3.6)

tro_plot <- trophic_data %>% ggplot(aes(x=trophic_binary, y=log(sexual_mean + 1))) + 
  plot_settings + scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") + #theme(axis.text.y = element_blank()) +
  geom_errorbar(data = tro_means, inherit.aes = FALSE,
                aes(x=trophic_binary, ymin = log_sex_mean - log_sex_mean_se, ymax = log_sex_mean + log_sex_mean_se),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_means, inherit.aes = FALSE,
             aes(x=trophic_binary, y = log_sex_mean), size = 3.6)



###############################################################################
                  #### Put it all together ####

# Group all the side plots together.
side_plots <- ggarrange(terr_plot + ylab("Sexual selection"), 
                        mig_plot,# + rremove("y.text"), 
                        devo_plot + ylab("Sexual selection"), 
                        seasonal_plot,# + rremove("y.text"), 
                        nrow = 2, ncol = 2, align = "h",
                        labels = c("b", "c", "d", "e"), 
                        font.label = list(size = 24),
                        hjust = c(-3.8, -1.7, -3.8, -1.7),
                        widths = c(1.2,1))

side_plots <- ggarrange(tro_plot + ylab("Sexual selection"), 
                        mig_plot,# + rremove("y.text"), 
                        terr_plot + ylab("Sexual selection"), 
                        seasonal_plot,# + rremove("y.text"), 
                        nrow = 2, ncol = 2, align = "h",
                        labels = c("b", "c", "d", "e"), 
                        font.label = list(size = 24),
                        hjust = c(-3.8, -1.7, -3.8, -1.7),
                        widths = c(1.2,1))

side_plots <- ggarrange(seasonal_plot + ylab("Sexual selection"), 
                        mig_plot,# + rremove("y.text"), 
                        tro_plot + ylab("Sexual selection"), 
                        terr_plot,# + rremove("y.text"), 
                        nrow = 2, ncol = 2, align = "h",
                        labels = c("b", "c", "d", "e"), 
                        font.label = list(size = 24),
                        hjust = c(-3.8, -1.7, -3.8, -1.7),
                        widths = c(1.2,1))

# Group together side plots and tree.
both_plots <- ggarrange(figure_tree, side_plots, widths = c(3,2), nrow = 1,
                        labels = c("a", ""), 
                        font.label = list(size = 28))

# Export.
ggsave("Plots/Trees/tree_and_plots.png", dpi = 900, width = 17, height = 10)



#################################################################################
                  ##### Add phylo lines for plots  ######


hist(sqrt(var_data$sexual_mean), breaks = 10)
hist(var_data$sexual_mean, breaks = 100)

hist(log(var_data$sexual_mean+1), breaks = 10)
library(phylolm)


seasonal_model <- phylolm(sqrt(sexual_mean) ~ seasonal_mean, data = var_data, phy = family_tree@phylo)
devo_model <- phylolm(sqrt(sexual_mean) ~ devo_mean, data = var_data, phy = family_tree@phylo)

hist(seasonal_model$residuals)
library(car)
qqPlot(seasonal_model$residuals)

library(brms)



sesasonal_brms <- brm(sexual_mean+0.0000001  ~ seasonal_mean + (1|gr(tree_tip, cov=A)), 
    data = var_data, 
    data2 = list(A=model_covar),
    family = exponential(), 
    prior = normal_priors,
    iter = 1000, chains = 2, cores = 2, init = 0,
    backend = "cmdstanr")

sqrt_hurdle_lognormal <- brm(log(sexual_mean+1) ~ seasonal_mean_z + (1|gr(tree_tip, cov=A)), 
                      data = var_data, 
                      data2 = list(A=model_covar),
                      family = hurdle_lognormal(), 
                      prior = normal_priors,
                      iter = 2000, chains = 2, cores = 2, init = 0,
                      backend = "cmdstanr")
hurdle_gamma <- brm(sexual_mean ~ seasonal_mean_z + (1|gr(tree_tip, cov=A)), 
                      data = var_data, 
                      data2 = list(A=model_covar),
                      family = hurdle_gamma(), 
                      prior = normal_priors,
                      iter = 2000, chains = 2, cores = 2, init = 0,
                      backend = "cmdstanr")







sesasonal_brms <- brm(sexual_mean ~ seasonal_mean + (1|gr(tree_tip, cov=A)), 
                      data = var_data, 
                      data2 = list(A=model_covar),
                      family = hurdle_lognormal(), 
                      prior = normal_priors,
                      iter = 1000, chains = 2, cores = 2, init = 0,
                      backend = "cmdstanr")
sesasonal_brms <- brm(sexual_mean ~ seasonal_mean + (1|gr(tree_tip, cov=A)), 
                      data = var_data, 
                      data2 = list(A=model_covar),
                      family = hurdle_lognormal(), 
                      prior = normal_priors,
                      iter = 1000, chains = 2, cores = 2, init = 0,
                      backend = "cmdstanr")

summary(sqrt_hurdle_lognormal)
summary(hurdle_gamma)
summary(devo_lognormal)
summary(seas_lognormal)

plot(conditional_effects(devo_lognormal), points = TRUE)

conditional_effects(devo_lognormal, line_args = list(color = "black"))

seasonal_preds <- conditional_effects(seas_lognormal, line_args = list(color = "black"))[[1]] 
devo_preds <- conditional_effects(devo_lognormal, line_args = list(color = "black"))[[1]] 





conditional_effects(seas_lognormal, line_args = list(color = "black")) +
  geom_point(aes(group=higher_clade, colour = higher_clade, 
                 size = sqrt(clade_sum), alpha = sqrt(clade_sum)))
  # Create settings for scatter plots.
  point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                        size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                         labs(x = "", y = NULL,  colour = "Clade"),
                         scale_colour_manual(values = clade_colours),
                         scale_y_continuous(breaks = c(0, 1, 2)),
                         theme_classic(base_size = 20),
                         theme(text = element_text(face = "bold"),
                               legend.position = "none",
                               axis.line = element_line(size = 1),
                               axis.ticks = element_line(size = 1),
                               axis.title = element_text(size = rel(0.8)))#, 
                         #geom_smooth(method = "lm", se = FALSE, linetype = "dashed", colour = "black")
                         )

# Make some nice plots.
seasonal_plot <- var_data %>% ggplot(aes(x=seasonal_mean, y=sqrt(sexual_mean))) + point_settings +
  xlab("Seasonality") + theme(axis.text.y = element_blank())




seas_lognormal
data_term <- conditional_effects(devo_lognormal, plot = FALSE)
data_term$devo_mean_z

conditional_effects(devo_lognormal)
?conditional_effects
library(bayesplot)
pp_check(sqrt_hurdle_lognormal) + xlim(0,10)
pp_check(hurdle_gamma) + xlim(0,10)
pp_check(devo_lognormal) + xlim(0,10)

# Run brms models.
brms_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  prior = normal_priors,
  iter = 10000,
  warmup = 5000,
  chains = 2,
  thin = 10,
  cores = 20,
  init = 0,
  file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  threads = threading(10),
)


### Extra code ####


# Create a function to group by different traits.
average_clade <- function(varible){
  var_data <- model_data %>% group_by(higher_clade, !!! syms(varible)) %>% 
    summarise(sexual_mean = mean(sexual_score),
              clade_size = length(sexual_score))
  double_clades <- var_data %>% count(higher_clade) %>% filter(n == 2) %>% 
    pull(higher_clade)
  var_data %<>% filter(higher_clade %in% double_clades)
  var_data %>% group_by(higher_clade) %>% summarise(clade_sum = sum(clade_size)) %>% right_join(var_data)
}

teritory_data <- average_clade("territory_binary")
migration_data <- average_clade("migration_binary")
trophic_data <- average_clade("trophic_binary")
devo_data <- average_clade("wang_edited")


#teritory_data %<>% group_by(higher_clade) %>% summarise(clade_sum = sum(clade_size)) %>% right_join(teritory_data)

plot_settings <- list(geom_line(aes(group=higher_clade, colour = higher_clade, alpha = clade_sum), size=1), 
                           #geom_point(aes(group=higher_clade, colour = higher_clade), size=4),
                           geom_point(aes(group=higher_clade, colour = higher_clade), size = 2),
                           labs(x = "", y = NULL,  colour = "Clade"),
                           scale_colour_manual(values = clade_colours),
                           #coord_cartesian(ylim = c(0,3), clip = "off"),
                           scale_y_continuous(breaks = c(0, 1, 2, 3)),
                           #scale_x_discrete(labels = c("Low", "High")),
                           theme_classic(base_size = 20),
                           theme(text = element_text(face = "bold"),
                                 legend.position = "none",
                                 axis.line = element_line(size = 1),
                                 axis.ticks = element_line(size = 1)))

teritory_data %>% ggplot(aes(x=territory_binary, y=sqrt(sexual_mean))) + plot_settings +
  geom_boxplot(col = "grey",fill = NA, outlier.shape = NA, width = 0.1, position = position_nudge(c(-0.1,0.1)))
migration_data %>% ggplot(aes(x=migration_binary, y=sqrt(sexual_mean))) + plot_settings +
  geom_boxplot(col = "grey",fill = NA, outlier.shape = NA, width = 0.1, position = position_nudge(c(-0.1,0.1)))
trophic_data %>% ggplot(aes(x=trophic_binary, y=sqrt(sexual_mean))) + plot_settings +
  geom_boxplot(col = "grey",fill = NA, outlier.shape = NA, width = 0.1, position = position_nudge(c(-0.1,0.1)))
devo_data %>% ggplot(aes(x=wang_edited, y=sqrt(sexual_mean))) + plot_settings +
  geom_boxplot(col = "grey",fill = NA, outlier.shape = NA, width = 0.1, position = position_nudge(c(-0.1,0.1)))

mean(model_data$chick_pc1)
model_data$devo_mode_full


  var_data <- model_data %>% group_by(family) %>% 
    summarise(higher_clade = first(higher_clade),
              tree_tip = first(tree_tip),
      sexual_mean = mean(sexual_score),
              clade_sum = length(sexual_score),
              seasonal_mean = mean(temp_log),
              devo_mean = mean(chick_pc1),
      trophic_mean = mean(as.numeric(as.factor(trophic_binary)))-1,
      migration_mean  = mean(as.numeric(as.factor(migration_binary)))-1,
      territory_mean  = mean(as.numeric(as.factor(territory_binary)))-1,
      
                )

point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                      size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                         labs(x = "", y = NULL,  colour = "Clade"),
                         scale_colour_manual(values = clade_colours),
                         #coord_cartesian(ylim = c(0,3), clip = "off") +
                         #scale_y_continuous(breaks = c(0, 1, 2, 3)) +
                         #scale_x_discrete(labels = c("Low", "High")),
                         theme_classic(base_size = 20),
                         theme(text = element_text(face = "bold"),
                               legend.position = "none",
                               axis.line = element_line(size = 1),
                               axis.ticks = element_line(size = 1)), 
                         geom_smooth(method = "lm", se = FALSE, linetype = "dashed", colour = "black"))


# library(phylolm)
# 
# row.names(var_data) <- var_data$tree_tip
# 
# phylolm(sexual_mean ~ seasonal_mean, phy = order_tree@phylo, data = var_data) %>% summary()

seasonal_plot <- var_data %>% ggplot(aes(x=seasonal_mean, y=sqrt(sexual_mean))) + point_settings
devo_plot <- var_data %>% ggplot(aes(x=devo_mean, y=sqrt(sexual_mean))) + point_settings
var_data %>% ggplot(aes(x=trophic_mean, y=sqrt(sexual_mean))) + point_settings
var_data %>% ggplot(aes(x=migration_mean, y=sqrt(sexual_mean))) + point_settings
var_data %>% ggplot(aes(x=territory_mean, y=sqrt(sexual_mean))) + point_settings

var_data %>% ggplot(aes(x=migration_mean, y=sqrt(sexual_mean))) + point_settings





# Create a function to group by different traits.
average_order <- function(varible){
  var_data <- model_data %>% group_by(order, !!! syms(varible)) %>% 
    summarise(higher_clade = first(higher_clade),
              sexual_mean = mean(sexual_score),
              clade_size = length(sexual_score))
  double_clades <- var_data %>% count(order) %>% filter(n == 2) %>% 
    pull(order)
  var_data %<>% filter(order %in% double_clades)
  var_data %>% group_by(order) %>% summarise(clade_sum = sum(clade_size)) %>% right_join(var_data)
}


teritory_data <- average_order("territory_binary")
migration_data <- average_order("migration_binary")
trophic_data <- average_order("trophic_binary")

plot_settings <- list(geom_line(aes(group=order, colour = higher_clade, alpha = sqrt(clade_sum))),# size=0.5), 
                      #geom_point(aes(group=higher_clade, colour = higher_clade), size=4),
                      geom_point(aes(group=order, colour = higher_clade)),
                      geom_boxplot(col = "grey",fill = NA, outlier.shape = NA, width = 0.1, position = position_nudge(c(-0.1,0.1))),
                      labs(x = "", y = NULL,  colour = "Clade"),
                      scale_colour_manual(values = clade_colours),
                     # coord_cartesian(ylim = c(0,4), clip = "off"),
                      scale_y_continuous(breaks = c(0, 1, 2, 3)),
                      #scale_x_discrete(labels = c("Low", "High")),
                      theme_classic(base_size = 20),
                      theme(text = element_text(face = "bold"),
                            legend.position = "none",
                            axis.line = element_line(size = 1),
                            axis.ticks = element_line(size = 1)))



teritory_data %>% ggplot(aes(x=territory_binary, y=sqrt(sexual_mean))) + plot_settings
migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + plot_settings
trophic_data %>% ggplot(aes(x=trophic_binary, y=sexual_mean)) + plot_settings








average_family <- function(varible){
  var_data <- model_data %>% group_by(family, !!! syms(varible)) %>% 
    summarise(higher_clade = first(higher_clade),
              sexual_mean = mean(sexual_score),
              clade_size = length(sexual_score))
 # double_clades <- var_data %>% count(family) %>% filter(n == 2) %>% 
#    pull(family)
  #var_data %<>% filter(family %in% double_clades)
  var_data %>% group_by(family) %>% summarise(clade_sum = sum(clade_size)) %>% right_join(var_data)
}

teritory_data <- average_family("territory_binary")
migration_data <- average_family("migration_binary")
trophic_data <- average_family("trophic_binary")

plot_settings <- list(geom_line(aes(group=family, colour = higher_clade, alpha = sqrt(clade_sum))),# size=0.5), 
                      #geom_point(aes(group=higher_clade, colour = higher_clade), size=4),
                      geom_boxplot(col = "grey",fill = NA, outlier.shape = NA, width = 0.1, position = position_nudge(c(-0.1,0.1))),
                      geom_point(aes(group=family, colour = higher_clade)),
                      labs(x = "", y = NULL,  colour = "Clade"),
                      scale_colour_manual(values = clade_colours),
                      coord_cartesian(ylim = c(0,4), clip = "off"),
                      scale_y_continuous(breaks = c(0, 1, 2, 3)),
                      #scale_x_discrete(labels = c("Low", "High")),
                      theme_classic(base_size = 20),
                      theme(text = element_text(face = "bold"),
                            legend.position = "none",
                            axis.line = element_line(size = 1),
                            axis.ticks = element_line(size = 1)))

plot_settings <- list(#geom_line(aes(group=family, colour = higher_clade, alpha = sqrt(clade_sum))),# size=0.5), 
  #geom_point(aes(group=higher_clade, colour = higher_clade), size=4),
  geom_point(aes(group=family, colour = higher_clade, size = sqrt(clade_size), 
                 alpha = sqrt(clade_sum)), position = position_jitter(0.2)),
  geom_boxplot(col = "black",fill = NA, outlier.shape = NA, width = 0.1),
  labs(x = "", y = NULL,  colour = "Clade"),
  scale_colour_manual(values = clade_colours),
  #coord_cartesian(ylim = c(0,4), clip = "off"),
  scale_y_continuous(breaks = c(0, 1, 2, 3)),
  #scale_x_discrete(labels = c("Low", "High")),
  theme_classic(base_size = 20),
  theme(text = element_text(face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1)))



terr_plot <- teritory_data %>% ggplot(aes(x=territory_binary, y=sqrt(sexual_mean))) + plot_settings
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sqrt(sexual_mean))) + plot_settings


terr_means <- teritory_data %>% group_by(territory_binary) %>% 
  summarise(sex_mean = mean(sexual_mean),
            sex_se = sd(sexual_mean)/sqrt(length(sexual_mean)))

mig_means <- migration_data %>% group_by(migration_binary) %>% 
  summarise(sex_mean = mean(sexual_mean),
            sex_se = sd(sexual_mean)/sqrt(length(sexual_mean)))

teritory_data$log_sex <- log(teritory_data$sexual_mean + 1)
migration_data$log_sex <- log(migration_data$sexual_mean + 1)

terr_means <- summarySE(teritory_data, measurevar="log_sex", groupvars=c("territory_binary"))
mig_means <- summarySE(migration_data, measurevar="log_sex", groupvars=c("migration_binary"))

teritory_data$log_sex <- log(teritory_data$sexual_mean + 1)

terr_plot +
  geom_errorbar(data = terr_means, inherit.aes = FALSE,
                aes(x=territory_binary, ymin = log_sex - se*2, ymax = log_sex + se*2), 
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.1) +
  geom_point(data = terr_means, inherit.aes = FALSE,
             aes(x=territory_binary, y = log_sex))

mig_plot +
  geom_errorbar(data = mig_means, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = log_sex - se*2, ymax = log_sex + se*2), 
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.1) +
  geom_point(data = mig_means, inherit.aes = FALSE,
             aes(x=migration_binary, y = log_sex))

teritory_data %>% ggplot(aes(x=territory_binary, y=sqrt(sexual_mean))) + plot_settings + 
  geom_errorbar(data = terr_means, inherit.aes = FALSE,
                aes(x=territory_binary, ymin = sqrt(sex_mean) - sex_se*2, ymax = sqrt(sex_mean) + sex_se*2), 
              position = position_dodge(width = 1), show.legend = FALSE, width = 0.1) +
  geom_point(data = terr_means, inherit.aes = FALSE,
             aes(x=territory_binary, y = sqrt(sex_mean)))


teritory_data %>% ggplot(aes(x=territory_binary, y=sexual_mean)) + plot_settings + 
  geom_errorbar(data = terr_means, inherit.aes = FALSE,
                aes(x=territory_binary, ymin = sex_mean - sex_se*2, ymax = sex_mean + sex_se*2), 
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.1) +
  geom_point(data = terr_means, inherit.aes = FALSE,
             aes(x=territory_binary, y = sex_mean))

migration_data %>% ggplot(aes(x=migration_binary, y=sqrt(sexual_mean))) + plot_settings + 
  geom_errorbar(data = mig_means, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = sex_mean - sex_se*2, ymax = sex_mean + sex_se*2), 
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.1) +
  geom_point(data = mig_means, inherit.aes = FALSE,
             aes(x=migration_binary, y = sex_mean))


trophic_data %>% ggplot(aes(x=trophic_binary, y=sqrt(sexual_mean))) + plot_settings

side_plots <- ggarrange(terr_plot, mig_plot, devo_plot, seasonal_plot, ncol = 1)


both_plots <- ggarrange(figure_tree, side_plots, widths = c(4,1), nrow = 1)

side_plots_1 <- ggarrange(terr_plot, devo_plot, ncol = 1, align = "hv")
side_plots_2 <- ggarrange(mig_plot, seasonal_plot, ncol = 1, align = "hv")

both_plots <- ggarrange(figure_tree, side_plots_1, side_plots_2, widths = c(3,1,1), nrow = 1)


ggsave("Plots/Trees/tree_and_plots.png", dpi = 900, width = 15, height = 10)

library(gridExtra)
  png("Plots/Trees/tree_and_plots.png", res = 900, width = 15000, height = 10000)
arrangeGrob(figure_tree,side_plots_1,#layout_matrix = rbind(c(1),c(2)),
                     widths=c(4,1),respect=FALSE) %>% plot()
dev.off()
plot(final)

figure_tree

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
