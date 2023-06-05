###############################################################################
                      ##### Preparing trait data #####
###############################################################################

# This script pre-processes the trait data so it can be used for phylogenetic
# models. Currently it's a edited version of the dichromatism project. It will 
# need extra work done, but has been done quickly to allow mapping and tree plots.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(tictoc)
library(stringr)
library(caper)
library(dplyr)
library(janitor)
library(ggpubr)
library(skimr)

# Read in the functions.
source("Code/functions.R")


###############################################################################
                        #### Read in the data #####


# Read in the life history traits.
bird_traits <- read.csv("Data/GBD_BirdTreeTaxo_ecotraits_2022_November_7th.csv") %>% clean_names()

# Select the traits we want to model.
phylo_data <- bird_traits %>% dplyr::select(
  # Taxonomy.
  bird_tree_name, order, family,
  # Parental Care.
  incubation_roles, incubation_certainty,
  # Sexual selection.                                   
  sexual_score, sexual_certainty,
  # Nest visibility.
  nest_placement, nest_certainty,
  # Social Selection.
  stability, territory, lifehistory_uncertainty, 
  # Extra hypotheses.
  cooperative, migration, trophic_level,
  # Extra controls.
  habitat, mass, centroid_latitude, range_size)

# Check the data. Still missing traits. Will drop these species now anyway.
skim(phylo_data)
phylo_data %<>% na.omit()

bio_data <- read.csv("Data/Extracts/jetz_bio1.csv")
for (x in list.files("Data/Extracts/", full.names = TRUE)[2:13]){
  bio_data <- left_join(bio_data, read.csv(x))
}




###############################################################################
            #### Create social selection metric #####


# Count both groups.
phylo_data %>% dplyr::count(territory)
phylo_data %>% dplyr::count(stability)

# Create new social selection metric.
phylo_data$social_selection <- NA
phylo_data[phylo_data$territory=="Strong" | phylo_data$stability == "Strong","social_selection"] <- "Strong"
phylo_data[is.na(phylo_data$social_selection), "social_selection"] <- "Weak"

# Check new groups.
phylo_data %>% dplyr::count(social_selection)



###############################################################################
                #### Binary model variables #####


### Sexual selection ###

# Count sexual selection
phylo_data %>% dplyr::count(sexual_score)
phylo_data %>% filter_high() %>% dplyr::count(sexual_score)

# Create a binary category for sexual selection.
phylo_data$sexual_binary <- as.character(phylo_data$sexual_score)
phylo_data[phylo_data$sexual_binary=="0","sexual_binary"] <- "Weak"
phylo_data[phylo_data$sexual_binary=="1","sexual_binary"] <- "Weak"
phylo_data[phylo_data$sexual_binary=="2","sexual_binary"] <- "Weak"
phylo_data[phylo_data$sexual_binary=="3","sexual_binary"] <- "Strong"
phylo_data[phylo_data$sexual_binary=="4","sexual_binary"] <- "Strong"

# Count again.
phylo_data %>% dplyr::count(sexual_binary)
phylo_data %>% filter_high() %>% dplyr::count(sexual_binary)

## Territory ##

# Numerical territory.
phylo_data$territory_num <- phylo_data$territory
phylo_data[phylo_data$territory_num=="None","territory_num"] <- "0"
phylo_data[phylo_data$territory_num=="Weak","territory_num"] <- "1"
phylo_data[phylo_data$territory_num=="Strong","territory_num"] <- "2"
phylo_data$territory_num %<>% as.numeric()

# Group weak and none together for binary.
phylo_data$territory_binary <- phylo_data$territory
phylo_data[phylo_data$territory_binary=="None","territory_binary"] <- "No territory"
phylo_data[phylo_data$territory_binary=="Weak","territory_binary"] <- "Territory"
phylo_data[phylo_data$territory_binary=="Strong","territory_binary"] <- "Territory"
phylo_data$territory_binary <- as.factor(phylo_data$territory_binary)

# Try year round territoriality.
phylo_data$territory_year_round <- phylo_data$territory
phylo_data[phylo_data$territory_year_round=="None","territory_year_round"] <- "No territory"
phylo_data[phylo_data$territory_year_round=="Weak","territory_year_round"] <- "No territory"
phylo_data[phylo_data$territory_year_round=="Strong","territory_year_round"] <- "Year Round Territory"
phylo_data$territory_year_round <- as.factor(phylo_data$territory_year_round)

### Habitat ###

# Count data.
phylo_data %>% dplyr::count(habitat)
phylo_data %<>% filter(habitat != "DD")

# Create a binary variable, grouping open and semi open together.
phylo_data$habitat_binary <- phylo_data$habitat
phylo_data[phylo_data$habitat_binary=="Semi-Open","habitat_binary"] <- "Dense"

phylo_data$habitat_num <- phylo_data$habitat
phylo_data[phylo_data$habitat_num=="Open","habitat_num"] <- "0"
phylo_data[phylo_data$habitat_num=="Semi-Open","habitat_num"] <- "1"
phylo_data[phylo_data$habitat_num=="Dense","habitat_num"] <- "2"
phylo_data$habitat_num %<>% as.numeric()


### Migration ###

# Count the data.
phylo_data %>% dplyr::count(migration)
phylo_data %<>% filter(migration != "DD")

# Create a binary migration variable, grouping none and weak/seasonal together.
phylo_data$migration_binary <- phylo_data$migration
phylo_data[phylo_data$migration_binary=="None","migration_binary"] <- "Weak"

phylo_data$migration_num <- phylo_data$migration
phylo_data[phylo_data$migration_num=="None","migration_num"] <- "0"
phylo_data[phylo_data$migration_num=="Weak","migration_num"] <- "1"
phylo_data[phylo_data$migration_num=="Strong","migration_num"] <- "2"
phylo_data$migration_num %<>% as.numeric()


### Trophic Level ###

# Count the data.
phylo_data %>% dplyr::count(trophic_level)

# Create a binary trophic factor, grouping into primary and secondary consumer.
phylo_data$trophic_binary <- phylo_data$trophic_level
phylo_data[phylo_data$trophic_level=="Carnivore","trophic_binary"] <- "Secondary"
phylo_data[phylo_data$trophic_level=="Omnivore","trophic_binary"] <- "Secondary"
phylo_data[phylo_data$trophic_level=="Scavenger","trophic_binary"] <- "Secondary"
phylo_data[phylo_data$trophic_level=="Herbivore","trophic_binary"] <- "Primary"

# Count the groupings.
phylo_data %>% dplyr::count(trophic_binary)


# Predatory
phylo_data$trophic_carn <- phylo_data$trophic_level
phylo_data[phylo_data$trophic_level=="Herbivore","trophic_carn"] <- "Non-Carnivore"
phylo_data[phylo_data$trophic_level=="Omnivore","trophic_carn"] <- "Non-Carnivore"
phylo_data[phylo_data$trophic_level=="Scavenger","trophic_carn"] <- "Non-Carnivore"
phylo_data$trophic_carn <- as.factor(phylo_data$trophic_carn)

# Trophic numbered.
phylo_data$trophic_num <- phylo_data$trophic_level
phylo_data[phylo_data$trophic_level=="Herbivore","trophic_num"] <- "0"
phylo_data[phylo_data$trophic_level=="Omnivore","trophic_num"] <- "1"
phylo_data[phylo_data$trophic_level=="Carnivore","trophic_num"] <- "2"
phylo_data[phylo_data$trophic_level=="Scavenger","trophic_num"] <- "2"
phylo_data$trophic_num %<>% as.numeric()


### Nest placement ###

# Count the data.
phylo_data %>% dplyr::count(nest_placement)

# Create an extra nest placement variable to test alongside nest predation.
phylo_data$placement_binary <- phylo_data$nest_placement
phylo_data[phylo_data$placement_binary=="Exposed Elevated","placement_binary"] <- "Exposed"
phylo_data[phylo_data$placement_binary=="Exposed Ground","placement_binary"] <- "Exposed"
phylo_data[phylo_data$placement_binary=="Brood Parasite","placement_binary"] <- "Concealed"
phylo_data[phylo_data$placement_binary=="Cavity","placement_binary"] <- "Concealed"

# Count the groupings.
phylo_data %>% dplyr::count(placement_binary)
phylo_data %>% filter_high() %>% dplyr::count(placement_binary)


### Incubation ###

# Count data.
phylo_data %>% dplyr::count(incubation_roles)

# Create a binary variable of incubation role to test alongside nest predation.
phylo_data$incubation_binary <- phylo_data$incubation_roles
phylo_data[phylo_data$incubation_binary=="Brood Parasite","incubation_binary"] <- "Both Sexes"
phylo_data[phylo_data$incubation_binary=="Megapode Mound","incubation_binary"] <- "Both Sexes"
phylo_data[phylo_data$incubation_binary=="Biparental","incubation_binary"] <- "Both Sexes"
phylo_data[phylo_data$incubation_binary=="Female Only","incubation_binary"] <- "Single Sex"
phylo_data[phylo_data$incubation_binary=="Male Only","incubation_binary"] <- "Single Sex"

# Count the groupings. # Still pretty even.
phylo_data %>% dplyr::count(incubation_binary)
phylo_data %>% filter_high() %>% dplyr::count(incubation_binary)



# ###############################################################################
#        #### Prepare response variables for dichromatism #####
# 
# 
# # Count dichro total.
# phylo_data %>% count(dichro_total)
# phylo_data %<>% filter(dichro_total != "DD")
# 
# # Create binary variables.
# phylo_data$head_bin <- phylo_data$head_dichro
# phylo_data$upper_bin <- phylo_data$upper_dichro
# phylo_data$under_bin <- phylo_data$under_dichro
# phylo_data$tail_bin <- phylo_data$tail_dichro
# phylo_data$wing_bin <- phylo_data$wing_dichro
# 
# # Count the scores. Can see very few intermediate.
# phylo_data %>% dplyr::count(head_bin)
# phylo_data %>% dplyr::count(upper_bin)
# phylo_data %>% dplyr::count(under_bin)
# phylo_data %>% dplyr::count(tail_bin)
# phylo_data %>% dplyr::count(wing_bin)
# 
# # Convert 2s to 1 for binary score.
# phylo_data[phylo_data$head_bin=="2","head_bin"] <- "1"
# phylo_data[phylo_data$upper_bin=="2","upper_bin"] <- "1"
# phylo_data[phylo_data$under_bin=="2","under_bin"] <- "1"
# phylo_data[phylo_data$tail_bin=="2","tail_bin"] <- "1"
# phylo_data[phylo_data$wing_bin=="2","wing_bin"] <- "1"
# 
# # Convert to numbers.
# phylo_data$head_bin %<>% as.numeric()
# phylo_data$upper_bin %<>% as.numeric()
# phylo_data$under_bin %<>% as.numeric()
# phylo_data$tail_bin %<>% as.numeric()
# phylo_data$wing_bin %<>% as.numeric()
# 
# # Look at the groupings.
# phylo_data %>% dplyr::count(head_bin)
# phylo_data %>% dplyr::count(upper_bin)
# phylo_data %>% dplyr::count(under_bin)
# phylo_data %>% dplyr::count(tail_bin)
# phylo_data %>% dplyr::count(wing_bin)
# 
# # For highest certainty.
# phylo_data %>% filter_high() %>% dplyr::count(head_bin)
# phylo_data %>% filter_high() %>% dplyr::count(upper_bin)
# phylo_data %>% filter_high() %>% dplyr::count(under_bin)
# phylo_data %>% filter_high() %>% dplyr::count(tail_bin)
# phylo_data %>% filter_high() %>% dplyr::count(wing_bin)
# 
# # Change dichro values to numeric.
# phylo_data$head_dichro %<>% as.numeric()
# phylo_data$upper_dichro %<>% as.numeric()
# phylo_data$under_dichro %<>% as.numeric()
# phylo_data$tail_dichro %<>% as.numeric()
# phylo_data$wing_dichro %<>% as.numeric()
# phylo_data$dichro_total %<>% as.numeric()
# 
# # Create a total dimorphism binary variable.
# phylo_data$dichro_bin <- as.numeric(phylo_data$dichro_total)
# phylo_data$dichro_bin[phylo_data$dichro_bin > 0] <- 1
# 
# # Create hidden and visible response variables.
# phylo_data %<>% mutate(
#   dichro_vis = upper_dichro + head_dichro + tail_dichro,
#   dichro_hid = under_dichro + wing_dichro,
#   
#   dichro_vis_bin = upper_bin + head_bin + tail_bin,
#   dichro_hid_bin = under_bin + wing_bin
# )
# 
# # Rescale as binary.
# phylo_data$dichro_vis_bin[phylo_data$dichro_vis_bin > 0] <- 1
# phylo_data$dichro_hid_bin[phylo_data$dichro_hid_bin > 0] <- 1
# 



###############################################################################
                 ##### Transform continuous variables #####


# Plot the densities of the continous variables.
(mass_plot <- ggplot(phylo_data, aes(x = mass)) + geom_density(cex=1) + theme_classic())
(lat_plot <- ggplot(phylo_data, aes(x = abs(centroid_latitude))) + geom_density(cex=1) + theme_classic())
(range_plot <- ggplot(phylo_data, aes(x = range_size)) + geom_density(cex=1) + theme_classic())

# Transform variables by taking logs.
phylo_data$body_mass_log <- log(phylo_data$mass)
phylo_data$centroid_sqrt <- sqrt(abs(phylo_data$centroid_latitude))
phylo_data$range_log <- log(phylo_data$range_size)

# Look at the new plots.
(log_mass_plot <- ggplot(phylo_data, aes(x = body_mass_log)) + geom_density(cex=1) + theme_classic())
(sqrt_lat_plot <- ggplot(phylo_data, aes(x = centroid_sqrt)) + geom_density(cex=1) + theme_classic())
(log_range_plot <- ggplot(phylo_data, aes(x = range_log)) + geom_density(cex=1) + theme_classic())


# Check together
ggarrange(mass_plot, log_mass_plot)
ggarrange(lat_plot, sqrt_lat_plot)
ggarrange(range_plot, log_range_plot)


###############################################################################
                 ##### Plot the predictors #####

# Quick function wrap for ggsave.
save_plot <- function(pathway){
  ggsave(pathway, width = 5, height = 5)
}

# # Plot incubation roles.
# sex_meanplot("incubation_roles")
# save_plot("Plots/Predictors/incubation_roles.tiff")
# 
# # Plot nest placement.
# sex_meanplot("nest_placement")
# save_plot("Plots/Predictors/nest_placement.tiff")

# Territory.
sex_meanplot("territory")
save_plot("Plots/Predictors/territory.tiff")

# Habitat
sex_meanplot("habitat")
save_plot("Plots/Predictors/habitat.tiff")

# Binary grouping open and semi open
sex_meanplot("habitat_binary")
save_plot("Plots/Predictors/habitat_bi.tiff")

# Stability.
sex_meanplot("stability")
save_plot("Plots/Predictors/stability.tiff")

# Migration.
sex_meanplot("migration")
save_plot("Plots/Predictors/migration.tiff")

# Migration Binary.
sex_meanplot("migration_binary")
save_plot("Plots/Predictors/migration_binary.tiff")

# Cooperative.
sex_meanplot("cooperative")
save_plot("Plots/Predictors/cooperative.tiff")

# Social selection.
sex_meanplot("social_selection")
save_plot("Plots/Predictors/social_selection.tiff")

# Diet.
sex_meanplot("trophic_level")
save_plot("Plots/Predictors/trophic_level.tiff")

# Diet.
sex_meanplot("trophic_binary")
save_plot("Plots/Predictors/trophic_binary.tiff")

# Body mass.
ggplot(phylo_data, aes(y=body_mass_log, x=sexual_score, group = sexual_score)) + 
    geom_boxplot(alpha=0.1) + ylab("Body Mass (Logged)") + theme_classic()
save_plot("Plots/Predictors/mass_log.tiff")

ggplot(phylo_data, aes(y=body_mass_log, x=sexual_binary, group = sexual_binary)) + 
  geom_boxplot(alpha=0.1) + ylab("Body Mass (Logged)") + theme_classic()
save_plot("Plots/Predictors/mass_log.tiff")

# Centroid.
ggplot(phylo_data, aes(y=centroid_sqrt, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Latitude (Sqrt)") + theme_classic()
save_plot("Plots/Predictors/centroid.tiff")


###############################################################################
##### Plot bio variables #####

phylo_data <- left_join(phylo_data, bio_data)


# Centroid.
ggplot(phylo_data, aes(y=npp, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Latitude (Sqrt)") + theme_classic()


###############################################################################
                 ##### Export data ######


# Export the data.
write.csv(phylo_data, "Data/sexual_traits.csv", row.names = FALSE)


