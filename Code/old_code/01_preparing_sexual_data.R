###############################################################################
                    ##### Preparing trait data #####
###############################################################################

# This script pre-processes the trait data so it can be used for phylogenetic
# models. It also reads in the generation length and bioclim data.

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
bird_traits <- read.csv("Data/ecotraits_09_05_without_invalids.csv") %>% 
  clean_names()

# Select the traits we want to model.
phylo_data <- bird_traits %>% dplyr::select(
  # Taxonomy.
  birdtree_name, order, family,
  # Sexual selection.                                   
  sexual_score, sexual_certainty,
  # Social Selection.
  territory, lifehistory_uncertainty, 
  # Extra hypotheses
  migration, trophic_level, trophic_niche,
  # Extra controls.
  mass, complete_latitude, gen_length, npp, bio4)

# Check the data. Still missing traits. Will drop these species now anyway.
skim(phylo_data)

# Remove NA species.
phylo_data %<>% tidyr::drop_na(bio4)

# Read in developmental data.
developmental_data <- read.csv("Data/developmental_data.csv") %>% 
  clean_names() %>% dplyr::select(species_cooney_taxonomy, chick_pc1, chick_pc_ascores, hatchling_pc1, hatchling_pc_ascores)
colnames(developmental_data)[1] <- "birdtree_name"
developmental_data$birdtree_name %<>% gsub(pattern = "_", replacement =  " ")

# Read in the population density data.
density_data <- read.csv("Data/population_density.csv") %>% clean_names()
colnames(density_data)[3] <- "birdtree_name"
density_data %<>% dplyr::select(-c(order, family))
density_data$birdtree_name %<>% gsub(pattern = "_", replacement =  " ")

# Read in the greisser data.
greisser_data <- read.csv("Data/greisser_data.csv") %>% clean_names()
colnames(greisser_data)[1] <- "birdtree_name"
greisser_data$birdtree_name %<>% gsub(pattern = "_", replacement =  " ")
greisser_data %<>% dplyr::select(birdtree_name, devo_mode, devo_mode_full, time_fed)

# Read in the wang data.
wang_data <- read.csv("Data/Developmental/wang_and_kimball.csv") %>% clean_names() %>% dplyr::select(-c(cb))
wang_data$birdtree_name %<>% gsub(pattern = "_", replacement =  " ")
colnames(wang_data)[2] <- "devo_mode_wang"
wang_data[wang_data$devo_mode_wang == "1","devo_mode_wang"] <- "altricial"
wang_data[wang_data$devo_mode_wang == "0","devo_mode_wang"] <- "precocial"

# Read in the cooney data.
cooney_data <- read.csv("Data/cooney_developmental_data.csv") %>% clean_names() %>% dplyr::select(binomial, dev_mod) %>% tidyr::drop_na(dev_mod)
colnames(cooney_data) <- c("birdtree_name", "devo_mode_cooney")
cooney_data$birdtree_name %<>% gsub(pattern = "_", replacement =  " ")
cooney_data$devo_binary_cooney <- cooney_data$devo_mode_cooney
cooney_data$devo_binary_cooney[cooney_data$devo_binary_cooney == "semi-precocial"] <- "precocial"

# Read in the remes data.
remes_data <- read.csv("Data/Developmental/remes_pnas_data.csv") %>% clean_names() %>% dplyr::select(species_name, chick_dev)
colnames(remes_data) <- c("birdtree_name", "devo_mode_remes")
remes_data$birdtree_name %<>% gsub(pattern = "_", replacement =  " ")
remes_data$devo_mode_remes %<>% as.character()
remes_data$devo_mode_remes[remes_data$devo_mode_remes == "1"] <- "altricial"
remes_data$devo_mode_remes[remes_data$devo_mode_remes == "2"] <- "precocial"


# Add the developmental and density data.
phylo_data %<>% left_join(developmental_data) %>% left_join(density_data) %>% left_join(greisser_data) %>% left_join(wang_data) %>% left_join(cooney_data) %>% left_join(remes_data)

# Check extra species that are missing.
skim(phylo_data)
#phylo_data %<>% na.omit()

# Check correlations between developmental datasets.
dev_corr <- phylo_data %>% dplyr::select(birdtree_name, order, devo_mode, devo_mode_full, devo_mode_wang, devo_binary_cooney, devo_mode_remes)

dev_corr$devo_mode %<>% as.factor() %>% as.numeric()
dev_corr$devo_binary_cooney %<>% as.factor() %>% as.numeric()
dev_corr$devo_mode_wang %<>% as.factor() %>% as.numeric()
dev_corr$devo_mode_remes %<>% as.factor() %>% as.numeric()
GGally::ggcorr(dev_corr, label = TRUE)

# Change these families to be correct altrcial / precocial families.
phylo_data$wang_edited <- phylo_data$devo_mode_wang
phylo_data$wang_edited[phylo_data$family %in% c("Ciconiidae", "Ardeidae", "Threskiornithidae", "Musophagidae", "Dromadidae", "Eurypygidae", "Heliornithidae")] <- "altricial"
phylo_data$wang_edited[phylo_data$family %in% c("Alcidae", "Laridae", "Caprimulgidae", "Stercorariidae","Diomedeidae", "Gaviidae")] <- "precocial"


sex_meanplot("devo_mode")
sex_meanplot("devo_binary_cooney")
sex_meanplot("devo_mode_wang")
sex_meanplot("wang_edited")


###############################################################################
                #### Binary model variables #####


# Group weak and none together for binary territory.
phylo_data$territory_binary <- phylo_data$territory
phylo_data[phylo_data$territory_binary=="None","territory_binary"] <- "No territory"
phylo_data[phylo_data$territory_binary=="Weak","territory_binary"] <- "Territory"
phylo_data[phylo_data$territory_binary=="Strong","territory_binary"] <- "Territory"
phylo_data$territory_binary <- as.factor(phylo_data$territory_binary)

# Create a binary migration variable, grouping none and weak/seasonal together.
phylo_data$migration_binary <- phylo_data$migration
phylo_data[phylo_data$migration_binary=="None","migration_binary"] <- "Weak"

# Create a binary trophic factor, grouping into primary and secondary consumer.
phylo_data$trophic_binary <- phylo_data$trophic_level
phylo_data[phylo_data$trophic_level=="Carnivore","trophic_binary"] <- "Secondary"
phylo_data[phylo_data$trophic_level=="Omnivore","trophic_binary"] <- "Secondary"
phylo_data[phylo_data$trophic_level=="Scavenger","trophic_binary"] <- "Secondary"
phylo_data[phylo_data$trophic_level=="Herbivore","trophic_binary"] <- "Primary"

# Create a binary trophic factor, grouping into primary and secondary consumer.
phylo_data[is.na(phylo_data$trophic_niche), "trophic_niche"] <- "Omnivore"
phylo_data$diet_binary <- NA
phylo_data[phylo_data$trophic_niche=="Nectarivore","diet_binary"] <- "Frug-nect"
phylo_data[phylo_data$trophic_niche=="Frugivore","diet_binary"] <- "Frug-nect"
phylo_data[is.na(phylo_data$diet_binary), "diet_binary"] <- "Non-frug-nect"

# Create binary sexual selection score for conservative analysis.
phylo_data$sexual_binary <- phylo_data$sexual_score
phylo_data[phylo_data$sexual_binary< 3,"sexual_binary"] <- 0
phylo_data[phylo_data$sexual_binary> 0,"sexual_binary"] <- 1

###############################################################################
              ##### Transform continuous variables #####


# Plot the densities of the continous variables.
(mass_plot <- ggplot(phylo_data, aes(x = mass)) + geom_density(cex=1) + theme_classic())
(lat_plot <- ggplot(phylo_data, aes(x = abs(complete_latitude))) + geom_density(cex=1) + theme_classic())
(temp_plot <- ggplot(phylo_data, aes(x = bio4)) + geom_density(cex=1) + theme_classic())
#(rain_plot <- ggplot(phylo_data, aes(x = bio15)) + geom_density(cex=1) + theme_classic())
(npp_plot <- ggplot(phylo_data, aes(x = npp)) + geom_density(cex=1) + theme_classic())
(gen_plot <- ggplot(phylo_data, aes(x = gen_length)) + geom_density(cex=1) + theme_classic())
#(survival_plot <- ggplot(phylo_data, aes(x = adult_survival)) + geom_density(cex=1) + theme_classic())
#(breeding_plot <- ggplot(phylo_data, aes(x = age_at_first_breeding)) + geom_density(cex=1) + theme_classic())
#(longevity_plot <- ggplot(phylo_data, aes(x = maximum_longevity)) + geom_density(cex=1) + theme_classic())
(chick_plot <- ggplot(phylo_data %>% filter(chick_pc_ascores == "RealData"), aes(x = chick_pc1)) + geom_density(cex=1) + theme_classic())
(hatch_plot <- ggplot(phylo_data %>% filter(hatchling_pc_ascores == "RealData"), aes(x = hatchling_pc1)) + geom_density(cex=1) + theme_classic())
(hatch_plot <- ggplot(phylo_data, aes(x = hatchling_pc1)) + geom_density(cex=1) + theme_classic())
(dens_plot <- ggplot(phylo_data, aes(x = predicted_density)) + geom_density(cex=1) + theme_classic())
(fed_plot <- ggplot(phylo_data, aes(x = time_fed)) + geom_density(cex=1) + theme_classic())



# Transform variables by taking logs.
phylo_data$body_mass_log <- log(phylo_data$mass)
phylo_data$centroid_sqrt <- sqrt(abs(phylo_data$complete_latitude))
phylo_data$dens_sqrt <- sqrt(phylo_data$predicted_density)
#phylo_data$range_log <- log(phylo_data$range_size)
phylo_data$temp_log <- log(phylo_data$bio4)

# This transformation works but it's a bit sketchy. Should assess both.
phylo_data$npp_sqrt <- sqrt(max(phylo_data$npp, na.rm = TRUE) - phylo_data$npp)*-1
phylo_data$gen_log <- log(phylo_data$gen_length)

phylo_data$fed_sqrt <- sqrt(phylo_data$time_fed)



# Regress seasonality against productivity.
seasonal_linear <- lm(phylo_data$temp_log ~ phylo_data$npp_sqrt)
plot(seasonal_linear)

# Doesn't work when you don't remove NAs.
#phylo_data$new_seasonal <- seasonal_linear$residuals
# cor(phylo_data$new_seasonal, phylo_data$temp_log)
# cor(phylo_data$new_seasonal, phylo_data$npp_sqrt)
# hist(phylo_data$new_seasonal)

#phylo_data$breed_log <- log(phylo_data$age_at_first_breeding)
#phylo_data$long_log <- log(phylo_data$maximum_longevity)

# Look at the new plots.
(log_mass_plot <- ggplot(phylo_data, aes(x = body_mass_log)) + geom_density(cex=1) + theme_classic())
(sqrt_lat_plot <- ggplot(phylo_data, aes(x = centroid_sqrt)) + geom_density(cex=1) + theme_classic())
#(log_range_plot <- ggplot(phylo_data, aes(x = range_log)) + geom_density(cex=1) + theme_classic())
(log_temp_plot <- ggplot(phylo_data, aes(x = temp_log)) + geom_density(cex=1) + theme_classic())
(sqrt_npp_plot <- ggplot(phylo_data, aes(x = npp_sqrt)) + geom_density(cex=1) + theme_classic())
(log_gen_plot <- ggplot(phylo_data, aes(x = gen_log)) + geom_density(cex=1) + theme_classic())
#(log_breeding_plot <- ggplot(phylo_data, aes(x = breed_log)) + geom_density(cex=1) + theme_classic())
#(log_longevity_plot <- ggplot(phylo_data, aes(x = long_log)) + geom_density(cex=1) + theme_classic())
(sqrt_dens_plot <- ggplot(phylo_data, aes(x = dens_sqrt)) + geom_density(cex=1) + theme_classic())
(sqrt_fed_plot <- ggplot(phylo_data, aes(x = fed_sqrt)) + geom_density(cex=1) + theme_classic())


# This data is pretty horrible so i guess don't transform.

phylo_data$hatch_sqrt <- sqrt(phylo_data$hatchling_pc1 + abs(min(phylo_data$hatchling_pc1)))
phylo_data$chick_sqrt <- sqrt(phylo_data$chick_pc1 + abs(min(phylo_data$chick_pc1)))
(log_hatch_plot <- ggplot(phylo_data, aes(x = hatch_sqrt)) + geom_density(cex=1) + theme_classic())
(sqrt_chick_plot <- ggplot(phylo_data, aes(x = chick_sqrt)) + geom_density(cex=1) + theme_classic())


# Check together
ggarrange(mass_plot, log_mass_plot)
ggarrange(lat_plot, sqrt_lat_plot)
#ggarrange(range_plot, log_range_plot)
ggarrange(temp_plot, log_temp_plot)
ggarrange(npp_plot, sqrt_npp_plot)
ggarrange(gen_plot, log_gen_plot)
#ggarrange(breeding_plot, log_breeding_plot)
#ggarrange(longevity_plot, log_longevity_plot)
ggarrange(dens_plot, sqrt_dens_plot)
ggarrange(chick_plot, sqrt_chick_plot)
cor(phylo_data$npp, phylo_data$npp_sqrt)



###############################################################################
            ##### Create variables used in plotting maps #####

# Grouping 0-2 as just monogamous.
phylo_data$sexual_tri <- phylo_data$sexual_score 
phylo_data$sexual_tri[phylo_data$sexual_tri < 3] <- 0
phylo_data$sexual_tri[phylo_data$sexual_tri == 3] <- 1
phylo_data$sexual_tri[phylo_data$sexual_tri == 4] <- 2

# Binary certainty map. # Don't really need now.
phylo_data$cert_dummy <- phylo_data$sexual_certainty
phylo_data$cert_dummy[phylo_data$cert_dummy > 1] <- 0

# Combine certainty for 1 and 2.
phylo_data$cert_dummy_2 <- phylo_data$sexual_certainty
phylo_data$cert_dummy_2[phylo_data$cert_dummy_2 > 2] <- 0
phylo_data$cert_dummy_2[phylo_data$cert_dummy_2 > 0] <- 1

# Make a reverse certainty (actually very useful.)
phylo_data$cert_reverse <- phylo_data$sexual_certainty
phylo_data$cert_reverse <- phylo_data$cert_reverse*-1
phylo_data$cert_reverse <- phylo_data$cert_reverse+5

# Make dummy territoriality metrics.
phylo_data$terr_dummy <- phylo_data$territory
phylo_data$terr_dummy[phylo_data$terr_dummy == "None"] <- "0"
phylo_data$terr_dummy[phylo_data$terr_dummy == "Weak"] <- "1"
phylo_data$terr_dummy[phylo_data$terr_dummy == "Strong"] <- "1"
phylo_data$terr_dummy %<>% as.numeric()

phylo_data$year_terr_dummy <- phylo_data$territory
phylo_data$year_terr_dummy[phylo_data$year_terr_dummy == "None"] <- "0"
phylo_data$year_terr_dummy[phylo_data$year_terr_dummy == "Weak"] <- "0"
phylo_data$year_terr_dummy[phylo_data$year_terr_dummy == "Strong"] <- "1"
phylo_data$year_terr_dummy %<>% as.numeric()

# Create a dummy herbivore variable to group by.
herbivores <- c("Herbivore aquatic", "Herbivore terrestrial")
phylo_data$herbivore_niche <- NA
phylo_data$herbivore_niche[phylo_data$trophic_niche %in% herbivores] <- "Herbivore"
phylo_data$herbivore_niche[is.na(phylo_data$herbivore_niche)] <- "Non-Herbivore"

###############################################################################
                 ##### Plot the predictors #####

# Quick function wrap for ggsave.
save_plot <- function(pathway){
  ggsave(pathway, width = 5, height = 5)
}

# Territory.
sex_meanplot("territory")
save_plot("Plots/Predictors/territory.tiff")

grouped_data <- phylo_data %>% 
  group_by(territory, diet_binary) %>% 
  summarise(trait = first(territory),
            sex_mean = mean(sexual_score),
            sex_sd = sd(sexual_score),
            sex_se = sd(sexual_score)/sqrt(length(sexual_score)))

ggplot(grouped_data, aes(x = trait, y = sex_mean, col = diet_binary, group = diet_binary)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se*1.96, ymax = sex_mean + sex_se*1.96), position = position_dodge(width = 1)) + 
  geom_point(position = position_dodge(width = 1)) + ylab("Sexual Selection") + xlab("Territory") + theme_classic()


# Migration.
sex_meanplot("migration")
save_plot("Plots/Predictors/migration.tiff")

# Migration Binary.
sex_meanplot("migration_binary")
save_plot("Plots/Predictors/migration_binary.tiff")

# Diet.
sex_meanplot("trophic_level")
save_plot("Plots/Predictors/trophic_level.tiff")
# Diet.
sex_meanplot("trophic_binary")
save_plot("Plots/Predictors/trophic_binary.tiff")

sex_meanplot("diet_binary")
save_plot("Plots/Predictors/diet_binary.tiff")


# Developmental mode (Greisser)
sex_meanplot("devo_mode")

phylo_data %>% tidyr::drop_na(devo_mode) %>% ggplot(aes(y=body_mass_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Body Mass (Logged)") + theme_classic() + facet_grid(~ devo_mode)
save_plot("Plots/Predictors/developmental_mode_body_mass.tiff")

phylo_data %>% tidyr::drop_na(devo_mode) %>% ggplot(aes(y=gen_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Generation length (Logged)") + theme_classic() + facet_grid(~ devo_mode)

# Developmental mode (Wang)
sex_meanplot("devo_mode_wang")
save_plot("Plots/Predictors/devo_mode_wang.tiff")

phylo_data %>% ggplot(aes(y=body_mass_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + geom_smooth(method = "lm") + ylab("Body Mass (Logged)") + theme_classic() + facet_grid(~ devo_mode_wang)
save_plot("Plots/Predictors/devo_wang_mode_body_mass.tiff")

phylo_data %>% ggplot(aes(y=gen_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + geom_smooth(method = "lm") + ylab("Generation Length (Logged)") + theme_classic() + facet_grid(~ devo_mode_wang)

# Time fed.
ggplot(phylo_data, aes(y=fed_sqrt, x=gen_log, col = sexual_score)) + 
  geom_point(alpha=0.1) + ylab("Time Fed") + theme_classic()

ggplot(phylo_data, aes(y=sqrt(time_fed), x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Time Fed") + theme_classic()

# Body mass.
ggplot(phylo_data, aes(y=body_mass_log, x=sexual_score, group = sexual_score)) + 
    geom_boxplot(alpha=0.1) + ylab("Body Mass (Logged)") + theme_classic() + facet_grid(~ devo_mode)
save_plot("Plots/Predictors/mass_log.tiff")

# Centroid.
ggplot(phylo_data, aes(y=centroid_sqrt, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Latitude (Sqrt)") + theme_classic()
save_plot("Plots/Predictors/centroid.tiff")

ggplot(phylo_data, aes(y=temp_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Temperature Stability") + theme_classic()

ggplot(phylo_data, aes(y=bio15, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Rainfall Stability") + theme_classic()

ggplot(phylo_data, aes(y=npp_sqrt, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Net Primary Productivity") + theme_classic()

ggplot(phylo_data, aes(y=gen_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Generation Length") + theme_classic()

ggplot(phylo_data, aes(y=adult_survival, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Adult Survival") + theme_classic()

ggplot(phylo_data, aes(y=breed_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Age at First Breeding") + theme_classic()

ggplot(phylo_data, aes(y=long_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Max Longevity") + theme_classic()

# Chick PC1
ggplot(phylo_data, aes(y=chick_pc1, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Chick PC1") + theme_classic()

ggplot(phylo_data, aes(y=chick_sqrt, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Chick PC1") + theme_classic()

ggplot(phylo_data, aes(y=chick_pc1, x=sexual_binary, group = sexual_binary)) + 
  geom_boxplot(alpha=0.1) + ylab("Chick PC1") + theme_classic()

ggplot(phylo_data, aes(y=chick_sqrt, x=sexual_binary, group = sexual_binary)) + 
  geom_boxplot(alpha=0.1) + ylab("Chick PC1") + theme_classic()


# Try binning 
seas_range <- seq(from = min(phylo_data$chick_pc1), to = max(phylo_data$chick_pc1), length.out = 20)
seas_labels <- seas_range[-1]
phylo_data$binned_chick <- cut(phylo_data$chick_pc1, breaks=seas_range, labels = seas_labels)

# Group the data by trophic level and bins.
chick_data <- phylo_data %>% 
  group_by(binned_chick) %>% 
  summarise(trait = first(binned_chick),
            sex_mean = mean(sexual_binary),
            sex_sd = sd(sexual_binary),
            sex_se = sd(sexual_binary)/sqrt(length(sexual_binary))) %>% na.omit()

# Plot chick developmental data.
chick_data %>% 
  ggplot(aes(x = trait, y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se)
                , show.legend = FALSE) + 
  geom_point() + 
  geom_smooth(se = FALSE, method = "gam") + 
  ylab("Sexual selection") + 
  xlab("Developmental mode") + theme_classic() + 
  theme(legend.position = c(0.1,0.9),
        text = element_text(face = "bold"))
chick_plot

plot(jitter(phylo_data$sexual_score) ~ phylo_data$chick_sqrt)

lm(phylo_data$sexual_score ~ phylo_data$chick_sqrt) %>% summary()

ggplot(phylo_data, aes(y=hatchling_pc1, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Hatch PC1") + theme_classic()

# Population density.
ggplot(phylo_data, aes(y=dens_sqrt, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("Population Density") + theme_classic()


ggplot(phylo_data, aes(y=new_seasonal, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("New Seasonal") + theme_classic()

ggplot(phylo_data, aes(y=temp_log, x=sexual_score, group = sexual_score)) + 
  geom_boxplot(alpha=0.1) + ylab("New Seasonal") + theme_classic()


###############################################################################
                 ##### Export data ######


# Export the data.
write.csv(phylo_data, "Data/sexual_traits.csv", row.names = FALSE)





################################################################################
                            #### END #####
################################################################################


# Showing territoriality interactions.
  grouped_data <- phylo_data %>% 
    group_by(territory_binary, trophic_binary) %>% 
    summarise(sex_mean = mean(sexual_score),
              sex_sd = sd(sexual_score),
              sex_se = sd(sexual_score)/sqrt(length(sexual_score)))
  
  grouped_data$label <- c("No territory x ")

  ggplot(grouped_data, aes(x = trophic_binary , y = sex_mean, 
                           col = territory_binary, group = territory_binary)) +
    geom_errorbar(aes(ymin = sex_mean - sex_se*1.96, ymax = sex_mean + sex_se*1.96), 
                  position = position_dodge(width = 1), linewidth = 1, width = 0.2) + 
    geom_point(position = position_dodge(width = 1)) + 
    ylab("Sexual Selection") + xlab("Territory x Diet") + theme_classic()

  
  
  
  grouped_data <- phylo_data %>% 
    group_by(territory_binary, migration_binary) %>% 
    summarise(sex_mean = mean(sexual_score),
              sex_sd = sd(sexual_score),
              sex_se = sd(sexual_score)/sqrt(length(sexual_score)))
  
  ggplot(grouped_data, aes(x = migration_binary , y = sex_mean, 
                           col = territory_binary, group = territory_binary)) +
    geom_errorbar(aes(ymin = sex_mean - sex_se*1.96, ymax = sex_mean + sex_se*1.96), 
                  position = position_dodge(width = 1), linewidth = 1, width = 0.2) + 
    geom_point(position = position_dodge(width = 1)) + 
    ylab("Sexual Selection") + xlab("Territory x Migration") + theme_classic()
  