###############################################################################
                   # Making sexual selection maps #
###############################################################################

# This script makes maps of sexual selection scores, accounting for diet.

# Load packages. 
library(sf)
library(fasterize)
library(magrittr)
library(dplyr)
library(raster)
library(tidyr)
library(ggplot2)
library(stringr)
library(maptools)
library(ggpubr)

# Clear environment. 
rm(list=ls())

# Functions.
source("Code/functions.R")

################################################################################
                         #### Data ####


# Load in the range maps.
load("../../Phd_data/Spatial/clean_jetz_ranges.RData")

# Load in the trait data.
phylo_data <- read.csv("Data/sexual_traits.csv")

# Join up the data.
sex_ranges <- left_join(phylo_data, new_jetz_ranges,  by = c("birdtree_name" = "SCINAME"))
sex_ranges <- st_as_sf(sex_ranges)
rm(new_jetz_ranges)



################################################################################
                #### Species richness raster ####


# Start by creating an empty raster stack to store our data in.
raster_template <- raster(ncols=8640, nrows=3600, ymn = -60)
#raster_template <- raster(ncols=17280, nrows=7200, ymn = -60)

# Create a map of species richness by summing overlapping ranges.
species_raster <- fasterize(sex_ranges, raster_template, fun = "sum")

# Mask for land.
data(wrld_simpl)
land <- wrld_simpl %>% st_as_sf() %>% fasterize(raster_template)
species_raster <- mask(species_raster, land)

# Remove less than 10 species.
species_raster[which(getValues(species_raster) < 10)] <- NA

# Plot the new map.
plot(species_raster, col=heat.colors(50))



################################################################################
                  #### Make rasters of life history traits ####


# Create sexual selection ranges with trait data.
sex_raster <- average_raster("sexual_score")

# Plot the raster.
plot(sex_raster)

# Make sexual selection certainty dummy score.
sex_ranges$cert_dummy <- sex_ranges$sexual_certainty
sex_ranges$cert_dummy[sex_ranges$cert_dummy > 1] <- 0
cert_raster <- average_raster("cert_dummy")

# Make inverse data certainty score for plotting. 
sex_ranges$cert_inverse <- 1/sex_ranges$sexual_certainty

# Make another certainty map averaging all certainty.
cert_raster_2 <- average_raster("cert_inverse")

# Combine certainty for 1 and 2.
sex_ranges$cert_dummy_2 <- sex_ranges$sexual_certainty
sex_ranges$cert_dummy_2[sex_ranges$cert_dummy_2 > 2] <- 0
sex_ranges$cert_dummy_2[sex_ranges$cert_dummy_2 > 0] <- 1
cert_raster_3 <- average_raster("cert_dummy_2")

# Make dummy territoriality maps.
sex_ranges$terr_dummy <- sex_ranges$territory
sex_ranges$terr_dummy[sex_ranges$terr_dummy == "None"] <- "0"
sex_ranges$terr_dummy[sex_ranges$terr_dummy == "Weak"] <- "1"
sex_ranges$terr_dummy[sex_ranges$terr_dummy == "Strong"] <- "1"
sex_ranges$terr_dummy %<>% as.numeric()

sex_ranges$year_terr_dummy <- sex_ranges$territory
sex_ranges$year_terr_dummy[sex_ranges$year_terr_dummy == "None"] <- "0"
sex_ranges$year_terr_dummy[sex_ranges$year_terr_dummy == "Weak"] <- "0"
sex_ranges$year_terr_dummy[sex_ranges$year_terr_dummy == "Strong"] <- "1"
sex_ranges$year_terr_dummy %<>% as.numeric()


################################################################################
                    #### Prepare climatic maps ####

# Read in bioclim data.
temp_seas <- raster("Data/Chelsa/CHELSA_bio4_1981-2010_V.2.1.tif")
npp <- raster("Data/Chelsa/CHELSA_npp_1981-2010_V.2.1.tif")

library(terra)
# Resample so they're same extent / res
temp_seas <- resample(rast(temp_seas), rast(species_raster))
npp <- resample(rast(npp), rast(species_raster))

# Mask.
temp_seas <- mask(raster(temp_seas), species_raster)
npp <- mask(raster(npp), species_raster)


################################################################################
         ###### Maps of primary vs secondary consumers #########

# Filter data for primary and secondary consumers.
primary_ranges <- sex_ranges %>% filter(trophic_binary == "Primary")
secondary_ranges <- sex_ranges %>% filter(trophic_binary == "Secondary")

# Create species richness rasters.
p_species_raster <- fasterize(primary_ranges, raster_template, fun = "sum")
s_species_raster <- fasterize(secondary_ranges, raster_template, fun = "sum")

# NA areas with less than 10 species.
p_species_raster[which(getValues(p_species_raster) < 10)] <- NA
s_species_raster[which(getValues(s_species_raster) < 10)] <- NA

# Mask for land.
p_species_raster <- mask(p_species_raster, land)
s_species_raster <- mask(s_species_raster, land)

# Create average rasters.
p_sex_raster <- average_raster("sexual_score", primary_ranges, p_species_raster)
s_sex_raster <- average_raster("sexual_score", secondary_ranges, s_species_raster)

# For territoriality as well.
p_terr_raster <- average_raster("terr_dummy", primary_ranges, p_species_raster)
s_terr_raster <- average_raster("terr_dummy", secondary_ranges, s_species_raster)

p_year_terr_raster <- average_raster("year_terr_dummy", primary_ranges, p_species_raster)
s_year_terr_raster <- average_raster("year_terr_dummy", secondary_ranges, s_species_raster)


#rm(sex_ranges)
# Delete secondary ranges now.
rm(primary_ranges)
rm(secondary_ranges)

################################################################################
            ###### Maps of high certainty species #########

# Filter for high certainty species.
high_ranges <- sex_ranges %>% filter(sexual_certainty == 1)
med_ranges <- sex_ranges %>% filter(sexual_certainty < 3)
low_ranges <- sex_ranges %>% filter(sexual_certainty > 2)

# Create species richness rasters.
h_species_raster <- fasterize(high_ranges, raster_template, fun = "sum")
m_species_raster <- fasterize(med_ranges, raster_template, fun = "sum")
l_species_raster <- fasterize(low_ranges, raster_template, fun = "sum")

# NA areas with less than 10 species.
h_species_raster[which(getValues(h_species_raster) < 10)] <- NA
m_species_raster[which(getValues(m_species_raster) < 10)] <- NA
l_species_raster[which(getValues(l_species_raster) < 10)] <- NA

# Mask for land.
h_species_raster <- mask(h_species_raster, land)
m_species_raster <- mask(m_species_raster, land)
l_species_raster <- mask(l_species_raster, land)

# Create average rasters.
h_sex_raster <- average_raster("sexual_score", high_ranges, h_species_raster)
m_sex_raster <- average_raster("sexual_score", med_ranges, m_species_raster)
l_sex_raster <- average_raster("sexual_score", low_ranges, l_species_raster)

# Delete sf dataframes now to free up memory.
rm(sex_ranges)
rm(high_ranges)
rm(med_ranges)
rm(low_ranges)
gc()


###############################################################################
                   #### Plot maps with ggplot ####


# Create data frame of land to use to crop ggplot maps.
land_data <- as.data.frame(land, xy=TRUE)

# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

## Sexual selection and data certainty ##

# Full sexual selection map.
sex_plot <- ggplot_ylim_raster(sex_raster, 6, "Sexual selection")
ggsave("Plots/Maps/sexual_selection_map.png", height = 7, width = 15, dpi = 1000)

# Dichotomised data certainty map.
cert_plot <- ggplot_ylim_raster(cert_raster, 6, "Data certainty")
ggsave("Plots/Maps/sexual_certainty_map.png", height = 7, width = 15, dpi = 1000)

# Put them together.
ggarrange(sex_plot, cert_plot, labels = c("a", "b"), nrow = 2, font.label = list(size = 24))
ggsave("Plots/Maps/scores_and_cert_plot.png", height = 13, width = 15, dpi = 1000)


## Various data certainties ##

# Try summing all certainty scores instead of dichotomised.
cert_plot_2 <- ggplot_ylim_raster(cert_raster_2, 6, "Data certainty")# +
ggsave("Plots/Maps/one_and_two_certainty_map.png", height = 7, width = 15, dpi = 1000)

# Dichotomised data certainty map.
cert_plot_3 <- ggplot_ylim_raster(cert_raster_3, 6, "Data certainty")
ggsave("Plots/Maps/sexual_certainty_map.png", height = 7, width = 15, dpi = 1000)

# Show high certainty species richness.
h_cert_sr_plot <- ggplot_ylim_raster(h_species_raster, 6, "Species Richness")# +
ggsave("Plots/Maps/high_cert_SR_map.png", height = 7, width = 15, dpi = 1000)

# Show high certainty species richness.
h_sex_plot <- ggplot_ylim_raster(h_sex_raster, 6, "Sexual selection")# +
ggsave("Plots/Maps/high_cert_sex_map.png", height = 7, width = 15, dpi = 1000)

m_sex_plot <- ggplot_ylim_raster(m_sex_raster, 6, "Sexual selection")# +

l_cert_sr_plot <- ggplot_ylim_raster(l_species_raster, 6, "Species Richness")# +
ggsave("Plots/Maps/low_cert_SR_map.png", height = 7, width = 15, dpi = 1000)


## Trophic level sexual selection ##

# Do primary and secondary consumers.
p_sex_plot <- ggplot_ylim_raster(p_sex_raster, 6, "Sexual selection")
ggsave("Plots/Maps/p_sexual_selection_map.png", height = 7, width = 15, dpi = 1000)

# Make a ggplot raster with s_sex_raster
s_sex_plot <- ggplot_ylim_raster(s_sex_raster, 6, "Sexual selection")
ggsave("Plots/Maps/s_sexual_selection_map.png", height = 7, width = 15, dpi = 1000)

diet_plots <- ggarrange(p_sex_plot, s_sex_plot, labels = c("a", "b"), nrow = 2, font.label = list(size = 28))
ggsave("Plots/Maps/pre_and_sec_plot.png", height = 13, width = 15, dpi = 1000)

sex_and_diet <- ggarrange(sex_plot, diet_plots, ncol = 2, widths = c(3,1))
ggsave("Plots/Maps/sex_and_diet_map.png", height = 7, width = 15, dpi = 1000)


## Territory maps ##

# Do primary and secondary consumers.
p_terr_plot <- ggplot_ylim_raster(p_terr_raster, 6, "")  +
  annotate("text", x = 20, y = -48, label  = "Primary consumers", size = 8, fontface = 2)
s_terr_plot <- ggplot_ylim_raster(s_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Secondary consumers", size = 8, fontface = 2)
p_year_terr_plot <- ggplot_terr_raster(p_year_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Primary consumers", size = 8, fontface = 2)
#hist(p_year_terr_raster, breaks = 1000)
#p_year_terr_raster_2 <- p_year_terr_raster
#values(p_year_terr_raster)[values(p_year_terr_raster) > 0.3 & !is.na(values(p_year_terr_raster))] <- NA
#p_year_terr_plot <- ggplot_cont_raster(log(p_year_terr_raster+1), 3, "") +
#  annotate("text", x = 20, y = -48, label  = "Primary consumers\n Year-Round Territoriality", size = 8, fontface = 2)
s_year_terr_plot <- ggplot_ylim_raster(s_year_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Secondary consumers", size = 8, fontface = 2)

# all_terrs <- ggarrange(p_terr_plot, s_terr_plot, p_year_terr_plot +
#                          annotate("text", x = 20, y = -48, label  = "Primary consumers\n Year-Round Territoriality", size = 8, fontface = 2), s_year_terr_plot, 
#           labels = c("a", "b", "c", "d"), nrow = 4, ncol = 1, font.label = list(size = 28), 
#           align = "hv")

#ggsave("Plots/Maps/test_terr.png", height = 35, width = 15, dpi = 1000)



## Seasonal ##

# Quick seasonality map.
temp_plot <- ggplot_legend_raster(temp_seas, 10, "Temperature\nseasonality")
ggsave("Plots/Maps/seasonality_map.png", height = 7, width = 15, dpi = 1000)


###############################################################################
                    #### Bin latitude ####

# Create 5 degree bins for latitude.
bin_range <- seq(0, 80, 5) 
bin_labels <- seq(0, 75, 5) 
phylo_data$binned_lat <- cut(abs(phylo_data$complete_latitude), breaks=bin_range, labels = bin_labels, include.lowest = TRUE)

# Make sexual selection certainty dummy score.
phylo_data$cert_dummy <- phylo_data$sexual_certainty
phylo_data$cert_dummy[phylo_data$cert_dummy > 1] <- 0

# Combine certainty for 1 and 2.
phylo_data$cert_dummy_2 <- phylo_data$sexual_certainty
phylo_data$cert_dummy_2[phylo_data$cert_dummy_2 > 2] <- 0
phylo_data$cert_dummy_2[phylo_data$cert_dummy_2 > 0] <- 1


# Make dummy territoriality maps.
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




# Group sex scores and certainty by lat bins.
lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait = first(binned_lat),
            sex_mean = mean(sexual_score),
            sex_sd = sd(sexual_score),
            sex_se = sd(sexual_score)/sqrt(length(sexual_score)),
            cert_mean = mean(cert_dummy),
            cert_sd = sd(cert_dummy),
            cert_se = sd(cert_dummy)/sqrt(length(cert_dummy)),
            cert_mean_2 = mean(cert_dummy_2),
            cert_sd_2 = sd(cert_dummy_2),
            cert_se_2 = sd(cert_dummy_2)/sqrt(length(cert_dummy_2))) %>% na.omit()

# Group the data by trophic level and bins.
diet_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% 
  summarise(trait = first(binned_lat),
            sex_mean = mean(sexual_score),
            sex_sd = sd(sexual_score),
            sex_se = sd(sexual_score)/sqrt(length(sexual_score))) %>% na.omit()

hi_lat_data <- phylo_data %>% filter(binned_lat != 75 & sexual_certainty == 1) %>% 
  group_by(binned_lat) %>% 
  summarise(trait = first(binned_lat),
            sex_mean = mean(sexual_score),
            sex_sd = sd(sexual_score),
            sex_se = sd(sexual_score)/sqrt(length(sexual_score)),
            cert_mean = mean(cert_dummy),
            cert_sd = sd(cert_dummy),
            cert_se = sd(cert_dummy)/sqrt(length(cert_dummy))) %>% na.omit()

med_lat_data <- phylo_data %>% filter(binned_lat != 75 & sexual_certainty < 3) %>% 
  group_by(binned_lat) %>% 
  summarise(trait = first(binned_lat),
            sex_mean = mean(sexual_score),
            sex_sd = sd(sexual_score),
            sex_se = sd(sexual_score)/sqrt(length(sexual_score)),
            cert_mean = mean(cert_dummy),
            cert_sd = sd(cert_dummy),
            cert_se = sd(cert_dummy)/sqrt(length(cert_dummy))) %>% na.omit()


# Group the data by trophic level territory and latitude.
terr_diet_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% 
  summarise(trait = first(binned_lat),
            sex_mean = mean(terr_dummy),
            sex_sd = sd(terr_dummy),
            sex_se = sd(terr_dummy)/sqrt(length(terr_dummy))) %>% na.omit()

year_terr_diet_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% 
  summarise(trait = first(binned_lat),
            sex_mean = mean(year_terr_dummy),
            sex_sd = sd(year_terr_dummy),
            sex_se = sd(year_terr_dummy)/sqrt(length(year_terr_dummy))) %>% na.omit()



###############################################################################
                      #### Add side plots for diet ####

library(ggnewscale)
# Create palettes for plotting.
light_colours <- c( "#77AD78","#7494EA", "#C98986", "#D8C77B")[2:3]
dark_colours <- c("#214F4B", "#05299E", "#8D0801","#A88A05")[2:3]


# Plot latitude for primary consumers.
pri_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
  ylab("Sexual selection") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", 
        text = element_text(face = "bold"))

# Plot latitude for secondary consumers.
sec_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) +
  ylab("Sexual selection") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", 
        text = element_text(face = "bold"))

p_sex_plot_2 <- p_sex_plot + annotate("text", x = 20, y = -48, label  = "Primary consumers", size = 8, fontface = 2)
s_sex_plot_2 <- s_sex_plot + annotate("text", x = 30, y = -48, label  = "Secondary consumers", size = 8, fontface = 2)

# Assemble them.
diet_lat_plots <- ggarrange(p_sex_plot_2, pri_lat_plot, s_sex_plot_2, sec_lat_plot, 
                           labels = c("a", "b", "c", "d"), nrow = 2, ncol = 2, 
                           widths = c(3,1.535), font.label = list(size = 28), 
                           align = "hv", hjust = c(-4,-3.5,-4,-3.5)) + theme(plot.margin = margin(l = -1.5, unit = "cm"))

# This is the one used in the figures document.
ggsave("Plots/Maps/diet_and_lat_plot.png", height = 10, width = 15, dpi = 1000)


###############################################################################
         #### Add side plots for sexual selection ####

# Sexual selection latitudinal gradient.
sex_lat_plot <- lat_data %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) +
  ylab("Sexual selection") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))

# Data certainty latitudinal gradient.
cert_lat_plot <- lat_data %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = cert_mean)) +
  geom_errorbar(aes(ymin = cert_mean - cert_se, ymax = cert_mean + cert_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col = "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) +
  ylab("Data certainty") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))

cert_lat_plot_2 <- lat_data %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = cert_mean_2)) +
  geom_errorbar(aes(ymin = cert_mean_2 - cert_se_2, ymax = cert_mean_2 + cert_se_2), 
                position = position_dodge(width = 1), show.legend = FALSE, col = "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) +
  ylab("Data certainty") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))


# Add the annotations.
sex_plot <- ggplot_ylim_raster(sex_raster, 6, "") + annotate("text", x = 20, y = -48, label  = "Sexual selection", size = 8, fontface = 2)
cert_plot <- ggplot_ylim_raster(cert_raster, 6, "") + annotate("text", x = 30, y = -48, label  = "Data certainty", size = 8, fontface = 2)
cert_plot_2 <- ggplot_ylim_raster(cert_raster_3, 6, "") + annotate("text", x = 30, y = -48, label  = "Data certainty", size = 8, fontface = 2)

# Assemble them.
sex_lat_plots <- ggarrange(sex_plot, sex_lat_plot, cert_plot, cert_lat_plot, 
                            labels = c("a", "b", "c", "d"), nrow = 2, ncol = 2, 
                            widths = c(3,1.535), font.label = list(size = 28), 
                           align = "hv", hjust = c(-4,-3.5,-4,-3.5))  + theme(plot.margin = margin(l = -1.5, unit = "cm"))

# Sex and latitude plot.
ggsave("Plots/Maps/sex_and_lat_plot.png", height = 10, width = 15, dpi = 1000)

# Assemble them.
sex_lat_plots_2 <- ggarrange(sex_plot, sex_lat_plot, cert_plot_2, cert_lat_plot_2, 
                           labels = c("a", "b", "c", "d"), nrow = 2, ncol = 2, 
                           widths = c(3,1.535), font.label = list(size = 28), 
                           align = "hv", hjust = c(-4,-3.5,-4,-3.5))  + theme(plot.margin = margin(l = -1.5, unit = "cm"))

ggsave("Plots/Maps/sex_and_lat_plot_2.png", height = 10, width = 15, dpi = 1000)

###############################################################################
     #### Add side plots for high and medium data certainty ####

# Sexual selection latitudinal gradient for high data.
hi_sex_lat_plot <- hi_lat_data %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) +
  ylab("Sexual selection") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))

med_sex_lat_plot <- med_lat_data %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) +
  ylab("Sexual selection") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))


# Add the annotations.
h_sex_plot_2 <- h_sex_plot + annotate("text", x = 20, y = -48, label  = "Certainty 1", size = 8, fontface = 2)
m_sex_plot_2 <- m_sex_plot + annotate("text", x = 30, y = -48, label  = "Certainty 1 & 2", size = 8, fontface = 2)

# Assemble them.
sex_lat_plots <- ggarrange(h_sex_plot_2, hi_sex_lat_plot, m_sex_plot_2, med_sex_lat_plot, 
                           labels = c("a", "b", "c", "d"), nrow = 2, ncol = 2, 
                           widths = c(3,1.535), font.label = list(size = 28), 
                           align = "hv", hjust = c(-4,-3.5,-4,-3.5)) + theme(plot.margin = margin(l = -1.5, unit = "cm"))
ggsave("Plots/Maps/figure_s1.png", height = 10, width = 15, dpi = 1000)



###############################################################################
            #### Add side plots for territoriality ####


pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% lat_side_plot(ylabel = "Proportion of species")
pri_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% lat_side_plot(ylabel = "Proportion of species")
sec_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% lat_side_plot(ylabel = "Proportion of species")
sec_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% lat_side_plot(ylabel = "Proportion of species")


# Sexual selection latitudinal gradient for high data.
pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) + 
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  ylab("Proportion of species") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))

pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% lat_side_plot(ylabel = "Proportion of species")
pri_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% lat_side_plot(ylabel = "Proportion of species")

terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% lm(formula = sex_mean ~ as.numeric(binned_lat))
year_terr_diet_lat_data


  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) + ylim(c(0,1.1)) +
  ylab("Year-round territoriality") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))


sec_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) + ylim(c(0,1.1)) +
  ylab("Territoriality") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))


sec_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  ggplot(aes(x = as.numeric(as.character(trait)), y = sex_mean)) +
  geom_errorbar(aes(ymin = sex_mean - sex_se, ymax = sex_mean + sex_se), 
                position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  geom_point(position = position_dodge(width = 1), col = "black") + 
  geom_smooth(se = FALSE, method = "glm", col = "black") + 
  scale_x_continuous(breaks = seq(from = 5, to = 75, by = 10)) + ylim(c(0,1.1)) +
  ylab("Year-round territoriality") + 
  xlab("Latitude") + theme_classic(base_size = 16) + 
  theme(legend.position = "none", # c(0.1,0.9),
        text = element_text(face = "bold"))





terr_lat_plots <- ggarrange(p_terr_plot, pri_terr_lat_plot, p_year_terr_plot, pri_yearterr_lat_plot,
                            s_terr_plot, sec_terr_lat_plot, s_year_terr_plot, sec_yearterr_lat_plot,
                           labels = c("a", "b", "c", "d", "e", "f", "g", "h"), nrow = 4, ncol = 2, 
                           widths = c(3,1.535), font.label = list(size = 28), 
                           align = "hv", hjust = c(-5,-4.5,-5,-4.5,-5,-9,-4.6,-4.5)) + 
  theme(plot.margin = margin(l = -1.5, unit = "cm"))
ggsave("Plots/Maps/territory_maps.png", height = 20, width = 15, dpi = 900)



###############################################################################
    #### Plot Seasonality without decimal places in legend ####

# Get the limits for breaks, using quantiles so each bin has an equal number of cells.
breaks <- quantile(values(temp_seas), seq(0,1,length.out=nbins+1), na.rm=TRUE)

# Cut the cell values into bins using the breaks.
cuts <- cut(values(temp_seas), breaks = breaks, include.lowest = TRUE)

# Replace cell values with bins.
temp_seas@data@values <- as.numeric(cuts)

# Create raster data.
temp_data <- as.data.frame(temp_seas, xy=TRUE)
colnames(temp_data) <- c("long", "lat", "values")
temp_data$values %<>% as.factor()

# Get rid of areas without land.
temp_data$land <- land_data$layer
temp_data %<>% drop_na(land)

# Create a label of the bins for plotting. Unicode is for en-dash.
labels <- as.character(round(breaks, 0))
labels <- paste0(labels[1:10], " \u2013 ", labels[2:11])

# Plot with ggplot.
temp_plot <- ggplot() + xlim(-180, 180) +  
  geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=temp_data) +
  scale_fill_manual(values = pal, breaks = 1:10, labels = labels, na.value = "lightgrey") +     #na.value = "grey"
  guides(fill = guide_legend("Temperature\nseasonality", byrow = TRUE)) +
  scale_x_discrete(expand=c(0,1))+
  scale_y_discrete(expand=c(0,1))+
  theme_classic(base_size = 18) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.position = c(0.075, 0.3),
        legend.key.height = unit(0.2, 'cm'), legend.spacing.y = unit(0.3, 'cm'),
        legend.title = element_text(face = "bold"),
        plot.margin = margin(t = -20, l = -20, b = -50, r = 5)) +
  ylab("") + xlab("") + coord_fixed()

# Put them together.
ggarrange(sex_plot, temp_plot, labels = c("a", "b"), nrow = 2, font.label = list(size = 24))
ggsave("Plots/Maps/scores_and_temp_plot.png", height = 13, width = 15, dpi = 1000)


###############################################################################
              #### Plot maps with latitude on side ####

# Create raster data.
p_raster_data <- as.data.frame(bin_raster(p_sex_raster), xy=TRUE)
colnames(p_raster_data) <- c("long", "lat", "values")

# Get rid of areas without land.
p_raster_data$land <- land_data$layer
p_raster_data %<>% drop_na(land)

# Do same for secondaries.
s_raster_data <- as.data.frame(bin_raster(s_sex_raster), xy=TRUE)
colnames(s_raster_data) <- c("long", "lat", "values")
s_raster_data$land <- land_data$layer
s_raster_data %<>% drop_na(land)


# Create 1 degree bins for latitude.
bin_range <- seq(-89.5, 89.5, 1) 
bin_labels <- seq(-89, 89, 1) 

# Bin the latitude.
p_raster_data$binned_lat <- cut(p_raster_data$lat, breaks=bin_range, labels = bin_labels)
p_raster_data$binned_lat %<>% as.character() %>% as.numeric()

s_raster_data$binned_lat <- cut(s_raster_data$lat, breaks=bin_range, labels = bin_labels)
s_raster_data$binned_lat %<>% as.character() %>% as.numeric()


# Group by 1 degree and take average.
p_lat_summarised <- p_raster_data %>% 
  group_by(binned_lat) %>% 
  summarize(sum = sum(values, na.rm = TRUE),
            mean = mean(values, na.rm = TRUE))

s_lat_summarised <- s_raster_data %>% 
  group_by(binned_lat) %>% 
  summarize(sum = sum(values, na.rm = TRUE),
            mean = mean(values, na.rm = TRUE))

wes_pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

(p_lat_plot <- p_lat_summarised %>% 
  group_by(binned_lat) %>% 
  summarize(sum = sum(sum, na.rm = TRUE),
            mean = mean(mean, na.rm = TRUE)) %>% 
  ggplot(aes(x = binned_lat, y = mean, colour = mean, fill=mean)) + 
  geom_bar(stat = "identity", position = position_dodge())  + 
  scale_colour_gradientn(colours = wes_pal) +
  scale_fill_gradientn(colours = wes_pal) +
  #expand_limits(x = c(-60,90)) +
  coord_flip() + theme_void() + theme(legend.position = "none"))

(s_lat_plot <- s_lat_summarised %>% 
    group_by(binned_lat) %>% 
    summarize(sum = sum(sum, na.rm = TRUE),
              mean = mean(mean, na.rm = TRUE)) %>% 
    ggplot(aes(x = binned_lat, y = mean, colour = mean, fill=mean)) + 
    geom_bar(stat = "identity", position = position_dodge())  + 
    scale_colour_gradientn(colours = wes_pal) +
    scale_fill_gradientn(colours = wes_pal) +
    #expand_limits(x = c(-60,90)) +
    coord_flip() + theme_void() + theme(legend.position = "none"))

plots_with_lat <- ggarrange(p_sex_plot, p_lat_plot, s_sex_plot, s_lat_plot, 
                        labels = c("a", "", "b", ""), nrow = 2, ncol=2, widths = c(15,1), 
                        font.label = list(size = 28), align = "v")
ggsave("Plots/Maps/lat_test.png", height = 13, width = 15, dpi = 1000)


p_raster_data$values %<>% as.factor()

# Create a label of the bins for plotting. Unicode is for en-dash.
labels <- as.character(format(round(breaks, 2), nsmall = 2))
labels <- paste0(labels[1:10], " \u2013 ", labels[2:11])

# Plot raster without margin edits.
sex_arrange_plot <- ggplot() +
  xlim(-180, 180) +  
  geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
  scale_fill_manual(values = pal, breaks = 1:10, labels = labels, na.value = "lightgrey") +     #na.value = "grey"
  guides(fill = guide_legend(variable, byrow = TRUE)) +
  scale_x_discrete(expand=c(0,1))+
  scale_y_discrete(expand=c(0,1))+
  theme_classic(base_size = 18) + theme(axis.text = element_blank(),
                                        axis.ticks = element_blank(),
                                        axis.line = element_blank(),
                                        legend.position = c(0.075, 0.3),
                                        legend.key.height = unit(0.2, 'cm'),
                                        legend.spacing.y = unit(0.3, 'cm'),
                                        legend.title = element_text(face = "bold")) +
  ylab("") + 
  xlab("") + 
  coord_fixed()

# Arrange together.
library(gridExtra)
png("Plots/Maps/test_lat.png", width = 15000, height = 7000)
grid.arrange(NULL, NULL, sex_arrange_plot, sex_lat_plot, NULL, NULL, nrow = 3,
             ncol = 2, widths = c(10,1), heights = c(0.1, 1, 0.1), respect = FALSE, clip = "on")
dev.off()




###############################################################################
      #### Make a histogram showing the colour scale ####

# Change scipen so that it plots 
options(scipen=999)

# Make the histogram.
sex_hist <- ggplot_colour_hist(sex_raster, "Sexual selection", scaled = FALSE,
                                x_axis_breaks = seq(0, 1.5, 0.5), x_axis_lim = c(0,1.5))
# Put together with map. Kept the code in for now for addling plot label.
sex_plot <- ggplot_raster(sex_raster) +
  annotation_custom(ggplotGrob(sex_hist), xmin = -180, xmax = -100, ymin = -60, ymax = 0)
  #annotate("text", x = -170, y = 75, label  = "a", size = 70, fontface = 2)
ggsave("Plots/Maps/sexual_selection_hist.png", height = 7, width = 15, dpi = 1000)

# Data certainty.
cert_hist <- ggplot_colour_hist(cert_raster, "Average Data Certainty", scaled = FALSE,
                                x_axis_breaks = seq(0, 1, 0.5), x_axis_lim = c(0,1))
# Put together with map. Kept the code in for now for addling plot label.
cert_plot <- ggplot_raster(cert_raster) +
  annotation_custom(ggplotGrob(cert_hist), xmin = -180, xmax = -100, ymin = -60, ymax = 0)
#annotate("text", x = -170, y = 75, label  = "a", size = 70, fontface = 2)

ggsave("Plots/Maps/sexual_certainty_hist.png", height = 7, width = 15, dpi = 1000)

# Seasonality.
temp_hist <- ggplot_colour_hist(temp_seas, "Temperature Seasonality", scaled = FALSE, 
                                    x_axis_breaks = seq(0, 20000, 10000), x_axis_lim = c(0,21600))
temp_plot <- ggplot_raster(temp_seas) +
  annotation_custom(ggplotGrob(temp_hist), xmin = -180, xmax = -100, ymin = -60, ymax = 0)
ggsave("Plots/Maps/temp_hist.png", height = 7, width = 15, dpi = 1000)


ggarrange(sex_plot, cert_plot, labels = c("a", "b"), nrow = 2, label.x = 1)
ggsave("Plots/Maps/sex_scores_plot.png", height = 7, width = 15, dpi = 1000)






# Primary consumers.
(p_sex_hist <- ggplot_colour_hist(p_sex_raster, "Sexual selection: Primary consumers", 
                                  scaled = FALSE,
                                x_axis_breaks = seq(0, 1, 0.5), x_axis_lim = c(0,1)))
p_sex_plot <- ggplot_raster(p_sex_raster) +
  annotation_custom(ggplotGrob(p_sex_hist), xmin = -180, xmax = -100, ymin = -60, ymax = 0)
ggsave("Plots/Maps/pri_sexual_selection_hist.png", height = 7, width = 15, dpi = 1000)

# Secondary.
(s_sex_hist <- ggplot_colour_hist(s_sex_raster, "Sexual selection: Secondary consumers", 
                                  scaled = FALSE,
                                  x_axis_breaks = seq(0, 1, 0.5), x_axis_lim = c(0,1)))
s_sex_plot <- ggplot_raster(s_sex_raster) +
  annotation_custom(ggplotGrob(s_sex_hist), xmin = -180, xmax = -100, ymin = -60, ymax = 0)
ggsave("Plots/Maps/sec_sexual_selection_hist.png", height = 7, width = 15, dpi = 1000)


sex_raster

plot(s_sex_raster)
plot(p_sex_raster)



function(){}


  # Get the limits for breaks, using quantiles so each bin has an equal number of cells.
  breaks <- quantile(values(s_sex_raster), seq(0,1,length.out=nbins+1), na.rm=TRUE)
  
  # Cut the cell values into bins using the breaks.
  cuts <- cut(values(s_sex_raster), breaks = breaks, include.lowest = TRUE)
  
  # Replace cell values with bins.
  s_sex_raster@data@values <- as.numeric(cuts)
  
  # Create raster data.
  raster_data <- as.data.frame(s_sex_raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  raster_data$values %<>% as.factor()
  
  # Get rid of areas without land.
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  
  # Create a label of the bins for plotting. Unicode is for en-dash.
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:10], " \u2013 ", labels[2:11])
  #labels <- c("", labels)
  # Plot with ggplot.
plot_test <- ggplot() +
    
    # Clip to world extent.
    xlim(-180, 180) +  
    
    # Add the raster data.
    geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_manual(values = pal, breaks = 1:10, labels = labels, na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend("Sexual selection", byrow = TRUE)) +
    
    # Make map closer to the edge.
    scale_x_discrete(expand=c(0,1))+
    scale_y_discrete(expand=c(0,1))+
    # Theme stuff.
    theme_classic(base_size = 18) + theme(axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            legend.position = c(0.075, 0.3),
                            legend.key.height = unit(0.2, 'cm'),
                            legend.spacing.y = unit(0.3, 'cm'),
                            legend.title = element_text(face = "bold"),
                            plot.margin = margin(t = -20, l = -20, b = -50, r = 5)     #unit(c(0,0,0,0), "null")
    ) +
    ylab("") + 
    xlab("") + 
    coord_fixed()

# 
# theme(legend.position = c(0.175, 0.85),
#       axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
#       axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
#       axis.title.x=element_text(size=rel(0.7), face = "bold"),
#       legend.text = element_text(size=rel(0.8), face = "bold"), 
#       legend.key.height = unit(0.2, 'cm'),
#       legend.spacing.y = unit(0.2, 'cm'),
#       # legend.margin = margin(b=7, t = 2, l = 2, r = 4),
#       # legend.background = element_rect(colour = 'black', 
#       #                                  fill = 'white', 
#       #                                  linetype='dashed', 
#       #                                  linewidth = 0.1)
# )
# 
# 

ggsave("Plots/Maps/sec_sexual_selection.png", height = 6, width = 15, dpi = 1000)
  











                             ##### End ######
########################################################################################


# Make a basic function for making maps with ggplot.
ggplot_raster <- function(raster, nbins = 10, variable = ""){
  
  # Get the limits for breaks, using quantiles so each bin has an equal number of cells.
  breaks <- quantile(values(raster), seq(0,1,length.out=nbins+1), na.rm=TRUE)
  
  # Cut the cell values into bins using the breaks.
  cuts <- cut(raster@data@values, breaks = breaks, include.lowest = TRUE)
  
  # Replace cell values with bins.
  raster@data@values <- as.numeric(cuts)
  
  # Create raster data.
  raster_data <- as.data.frame(raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  raster_data$values %<>% as.factor()
  
  # Get rid of areas without land.
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  
  # Create a label of the bins for plotting. Unicode is for en-dash.
  labels <- rev(as.character(round(breaks, 2)))
  labels <- paste0(labels[2:11], " \u2013 ", labels[1:10])
  
  # Plot with ggplot.
  ggplot() +
    
    # Clip to world extent.
    xlim(-180, 180) +  
    
    # Add the raster data.
    geom_raster(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_manual(values = pal, labels = labels, na.value = "grey") +
    guides(fill = guide_legend(variable)) +
    
    # Make map closer to the edge.
    scale_x_discrete(expand=c(0,1))+
    scale_y_discrete(expand=c(0,1))+
    # Theme stuff.
    theme_classic() + theme(axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"
                            #plot.margin = unit(c(0,0,0,0), "null"),
                            #panel.margin = unit(c(0, 0, 0, 0), "null")
    ) +
    ylab("") + 
    xlab("") + 
    coord_fixed()
}













# First convert dichromatism raster into xy data.
sex_data <- as.data.frame(sex_raster, xy=TRUE)
colnames(sex_data) <- c("long", "lat", "values")

# Get rid of areas without land.
sex_data$land <- land_data$layer
sex_data %<>% drop_na(land)

# Scale the values so that zero is near the middle and looks nicer.
sex_data$scaled_values <- scale(sex_data$values)

# Create breaks for the new scaled values.
breaks <- quantile(sex_data$scaled_values, seq(0,1,length.out=(nbins+1)), na.rm=TRUE)
breaks <- quantile(log(sex_data$values + 1), seq(0,1,length.out=(nbins+1)), na.rm=TRUE)
breaks <- quantile(sex_data$values, seq(0,1,length.out=(nbins+1)), na.rm=TRUE)


library(ggpubr)

# Create a plot that could be used as a legend.
(colour_plot <- ggplot(sex_data,
                       aes(x = values, 
                           # Get fill to be discrete version of x axis, using cuts from quantile breaks.
                           fill = cut(..x.., breaks = breaks, 
                                      include.lowest = TRUE, na.rm = TRUE))) + 
    # Set the number of breaks. For some reason varying this number seems to make certain
    # bins on the outskirts of the plot go grey coloured. Not sure why.
    geom_histogram(bins = 150) +
    
    # Colour scale.
    scale_fill_manual(values = pal, labels = labels) +
    
    # X axis label.
    labs(x = "Dichromatism", y = NULL, fill = NULL, colour = NULL) + 
    
    # Standard theme stuff.
    theme_pubr(base_size = 5, legend = "none") + 
    theme(axis.text.x=element_text(size=rel(0.5), face = "bold"), 
          axis.text.y=element_text(size=rel(0.5), face = "bold"),
          axis.title.x=element_text(size=rel(0.6), face = "bold")))

# # New device.
# dev.new(height = 10, width = 20, units = "mm", res = 1500, noRStudioGD = TRUE)

# Remove the original legend and add the new plot.
ggplot_raster(human_raster, 10, "Dichromatism") + rremove("legend") +
  annotation_custom(ggplotGrob(colour_plot), xmin = -180, xmax = -90, ymin = -60, ymax = 0)

# Export the plot.
ggsave("Plots/Maps/human_dichromatism_hist.png", height = 3, width = 6, dpi = 5000)





# First convert dichromatism raster into xy data.
avian_data <- as.data.frame(avian_raster, xy=TRUE)
colnames(avian_data) <- c("long", "lat", "values")

# Get rid of areas without land.
avian_data$land <- land_data$layer
avian_data %<>% drop_na(land)

# Create breaks for the new scaled values.
breaks <- quantile(log(avian_data$values), seq(0,1,length.out=(nbins+1)), na.rm=TRUE)

(avian_colour_plot <- ggplot(avian_data,
                       aes(x = log(values), 
                           # Get fill to be discrete version of x axis, using cuts from quantile breaks.
                           fill = cut(..x.., breaks = breaks, 
                                      include.lowest = TRUE, na.rm = TRUE))) + 
    # Set the number of breaks. For some reason varying this number seems to make certain
    # bins on the outskirts of the plot go grey coloured. Not sure why.
    geom_histogram(bins = 150) +
    
    # Colour scale.
    scale_fill_manual(values = pal, labels = labels) +
    
    # X axis label.
    labs(x = "Dichromatism", y = NULL, fill = NULL, colour = NULL) + 
    
    # Standard theme stuff.
    theme_pubr(base_size = 5, legend = "none") + 
    theme(axis.text.x=element_text(size=rel(0.5), face = "bold"), 
          axis.text.y=element_text(size=rel(0.5), face = "bold"),
          axis.title.x=element_text(size=rel(0.6), face = "bold")))


# Remove the original legend and add the new plot.
ggplot_raster(avian_raster, 10, "Dichromatism") + rremove("legend") +
  annotation_custom(ggplotGrob(avian_colour_plot), xmin = -180, xmax = -90, ymin = -60, ymax = 0)

# Export the plot.
ggsave("Plots/Maps/avian_dichromatism_hist.png", height = 3, width = 6, dpi = 5000)



.ggplot(sex_ranges, aes(x = trophic_dummy)) + geom_density() 

hist(log(values(trophic_raster)+1))
hist(log(values(human_raster)+1))

min(na.omit(values(human_raster))+1)

cor(sqrt(na.omit(values(trophic_raster))), sqrt(na.omit(values(human_raster))))


toms_colours <- c("dodgerblue4", "dodgerblue3", 
  "dodgerblue2","dodgerblue1", 
  "skyblue", "tomato", 
  "red1","firebrick")

chris_colours <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00') # # colorRampPalette(GetColors(nbins, scheme = "ocean"))

my_colours <- c("skyblue", "red")

# Convert the raster into a raster dataframe.
raster_data <- as.data.frame(human_raster, xy=TRUE) %>% drop_na()
colnames(raster_data) <- c("long", "lat", "dichromatism")

# Plot with ggplot.
dicrho_plot <- ggplot() +
  borders(ylim = c(-60,90), fill = "grey90", colour = "grey90") +
  xlim(-180, 180) + 
  geom_tile(aes(x=long, y=lat, fill= log(dichromatism+1)), data=raster_data) +
  
  # Here we add a name to the legend, and set manual colours for either end of a gradient.
  scale_fill_gradientn(name = "Average\nDichromatism", colors = my_colours) +
  
  # scale_fill_gradientn(name = paste0("Average", '\n', "Dichromatism"),
  #                      colors = c("dodgerblue4", "dodgerblue3", "dodgerblue2","dodgerblue1", 
  #                                 "skyblue", "tomato", "red1","firebrick"),
  #                      na.value = "grey80",
                       #breaks = c(-3.2, -3.6, -4.0)
  #) +
  
  # You should be getting used to this code!
  #ggtitle("Accipitridae Species Richness Heat Map") + 
  theme_classic() +
  ylab("") + 
  xlab("") + 
  coord_fixed()

# New device.
dev.new(height = 10, width = 20, units = "mm", res = 1500, noRStudioGD = TRUE)

# Return the plot so we can view it.
dicrho_plot


trophic_data <- as.data.frame(trophic_raster, xy=TRUE) %>% drop_na()
colnames(trophic_data) <- c("long", "lat", "herbivory")

herbivory_plot <- ggplot() +
  borders(ylim = c(-60,90), fill = "grey90", colour = "grey90") +
  xlim(-180, 180) + 
  geom_tile(aes(x=long, y=lat, fill= log(herbivory+1)), data=trophic_data) +
  
  # Here we add a name to the legend, and set manual colours for either end of a gradient.
  scale_fill_gradientn(name = "Percent\nHerbivory", colors = my_colours) +
  #scale_fill_gradientn(name = "Average Dichromatism", colors = c("navy", "firebrick")) +
  
  # scale_fill_gradientn(name = paste0("Average", '\n', "Dichromatism"),
  #                      colors = c("dodgerblue4", "dodgerblue3", "dodgerblue2","dodgerblue1", 
  #                                 "skyblue", "tomato", "red1","firebrick"),
  #                      na.value = "grey80",
  #breaks = c(-3.2, -3.6, -4.0)
  #) +
  
  # You should be getting used to this code!
  #ggtitle("Accipitridae Species Richness Heat Map") + 
theme_classic() +
  ylab("Latitude") + 
  xlab("Longitude") + 
  coord_fixed()


herbivory_plot

library(ggpubr)
tiff("range_maps.tiff", width = 7000, height = 4000, res = 400, compression = "lzw+p")
ggarrange(dicrho_plot, herbivory_plot, ncol = 1)
dev.off()

tiff("ggplot_sqrt_range_maps.tiff", width = 7000, height = 4000, res = 400, compression = "lzw+p")
ggarrange(dicrho_plot, herbivory_plot, ncol = 1)
dev.off()


tiff("log_range_maps.tiff", width = 7000, height = 4000, res = 400, compression = "lzw+p")
ggarrange(dicrho_plot, herbivory_plot, ncol = 1)
dev.off()

ggsave("range_maps.tiff", dpi = 600, width = 5, height = 3)

tiff("sqrt_range_maps.tiff", width = 7000, height = 4000, res = 400, compression = "lzw+p")
par(mfrow=c(2,1))
plot(sqrt(human_raster))
plot(sqrt(trophic_raster))
dev.off()
## Tom's version ##


human_raster[which(is.na(getValues(land)))] <- NA
all_rasdf <- as.data.frame(human_raster, xy = TRUE) %>% drop_na()
all_rasdf$layer[all_rasdf$layer == 0] <- NA
colnames(all_rasdf) <- c("long", "lat", "dichromatism")



# Plot with ggplot.
dicrho_plot <- ggplot() +
  #borders(ylim = c(-60,90), fill = "grey90", colour = "grey90") +
  xlim(-180, 180) + 
  geom_tile(aes(x=long, y=lat, fill= dichromatism), data=all_rasdf) +
  
  # Here we add a name to the legend, and set manual colours for either end of a gradient.
  #scale_fill_gradientn(name = "Average Dichromatism", colors = c("skyblue", "red")) +
  #scale_fill_gradientn(name = "Average Dichromatism", colors = c("navy", "firebrick")) +
  
  scale_fill_gradientn(name = paste0("Average", '\n', "Dichromatism"),
                       colors = c("dodgerblue4", "dodgerblue3", "dodgerblue2","dodgerblue1", 
                                  "skyblue", "tomato", "red1","firebrick"),
                       na.value = "grey80",
                       #breaks = c(-3.2, -3.6, -4.0)
                       ) +
  
  # You should be getting used to this code!
  #ggtitle("Accipitridae Species Richness Heat Map") + 
  theme_classic() +
  ylab("Latitude") + 
  xlab("Longitude") + 
  coord_fixed()

# Return the plot so we can view it.
dicrho_plot




#
# library(tidyr)
# library(ggplot2)
# library(raster)
# # Convert the raster into a raster dataframe. This will be coordinates of the
# # raster pixels (cols x and y) and the value of the raster pixels.
# raster_data <- as.data.frame(range_raster, xy=TRUE)
#
# # Remove NAs with no information (like parts of the sea)
# raster_data <- na.omit(raster_data)
#
# # Change the column names to something sensible.
# colnames(raster_data) <- c("long", "lat", "richness")
#
# # Create a plot with ggplot (the plus signs at the end of a line carry over to the next line).
# range_plot <- ggplot(raster_data) +
#
#   # borders imports all the country outlines onto the map.
#   # colour changes the colour of the outlines.
#   # fill changes the colour of the insides of the countries.
#   # this will grey out any terrestrial area which isn't part of a range.
#   borders(ylim = c(-60,90), fill = "grey90", colour = "grey90") +
#
#   # Borders() xlim is -160/200 to catch the edge of russia.
#   # We need to reset the xlim to -180/180 to fit our raster_stack.
#   xlim(-180, 180) +
#
#   # Add the range information on top.
#   geom_tile(aes(x=long, y=lat, fill= richness)) +
#
#   # Add a colour blind friendly scale.
#   scale_fill_viridis_c() +
#
#   # Add title.
#   ggtitle("Accipitidae Species Richness") +
#
#   # Add the classic theme (things like gridlines, font etc.)
#   theme_classic() +
#
#   # Add axes labels.
#   ylab("Latitude") +
#   xlab("Longitude") +
#
#   # coord_fixed() makes ggplot keep our aspect ratio the same.
#   coord_fixed()
#
# # Change the size of the plot window and return the plot.
# #options(repr.plot.width=15, repr.plot.height=10)
# range_plot
#
# library(magrittr)
# avonet_data %<>% filter(jetz_name %in% new_jetz_ranges$SCINAME)
#
# # And lets add a column to our data for storing if it's a small or large range.
# avonet_data$range_large <- NA
#
# # We'll use a basic loop that goes from 1 to 237.
# row_numbers <- 1:nrow(avonet_data)
#
# # The curly brackets show the beginning and the end of the loop.
# for (x in row_numbers){
#
#   # Pull out the range size we want for each iteration (x) of the loop.
#   range <- avonet_data$range_size[x]
#
#   # Calculate if it's small range or a large range.
#   if (range > 1000000){
#     range_large <- 1
#   } else {
#     range_large <- 0
#   }
#
#   # Lastly we want to add our new value to the dataframe.
#   avonet_data$range_large[x] <- range_large
# }
#
#
#
# # Load fasterize package.
# library(fasterize)
#
# # Combine the two datasets into one object so we have range size info and the polygons together.
# # This turns our sf object into a normal dataframe.
# Accip_all <- left_join(avonet_data, new_jetz_ranges, by = c("jetz_name" = "SCINAME"))
#
# # Start by creating an empty raster stack to store our data in.
# raster_template <- raster(ncols=2160, nrows = 900, ymn = -60)
#
# # 'fasterize' needs objects to be an sf class so we'll convert it back.
# Accip_all <- st_sf(Accip_all)
#
# # Use the fasterize function with the raster template. We want to use the
# # range_large field, and the function min takes the smallest value when they overlap.
# # (so small ranges are shown on top of large ranges)
# range_raster <- fasterize(Accip_all, raster_template, field = "range_size", fun = "sum")
# ?fasterize
# # Plot the new map.
# plot(range_raster, col=rainbow(2))
#
#
# library(tidyr)
# library(ggplot2)
#
# # Convert the raster into a raster dataframe.
# raster_data <- as.data.frame(range_raster, xy=TRUE) %>% drop_na()
# colnames(raster_data) <- c("long", "lat", "index")
#
# # Add labels for the range sizes so that ggplot colours them as discrete, rather than a continuous number.
# raster_data$ranges[raster_data$index == 0] <- "Small"
# raster_data$ranges[raster_data$index == 1] <- "Large"
#
# # We can then plot this in ggplot. We have to first create the colour scheme for our map.
# myColors <- c("grey80", "red")
#
# # Assign names to these colors that correspond to each range size.
# names(myColors) <- unique(raster_data$ranges)
#
# # Create the colour scale.
# colScale <- scale_fill_manual(name = "Range Sizes", values = myColors)
#
# # Create a plot with ggplot (the plus signs at the end of a line carry over to the next line).
# range_plot <- ggplot() +
#   # borders imports all the country outlines onto the map.
#   # colour changes the colour of the outlines,
#   # fill changes the colour of the insides of the countries.
#   # This will grey out any terrestrial area which isn't part of a range.
#   borders(ylim = c(-60,90), fill = "grey90", colour = "grey90") +
#
#   # Borders() xlim is -160/200 to catch the edge of Russia. We need to reset the
#   # xlim to -180/180 to fit our raster_stack.
#   xlim(-180, 180) +
#
#   # Add the range information on top.
#   geom_tile(aes(x=long, y=lat, fill= ranges), data=raster_data) +
#   # Add colours.
#   colScale +
#   # Add title.
#   ggtitle("Small range sizes in the Accipitidae") +
#   # Add the classic theme (things like gridlines, font etc.)
#   theme_classic() +
#   # Add axes labels.
#   ylab("Latitude") +
#   xlab("Longitude") +
#   # coord_fixed() makes ggplot keep our aspect ratio the same.
#   coord_fixed()
#
# # Return the plot so we can view it.
# options(repr.plot.width=15, repr.plot.height=10)
# range_plot
#
#
#
#
#
#
# plot(range_raster)
#
# test_raster <- range_raster / species_raster
#
# plot(test_raster)
#
#
# # raster pixels (cols x and y) and the value of the raster pixels.
# raster_data <- as.data.frame(test_raster, xy=TRUE)
#
# # Remove NAs with no information (like parts of the sea)
# raster_data <- na.omit(raster_data)
#
# # Change the column names to something sensible.
# colnames(raster_data) <- c("long", "lat", "richness")
#
# # Create a plot with ggplot (the plus signs at the end of a line carry over to the next line).
# range_plot <- ggplot(raster_data) +
#
#   # borders imports all the country outlines onto the map.
#   # colour changes the colour of the outlines.
#   # fill changes the colour of the insides of the countries.
#   # this will grey out any terrestrial area which isn't part of a range.
#   borders(ylim = c(-60,90), fill = "grey90", colour = "grey90") +
#
#   # Borders() xlim is -160/200 to catch the edge of russia.
#   # We need to reset the xlim to -180/180 to fit our raster_stack.
#   xlim(-180, 180) +
#
#   # Add the range information on top.
#   geom_tile(aes(x=long, y=lat, fill= richness)) +
#
#   # Add a colour blind friendly scale.
#   scale_fill_viridis_c() +
#
#   # Add title.
#   ggtitle("Accipitidae Species Richness") +
#
#   # Add the classic theme (things like gridlines, font etc.)
#   theme_classic() +
#
#   # Add axes labels.
#   ylab("Latitude") +
#   xlab("Longitude") +
#
#   # coord_fixed() makes ggplot keep our aspect ratio the same.
#   coord_fixed()
#
# # Change the size of the plot window and return the plot.
# #options(repr.plot.width=15, repr.plot.height=10)
# range_plot
#
#
