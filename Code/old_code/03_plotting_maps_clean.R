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
# temp_seas <- raster("Data/Chelsa/CHELSA_bio4_1981-2010_V.2.1.tif")
# npp <- raster("Data/Chelsa/CHELSA_npp_1981-2010_V.2.1.tif")
# 
# library(terra)
# # Resample so they're same extent / res
# temp_seas <- resample(rast(temp_seas), rast(species_raster))
# npp <- resample(rast(npp), rast(species_raster))
# 
# # Mask.
# temp_seas <- mask(raster(temp_seas), species_raster)
# npp <- mask(raster(npp), species_raster)


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
                        #### Bin latitude ####

# Create 5 degree bins for latitude.
bin_range <- seq(0, 80, 5) 
bin_labels <- seq(0, 75, 5) 
phylo_data$binned_lat <- cut(abs(phylo_data$complete_latitude), breaks=bin_range, labels = bin_labels, include.lowest = TRUE)
phylo_data$binned_lat %<>% as.character() %>% as.numeric()

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
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se,
            cert_mean = mean(cert_dummy),
            cert_sd = sd(cert_dummy),
            cert_se = sd(cert_dummy)/sqrt(length(cert_dummy)),
            cert_mean_2 = mean(cert_dummy_2),
            cert_sd_2 = sd(cert_dummy_2),
            cert_se_2 = sd(cert_dummy_2)/sqrt(length(cert_dummy_2))) %>% na.omit()

cert_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(cert_dummy),
            trait_sd = sd(cert_dummy),
            trait_se = sd(cert_dummy)/sqrt(length(cert_dummy)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()


# Group the data by trophic level and bins.
diet_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()

hi_lat_data <- phylo_data %>% filter(binned_lat != 75 & sexual_certainty == 1) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se,
            cert_mean = mean(cert_dummy),
            cert_sd = sd(cert_dummy),
            cert_se = sd(cert_dummy)/sqrt(length(cert_dummy))) %>% na.omit()

med_lat_data <- phylo_data %>% filter(binned_lat != 75 & sexual_certainty < 3) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se,
            cert_mean = mean(cert_dummy),
            cert_sd = sd(cert_dummy),
            cert_se = sd(cert_dummy)/sqrt(length(cert_dummy))) %>% na.omit()


# Group the data by trophic level territory and latitude.
terr_diet_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% 
  summarise(trait_mean = mean(terr_dummy),
            trait_sd = sd(terr_dummy),
            trait_se = sd(terr_dummy)/sqrt(length(terr_dummy)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()
terr_diet_lat_data$trait_min[terr_diet_lat_data$trait_min < 0] <- 0

year_terr_diet_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% 
  summarise(trait_mean = mean(year_terr_dummy),
            trait_sd = sd(year_terr_dummy),
            trait_se = sd(year_terr_dummy)/sqrt(length(year_terr_dummy)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()



###############################################################################
                 #### Plot sexual selection figure 2 ####


# Create data frame of land to use to crop ggplot maps.
land_data <- as.data.frame(land, xy=TRUE)

# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

# Blank plot.
blank_plot <- ggplot() + geom_blank() + theme_classic() +
  theme(axis.line = element_blank())

# Full sexual selection map.
sex_plot <- ggplot_raster(sex_raster, 6, "Sexual selection") + 
  annotate("text", x = 20, y = -48, label  = "All birds", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "a", size = 12, fontface = 2)
  
# Dichotomised data certainty map.
cert_plot <- ggplot_raster(cert_raster, 6, "Data certainty") + 
  annotate("text", x = 20, y = -48, label  = "Data certainty", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "c", size = 12, fontface = 2)

# Latiduinal relationship with sex.
sex_lat_plot <- lat_data %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1), 
                ybreaks =  c(0,0.5,1.0), lab_ypos = 0.2) + 
  annotate("text", x = 0, y = 1, label = "b", size = 12, fontface = 2)  + rremove("xlab") + rremove("x.text")

# Lat relationship with data certainty.
cert_lat_plot <- cert_lat_data %>% 
  lat_side_plot(ylabel = "Data certainty", ylimits = c(0,1), 
                ybreaks =  c(0,0.5,1.0), lab_ypos = 0.2) + 
  annotate("text", x = 0, y = 1, label = "d", size = 12, fontface = 2)

# Arrange side plots.
sex_side_plots <- ggarrange(sex_lat_plot + rremove("xlab") + rremove("x.text"), 
                            cert_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
# Arrange maps.
sex_maps <- ggarrange(blank_plot, sex_plot, cert_plot, blank_plot, nrow = 4, ncol = 1,
                       heights = c(0.1,1,1,0.1))

# Put all together.
sex_both_plots <- ggarrange(sex_maps, sex_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Maps/figure_2.png", height = 10, width = 15, dpi = 600)

library(patchwork)

layout <- "
AAAAAAAABBBB
AAAAAAAABBBB
AAAAAAAABBBB
AAAAAAAABBBB
CCCCCCCCDDDD
CCCCCCCCDDDD
CCCCCCCCDDDD
CCCCCCCCDDDD
"
figure_2 <- wrap_plots(A = sex_plot, B = sex_lat_plot, C = cert_plot, D = cert_lat_plot, design = layout) 
ggsave("Plots/Maps/figure_2_patch.png", height = 10, width = 15, dpi = 600)


###############################################################################
               #### Plot trophic level figure 4 ####

# Plot sexual selection maps for both trophic levels.
p_sex_plot <- ggplot_raster(p_sex_raster, 6, "Sexual selection") + 
  annotate("text", x = 20, y = -48, label  = "Primary consumers", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "a", size = 12, fontface = 2)
s_sex_plot <- ggplot_raster(s_sex_raster, 6, "Sexual selection") + 
  annotate("text", x = 20, y = -48, label  = "Secondary consumers", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "c", size = 12, fontface = 2)

# Lat relationship for primary & secondary consumers.
pri_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
                ybreaks =  c(0,1.0, 2.0), lab_ypos = 0.1) + 
  annotate("text", x = 0, y = 2, label = "b", size = 12, fontface = 2)
sec_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1), 
                ybreaks =  c(0,0.5,1.0), lab_ypos = 0.1) + 
  annotate("text", x = 0, y = 1, label = "d", size = 12, fontface = 2)

# Arrange maps together.
diet_side_plots <- ggarrange(pri_lat_plot + rremove("xlab") + rremove("x.text"), 
                            sec_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
diet_maps <- ggarrange(blank_plot, p_sex_plot, s_sex_plot, blank_plot, nrow = 4, ncol = 1,
                      heights = c(0.1,1,1,0.1))
diet_both_plots <- ggarrange(diet_maps, diet_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Maps/figure_4.png", height = 10, width = 15, dpi = 600)


###############################################################################
          #### Plot data certainty supplementary figure ####

# Show high certainty sexual selection maps.
h_sex_plot <- ggplot_raster(h_sex_raster, 6, "Sexual selection") +
  annotate("text", x = 20, y = -48, label  = "Certainty 1", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "a", size = 12, fontface = 2)
m_sex_plot <- ggplot_raster(m_sex_raster, 6, "Sexual selection") + 
  annotate("text", x = 30, y = -48, label  = "Certainty 1 & 2", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "c", size = 12, fontface = 2)
# And side plots.
hi_sex_lat_plot <- hi_lat_data  %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
                ybreaks =  c(0,1.0, 2.0), lab_ypos = 0.1) + 
  annotate("text", x = 0, y = 2, label = "b", size = 12, fontface = 2)
med_sex_lat_plot <- med_lat_data  %>% 
    lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1), 
                  ybreaks =  c(0,0.5,1.0), lab_ypos = 0.1) + 
    annotate("text", x = 0, y = 1, label = "d", size = 12, fontface = 2)

# Arrange maps together.
cert_side_plots <- ggarrange(hi_sex_lat_plot + rremove("xlab") + rremove("x.text"), 
                             med_sex_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
cert_maps <- ggarrange(blank_plot, h_sex_plot, m_sex_plot, blank_plot, nrow = 4, ncol = 1,
                       heights = c(0.1,1,1,0.1))
cert_both_plots <- ggarrange(cert_maps, cert_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Maps/figure_S1.png", height = 10, width = 15, dpi = 600)

# Species richness of 3 and 4 cert species.
l_cert_sr_plot <- ggplot_raster(l_species_raster, 6, "Species Richness")# +
ggsave("Plots/Maps/low_cert_SR_map.png", height = 7, width = 15, dpi = 1000)


###############################################################################
        #### Plot latitudinal gradient in territoriality ####

# Territoriality maps.
p_terr_plot <- ggplot_raster(p_terr_raster, 6, "")  +
  annotate("text", x = 20, y = -48, label  = "Primary consumers", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "a", size = 12, fontface = 2)
s_terr_plot <- ggplot_raster(s_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Secondary consumers", size = 8, fontface = 2)+
  annotate("text", x =  -170, y = 80, label  = "c", size = 12, fontface = 2)
p_year_terr_plot <- ggplot_terr_raster(p_year_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Primary consumers", size = 8, fontface = 2)+
  annotate("text", x =  -170, y = 80, label  = "e", size = 12, fontface = 2)
s_year_terr_plot <- ggplot_raster(s_year_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Secondary consumers", size = 8, fontface = 2)+
  annotate("text", x =  -170, y = 80, label  = "g", size = 12, fontface = 2)

# Territoriality side plots.
pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Proportion of species") + 
  annotate("text", x = 0, y = 1.05, label = "b", size = 10, fontface = 2)
pri_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Proportion of species") + 
  annotate("text", x = 0, y = 1.05, label = "d", size = 10, fontface = 2)
sec_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Proportion of species") + 
  annotate("text", x = 0, y = 1.05, label = "f", size = 10, fontface = 2)
sec_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Proportion of species") + 
  annotate("text", x = 0, y = 1.05, label = "h", size = 10, fontface = 2)

# Arrange side plots and maps.
terr_side_plots <- ggarrange(pri_terr_lat_plot + rremove("xlab") + rremove("x.text"), 
                             pri_yearterr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_terr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_yearterr_lat_plot,
                             nrow = 4, ncol = 1, heights = c(1,1,1,1.2))
terr_maps <- ggarrange(blank_plot, p_terr_plot,p_year_terr_plot, 
                       s_terr_plot, s_terr_plot, blank_plot,
                       nrow = 6, ncol = 1, heights = c(0.1,1,1,1,1,0.1))
# Arrange together.
terr_both_plots <- ggarrange(terr_maps, terr_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))
ggsave("Plots/Maps/territory_maps.png", height = 20, width = 15, dpi = 900)






## Seasonal ##

# Quick seasonality map.
temp_plot <- ggplot_legend_raster(temp_seas, 10, "Temperature\nseasonality")
ggsave("Plots/Maps/seasonality_map.png", height = 7, width = 15, dpi = 1000)



###############################################################################
                      #### Add side plots for diet ####

library(ggnewscale)
# Create palettes for plotting.
light_colours <- c( "#77AD78","#7494EA", "#C98986", "#D8C77B")[2:3]
dark_colours <- c("#214F4B", "#05299E", "#8D0801","#A88A05")[2:3]





###############################################################################
            #### Add side plots for territoriality ####

# Territoriality maps.
pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Proportion of species") + 
  annotate("text", x = 0, y = 1.05, label = "b", size = 10, fontface = 2)
pri_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Proportion of species") + 
  annotate("text", x = 0, y = 1.05, label = "d", size = 10, fontface = 2)
sec_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Proportion of species") + 
  annotate("text", x = 0, y = 1.05, label = "f", size = 10, fontface = 2)
sec_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Proportion of species") + 
  annotate("text", x = 0, y = 1.05, label = "h", size = 10, fontface = 2)


terr_side_plots <- ggarrange(pri_terr_lat_plot + rremove("xlab") + rremove("x.text"), 
                             pri_yearterr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_terr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_yearterr_lat_plot,
                             #labels = c("b", "d", "f", "h"), 
                             nrow = 4, ncol = 1,
                             heights = c(1,1,1,1.2) #vjust = c(0.01),
                             #hjust = c(-4.6,-4.6,-8.75,-4.75), font.label = list(size = 32)
                             )
blank_plot <- ggplot() + geom_blank() + theme_classic() +
  theme(axis.line = element_blank())

terr_maps <- ggarrange(blank_plot, p_terr_plot,p_year_terr_plot, s_terr_plot, s_terr_plot, blank_plot,
                       labels = c("", "a", "c", "e", "g", ""), nrow = 6, ncol = 1,
                       heights = c(0.1,1,1,1,1,0.1), vjust = c(0.025),
                       hjust = c(-2.8), font.label = list(size = 32))

terr_both_plots <- ggarrange(terr_maps, terr_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# terr_lat_plots <- ggarrange(p_terr_plot, pri_terr_lat_plot + rremove("xlab"), 
#                             p_year_terr_plot, pri_yearterr_lat_plot + rremove("xlab"),
#                             s_terr_plot, sec_terr_lat_plot + rremove("xlab"), 
#                             s_year_terr_plot, sec_yearterr_lat_plot,
#                            labels = c("a", "b", "c", "d", "e", "f", "g", "h"), nrow = 4, ncol = 2, 
#                            widths = c(3,1.535), font.label = list(size = 30), 
#                            align = "hv", hjust = c(-4,-4.5,-4,-4.5,-4,-9,-4,-4.5)) + 
#   theme(plot.margin = margin(l = -1.5, unit = "cm"))
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


## ## ## ## End ## ## ## ## 


