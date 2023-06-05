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

# Load in the trait da
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
#land <- mask(land, species_raster)          ## Could try this to remove grey areas by coastline.

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

sex_ranges$sexual_tri <- sex_ranges$sexual_score 
sex_ranges$sexual_tri[sex_ranges$sexual_tri < 3] <- 0
sex_ranges$sexual_tri[sex_ranges$sexual_tri == 3] <- 1
sex_ranges$sexual_tri[sex_ranges$sexual_tri == 4] <- 2

sex_tri_raster <- average_raster("sexual_tri")

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

# Make a reverse
sex_ranges$cert_reverse <- sex_ranges$sexual_certainty
sex_ranges$cert_reverse <- sex_ranges$cert_reverse*-1
sex_ranges$cert_reverse <- sex_ranges$cert_reverse+5
cert_raster_4 <- average_raster("cert_reverse")

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

#Read in bioclim data.
temp_seas <- raster("Data/Chelsa/CHELSA_bio4_1981-2010_V.2.1.tif")
npp <- raster("Data/Chelsa/CHELSA_npp_1981-2010_V.2.1.tif")
rainfall <- raster("Data/Chelsa/CHELSA_bio15_1981-2010_V.2.1.tif")

library(terra)
# Resample so they're same extent / res
temp_seas <- resample(rast(temp_seas), rast(species_raster))
npp <- resample(rast(npp), rast(species_raster))
rainfall <- resample(rast(rainfall), rast(species_raster))

# Mask.
temp_seas <- mask(raster(temp_seas), species_raster)
npp <- mask(raster(npp), species_raster)
rainfall <- mask(raster(rainfall), species_raster)
gc()

################################################################################
         ###### Maps of primary vs secondary consumers #########

# Primary and secondary sexual selection.
p_sex_raster <- sex_ranges %>% filter(trophic_binary == "Primary") %>% average_raster_2()
s_sex_raster <- sex_ranges %>% filter(trophic_binary == "Secondary") %>% average_raster_2()

# Primary and secondary sexual selection.
p_terr_raster <- sex_ranges %>% filter(trophic_binary == "Primary") %>% average_raster_2(var_name = "terr_dummy")
s_terr_raster <- sex_ranges %>% filter(trophic_binary == "Secondary") %>% average_raster_2(var_name = "terr_dummy")

# Primary and secondary sexual selection.
p_year_terr_raster <- sex_ranges %>% filter(trophic_binary == "Primary") %>% average_raster_2(var_name = "year_terr_dummy")
s_year_terr_raster <- sex_ranges %>% filter(trophic_binary == "Secondary") %>% average_raster_2(var_name = "year_terr_dummy")

gc()

################################################################################
              ###### Maps of primary consumer niches #########


fruit_sex_raster <- sex_ranges %>% filter(trophic_niche == "Frugivore") %>% average_raster_2()
nect_sex_raster <- sex_ranges %>% filter(trophic_niche == "Nectarivore") %>% average_raster_2()
gran_sex_raster <- sex_ranges %>% filter(trophic_niche == "Granivore") %>% average_raster_2()
herbivores <- c("Herbivore aquatic", "Herbivore terrestrial")
herb_sex_raster <- sex_ranges %>% filter(trophic_niche %in% herbivores) %>% average_raster_2()
gc()

# Without removing species.
nect_sex_raster <- sex_ranges %>% filter(trophic_niche == "Nectarivore") %>% average_raster_3()

plot(nect_sex_plot)

################################################################################
            ###### Maps of high certainty species #########


h_sex_raster <- sex_ranges %>% filter(sexual_certainty == 1) %>% average_raster_2()
m_sex_raster <- sex_ranges %>% filter(sexual_certainty < 3) %>% average_raster_2()
l_sex_raster <- sex_ranges %>% filter(sexual_certainty > 2) %>% average_raster_2()
gc()


################################################################################

# Pull out those socially monogamous species.
mono_ranges <- sex_ranges %>% filter(sexual_score < 3)

# Create species richness rasters.
mono_species_raster <- fasterize(mono_ranges, raster_template, fun = "sum")

# NA areas with less than 10 species.
mono_species_raster[which(getValues(mono_species_raster) < 10)] <- NA

# Mask for land.
mono_species_raster <- mask(mono_species_raster, land)

# Create average rasters.
mono_sex_raster <- average_raster("sexual_score", mono_ranges, mono_species_raster)

# Make a binary sexual selection map.
sex_bin_raster <- average_raster("sexual_binary")

# Plot new rasters.
plot(mono_sex_raster)
plot(sex_bin_raster)

# Remove big objects.
rm(mono_ranges)
rm(sex_ranges)
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

# Make inverse data certainty score for plotting. 
phylo_data$cert_inverse <- 1/phylo_data$sexual_certainty

# Make a reverse
phylo_data$cert_reverse <- phylo_data$sexual_certainty
phylo_data$cert_reverse <- phylo_data$cert_reverse*-1
phylo_data$cert_reverse <- phylo_data$cert_reverse+5

# Create tri sex score.
phylo_data$sexual_tri <- phylo_data$sexual_score 
phylo_data$sexual_tri[phylo_data$sexual_tri < 3] <- 0
phylo_data$sexual_tri[phylo_data$sexual_tri == 3] <- 1
phylo_data$sexual_tri[phylo_data$sexual_tri == 4] <- 2

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

# Create a dummy herbivore variable to group by.
phylo_data$herbivore_niche <- NA
phylo_data$herbivore_niche[phylo_data$trophic_niche %in% herbivores] <- "Herbivore"
phylo_data$herbivore_niche[is.na(phylo_data$herbivore_niche)] <- "Non-Herbivore"

# Group sex scores and certainty by lat bins.
lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()
cert_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(cert_dummy),
            trait_sd = sd(cert_dummy),
            trait_se = sd(cert_dummy)/sqrt(length(cert_dummy)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()

cert_lat_data_2 <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(cert_inverse),
            trait_sd = sd(cert_inverse),
            trait_se = sd(cert_inverse)/sqrt(length(cert_inverse)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()

cert_lat_data_4 <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(cert_reverse),
            trait_sd = sd(cert_reverse),
            trait_se = sd(cert_reverse)/sqrt(length(cert_reverse)),
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


# Group the data by trophic niche and bins.
niche_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_niche) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()

# Group the data by herbivores or not and bins.
herb_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, herbivore_niche) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()


# Do same for high and medium certainty.
hi_lat_data <- phylo_data %>% filter(binned_lat != 75 & sexual_certainty == 1) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()
med_lat_data <- phylo_data %>% filter(binned_lat != 75 & sexual_certainty < 3) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()

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

## Other metrics of sexual selection ##
mono_lat_data <- phylo_data %>% filter(sexual_score < 3) %>%  filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(sexual_score),
            trait_sd = sd(sexual_score),
            trait_se = sd(sexual_score)/sqrt(length(sexual_score)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()

sexbin_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(sexual_binary),
            trait_sd = sd(sexual_binary),
            trait_se = sd(sexual_binary)/sqrt(length(sexual_binary)),
            trait_max = trait_mean + trait_se,
            trait_min = trait_mean - trait_se) %>% na.omit()


sextri_lat_data <- phylo_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat) %>% 
  summarise(trait_mean = mean(sexual_tri),
            trait_sd = sd(sexual_tri),
            trait_se = sd(sexual_tri)/sqrt(length(sexual_tri)),
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


###############################################################################
         #### Testing normal certainty scores instead of dichotomised ####

# Dichotomised data certainty map.
cert_plot_4 <- ggplot_raster(cert_raster_4, 6, "Data certainty") + 
  annotate("text", x = 20, y = -48, label  = "Data certainty", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "c", size = 12, fontface = 2)

# Lat relationship with data certainty.
cert_lat_plot_4 <- cert_lat_data_4 %>% 
  lat_side_plot(ylabel = "Data certainty", ylimits = c(0,4), 
                ybreaks =  c(0,1,2,3,4), lab_ypos = 0.8) + 
  annotate("text", x = 0, y = 4, label = "d", size = 12, fontface = 2)

# Arrange side plots.
sex_side_plots_2 <- ggarrange(sex_lat_plot + rremove("xlab") + rremove("x.text"), 
                              cert_lat_plot_4, nrow = 2, ncol = 1, heights = c(1,1.2))
# Arrange maps.
sex_maps_2 <- ggarrange(blank_plot, sex_plot, cert_plot_4, blank_plot, nrow = 4, ncol = 1,
                      heights = c(0.1,1,1,0.1))

# Put all together.
sex_both_plots_2 <- ggarrange(sex_maps_2, sex_side_plots_2, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Maps/alt_figure_2.png", height = 10, width = 15, dpi = 600)






###############################################################################
               #### Plot trophic level figure 4 ####


# Plot sexual selection maps for both trophic levels.
p_sex_plot <- panel_ggplot_raster(p_sex_raster, plot_title = "Primary consumers", plot_label = "a")
s_sex_plot <- panel_ggplot_raster(s_sex_raster, plot_title = "Secondary consumers", plot_label = "c")

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
                ybreaks =  c(0,1.0, 2.0), lab_x_pos = 10, lab_ypos = 1.5) + 
  annotate("text", x = 0, y = 2, label = "b", size = 12, fontface = 2)
sec_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
                #ybreaks =  c(0,0.5,1.0), 
                ybreaks =  c(0,1.0,2.0), 
                lab_x_pos = 10, lab_ypos = 1.5) + 
  annotate("text", x = 0, y = 2, label = "d", size = 12, fontface = 2)

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
         #### Plot trophic niche maps for primary consumers ####


# Create maps of the four niches.
fruit_sex_plot <- panel_ggplot_raster(fruit_sex_raster, plot_title = "Frugivores", plot_label = "a")
fruit_lat_plot <-  niche_lat_data %>% filter(trophic_niche == "Frugivore") %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,4), 
                ybreaks =  c(0,2, 4), lab_x_pos = 10, lab_ypos = 3.5, plot_label = "b")

nect_sex_plot <- panel_ggplot_raster(nect_sex_raster, plot_title = "Nectarivores", plot_label = "c")
nect_lat_plot <-  niche_lat_data %>% filter(trophic_niche == "Nectarivore") %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,4), 
                ybreaks =  c(0,2, 4), lab_x_pos = 10, lab_ypos = 3.5, plot_label = "d")

gran_sex_plot <- panel_ggplot_raster(gran_sex_raster, plot_title = "Granivores", plot_label = "e")
gran_lat_plot <-  niche_lat_data %>% filter(trophic_niche == "Granivore") %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,4), 
                ybreaks =  c(0,2, 4), lab_x_pos = 10, lab_ypos = 3.5, plot_label = "f")

herb_sex_plot <- panel_ggplot_raster(herb_sex_raster, plot_title = "Herbivores", plot_label = "g")
herb_lat_plot <-  herb_lat_data %>% filter(herbivore_niche == "Herbivore") %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,4), 
                ybreaks =  c(0,2, 4), lab_x_pos = 10, lab_ypos = 3.5, plot_label = "h")
gc()

# Arrange side plots and maps.
niche_side_plots <- ggarrange(fruit_lat_plot + rremove("xlab") + rremove("x.text"), 
                              nect_lat_plot + rremove("xlab") + rremove("x.text"),
                              gran_lat_plot + rremove("xlab") + rremove("x.text"),
                              herb_lat_plot,
                             nrow = 4, ncol = 1, heights = c(1,1,1,1.2))
niche_maps <- ggarrange(blank_plot, fruit_sex_plot, nect_sex_plot, 
                        gran_sex_plot, herb_sex_plot, blank_plot,
                       nrow = 6, ncol = 1, heights = c(0.1,1,1,1,1,0.1))
# Arrange together.
niche_both_plots <- ggarrange(niche_maps, niche_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))
ggsave("Plots/Maps/niche_maps.png", height = 20, width = 15, dpi = 900)


# frug_both_plots <- ggarrange(fruit_sex_plot, fruit_lat_plot, ncol = 2, widths = c(3,1.535)) +
#   theme(plot.margin = margin(l = -0.3, unit = "cm"))
# 
# # Export figure.
# ggsave("Plots/Maps/frugivore_map.png", height = 5, width = 15, dpi = 600)



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
                ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.9) + 
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
p_year_terr_plot <- ggplot_terr_raster(p_year_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Primary consumers", size = 8, fontface = 2)+
  annotate("text", x =  -170, y = 80, label  = "c", size = 12, fontface = 2)
s_terr_plot <- ggplot_raster(s_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Secondary consumers", size = 8, fontface = 2)+
  annotate("text", x =  -170, y = 80, label  = "e", size = 12, fontface = 2)
s_year_terr_plot <- ggplot_raster(s_year_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = "Secondary consumers", size = 8, fontface = 2)+
  annotate("text", x =  -170, y = 80, label  = "g", size = 12, fontface = 2)

# Territoriality side plots.
pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Proportion of species", lab_ypos = 0.9) + 
  annotate("text", x = 0, y = 1.05, label = "b", size = 12, fontface = 2) + 
  annotate("text", x = 35, y = 1.05, label = "Territorial", size = 7, fontface = 2)
pri_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Proportion of species", lab_ypos = 0.9) + 
  annotate("text", x = 0, y = 1.05, label = "d", size = 12, fontface = 2) + 
  annotate("text", x = 35, y = 1.05, label = "Year-round territorial", size = 7, fontface = 2)
sec_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Proportion of species", lab_ypos = 0.9) + 
  annotate("text", x = 0, y = 1.05, label = "f", size = 12, fontface = 2) + 
  annotate("text", x = 35, y = 1.05, label = "Territorial", size = 7, fontface = 2)
sec_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Proportion of species", lab_ypos = 0.9) + 
  annotate("text", x = 0, y = 1.05, label = "h", size = 12, fontface = 2) + 
  annotate("text", x = 35, y = 1.05, label = "Year-round territorial", size = 7, fontface = 2)

# Arrange side plots and maps.
terr_side_plots <- ggarrange(pri_terr_lat_plot + rremove("xlab") + rremove("x.text"), 
                             pri_yearterr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_terr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_yearterr_lat_plot,
                             nrow = 4, ncol = 1, heights = c(1,1,1,1.2))
terr_maps <- ggarrange(blank_plot, p_terr_plot,p_year_terr_plot, 
                       s_terr_plot, s_year_terr_plot, blank_plot,
                       nrow = 6, ncol = 1, heights = c(0.1,1,1,1,1,0.1))
# Arrange together.
terr_both_plots <- ggarrange(terr_maps, terr_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))
ggsave("Plots/Maps/territory_maps.png", height = 20, width = 15, dpi = 900)

# Try different axis label.
pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Proporition territorial", lab_ypos = 0.9) + 
  annotate("text", x = 0, y = 1.05, label = "b", size = 12, fontface = 2)
pri_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Proportion year-round territorial", lab_ypos = 0.9) + 
  annotate("text", x = 0, y = 1.05, label = "d", size = 12, fontface = 2)
sec_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Proporition territorial", lab_ypos = 0.9) + 
  annotate("text", x = 0, y = 1.05, label = "f", size = 12, fontface = 2)
sec_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Proportion year-round territorial", lab_ypos = 0.9) + 
  annotate("text", x = 0, y = 1.05, label = "h", size = 12, fontface = 2)

# Arrange side plots and maps.
terr_side_plots <- ggarrange(pri_terr_lat_plot + rremove("xlab") + rremove("x.text"), 
                             pri_yearterr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_terr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_yearterr_lat_plot,
                             nrow = 4, ncol = 1, heights = c(1,1,1,1.2))
terr_both_plots <- ggarrange(terr_maps, terr_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))
ggsave("Plots/Maps/territory_maps_2.png", height = 20, width = 15, dpi = 900)

# Try different axis label.
pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Prop. terr.") + 
  annotate("text", x = 0, y = 1.05, label = "b", size = 12, fontface = 2)
pri_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  lat_side_plot(ylabel = "Prop. year-round terr.") + 
  annotate("text", x = 0, y = 1.05, label = "d", size = 12, fontface = 2)
sec_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Prop. terr.") + 
  annotate("text", x = 0, y = 1.05, label = "f", size = 12, fontface = 2)
sec_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  lat_side_plot(ylabel = "Prop. year-round terr.") + 
  annotate("text", x = 0, y = 1.05, label = "h", size = 12, fontface = 2)

ggsave("Plots/Maps/territory_maps_3.png", height = 20, width = 15, dpi = 900)

################################################################################
              #### Plot sexual selection alternatives #####

# Show high certainty sexual selection maps.
mono_sex_plot <- ggplot_raster(mono_sex_raster, 6, "Sexual selection") +
  annotate("text", x = 20, y = -48, label  = "EPP data", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "a", size = 12, fontface = 2)
bin_sex_plot <- ggplot_raster(sex_bin_raster, 6, "Sexual selection") + 
  annotate("text", x = 30, y = -48, label  = "Binary sex scores", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "c", size = 12, fontface = 2)
# And side plots.
mono_sex_lat_plot <- mono_lat_data  %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
                ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.9) + 
  annotate("text", x = 0, y = 2, label = "b", size = 12, fontface = 2)
bin_sex_lat_plot <- sexbin_lat_data  %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1), 
                ybreaks =  c(0,0.5,1.0), lab_ypos = 0.3) + 
  annotate("text", x = 0, y = 1, label = "d", size = 12, fontface = 2)

# Arrange maps together.
alt_side_plots <- ggarrange(mono_sex_lat_plot + rremove("xlab") + rremove("x.text"), 
                             bin_sex_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
alt_maps <- ggarrange(blank_plot, mono_sex_plot, bin_sex_plot, blank_plot, nrow = 4, ncol = 1,
                       heights = c(0.1,1,1,0.1))
alt_both_plots <- ggarrange(alt_maps, alt_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Maps/alternative_sexual_groupings.png", height = 10, width = 15, dpi = 600)


# Do maps grouping 0,1,2 together.
tri_sex_plot <- ggplot_raster(sex_tri_raster, 6, "Sexual selection") + 
  annotate("text", x = 30, y = -48, label  = "Tertiary sex scores", size = 8, fontface = 2) +
  annotate("text", x = -170, y = 80, label  = "c", size = 12, fontface = 2)
tri_sex_lat_plot <- sextri_lat_data  %>% 
  lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1), 
                ybreaks =  c(0,0.5,1.0), lab_ypos = 0.3) + 
  annotate("text", x = 0, y = 1, label = "d", size = 12, fontface = 2)

# Arrange maps together.
alt_side_plots_2 <- ggarrange(mono_sex_lat_plot + rremove("xlab") + rremove("x.text"), 
                              tri_sex_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
alt_maps_2 <- ggarrange(blank_plot, mono_sex_plot, tri_sex_plot, blank_plot, nrow = 4, ncol = 1,
                      heights = c(0.1,1,1,0.1))
alt_both_plots_2 <- ggarrange(alt_maps_2, alt_side_plots_2, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

ggsave("Plots/Maps/alternative_sexual_groupings_2.png", height = 10, width = 15, dpi = 600)


gc()




################################################################################ 
                      #### Test linear model fits #####

# Prepare datasets.
primary_data <- diet_lat_data %>% filter(trophic_binary == "Primary")
secondary_data <- diet_lat_data %>% filter(trophic_binary == "Secondary")


# Run models.
primary_lat_model <- lm(trait_mean ~ binned_lat, data = primary_data)
secondary_lat_model <- lm(trait_mean ~ binned_lat, data = secondary_data)



summary(primary_lat_model)
summary(secondary_lat_model)

plot(primary_lat_model)




## Seasonal ##
# 
# # Quick seasonality map.
temp_plot <- ggplot_raster(temp_seas, 10, "Temperature\nseasonality")
ggsave("Plots/Maps/seasonality_map.png", height = 7, width = 15, dpi = 1000)

# Rainfall stability map.
rain_plot <- ggplot_raster(rainfall, 10, "Rainfall\nstability")
ggsave("Plots/Maps/rainfall_map.png", height = 7, width = 15, dpi = 1000)




## ## ## ## End ## ## ## ## 

# Code for using patchwork to assemble figures.
library(patchwork)

# Create a text based layout.
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
# Add the figures using the layout to show.
figure_2 <- wrap_plots(A = sex_plot, B = sex_lat_plot, C = cert_plot, D = cert_lat_plot, design = layout) 
ggsave("Plots/Maps/figure_2_patch.png", height = 10, width = 15, dpi = 600)

