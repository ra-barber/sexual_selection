###############################################################################
                       # Making map for figure 1 #
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


###############################################################################
                             #### Data ####

# Read in the world data.
world_continents <- read_sf("Data/World_Continents/World_Continents.shp")

# Start by creating an empty raster stack to store our data in.
raster_template <- raster(ncols=8640, nrows=3600, ymn = -60)

# Filter world data for specific continents and turn into rasters.
east <- world_continents %>% filter(CONTINENT %in% c("Asia", "Oceania"))
east_raster <- east %>% fasterize(raster_template) 
west <- world_continents %>% filter(CONTINENT %in% c("North America", "South America"))
west_raster <- west %>% fasterize(raster_template) 

# Read in elevation raster.
elevation <- raster("Data/Chelsa/wc2.1_30s_elev.tif")
plot(elevation)

library(terra)
# Resample so they're same extent / res
resampled_elev <- resample(rast(elevation), rast(raster_template))

# Mask.
# temp_seas <- mask(raster(temp_seas), species_raster)
# npp <- mask(raster(npp), species_raster)
# rainfall <- mask(raster(rainfall), species_raster)


# Crop elevation to east and west rasters.
resampled_elev <- resample(elevation, east_raster)

# Bin elevation to better show variation.
elevation_binned <- bin_raster(resampled_elev, nbins = 100)

###############################################################################
                           #### Section 3 ####
# Create raster data.
east_data <- as.data.frame(east_raster, xy=TRUE)
colnames(east_data) <- c("long", "lat", "values")
#east_data$long[east_data$long < 0] <- (180 - abs(east_data$long[east_data$long < 0])) + 180

west_data <- as.data.frame(west_raster, xy=TRUE)
colnames(west_data) <- c("long", "lat", "values")

# Crop east and west maps to correct size.
east_cropped <- east_data %>% filter(long > 100 & long < 160  & lat > -5 & lat < 50)
west_cropped <- west_data %>% filter(long < -40 & long > -100 & lat > -5 & lat < 50)

ggplot() + geom_tile(aes(x=long, y=lat, fill= as.factor(values)), data=east_cropped) + 
  scale_fill_manual(values = "forestgreen", na.value = "white") + 
  scale_colour_manual(values = "forestgreen", na.value = "white") + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
  geom_hline(yintercept = 45, linetype = "dashed", size = 2) +
  theme_classic(base_size = 18) + coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0, l = 0, b = 0, r = 0))
ggsave("Plots/Maps/east_world.png", height = 10, width = 10, dpi = 600)

# Plot with ggplot.
ggplot() + geom_tile(aes(x=long, y=lat, fill= as.factor(values)), data=west_cropped) + 
  scale_fill_manual(values = "forestgreen", na.value = "white") + 
  scale_colour_manual(values = "forestgreen", na.value = "white") + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
  geom_hline(yintercept = 45, linetype = "dashed", size = 2) +
  theme_classic(base_size = 18) + coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0, l = 0, b = 0, r = 0))

ggsave("Plots/Maps/west_world.png", height = 10, width = 10, dpi = 600)

###############################################################################
               #### Try using elevation instead ####


elev_data <- as.data.frame(elevation_binned, xy=TRUE)
colnames(elev_data) <- c("long", "lat", "values")
west_elev_cropped <- elev_data %>% filter(long < -40 & long > -100 & lat > -5 & lat < 50)
east_elev_cropped <- elev_data %>% filter(long > 95 & long < 155  & lat > -5 & lat < 50)

"#99D19C"
"#5B8C5A"
"#214E34"  # Really liked this one.
"#034C3C"
"#4B644A"
"#1A351D"
"#5C7457"
"#40531B"
"#618B4A"
"#AFBC88"
scale_fill_gradient(low = "#5B8C5A", high = "#214E34", na.value = "white")
scale_fill_gradient(low = "#618B4A", high = "#214E34", na.value = "white")
scale_fill_gradient(low = "#AFBC88", high = "#214E34", na.value = "white")

"#228B22"
"#254117"
"forestgreen"
"#437C17"

ggplot() + geom_tile(aes(x=long, y=lat, fill= values), data=east_elev_cropped) + 
  scale_fill_gradient(low = "forestgreen", high = "#437C17", na.value = "white") +
  scale_colour_gradient(low = "#AFBC88", high = "#214E34", na.value = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
  geom_hline(yintercept = 45, linetype = "dashed", size = 2) +
  theme_classic(base_size = 18) + coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0, l = 0, b = 0, r = 0))

ggsave("Plots/Maps/elev_east_world.png", height = 10, width = 10, dpi = 600)
  
ggplot() + geom_tile(aes(x=long, y=lat, fill= values), data=west_elev_cropped) + 
  scale_fill_gradient(low = "forestgreen", high = "#214E34", na.value = "white") +
  scale_colour_gradient(low = "#AFBC88", high = "#214E34", na.value = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
  geom_hline(yintercept = 45, linetype = "dashed", size = 2) +
  theme_classic(base_size = 18) + coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0, l = 0, b = 0, r = 0))

ggsave("Plots/Maps/elev_west_world.png", height = 10, width = 10, dpi = 600)




###############################################################################
                   #### Redo plots with colour and fill ####


# Plot with ggplot.
col_east <- ggplot(east_cropped) + 
  geom_tile(aes(x=long, y=lat, fill= as.factor(values), 
                colour = as.factor(values)), size = 0.1, data=east_cropped) + 
  scale_fill_manual(values = "forestgreen", na.value = "white") + 
  scale_colour_manual(values = "forestgreen", na.value = "white") + 
  theme_classic(base_size = 18) + coord_fixed() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0, l = 0, b = 0, r = 0))


ggsave("Plots/Maps/test_east_world.png", height = 10, width = 10, dpi = 600)

###############################################################################
                           #### Section 5 ####

###############################################################################
                           #### Section 6 ####


###############################################################################
                           #### Section 7 ####


###############################################################################
                           #### Section 8 ####

# Look Rob, you've had your fun with the sectioning. 
# They'll be no more sectioning today.


###############################################################################
                             #### END ####
###############################################################################


###############################################################################
              #### All the stuff I'm afraid to delete ####


