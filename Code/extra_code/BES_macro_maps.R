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
gc()

# Functions.
source("Code/functions.R")

################################################################################
                           #### Data ####

# Load in the range maps.
load("../../Phd_data/Spatial/clean_jetz_ranges.RData")

# Load in the trait da
phylo_data <- read.csv("Data/sexual_traits.csv")

# Join up the data.
sex_ranges <- left_join(phylo_data, new_jetz_ranges,  
                        by = c("birdtree_name" = "SCINAME"))
sex_ranges <- st_as_sf(sex_ranges)
rm(new_jetz_ranges)


################################################################################
                  #### Species richness raster ####

# Start by creating an empty raster stack to store our data in.
raster_template <- raster(ncols=8640, nrows=3600, ymn = -60)
#raster_template <- raster(ncols=17280, nrows=7200, ymn = -60)

# Code for figuring out the resolution.

# library(terra)
# distance(cbind(0,-21), cbind(0,-21.04166667), lonlat=TRUE)

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
cert_raster <- average_raster("cert_reverse")  # Only one we actually need.

# # Primary and secondary sexual selection.
p_sex_raster <- sex_ranges %>% filter(trophic_binary == "Primary") %>% average_raster_2()
s_sex_raster <- sex_ranges %>% filter(trophic_binary == "Secondary") %>% average_raster_2()

fruit_sex_raster <- sex_ranges %>% filter(trophic_niche == "Frugivore") %>% average_raster_2()
invert_sex_raster <- sex_ranges %>% filter(trophic_niche == "Invertivore") %>% average_raster_2()
gran_sex_raster <- sex_ranges %>% filter(trophic_niche == "Granivore") %>% average_raster_2()

h_sex_raster <- sex_ranges %>% filter(sexual_certainty == 1) %>% average_raster_2()
m_sex_raster <- sex_ranges %>% filter(sexual_certainty < 3) %>% average_raster_2()

# Remove big objects.
rm(sex_ranges)
gc()

###############################################################################
               #### Plot sexual selection figure 2 ####

# Load in the side plots.
load("Plots/Presentation/presentation_lat_sideplots.Rdata")

# Create data frame of land to use to crop ggplot maps.
land_data <- as.data.frame(land, xy=TRUE)

# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

# Blank plot.
blank_plot <- ggplot() + geom_blank() + theme_classic() +
  theme(axis.line = element_blank())

# Full sexual selection map.
sex_plot <- panel_ggplot_raster(sex_raster, plot_title = "All birds (n = 9836)", plot_label = "")
sex_both_plots <- ggarrange(sex_plot, sex_lat_plot, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Presentation/all_birds_map.png", height = 5, width = 15, dpi = 600)

# Certainty.
cert_plot <- panel_ggplot_raster(cert_raster, plot_title = "Data certainty (n = 9836)", plot_label = "")
cert_both_plots <- ggarrange(cert_plot, cert_lat_plot, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))
ggsave("Plots/Presentation/cert_map.png", height = 5, width = 15, dpi = 600)



###############################################################################
            #### Plot trophic niche maps for new figure 4  ####


# Create maps of the four niches.
p_sex_plot <- panel_ggplot_raster(p_sex_raster, plot_title = expression(bold("1"^ry*" consumers (n = 2753)")), plot_label = "")
s_sex_plot <- panel_ggplot_raster(s_sex_raster, plot_title = expression(bold("2"^ry*" consumers (n = 7083)")), plot_label = "")

# Arrange maps together.
trophic_side_plots <- ggarrange(pri_lat_plot + rremove("xlab") + rremove("x.text"), 
                                sec_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
trophic_maps <- ggarrange(blank_plot, p_sex_plot, s_sex_plot, blank_plot, nrow = 4, ncol = 1,
                       heights = c(0.1,1,1,0.1))
trophic_both_plots <- ggarrange(trophic_maps, trophic_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

ggsave("Plots/Presentation/trophic_maps.png", height = 10, width = 15, dpi = 600)


fruit_sex_plot <- panel_ggplot_raster(fruit_sex_raster, plot_title = "Frugivores (n = 1025)", plot_label = "")
invert_sex_plot <- panel_ggplot_raster(invert_sex_raster, plot_title = "Invertivores (n = 4694)", plot_label = "")

niche_side_plots <- ggarrange(fruit_lat_plot + rremove("xlab") + rremove("x.text"), 
                              invert_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
niche_maps <- ggarrange(blank_plot, fruit_sex_plot, invert_sex_plot, blank_plot, nrow = 4, ncol = 1,
                          heights = c(0.1,1,1,0.1))
niche_both_plots <- ggarrange(niche_maps, niche_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

ggsave("Plots/Presentation/niche_maps.png", height = 10, width = 15, dpi = 600)



###############################################################################
          #### Plot data certainty supplementary figure ####

# Show high certainty sexual selection maps.
m_sex_plot <- panel_ggplot_raster(m_sex_raster, plot_title = "Certainty 3 & 4 (n = 7700)", plot_label = "")
h_sex_plot <- panel_ggplot_raster(h_sex_raster, plot_title = "Certainty 4 (n = 2851)", plot_label = "")

# Arrange maps together.
cert_side_plots <- ggarrange(med_sex_lat_plot + rremove("xlab") + rremove("x.text"), 
                             hi_sex_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
cert_maps <- ggarrange(blank_plot, m_sex_plot, h_sex_plot, blank_plot, nrow = 4, ncol = 1,
                       heights = c(0.1,1,1,0.1))
cert_both_plots <- ggarrange(cert_maps, cert_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Presentation/sensitivity_maps.png", height = 10, width = 15, dpi = 600)

phylo_data %>% count(trophic_niche)

gran_sex_plot <- panel_ggplot_raster(gran_sex_raster, plot_title = "Granivores (n = 653)", plot_label = "")
niche_maps <- ggarrange(blank_plot, fruit_sex_plot, gran_sex_plot, blank_plot, nrow = 4, ncol = 1,
                        heights = c(0.1,1,1,0.1))
niche_side_plots <- ggarrange(fruit_lat_plot + rremove("xlab") + rremove("x.text"), 
                              gran_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))

niche_both_plots <- ggarrange(niche_maps, niche_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

ggsave("Plots/Presentation/gran_niche_maps.png", height = 10, width = 15, dpi = 600)

