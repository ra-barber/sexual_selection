###############################################################################
                   # Making sexual selection maps #
###############################################################################

# This script makes maps of sexual selection scores, accounting for different 
# ecological groups and data certainties.

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

# Load in the trait data.
model_data <- read_ss_data()

# Join up the data.
sex_ranges <- left_join(model_data, new_jetz_ranges,  
                        by = c("scientific_name_bird_tree" = "SCINAME"))
sex_ranges <- st_as_sf(sex_ranges)
rm(new_jetz_ranges)

# Start by creating an empty raster stack to store our data in. 
# (This creates 5km res, which is only used for visualisation)
raster_template <- raster(ncols=8640, nrows=3600, ymn = -60)
# raster_template <- raster(ncols=8640, nrows=3600, ymn = -60, crs = behrmannCRS)
# 
# behrmannCRS <- CRS('+proj=cea +lat_ts=30')



# Code for figuring out the resolution.
# library(terra)
# distance(cbind(0,-21), cbind(0,-21.04166667), lonlat=TRUE)

# Create a raster of landmass to mask maps.
data(wrld_simpl)
land <- wrld_simpl %>% st_as_sf() %>% fasterize(raster_template)


################################################################################
            #### Make rasters for ecological partitions ####


# Figure 2.
sex_raster <- average_raster()
cert_raster <- average_raster(var_name = "data_certainty")

# Figure 3.
p_sex_raster <- sex_ranges %>% filter(trophic_level_binary == "Primary") %>% average_raster()
s_sex_raster <- sex_ranges %>% filter(trophic_level_binary == "Secondary") %>% average_raster()
fruit_sex_raster <- sex_ranges %>% filter(trophic_niche == "Frugivore") %>% average_raster()
invert_sex_raster <- sex_ranges %>% filter(trophic_niche == "Invertivore") %>% average_raster()

# Extended data figure 3.
h_sex_raster <- sex_ranges %>% filter(data_certainty == 4) %>% average_raster()
m_sex_raster <- sex_ranges %>% filter(data_certainty > 2) %>% average_raster()

# Extended data figure 4.
mig_sex_raster <- sex_ranges %>% filter(migration_binary == "Strong") %>% average_raster()
no_mig_sex_raster <- sex_ranges %>% filter(migration_binary == "Weak") %>% average_raster()

# Extended data figure 5.
terr_sex_raster <- sex_ranges %>% filter(territoriality_binary == "Territorial") %>% average_raster()
no_terr_sex_raster <- sex_ranges %>% filter(territoriality_binary == "Non-territorial") %>% average_raster()

# Extended data figure 8.
p_terr_raster <- sex_ranges %>% filter(trophic_level_binary == "Primary") %>% average_raster(var_name = "terr_dummy")
s_terr_raster <- sex_ranges %>% filter(trophic_level_binary == "Secondary") %>% average_raster(var_name = "terr_dummy")
p_year_terr_raster <- sex_ranges %>% filter(trophic_level_binary == "Primary") %>% average_raster(var_name = "year_terr_dummy")
s_year_terr_raster <- sex_ranges %>% filter(trophic_level_binary == "Secondary") %>% average_raster(var_name = "year_terr_dummy")

# Figure S2.
h_spec_raster <- sex_ranges %>% filter(data_certainty == 4) %>% spec_raster_func()

# Remove big objects.
rm(sex_ranges)
gc()


###############################################################################
               #### Plot sexual selection figure 2 ####


# Load in the side plots.
#load("Plots/Maps/latitudinal_sideplots.Rdata")
#load("Plots/Maps/relative_size_latitudinal_sideplots.Rdata")
load("Plots/Maps/purple_relative_remove_10_latitudinal_sideplots.Rdata")


# Create data frame of land to use to crop ggplot maps.
land_data <- as.data.frame(land, xy=TRUE)

# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

# Blank plot.
blank_plot <- ggplot() + geom_blank() + theme_classic() +
  theme(axis.line = element_blank())

# Full sexual selection map.
sex_plot <- panel_ggplot_raster(sex_raster, plot_title = "All birds", plot_label = "a")
cert_plot <- panel_ggplot_raster(cert_raster, plot_title = "Data certainty", plot_label = "c")

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
ggsave("Figures_cell_models/figure_2.png", height = 10, width = 15, dpi = 600, bg="white")
ggsave("Figures_cell_models/figure_2.pdf", height = 10, width = 15, dpi = 600, bg="white")
gc()

# Try removing excess objects to free up space.
rm(sex_plot, cert_plot, sex_maps, sex_side_plots, sex_both_plots)
#rm(sex_raster, cert_raster)


###############################################################################
                        #### Figure 3  ####


# Create maps of the four niches.
p_sex_plot <- panel_ggplot_raster(p_sex_raster, plot_title = expression(bold("1"^ry*" consumers (n = 2753)")), plot_label = "a")

fruit_sex_plot <- panel_ggplot_raster(fruit_sex_raster, plot_title = "Frugivores (n = 1025)", plot_label = "c") 

# fruit_sex_plot <- ggplot_terr_raster(fruit_sex_raster) +
#   annotate("text", x = 20, y = -48, label  = "Frugivores (n = 1025)", size = 8, fontface = 2)+
#   annotate("text", x =  -170, y = 80, label  = "c", size = 12, fontface = 2)

s_sex_plot <- panel_ggplot_raster(s_sex_raster, plot_title = expression(bold("2"^ry*" consumers (n = 7083)")), plot_label = "e")
invert_sex_plot <- panel_ggplot_raster(invert_sex_raster, plot_title = "Invertivores (n = 4694)", plot_label = "g")

# Arrange side plots and maps.
niche_side_plots <- ggarrange(pri_lat_plot + rremove("xlab") + rremove("x.text"), 
                              fruit_lat_plot + rremove("xlab") + rremove("x.text"),
                              sec_lat_plot + rremove("xlab") + rremove("x.text"),
                              invert_lat_plot,
                             nrow = 4, ncol = 1, heights = c(1,1,1,1.2))
niche_maps <- ggarrange(blank_plot, p_sex_plot, fruit_sex_plot, 
                        s_sex_plot, invert_sex_plot, blank_plot,
                       nrow = 6, ncol = 1, heights = c(0.1,1,1,1,1,0.1))
niche_both_plots <- ggarrange(niche_maps, niche_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export.
ggsave("Figures_cell_models/figure_3.png", height = 20, width = 15, dpi = 900, bg="white")
ggsave("Figures_cell_models/figure_3.pdf", height = 20, width = 15, dpi = 900, bg="white")


###############################################################################
                 #### Extended data figure 3 ####


# Show high certainty sexual selection maps.
m_sex_plot <- panel_ggplot_raster(m_sex_raster, plot_title = "Certainty 3 & 4 (n = 7592)", plot_label = "a")
h_sex_plot <- panel_ggplot_raster(h_sex_raster, plot_title = "Certainty 4 (n = 2851)", plot_label = "c")

# Arrange side plots and maps.
cert_side_plots <- ggarrange(med_sex_lat_plot + rremove("xlab") + rremove("x.text"), 
                             hi_sex_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
cert_maps <- ggarrange(blank_plot, m_sex_plot, h_sex_plot, blank_plot, nrow = 4, ncol = 1,
                       heights = c(0.1,1,1,0.1))
cert_both_plots <- ggarrange(cert_maps, cert_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Figures_cell_models/figure_ED_3.png", height = 10, width = 15, dpi = 600, bg="white")   
ggsave("Figures_cell_models/figure_ED_3.pdf", height = 10, width = 15, dpi = 600, bg="white")


###############################################################################
               #### Extended data figure 4 & 5 ####


# Show high certainty sexual selection maps.
mig_sex_plot <- panel_ggplot_raster(mig_sex_raster, plot_title = "Migrants (n = 901)", plot_label = "a")
no_mig_sex_plot <- panel_ggplot_raster(no_mig_sex_raster, plot_title = "Non-migrants (n = 8935)", plot_label = "c")

# Arrange maps together.
mig_side_plots <- ggarrange(mig_lat_plot + rremove("xlab") + rremove("x.text"), 
                            no_mig_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
mig_maps <- ggarrange(blank_plot, mig_sex_plot, no_mig_sex_plot, blank_plot, nrow = 4, ncol = 1,
                      heights = c(0.1,1,1,0.1))
mig_both_plots <- ggarrange(mig_maps, mig_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Figures_cell_models/figure_ED_4.png", height = 10, width = 15, dpi = 600, bg="white")
ggsave("Figures_cell_models/figure_ED_4.pdf", height = 10, width = 15, dpi = 600, bg="white")

# Show high certainty sexual selection maps.
terr_sex_plot <- panel_ggplot_raster(terr_sex_raster, plot_title = "Territorial (n = 7261)", plot_label = "a")
no_terr_sex_plot <- panel_ggplot_raster(no_terr_sex_raster, plot_title = "Non-territorial (n = 2575)", plot_label = "c")

# Arrange maps together.
terr_side_plots <- ggarrange(terr_lat_plot + rremove("xlab") + rremove("x.text"), 
                             no_terr_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
terr_maps <- ggarrange(blank_plot, terr_sex_plot, no_terr_sex_plot, blank_plot, nrow = 4, ncol = 1,
                       heights = c(0.1,1,1,0.1))
terr_both_plots <- ggarrange(terr_maps, terr_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Figures_cell_models/figure_ED_5.png", height = 10, width = 15, dpi = 600, bg="white")
ggsave("Figures_cell_models/figure_ED_5.pdf", height = 10, width = 15, dpi = 600, bg="white")


###############################################################################
                  #### Extended data figure 8 ####

# Territoriality maps.
p_terr_plot <- panel_ggplot_raster(p_terr_raster, plot_title = expression(bold("1"^ry*" consumers (n = 1379)")), plot_label = "a")
p_year_terr_plot <- ggplot_terr_raster(p_year_terr_raster, 6, "") +   # Has a custom function due to lack of variation meaning it needs custom breaks.
  annotate("text", x = 20, y = -48, label  = expression(bold("1"^ry*" consumers (n = 333)")), size = 8, fontface = 2)+
  annotate("text", x =  -170, y = 80, label  = "c", size = 12, fontface = 2)
s_terr_plot <- panel_ggplot_raster(s_terr_raster, plot_title = expression(bold("2"^ry*" consumers (n = 5882)")), plot_label = "e")
s_year_terr_plot <- panel_ggplot_raster(s_year_terr_raster, plot_title = expression(bold("2"^ry*" consumers (n = 2640)")), plot_label = "g")

# Arrange side plots and maps.
terr_side_plots <- ggarrange(pri_terr_lat_plot + rremove("xlab") + rremove("x.text"), 
                             pri_yearterr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_terr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_yearterr_lat_plot,
                             nrow = 4, ncol = 1, heights = c(1,1,1,1.2))
terr_maps <- ggarrange(blank_plot, p_terr_plot, p_year_terr_plot, 
                       s_terr_plot, s_year_terr_plot, blank_plot,
                       nrow = 6, ncol = 1, heights = c(0.1,1,1,1,1,0.1))
terr_both_plots <- ggarrange(terr_maps, terr_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export.
#ggsave("Figures_relative_error/figure_ED_8.png", height = 20, width = 15, dpi = 900, bg="white")
#ggsave("Figures_relative_error/figure_ED_8.pdf", height = 20, width = 15, dpi = 900, bg="white")

ggsave("Figures_cell_models/test_figure_ED_8.png", height = 20, width = 15, dpi = 600, bg="white")
ggsave("Figures_cell_models/test_figure_ED_8.pdf", height = 20, width = 15, dpi = 600, bg="white")

rm(sex_ranges)

###############################################################################
                   ##### Figure S2 ######


# Show high certainty sexual selection maps.
#m_spec_plot <- panel_ggplot_raster(m_spec_raster, plot_title = "Certainty 3 & 4 (n = 7592)", plot_label = "a")
h_spec_plot <- panel_ggplot_raster(h_spec_raster, plot_title = "Certainty 4 (n = 2851)", plot_label = "") +
  theme(plot.margin = margin(l = -0.3, r =0.3, unit = "cm"))

# Export figure.
ggsave("Figures_relative_error/figure_S2.png", height = 5, width = 10, dpi = 600)


################################################################################
                        ###### End ######
################################################################################