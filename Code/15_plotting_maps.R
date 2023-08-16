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

### should try using this function instead of the previous code to make a species richness map.
# spec_raster_func()



################################################################################
                #### Make rasters of life history traits ####

# Create sexual selection ranges with trait data.
sex_raster <- average_raster("sexual_score")
cert_raster <- average_raster("cert_reverse")  # Only one we actually need.
gc()


################################################################################
         ###### Maps of primary vs secondary consumers #########

# # Primary and secondary sexual selection.
p_sex_raster <- sex_ranges %>% filter(trophic_binary == "Primary") %>% average_raster_2()
s_sex_raster <- sex_ranges %>% filter(trophic_binary == "Secondary") %>% average_raster_2()

# # Primary and secondary sexual selection.
p_terr_raster <- sex_ranges %>% filter(trophic_binary == "Primary") %>% average_raster_2(var_name = "terr_dummy")
s_terr_raster <- sex_ranges %>% filter(trophic_binary == "Secondary") %>% average_raster_2(var_name = "terr_dummy")

# Primary and secondary sexual selection.
p_year_terr_raster <- sex_ranges %>% filter(trophic_binary == "Primary") %>% average_raster_2(var_name = "year_terr_dummy")
s_year_terr_raster <- sex_ranges %>% filter(trophic_binary == "Secondary") %>% average_raster_2(var_name = "year_terr_dummy")
gc()

################################################################################
              ###### Maps of primary consumer niches #########


fruit_sex_raster <- sex_ranges %>% filter(trophic_niche == "Frugivore") %>% average_raster_2()
invert_sex_raster <- sex_ranges %>% filter(trophic_niche == "Invertivore") %>% average_raster_2()
gc()



################################################################################
            ###### Maps of high certainty species #########


h_sex_raster <- sex_ranges %>% filter(sexual_certainty == 1) %>% average_raster_2()
m_sex_raster <- sex_ranges %>% filter(sexual_certainty < 3) %>% average_raster_2()
gc()


################################################################################
          #### Maps of species richness for high certainty ####

h_spec_raster <- sex_ranges %>% filter(sexual_certainty == 1) %>% spec_raster_func()
m_spec_raster <- sex_ranges %>% filter(sexual_certainty < 3) %>% spec_raster_func()

################################################################################
       ###### Maps of migration and territorialty #########


mig_sex_raster <- sex_ranges %>% filter(migration_binary == "Strong") %>% average_raster_2()
no_mig_sex_raster <- sex_ranges %>% filter(migration_binary == "Weak") %>% average_raster_2()

terr_sex_raster <- sex_ranges %>% filter(territory_binary == "Territory") %>% average_raster_2()
no_terr_sex_raster <- sex_ranges %>% filter(territory_binary == "No territory") %>% average_raster_2()
gc()

# Remove big objects.
rm(sex_ranges)
gc()

###############################################################################
               #### Plot sexual selection figure 2 ####

# Load in the side plots.
load("Plots/Maps/latitudinal_sideplots.Rdata")

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
ggsave("Plots/Maps/figure_2.png", height = 10, width = 15, dpi = 600)
ggsave("Plots/Maps/figure_2.pdf", height = 10, width = 15, dpi = 600)

# rm(sex_raster)
# rm(cert_raster)

###############################################################################
            #### Plot trophic niche maps for new figure 4  ####


# Create maps of the four niches.
p_sex_plot <- panel_ggplot_raster(p_sex_raster, plot_title = expression(bold("1"^ry*" consumers (n = 2753)")), plot_label = "a")
fruit_sex_plot <- panel_ggplot_raster(fruit_sex_raster, plot_title = "Frugivores (n = 1025)", plot_label = "c")
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
# Arrange together.
niche_both_plots <- ggarrange(niche_maps, niche_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))
ggsave("Plots/Maps/figure_3.png", height = 20, width = 15, dpi = 900)
ggsave("Plots/Maps/figure_3.pdf", height = 20, width = 15, dpi = 900)



###############################################################################
          #### Plot data certainty supplementary figure ####

# Show high certainty sexual selection maps.
m_sex_plot <- panel_ggplot_raster(m_sex_raster, plot_title = "Certainty 3 & 4 (n = 7592)", plot_label = "a")
h_sex_plot <- panel_ggplot_raster(h_sex_raster, plot_title = "Certainty 4 (n = 2851)", plot_label = "c")

# Arrange maps together.
cert_side_plots <- ggarrange(med_sex_lat_plot + rremove("xlab") + rremove("x.text"), 
                             hi_sex_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
cert_maps <- ggarrange(blank_plot, m_sex_plot, h_sex_plot, blank_plot, nrow = 4, ncol = 1,
                       heights = c(0.1,1,1,0.1))
cert_both_plots <- ggarrange(cert_maps, cert_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Maps/figure_ED_2.png", height = 10, width = 15, dpi = 600)
ggsave("Plots/Maps/figure_ED_2.pdf", height = 10, width = 15, dpi = 600)




###############################################################################
        #### Plot latitudinal gradient in territoriality ####

# Territoriality maps.
p_terr_plot <- panel_ggplot_raster(p_terr_raster, plot_title = expression(bold("1"^ry*" consumers")), plot_label = "a")

# Has a custom function due to lack of variation meaning it needs custom breaks.
p_year_terr_plot <- ggplot_terr_raster(p_year_terr_raster, 6, "") +
  annotate("text", x = 20, y = -48, label  = expression(bold("1"^ry*" consumers")), size = 8, fontface = 2)+
  annotate("text", x =  -170, y = 80, label  = "c", size = 12, fontface = 2)
s_terr_plot <- panel_ggplot_raster(s_terr_raster, plot_title = expression(bold("2"^ry*" consumers")), plot_label = "e")
s_year_terr_plot <- panel_ggplot_raster(s_year_terr_raster, plot_title = expression(bold("2"^ry*" consumers")), plot_label = "g")

# Arrange side plots and maps.
terr_side_plots <- ggarrange(pri_terr_lat_plot + rremove("xlab") + rremove("x.text"), 
                             pri_yearterr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_terr_lat_plot + rremove("xlab") + rremove("x.text"),
                             sec_yearterr_lat_plot,
                             nrow = 4, ncol = 1, heights = c(1,1,1,1.2))
terr_maps <- ggarrange(blank_plot, p_terr_plot, p_year_terr_plot, 
                       s_terr_plot, s_year_terr_plot, blank_plot,
                       nrow = 6, ncol = 1, heights = c(0.1,1,1,1,1,0.1))
# Arrange together.
terr_both_plots <- ggarrange(terr_maps, terr_side_plots, ncol = 2, widths = c(3,1.535)) +
  theme(plot.margin = margin(l = -0.3, unit = "cm"))
ggsave("Plots/Maps/figure_ED_7.png", height = 20, width = 15, dpi = 900)
ggsave("Plots/Maps/figure_ED_7.pdf", height = 20, width = 15, dpi = 900)



###############################################################################
    #### Plot lat gradient for migration and territory partitions ####


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
ggsave("Plots/Maps/figure_ED_3.png", height = 10, width = 15, dpi = 600)
ggsave("Plots/Maps/figure_ED_3.pdf", height = 10, width = 15, dpi = 600)

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
ggsave("Plots/Maps/figure_ED_4.png", height = 10, width = 15, dpi = 600)
ggsave("Plots/Maps/figure_ED_4.pdf", height = 10, width = 15, dpi = 600)


###############################################################################
                   ##### Species richness map ######



# Show high certainty sexual selection maps.
#m_spec_plot <- panel_ggplot_raster(m_spec_raster, plot_title = "Certainty 3 & 4 (n = 7592)", plot_label = "a")
h_spec_plot <- panel_ggplot_raster(h_spec_raster, plot_title = "Certainty 4 (n = 2851)", plot_label = "") +
  theme(plot.margin = margin(l = -0.3, r =0.3, unit = "cm"))

# Arrange maps together.
# cert_side_plots <- ggarrange(med_sex_lat_plot + rremove("xlab") + rremove("x.text"), 
#                              hi_sex_lat_plot, nrow = 2, ncol = 1, heights = c(1,1.2))
# cert_maps <- ggarrange(blank_plot, m_sex_plot, h_sex_plot, blank_plot, nrow = 4, ncol = 1,
#                        heights = c(0.1,1,1,0.1))
# cert_both_plots <- ggarrange(cert_maps, cert_side_plots, ncol = 2, widths = c(3,1.535)) +
#   theme(plot.margin = margin(l = -0.3, unit = "cm"))

# Export figure.
ggsave("Plots/Maps/figure_S2.png", height = 5, width = 10, dpi = 600)
















################################################################################
                     ##### Extra stuff ######



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


ranges <- sex_ranges %>% filter(diet_binary == "Frug-nect") 
# Function to make average raster using different species richness maps.
average_raster_2 <- function(ranges = sex_ranges, var_name = "sexual_score"){
  
  # Create the species richness raster.
  spec_raster <- fasterize(ranges, raster_template, fun = "sum")
  spec_raster[which(getValues(spec_raster) < 10)] <- NA
  spec_raster <- mask(spec_raster, land)
  
  # Create the average raster.
  var_raster <- fasterize(ranges, raster_template, fun = "sum", field = var_name)
  var_raster <- mask(var_raster, spec_raster)
  var_raster <- var_raster / spec_raster
  nbins = 6 
  variable = ""
  
   raster <- var_raster
  
  
  
  # Get the limits for breaks, using quantiles so each bin has an equal number of cells.
  breaks <- quantile(values(raster), seq(0,1,length.out=nbins+1), na.rm=TRUE)
  
  # Cut the cell values into bins using the breaks.
  cuts <- cut(values(raster), breaks = breaks, include.lowest = TRUE)
  
  # Replace cell values with bins.
  raster@data@values <- as.numeric(cuts)
  
  # Create raster data.
  raster_data <- as.data.frame(raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  raster_data$values %<>% as.factor()
  
  spec_data <- as.data.frame(spec_raster, xy=TRUE)
  colnames(spec_data) <- c("long", "lat", "n")
  
  hist(spec_data$n)
  raster_data <- left_join(raster_data, spec_data)
  # Get rid of areas without land.
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  
  # Create a label of the bins for plotting. Unicode is for en-dash.
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])
  
  # Plot with ggplot.
  ggplot() +
    
    # Add the raster data.
    geom_raster(aes(x=long, y=lat, fill= values, alpha = n), data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    scale_fill_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    # Make map closer to the edge.
    scale_y_continuous(limits = c(-57, 85), expand = expansion(), ) +
    scale_x_continuous(limits = c(-180, 180), expand = expansion()) +
    # Theme stuff.
    theme_classic(base_size = 18) + theme(axis.text = element_blank(),
                                          axis.ticks = element_blank(),
                                          axis.line = element_blank(),
                                          legend.position = c(0.13, 0.2),
                                          legend.key.height = unit(0.2, 'cm'),
                                          legend.spacing.y = unit(0.3, 'cm'),
                                          legend.title = element_blank(),
                                          legend.background = element_rect(fill = NA, colour = "lightgrey"),
                                          legend.margin = margin(t = 0, r = 0.4, b = 0.3, l = 0.4, unit = "cm"),
                                          plot.margin = margin(t = 0, l = 0, b = 0, r = 0)
                                          #  plot.margin = margin(t =0, l = -10, b = 0, r = 0.5, unit = "cm")#unit(c(0,0,0,0), "null")
    ) +
    ylab("") + 
    xlab("") 
}

















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

