


load("Data/Psittacidae_range_data.RData")

load("../../Teaching/BCB 2021/bcb_notebook_2021/corvid.RData")


library(sf)
library(fasterize)

clade_ranges %>% head()

bcb_ring_necked <- clade_ranges %>% filter(SCINAME =="Psittacula krameri")

# Start by creating an empty raster stack to store our data in.
raster_template <- raster(ncols=8640, nrows=3600, ymn = -60) #%>% crop(extent(-10,50,0,70))
#raster_template <- raster(ncols=17280, nrows=7200, ymn = -60)

# Create a map of species richness by summing overlapping ranges.
species_raster <- fasterize(new_corvid_ranges, raster_template, fun = "sum")

# Mask for land.
library(maptools)
data(wrld_simpl)
land <- wrld_simpl %>% st_as_sf() %>% fasterize(raster_template)
# species_raster <- mask(species_raster, land)
avonet <- read.csv("Data/sexual_traits.csv")

avonet %>% filter(family %in% c("Corvidae"))
sexy$family
library(raster)
parrot_raster <- fasterize(bcb_ring_necked, raster_template) %>% crop(extent(-10,50,0,70))

# Create data frame of land to use to crop ggplot maps.
land_data <- as.data.frame(land, xy=TRUE)


  # Create raster data.
  raster_data <- as.data.frame(species_raster, xy=TRUE)
  colnames(raster_data) <- c("long", "lat", "values")
  #raster_data$values %<>% as.factor()
  
  library(magrittr)
  library(tidyr)
  library(ggplot2)
  # Get rid of areas without land.
  raster_data$land <- land_data$layer
  raster_data %<>% drop_na(land)
  # Create a palette to match bin length.
  pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
  pal <- colorRampPalette(pal)(10)
  

  
  # Plot with ggplot.
  ggplot() +
    
    # Clip to world extent.
    xlim(-10, 50) +  
    
    # Add the raster data.
    geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
    
    # Specify colours, legend labels and legend title.
    #scale_fill_manual(values = pal, labels = labels, na.value = "grey") +
    scale_fill_viridis_c() +
    #guides(fill = guide_legend(variable)) +
    
    # Make map closer to the edge.
    scale_x_discrete(expand=c(0,1))+
    scale_y_discrete(expand=c(0,1))+
    # Theme stuff.
    theme_classic() + theme(axis.text = element_blank(),
                            axis.ticks = element_blank()
                            #plot.margin = unit(c(0,0,0,0), "null"),
                            #panel.margin = unit(c(0, 0, 0, 0), "null")
    ) +
    ylab("") + 
    xlab("") + 
    coord_fixed()



load("../../Phd_data/Spatial/jetz_ranges.RData")

jetz_ranges %>% head()


jetz_ring_necked <- jetz_ranges %>% filter(SCINAME =="Psittacula krameri")

filtered_jetz_ring_necked <- jetz_ring_necked %>% filter(PRESENCE %in% c(1,2, 3) & ORIGIN %in% c(1,2) & SEASONAL %in% c(1,2))


fil_parrot_raster <- fasterize(filtered_jetz_ring_necked, raster_template) %>% crop(extent(-10,50,0,70))


# Create raster data.
raster_data <- as.data.frame(fil_parrot_raster, xy=TRUE)
colnames(raster_data) <- c("long", "lat", "values")
raster_data$values %<>% as.factor()


library(tidyr)
# Get rid of areas without land.
raster_data$land <- land_data$layer
raster_data %<>% drop_na(land)
# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
pal <- colorRampPalette(pal)(10)



# Plot with ggplot.
ggplot() +
  
  # Clip to world extent.
  xlim(-10, 50) +  
  
  # Add the raster data.
  geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
  
  # Specify colours, legend labels and legend title.
  scale_fill_manual(values = pal, labels = labels, na.value = "grey") +
  #guides(fill = guide_legend(variable)) +
  
  # Make map closer to the edge.
  scale_x_discrete(expand=c(0,1))+
  scale_y_discrete(expand=c(0,1))+
  # Theme stuff.
  theme_classic() + theme(axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          legend.position = "none"
  ) +
  ylab("") + 
  xlab("") + 
  coord_fixed()


load("../../Phd_data/Spatial/clean_jetz_ranges.RData")

new_jetz_ring_necked <- new_jetz_ranges %>% filter(SCINAME =="Psittacula krameri")


new_parrot_raster <- fasterize(new_jetz_ring_necked, raster_template) %>% crop(extent(-10,50,0,70))


# Create raster data.
raster_data <- as.data.frame(new_parrot_raster, xy=TRUE)
colnames(raster_data) <- c("long", "lat", "values")
raster_data$values %<>% as.factor()


library(tidyr)
# Get rid of areas without land.
raster_data$land <- land_data$layer
raster_data %<>% drop_na(land)


# Plot with ggplot.
ggplot() +
  
  # Clip to world extent.
  xlim(-10, 50) +  
  
  # Add the raster data.
  geom_tile(aes(x=long, y=lat, fill= values), colour = NA, data=raster_data) +
  
  # Specify colours, legend labels and legend title.
  scale_fill_manual(values = pal, labels = labels, na.value = "grey") +
  #guides(fill = guide_legend(variable)) +
  
  # Make map closer to the edge.
  scale_x_discrete(expand=c(0,1))+
  scale_y_discrete(expand=c(0,1))+
  # Theme stuff.
  theme_classic() + theme(axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          legend.position = "none"
  ) +
  ylab("") + 
  xlab("") + 
  coord_fixed()

jetz_species <- read.csv("Data/GBD_BirdTreeTaxo_ecotraits_2022_November_7th.csv")

str(new_jetz_ranges)
library(dplyr)

jetz_species %<>% filter(Family == "Corvidae")

new_corvid_ranges <- new_jetz_ranges %>% filter(SCINAME %in% jetz_species$Bird_tree_name)

save(new_parrot_ranges, file = "Data/new_Psittacidae_ranges.RData")
