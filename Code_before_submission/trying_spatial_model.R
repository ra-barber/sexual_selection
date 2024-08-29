###############################################################################
                        # Running spatial models #
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

# Create the Behrmann's CRS.
behr <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
# 
# behr_crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
# 
# behr_proj <- "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
# 
# wgs84_proj <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

# Start by creating an empty raster stack to store our data in. 
# (This creates 5km res, which is only used for visualisation)
superhigh_res_template <- raster(ncols=8640, nrows=3600, ymn = -60)
low_res_template <- raster(ncols=86, nrows=36, ymn = -60, crs = behr)

low_res_template <- raster(ncols=36, nrows=18, ymn = -90, crs = behr)
low_res_template <- raster(ncols=360, nrows=180, ymn = -90, crs = behr)

low_res_template <- raster(ncols=180, nrows=90, ymn = -90, crs = behr)


low_res_template <- raster(ncols=86, nrows=36, ymn = -60)  # This is the only one that works with eco regions.
high_res_template <- raster(ncols=864, nrows=360, ymn = -60, crs = behr)

raster_template <- low_res_template

# Create a raster of landmass to mask maps.
data(wrld_simpl)
land <- wrld_simpl %>% st_as_sf() %>% fasterize(raster_template)



################################################################################
                ##### Projecting into behrmann #####

# Read in the blank raster.
blank_raster <- raster("Data/spatial/blank_map_1deg.tif")

# Define behrman's proj.
behr_proj <- "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

# Define normal long lat.
wgs84_proj <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

# See what sex ranges crs is.
crs(sex_ranges)

# Transform land.
land <- wrld_simpl %>% st_as_sf()
behr_land <- st_transform(land, crs = behr_proj)


# Transform bird ranges.
behr_ranges <- st_transform(sex_ranges, crs = behr_proj)

# Transform blank raster.
behr_blank <- projectRaster(blank_raster, crs = behr_proj)

# Fasterize land for a mask.
behr_land_raster <- fasterize(behr_land, behr_blank)


plot(behr_blank)

plot(behr_land_raster)

land <- behr_land_raster
raster_template <- behr_blank
# Blank raster (1degree)

#blank_raster <- projectRaster(blank_raster, crs=behr)
plot(raster_template)
plot(sex_raster)

spec_raster <- fasterize(behr_ranges, behr_blank, fun = "sum") # Rasterize
spec_raster[which(getValues(spec_raster) < 10)] <- NA # Remove low SR cells
spec_raster <- mask(spec_raster, behr_land_raster) # Mask for land.
return(spec_raster)

# Create an average raster for the variable and data in question.
var_raster <- fasterize(behr_ranges, behr_blank, fun = "sum", field = "sexual_selection")
var_raster <- mask(var_raster, spec_raster)
var_raster <- var_raster / spec_raster

sex_ranges <- behr_ranges

plot(var_raster)

sex_raster <- var_raster

raster_template <- behr_blank
land <- behr_land_raster
################################################################################
                  ##### Create rasters for each group ######


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


# Remove big objects.
#rm(sex_ranges)
gc()

################################################################################
          ##### Function for spatial auto-correlation ####

library(spdep)
library(spatialreg)

# Get the cell resolution
cellsize <- res(sex_raster)[[1]] *125 # (Convert from degree to km with buffer for rounding issues)
cellsize <- res(sex_raster)[[1]] /1000 # (Convert to km with buffer for rounding issues)

cellsize <- res(sex_raster)[[1]]

cellsize <- 200 #(km)
cellsize <- 96500
spatial_raster <- sex_raster

# Extract latutide from cells.
extract_lat <- function(x){
  nums <- as.numeric(x)
  abs(nums[2])
}

spatial_raster <- sex_raster 

# Function to convert raster in sf_data with neighbourhood weights.
create_spatial_data <- function(spatial_raster){
  
  # Create spatial feature dataframe.
  data_sf <- st_as_sf(as(spatial_raster, 'SpatialPixelsDataFrame'))
  
  # Use queen nearest neighbours.
  queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize)
  data_sf$card_queen <- card(queen)
  
  # Remove the cells with zero neighbours
  data_sf <- subset(data_sf, card_queen > 0)
  
  # Recalculate neighbours.
  queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize) # + 1)
  data_sf$card_queen <- card(queen)
  
  # Convert to weights.
  queen <- nb2listw(queen, style='W')
  
  # Add latitude.
  list_test <- lapply(data_sf$geometry, extract_lat)
  data_sf$lat <- as.numeric(list_test)
  data_sf$abs_lat <- abs(data_sf$lat)
  
  return(data_sf)
  
}

# Function to run spatial model.
run_spatial_model <- function(data_sf){
  
  # Calculate neighbourhood weights.
  queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize) # + 1)
  queen <- nb2listw(queen, style='W')
  
  # Run model.
  # errorsarlm(scale(layer) ~ scale(abs_lat), 
  #            data=data_sf, listw=queen)
  errorsarlm(layer ~ abs_lat, 
             data=data_sf, listw=queen)
}


# Create a function to plot cells against data.
plot_model_preds <- function(data, model){
  
  # extract the predictions from the model into the spatial data frame
  data$model_fit <- predict(model)
  
  # Plot average sexual selection against pixel latitude.
  plot(data$layer ~ data$abs_lat)
  
  # Create new data for predictions.
  new_data <- data.frame(abs_lat = seq(from = 0, to = round(max(data$abs_lat)), by = 1 ))
  
  # Extract model predictions.
  new_data$preds <- predict(model, newdata = new_data)[,1]
  
  # Add the line from model preds.
  lines(new_data$preds ~ new_data$abs_lat)
  
}




# Figure 2.
sex_sf <- create_spatial_data(sex_raster)
cert_sf <- create_spatial_data(cert_raster)

# Figure 3.
p_sex_sf <- create_spatial_data(p_sex_raster)
s_sex_sf <- create_spatial_data(s_sex_raster)
fruit_sex_sf <- create_spatial_data(fruit_sex_raster)
invert_sex_sf <- create_spatial_data(invert_sex_raster)


# Extended data figure 3.
h_sex_sf <- create_spatial_data(h_sex_raster)
m_sex_sf <- create_spatial_data(m_sex_raster)

# Extended data figure 4.
mig_sex_sf <- create_spatial_data(mig_sex_raster)
no_mig_sex_sf <- create_spatial_data(no_mig_sex_raster)

# Extended data figure 5.
terr_sex_sf <- create_spatial_data(terr_sex_raster)
no_terr_sex_sf <- create_spatial_data(no_terr_sex_raster)

# Extended data figure 8.
p_terr_sf <- create_spatial_data(p_terr_raster)
s_terr_sf <- create_spatial_data(s_terr_raster)
p_year_terr_sf <- create_spatial_data(p_year_terr_raster)
s_year_terr_sf <- create_spatial_data(s_year_terr_raster)


##########

library(tictoc)
tic()
sex_model <- run_spatial_model(sex_sf)
toc()
# Figure 2.
sex_model <- run_spatial_model(sex_sf)   # 17 minutes.
cert_model <- run_spatial_model(cert_sf)

# Figure 3.
p_sex_model <- run_spatial_model(p_sex_sf)
s_sex_model <- run_spatial_model(s_sex_sf)
fruit_sex_model <- run_spatial_model(fruit_sex_sf) #62 seconds.
invert_sex_model <- run_spatial_model(invert_sex_sf)


# Extended data figure 3.
h_sex_model <- run_spatial_model(h_sex_sf)
m_sex_model <- run_spatial_model(m_sex_sf)

# Extended data figure 4.
mig_sex_model <- run_spatial_model(mig_sex_sf)
no_mig_sex_model <- run_spatial_model(no_mig_sex_sf)

# Extended data figure 5.
terr_sex_model <- run_spatial_model(terr_sex_sf)
no_terr_sex_model <- run_spatial_model(no_terr_sex_sf)

# Extended data figure 8.
p_terr_model <- run_spatial_model(p_terr_sf)
s_terr_model <- run_spatial_model(s_terr_sf)
p_year_terr_model <- run_spatial_model(p_year_terr_sf)
s_year_terr_model <- run_spatial_model(s_year_terr_sf)


#######


# Figure 2.
summary(sex_model) 
summary(cert_model)

# Figure 3.
summary(p_sex_model) #Insignificant
summary(s_sex_model)
summary(fruit_sex_model)
summary(invert_sex_model)

# Extended data figure 3.
summary(h_sex_model)
summary(m_sex_model)

# Extended data figure 4.
summary(mig_sex_model)
summary(no_mig_sex_model)

# Extended data figure 5.
summary(terr_sex_model)
summary(no_terr_sex_model) #Insignificant

# Extended data figure 8.
summary(p_terr_model)
summary(s_terr_model)
summary(p_year_terr_model)
summary(s_year_terr_model)




#################################################################################
                    ####### Try making a fancier plot ######


# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
nbins <- 6

# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
sar_lat_side_plot <- function(data_set, ylabel = "", ylimits = c(0,1.2), ybreaks = c(0,0.5,1), 
                               lab_x_pos = 70, lab_ypos = 1.2, plot_label = "b", 
                               plot_model = sex_model, #plot_raster = sex_raster,
                               sex_score = TRUE, p_year_terr = FALSE, p_include = TRUE){
  
  # Extract model predictions with new data.
  new_data <- data.frame(abs_lat = seq(from = 0, to = round(max(data_set$abs_lat)), by = 1 ))
  new_data$preds <- predict(plot_model, newdata = new_data)[,1]
  
  # Estimate
  estimate <- last(coef(summary(plot_model))[,1])
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- last(coef(summary(plot_model))[,1])
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }

  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)

  # # Extract p-values using probability of direction two-tailed test.
  p_value <- last(coef(summary(plot_model))[,4])
  
  # Change p value to a string, using standard thresholds. 
  if (p_value < 0.001 ){
    p_value <- "p < 0.001"
  } else if (p_value < 0.01) {
    p_value <- "p < 0.01"
  } else if (p_value < 0.05) {
    p_value <- "p < 0.05"
  } else {
    p_value <- paste0("p = ", as.character(format(round(p_value, 2), nsmall = 2)))
  }
  
  # Create a label 
  stats_label <- paste0(estimate, "\n", p_value)
  
  if (p_year_terr){
    # Special case for when primary territory raster is skewed by too many zeroes.
    breaks <- quantile(data_set$layer[data_set$layer>0], seq(0,1,length.out=nbins), na.rm=TRUE)
    breaks <- c(0,breaks) 
    cuts <- cut(data_set$layer, breaks = breaks, include.lowest = TRUE)
  } else {
    # Create binned colour scale.
    breaks <- quantile(data_set$layer, seq(0,1,length.out=nbins+1), na.rm=TRUE)
    cuts <- cut(data_set$layer, breaks = breaks, include.lowest = TRUE)
  }
  data_set$binned_layer <- as.factor(as.numeric(cuts))
  
  # Create a label of the bins for plotting. (Only used in map legends)
  labels <- as.character(format(round(breaks, 2), nsmall = 2))
  labels <- paste0(labels[1:nbins], " \u2013 ", labels[2:(nbins+1)])
  
  ggplot(data_set, aes(x = abs_lat, y = layer)) +
    geom_point(aes(colour = binned_layer), size = 3, alpha = 0.9) + 
    scale_x_continuous(breaks = seq(from = 0, to = 80, by = 40)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 80), clip = 'off') +
    scale_size_continuous(range = c(1,6))+
    ylab(ylabel) +
    xlab("Latitude") + theme_classic(base_size = 25) + 
    
    scale_colour_manual(values = rev(colorRampPalette(pal)(nbins)), breaks = rev(1:nbins), labels = rev(labels), na.value = "lightgrey") +     #na.value = "grey"
    guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
    
    theme(legend.position = "none", 
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7) + 
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
    geom_line(data = new_data, aes(x = abs_lat, y = preds), linetype = "dashed", linewidth = 1)
}



# Figure 2.
sex_lat_plot <- sar_lat_side_plot(sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.2), 
                                  ybreaks =  c(0,0.5,1.0), lab_ypos = 0.12, plot_label = "b")

cert_lat_plot <- sar_lat_side_plot(cert_sf, ylabel = "Data certainty", ylimits = c(1,4), 
                   ybreaks =  c(1,2,3,4), lab_ypos = 1.3, plot_label = "d",
                   plot_model = cert_model)


# Primary and secondary consumers.
pri_lat_plot <- sar_lat_side_plot(p_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.5, 
                     plot_label = "b", plot_model = p_sex_model) 

sec_lat_plot <- sar_lat_side_plot(s_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                     ybreaks =  c(0,1.0,2.0), lab_ypos = 1.5,
                     plot_label = "f", plot_model = s_sex_model) 

# Trophic niche models.
fruit_lat_plot <- sar_lat_side_plot(fruit_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.5, 
                     plot_label = "d", plot_model = fruit_sex_model) 

# fruit_lat_plot <- sar_lat_side_plot(fruit_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
#                                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.5, 
#                                     plot_label = "d", plot_model = fruit_sex_model, p_year_terr = TRUE) 

# 
# fruit_lat_plot_3 <- sar_lat_side_plot(fruit_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
#                                       ybreaks =  c(0,1.0, 2.0), lab_ypos = 1.5, 
#                                       plot_label = "d", plot_model = fruit_sex_model, p_year_terr = TRUE) 

# over_10_plot <- fruit_lat_plot
# over_5_plot <- fruit_lat_plot_2
# over_0_plot <- fruit_lat_plot_3


invert_lat_plot <- sar_lat_side_plot(invert_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_x_pos = 20, lab_ypos = 1.5, 
                     plot_label = "h", plot_model = invert_sex_model) 


# Sexual lat gradient for high data certainty.
med_sex_lat_plot <- sar_lat_side_plot(m_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.2), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1.11, plot_label = "b",
                     plot_model = m_sex_model,
                     sex_score = TRUE)

hi_sex_lat_plot <- sar_lat_side_plot(h_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.8), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.67, plot_label = "d",
                     plot_model = h_sex_model,
                     sex_score = TRUE)

# Make the side plots.
mig_lat_plot <- sar_lat_side_plot(mig_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.5), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39,  plot_label = "b",
                     plot_model = mig_sex_model, 
                     sex_score = TRUE)

no_mig_lat_plot <- sar_lat_side_plot(no_mig_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.5), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "d",
                     plot_model = no_mig_sex_model, 
                     sex_score = TRUE)

terr_lat_plot <-  sar_lat_side_plot(terr_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.5), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "b",
                     plot_model = terr_sex_model,
                     sex_score = TRUE)

no_terr_lat_plot <- sar_lat_side_plot(no_terr_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.5), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.39, plot_label = "d",
                     plot_model = no_terr_sex_model, 
                     sex_score = TRUE)

# Primary terr.
pri_terr_lat_plot <- sar_lat_side_plot(p_terr_sf, ylabel = "Proportion of species", ylimits = c(0,1), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "b",
                     plot_model = p_terr_model,
                     sex_score = FALSE)
# Primary year terr.
pri_yearterr_lat_plot <-  sar_lat_side_plot(p_year_terr_sf, ylabel = "Proportion of species", ylimits = c(0,1), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "d",
                     plot_model = p_year_terr_model, p_year_terr = TRUE,
                     sex_score = FALSE)

# Secondary terr.
sec_terr_lat_plot <- sar_lat_side_plot(s_terr_sf, ylabel = "Proportion of species", ylimits = c(0,1), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "f",
                     plot_model = s_terr_model,
                     sex_score = FALSE)

# Secondary year terr.
sec_yearterr_lat_plot <-  sar_lat_side_plot(s_year_terr_sf, ylabel = "Proportion of species", ylimits = c(0,1), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 0.93, plot_label = "h",
                     plot_model = s_year_terr_model,
                     sex_score = FALSE)




# Export the plots.
save(list = ls(pattern =  "lat_plot"), file = "Plots/Maps/cell_model_latitudinal_sideplots.Rdata")
gc()





###############################################################################
                      ##### Look at model fit ######





plot_model_preds(sex_sf, sex_model)


# extract the predictions from the model into the spatial data frame
sex_sf$model_fit <- predict(sex_model)

data_sf$abs_lat <- abs(data_sf$lat)

plot(sex_sf$layer ~ sex_sf$abs_lat)

# Create 5 degree bins for latitude.
bin_range <- seq(0, 80, 5)
bin_labels <- seq(0, 75, 5) 
sex_sf$binned_lat <- cut(sex_sf$abs_lat, breaks=bin_range, labels = bin_labels, include.lowest = TRUE)
sex_sf$binned_lat %<>% as.character() %>% as.numeric()

# Create grouped data to reuse with different varaibles below.
grouped_lat_data <- sex_sf %>% filter(binned_lat != 75) %>% group_by(binned_lat)

# Function for grouping sexual selection by trait.
average_lat_bins <- function(grouped_data, predictor = "sexual_selection"){
  grouped_data %<>% 
    summarise(trait_mean = mean(!!! syms(predictor)),
              trait_sd = sd(!!! syms(predictor)),
              trait_se = sd(!!! syms(predictor))/sqrt(length(!!! syms(predictor))),
              trait_max = trait_mean + (trait_se*1.96),
              trait_min = trait_mean - (trait_se*1.96),
              trait_n = length(!!! syms(predictor))) %>% na.omit()
  grouped_data %>% filter(trait_n > 10)
}

# Group sex scores and certainty by lat bins.
lat_data <- grouped_lat_data %>% average_lat_bins(predictor = "layer")

fit_data <- grouped_lat_data %>% average_lat_bins(predictor = "sar_fit")


plot(lat_data$trait_mean ~ lat_data$binned_lat)
plot(fit_data$trait_mean ~ fit_data$binned_lat)

hist(log(data_sf$layer+ 1))

new_data <- data.frame(abs_lat = seq(from = 0, to = 70, by = 1 ))

new_data$preds <- predict(sex_model, newdata = new_data)[,1]

lines(new_data$preds ~ new_data$abs_lat)
















################################################################################
             ##### Caclulate spatial autocorrelation ####



data_spdf <- as(sex_raster, 'SpatialPixelsDataFrame')
summary(data_spdf)

data_sf <- st_as_sf(data_spdf)
print(data_sf)

# Calculate spatial autocorrelation.
hist(sex_raster)

# Get the cell resolution
cellsize <- res(sex_raster)[[1]] + 0.1


# Give dnearneigh the coordinates of the points and the distances to use
rook <- dnearneigh(data_sf, d1=0, d2=cellsize)
queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize)


# Store the neighbourhood cardinalities in data_sf
data_sf$card_rook <- card(rook)
data_sf$card_queen <- card(queen)
# Look at the count of rook and queen neighbours for each point
plot(data_sf[c('card_rook', 'card_queen')], key.pos=4)


# Remove the island cells with zero neighbours
data_sf <- subset(data_sf, card_rook > 0)

# Recalculate weights.
rook <- dnearneigh(data_sf, d1=0, d2=cellsize)# + 1)
queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize) # + 1)
data_sf$card_rook <- card(rook)
data_sf$card_queen <- card(queen)
knn <- knearneigh(data_sf, k=8)

# Convert to weights.
rook <- nb2listw(rook, style='W')
queen <- nb2listw(queen, style='W')
knn <- nb2listw(knn2nb(knn), style='W')

# Test for spatial autocorrelation.
moran_sexual_selection <- moran.test(data_sf$layer, rook)
print(moran_sexual_selection)

# Look on a local level.
local_moran_avr <- localmoran(data_sf$layer, rook)
data_sf$local_moran_avr <- local_moran_avr[, 'Ii'] # Observed Moran's I
plot(data_sf['local_moran_avr'], cex=0.6, pch=20)


data_sf$local_g_avian_richness <- localG(data_sf$layer, rook)
plot(data_sf['local_g_avian_richness'], cex=0.6, pch=20)


###############################################################################
                      ###### Model the data ######


# Add latutide.
extract_lat <- function(x){
  nums <- as.numeric(x)
  abs(nums[2])
}
list_test <- lapply(data_sf$geometry, extract_lat)
data_sf$lat <- as.numeric(list_test)


# Fit a simple linear model
simple_model <- lm(layer ~ lat, data=data_sf)
summary(simple_model)

# Fit a spatial autoregressive model: this is much slower and can take minutes to calculate
library(spatialreg)

sar_model <- errorsarlm(layer ~ lat, 
                        data=data_sf, listw=queen)

summary(sar_model)
hist(sar_model$residuals, breaks = 50)

################################################################################
                   ##### Look at model fit #####



# extract the predictions from the model into the spatial data frame
data_sf$simple_fit <- predict(simple_model)
data_sf$sar_fit <- predict(sar_model)

data_sf$abs_lat <- abs(data_sf$lat)

plot(data_sf$sar_fit ~ data_sf$abs_lat)

# Create 5 degree bins for latitude.
bin_range <- seq(0, 80, 5)
bin_labels <- seq(0, 75, 5) 
data_sf$binned_lat <- cut(data_sf$abs_lat, breaks=bin_range, labels = bin_labels, include.lowest = TRUE)
data_sf$binned_lat %<>% as.character() %>% as.numeric()

# Create grouped data to reuse with different varaibles below.
grouped_lat_data <- data_sf %>% filter(binned_lat != 75) %>% group_by(binned_lat)

# Function for grouping sexual selection by trait.
average_lat_bins <- function(grouped_data, predictor = "sexual_selection"){
  grouped_data %<>% 
    summarise(trait_mean = mean(!!! syms(predictor)),
              trait_sd = sd(!!! syms(predictor)),
              trait_se = sd(!!! syms(predictor))/sqrt(length(!!! syms(predictor))),
              trait_max = trait_mean + (trait_se*1.96),
              trait_min = trait_mean - (trait_se*1.96),
              trait_n = length(!!! syms(predictor))) %>% na.omit()
  grouped_data %>% filter(trait_n > 10)
}

# Group sex scores and certainty by lat bins.
lat_data <- grouped_lat_data %>% average_lat_bins(predictor = "layer")

fit_data <- grouped_lat_data %>% average_lat_bins(predictor = "sar_fit")


plot(lat_data$trait_mean ~ lat_data$binned_lat)
plot(fit_data$trait_mean ~ fit_data$binned_lat)

hist(log(data_sf$layer+ 1))

new_data <- data.frame(lat = seq(from = 0, to = 30, by = 1 ))

new_data$preds <- predict(sar_model, newdata = new_data, pred.type = "trend")[,1]
new_data$preds <- as.data.frame(predict(sar_model))[,1]

lines(new_data$preds ~ new_data$lat)




?predict.Sarlm



# Compare those two predictions with the data
plot(data_sf[c('layer','simple_fit','sar_fit')], 
     pal=function(n) hcl.colors(n, pal='Blue-Red'), key.pos=4, pch=19)


# extract the residuals from the model into the spatial data frame
data_sf$simple_resid <- residuals(simple_model)
data_sf$sar_resid <- residuals(sar_model)
plot(data_sf[c('simple_resid', 'sar_resid')], 
     pal=function(n) hcl.colors(n, pal='Blue-Red'), key.pos=4, pch=19)

# Correlograms
library(ncf)

# add the X and Y coordinates to the data frame
data_xy <- data.frame(st_coordinates(data_sf))
data_sf$x <- data_xy$X
data_sf$y <- data_xy$Y

# calculate a correlogram for avian richness: a slow process!
rich.correlog <- with(data_sf, correlog(x, y, layer, increment=cellsize, resamp=0))
plot(rich.correlog)























distance(sex_raster)

area(sex_raster)
# Code for figuring out the resolution.
#library(terra)
distance(cbind(0,-21), cbind(0,-21.04166667), lonlat=TRUE)
?distance
distance(raster_template)

plot(sex_raster)

# Should do it as a spatial pixel data frame in david's course.
sex_sf <- rasterToPoints(sex_raster, spatial = TRUE)

sex_sf <- st_as_sf(sex_sf)

library(spdep)

cellsize <- res(sex_raster)[[1]]

?dnearneigh

cell_km <- cellsize*111
rook <- dnearneigh(sex_sf, d1=0, d2=cell_km)

queen <- dnearneigh(sex_sf, d1=0, d2=cell_km*sqrt(2))

knn <- knearneigh(sex_sf, k=8)
# We have to look at the `nn` values in `knn` to see the matrix of neighbours
head(knn$nn, n=3)
head(rook)

#card(knn)


# Store the neighbourhood cardinalities in data_sf
sex_sf$card_rook <- card(rook)
sex_sf$card_queen <- card(queen)
# Look at the count of rook and queen neighbours for each point
plot(sex_sf[c('card_rook', 'card_queen')], key.pos=4)
