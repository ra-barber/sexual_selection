###############################################################################
                        # Running spatial models #
###############################################################################

# This script makes runs SAR models of sexual selection scores, accounting for 
# different ecological groups and data certainties.

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
library(spdep)
library(spatialreg)

# Clear environment. 
rm(list=ls())
gc()

# Functions.
source("Code/functions.R")
options(scipen = 999)


################################################################################
                              #### Data ####


# Load in the range maps.
load("Data/Spatial/clean_jetz_ranges.RData")

# Load in the trait data.
model_data <- read_ss_data()

# Join up the data.
sex_ranges <- left_join(model_data, new_jetz_ranges,  
                        by = c("scientific_name_bird_tree" = "SCINAME"))
sex_ranges <- st_as_sf(sex_ranges)
rm(new_jetz_ranges)


################################################################################
                ##### Projecting into behrmann #####


# Define behrman's proj.
behr_proj <- "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

# Define boundaries.
xmin = -18000000
xmax = 18000000
ymin = -6000000 # Could set as 60 degrees.
ymax = 8000000

# Create a blank raster template.
raster_template <- raster(ncols=180, nrows=70, xmn = xmin, xmx = xmax, 
                          ymn = ymin, ymx = ymax, crs = behr_proj)

# Transform bird ranges.
sex_ranges <- st_transform(sex_ranges, crs = behr_proj)

# Transform land.
data(wrld_simpl)
land <- wrld_simpl %>% st_as_sf()
land <- st_transform(land, crs = behr_proj)
land <- fasterize(land, raster_template)


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

# Supplementary figures.
h_sex_raster <- sex_ranges %>% filter(data_certainty == 4) %>% average_raster()
m_sex_raster <- sex_ranges %>% filter(data_certainty > 2) %>% average_raster()
mig_sex_raster <- sex_ranges %>% filter(migration_binary == "Strong") %>% average_raster()
no_mig_sex_raster <- sex_ranges %>% filter(migration_binary == "Weak") %>% average_raster()
terr_sex_raster <- sex_ranges %>% filter(territoriality_binary == "Territorial") %>% average_raster()
no_terr_sex_raster <- sex_ranges %>% filter(territoriality_binary == "Non-territorial") %>% average_raster()
p_terr_raster <- sex_ranges %>% filter(trophic_level_binary == "Primary") %>% average_raster(var_name = "terr_dummy")
s_terr_raster <- sex_ranges %>% filter(trophic_level_binary == "Secondary") %>% average_raster(var_name = "terr_dummy")
p_year_terr_raster <- sex_ranges %>% filter(trophic_level_binary == "Primary") %>% average_raster(var_name = "year_terr_dummy")
s_year_terr_raster <- sex_ranges %>% filter(trophic_level_binary == "Secondary") %>% average_raster(var_name = "year_terr_dummy")

##### High cert rasters #####

# Filter for high certainty species.
high_ranges <- sex_ranges %>% filter(data_certainty > 2)

# Diet rasters.
high_p_sex_raster <- high_ranges %>% filter(trophic_level_binary == "Primary") %>% average_raster()
high_s_sex_raster <- high_ranges %>% filter(trophic_level_binary == "Secondary") %>% average_raster()
high_fruit_sex_raster <- high_ranges %>% filter(trophic_niche == "Frugivore") %>% average_raster()
high_invert_sex_raster <- high_ranges %>% filter(trophic_niche == "Invertivore") %>% average_raster()

# Migration and territoriality.
high_mig_sex_raster <- high_ranges %>% filter(migration_binary == "Strong") %>% average_raster()
high_no_mig_sex_raster <- high_ranges %>% filter(migration_binary == "Weak") %>% average_raster()
high_terr_sex_raster <- high_ranges %>% filter(territoriality_binary == "Territorial") %>% average_raster()
high_no_terr_sex_raster <- high_ranges %>% filter(territoriality_binary == "Non-territorial") %>% average_raster()

# Remove big objects.
rm(sex_ranges)
save(list = ls(pattern =  "_raster"), file = "Results/Spatial/SAR_rasters.Rdata")
gc()


################################################################################
                    #### Create spatial data ####


# Get the cell resolution
cellsize <- res(sex_raster)[[1]]

# Figure 3.
sex_sf <- create_spatial_data(sex_raster)
cert_sf <- create_spatial_data(cert_raster)

# Figure 4.
p_sex_sf <- create_spatial_data(p_sex_raster)
s_sex_sf <- create_spatial_data(s_sex_raster)
fruit_sex_sf <- create_spatial_data(fruit_sex_raster)
invert_sex_sf <- create_spatial_data(invert_sex_raster)

# Supplementary figures.
h_sex_sf <- create_spatial_data(h_sex_raster)
m_sex_sf <- create_spatial_data(m_sex_raster)
mig_sex_sf <- create_spatial_data(mig_sex_raster)
no_mig_sex_sf <- create_spatial_data(no_mig_sex_raster)
terr_sex_sf <- create_spatial_data(terr_sex_raster)
no_terr_sex_sf <- create_spatial_data(no_terr_sex_raster)
p_terr_sf <- create_spatial_data(p_terr_raster)
s_terr_sf <- create_spatial_data(s_terr_raster)
p_year_terr_sf <- create_spatial_data(p_year_terr_raster)
s_year_terr_sf <- create_spatial_data(s_year_terr_raster)

# High certainty rasters.
high_p_sex_sf <- create_spatial_data(high_p_sex_raster)
high_s_sex_sf <- create_spatial_data(high_s_sex_raster)
high_fruit_sex_sf <- create_spatial_data(high_fruit_sex_raster)
high_invert_sex_sf <- create_spatial_data(high_invert_sex_raster)
high_mig_sex_sf <- create_spatial_data(high_mig_sex_raster)
high_no_mig_sex_sf <- create_spatial_data(high_no_mig_sex_raster)
high_terr_sex_sf <- create_spatial_data(high_terr_sex_raster)
high_no_terr_sex_sf <- create_spatial_data(high_no_terr_sex_raster)

# Save data for later.
save(list = ls(pattern =  "_sf"), file = "Results/Spatial/SAR_sfdata.Rdata")


################################################################################
                     ##### Run models ######


# Figure 3.
sex_model <- run_spatial_model(sex_sf)
cert_model <- run_spatial_model(cert_sf)

# Figure 4.
p_sex_model <- run_spatial_model(p_sex_sf) 
s_sex_model <- run_spatial_model(s_sex_sf)
fruit_sex_model <- run_spatial_model(fruit_sex_sf) 
invert_sex_model <- run_spatial_model(invert_sex_sf)

# Supplementary figures.
h_sex_model <- run_spatial_model(h_sex_sf) 
m_sex_model <- run_spatial_model(m_sex_sf) 
mig_sex_model <- run_spatial_model(mig_sex_sf) 
no_mig_sex_model <- run_spatial_model(no_mig_sex_sf)
terr_sex_model <- run_spatial_model(terr_sex_sf)
no_terr_sex_model <- run_spatial_model(no_terr_sex_sf)  
p_terr_model <- run_spatial_model(p_terr_sf)
s_terr_model <- run_spatial_model(s_terr_sf) s.
p_year_terr_model <- run_spatial_model(p_year_terr_sf)
s_year_terr_model <- run_spatial_model(s_year_terr_sf)

##### Scaled models for estimate comparisons ####

scaled_sex_model <- scaled_run_spatial_model(sex_sf)
scaled_cert_model <- scaled_run_spatial_model(cert_sf)
scaled_p_sex_model <- scaled_run_spatial_model(p_sex_sf)
scaled_s_sex_model <- scaled_run_spatial_model(s_sex_sf)
scaled_fruit_sex_model <- scaled_run_spatial_model(fruit_sex_sf) 
scaled_invert_sex_model <- scaled_run_spatial_model(invert_sex_sf)
scaled_h_sex_model <- scaled_run_spatial_model(h_sex_sf)
scaled_m_sex_model <- scaled_run_spatial_model(m_sex_sf)
scaled_mig_sex_model <- scaled_run_spatial_model(mig_sex_sf)
scaled_no_mig_sex_model <- scaled_run_spatial_model(no_mig_sex_sf)
scaled_terr_sex_model <- scaled_run_spatial_model(terr_sex_sf)
scaled_no_terr_sex_model <- scaled_run_spatial_model(no_terr_sex_sf)
scaled_p_terr_model <- scaled_run_spatial_model(p_terr_sf)
scaled_s_terr_model <- scaled_run_spatial_model(s_terr_sf)
scaled_p_year_terr_model <- scaled_run_spatial_model(p_year_terr_sf)
scaled_s_year_terr_model <- scaled_run_spatial_model(s_year_terr_sf)

#### High certainty ####
high_sex_model <- scaled_run_spatial_model(m_sex_sf)
high_p_sex_model <- scaled_run_spatial_model(high_p_sex_sf)
high_s_sex_model <- scaled_run_spatial_model(high_s_sex_sf)
high_fruit_sex_model <- scaled_run_spatial_model(high_fruit_sex_sf)
high_invert_sex_model <- scaled_run_spatial_model(high_invert_sex_sf)
high_mig_sex_model <- scaled_run_spatial_model(high_mig_sex_sf)
high_no_mig_sex_model <- scaled_run_spatial_model(high_no_mig_sex_sf)
high_terr_sex_model <- scaled_run_spatial_model(high_terr_sex_sf)
high_no_terr_sex_model <- scaled_run_spatial_model(high_no_terr_sex_sf)

# Save the models for quick editing later.
save(list = ls(pattern =  "_model"), file = "Results/SAR_200_models.Rdata")
load("Results/SAR_200_models.Rdata")


#################################################################################
               ####### Create side plots for figures ######


# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
nbins <- 6

# Figure 2.
sex_lat_plot <- sar_lat_side_plot(sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.2), 
                                  ybreaks =  c(0,0.5,1.0), lab_ypos = 0.12, plot_label = "b")

cert_lat_plot <- sar_lat_side_plot(cert_sf, ylabel = "Data certainty", ylimits = c(1,4), 
                   ybreaks =  c(1,2,3,4), lab_ypos = 1.3, plot_label = "d",
                   plot_model = cert_model, stats_model = scaled_cert_model)

# Primary and secondary consumers.
pri_lat_plot <- sar_lat_side_plot(p_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 2, 
                     plot_label = "b", plot_model = p_sex_model, stats_model = scaled_p_sex_model) 

sec_lat_plot <- sar_lat_side_plot(s_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                     ybreaks =  c(0,1.0,2.0), lab_ypos = 2,
                     plot_label = "f", plot_model = s_sex_model, stats_model = scaled_s_sex_model) 

# Trophic niche models.
fruit_lat_plot <- sar_lat_side_plot(fruit_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 2, 
                     plot_label = "d", plot_model = fruit_sex_model, stats_model = scaled_fruit_sex_model) 

invert_lat_plot <- sar_lat_side_plot(invert_sex_sf, ylabel = "Sexual selection", ylimits = c(-0.1,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 2, 
                     plot_label = "h", plot_model = invert_sex_model, stats_model = scaled_invert_sex_model) 

# Sexual lat gradient for high data certainty.
med_sex_lat_plot <- sar_lat_side_plot(m_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.3), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1.3, plot_label = "b",
                     plot_model = m_sex_model, stats_model = scaled_m_sex_model)

hi_sex_lat_plot <- sar_lat_side_plot(h_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.5), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.5, plot_label = "d",
                     plot_model = h_sex_model, stats_model = scaled_h_sex_model)

# Migration and territoriality.
mig_lat_plot <- sar_lat_side_plot(mig_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.5), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.5,  plot_label = "b",
                     plot_model = mig_sex_model, stats_model = scaled_mig_sex_model)

no_mig_lat_plot <- sar_lat_side_plot(no_mig_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.5), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.5, plot_label = "d",
                     plot_model = no_mig_sex_model, stats_model = scaled_no_mig_sex_model)

terr_lat_plot <-  sar_lat_side_plot(terr_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.6), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.6, plot_label = "b",
                     plot_model = terr_sex_model, stats_model = scaled_terr_sex_model)

no_terr_lat_plot <- sar_lat_side_plot(no_terr_sex_sf, ylabel = "Sexual selection", ylimits = c(0,1.6), 
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.6, plot_label = "d",
                     plot_model = no_terr_sex_model, stats_model = scaled_no_terr_sex_model)

# Primary terr.
pri_terr_lat_plot <- sar_lat_side_plot(p_terr_sf, ylabel = "Proportion of species", ylimits = c(0,1), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1, plot_label = "b",
                     plot_model = p_terr_model, stats_model = scaled_p_terr_model)
# Primary year terr.
pri_yearterr_lat_plot <-  sar_lat_side_plot(p_year_terr_sf, ylabel = "Proportion of species", ylimits = c(0,1), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1, plot_label = "d",
                     plot_model = p_year_terr_model, p_year_terr = TRUE, stats_model = scaled_p_year_terr_model)

# Secondary terr.
sec_terr_lat_plot <- sar_lat_side_plot(s_terr_sf, ylabel = "Proportion of species", ylimits = c(0,1), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1, plot_label = "f",
                     plot_model = s_terr_model, stats_model = scaled_s_terr_model)

# Secondary year terr.
sec_yearterr_lat_plot <-  sar_lat_side_plot(s_year_terr_sf, ylabel = "Proportion of species", ylimits = c(0,1), 
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1, plot_label = "h",
                     plot_model = s_year_terr_model, stats_model = scaled_s_year_terr_model)


# Export the plots.
save(list = ls(pattern =  "lat_plot"), file = "Plots/Maps/behr_200_latitudinal_sideplots.Rdata")
gc()


################################################################################
                ##### Extract model coefficients #####


# Model names for spatial models.
model_names <- c("all_birds", "primary", "fruit", "secondary", "invert", 
                 "migration", "no_migration", "territorial", "non_territorial")

# Blank dataframe.
all_estimates <- data.frame(model = model_names, estimate = 0, lower = 0, upper = 0)
high_estimates <- all_estimates

# Function for extracting coeffs.
extract_coeffs <- function(model){
  estimate <- last(coef(summary(model))[,1])
  se <- last(coef(summary(model))[,2])
  upper <- estimate + 1.96*se
  lower <-  estimate - 1.96*se
  c(estimate, lower, upper)
}

# Extract models. 
all_estimates[1,2:4] <- extract_coeffs(scaled_sex_model)
all_estimates[2,2:4] <- extract_coeffs(scaled_p_sex_model)
all_estimates[3,2:4] <- extract_coeffs(scaled_fruit_sex_model)
all_estimates[4,2:4] <- extract_coeffs(scaled_s_sex_model)
all_estimates[5,2:4] <- extract_coeffs(scaled_invert_sex_model)
all_estimates[6,2:4] <- extract_coeffs(scaled_mig_sex_model)
all_estimates[7,2:4] <- extract_coeffs(scaled_no_mig_sex_model)
all_estimates[8,2:4] <- extract_coeffs(scaled_terr_sex_model)
all_estimates[9,2:4] <- extract_coeffs(scaled_no_terr_sex_model)


# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(estimate, 2),
  intervals = paste0(round(lower, 2), ", ", 
                     round(upper, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export.
write.csv(all_estimates, "Results/Tables/all_spatial_regression.csv", row.names = FALSE)

# Reapeat for high certainty.
high_estimates[1,2:4] <- extract_coeffs(high_sex_model)
high_estimates[2,2:4] <- extract_coeffs(high_p_sex_model)
high_estimates[3,2:4] <- extract_coeffs(high_fruit_sex_model)
high_estimates[4,2:4] <- extract_coeffs(high_s_sex_model)
high_estimates[5,2:4] <- extract_coeffs(high_invert_sex_model)
high_estimates[6,2:4] <- extract_coeffs(high_mig_sex_model)
high_estimates[7,2:4] <- extract_coeffs(high_no_mig_sex_model)
high_estimates[8,2:4] <- extract_coeffs(high_terr_sex_model)
high_estimates[9,2:4] <- extract_coeffs(high_no_terr_sex_model)

# Paste together values for reporting in a table.
high_estimates %<>% mutate(
  round_est = round(estimate, 2),
  intervals = paste0(round(lower, 2), ", ", 
                     round(upper, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export.
write.csv(high_estimates, "Results/Tables/high_spatial_regression.csv", row.names = FALSE)

###############################################################################
                             #### END ####
###############################################################################