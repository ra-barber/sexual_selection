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
library(spdep)
library(spatialreg)

# Clear environment. 
rm(list=ls())
gc()

# Functions.
source("Code/functions.R")




# Read in models #

# Directory where models are saved.
first_half <- "Z:/home/sexual_selection/Results/Models/Spatial/"
list.files(first_half)

# Function to read in models.
read_lat_model <- function(model_name){
  model_file <- paste0(first_half, model_name, ".rds")
  readRDS(model_file)
}


#### Scaled models #####

# Figure 2.
scaled_all_allbirds_model <- read_lat_model("scaled_all_allbirds")
scaled_all_certainty_model <- read_lat_model("scaled_all_certainty")

# Diet.
scaled_all_primary_model <- read_lat_model("scaled_all_primary")
scaled_all_fruit_model <- read_lat_model("scaled_all_fruit")

scaled_all_secondary_model <- read_lat_model("scaled_all_secondary")
scaled_all_invert_model <- read_lat_model("scaled_all_invert")

# Migration and territoriality.
scaled_all_mig_model <- read_lat_model("scaled_all_mig")
scaled_all_nomig_model <- read_lat_model("scaled_all_nomig")

scaled_all_terr_model <- read_lat_model("scaled_all_terr")
scaled_all_noterr_model <- read_lat_model("scaled_all_noterr")

# Proportion of terr models.
scaled_all_primterr_model <- read_lat_model("scaled_all_primterr")
scaled_all_primyearterr_model <- read_lat_model("scaled_all_primyearterr")

scaled_all_secterr_model <- read_lat_model("scaled_all_secterr")
scaled_all_secyearterr_model <- read_lat_model("scaled_all_secyearterr")

# High certainty sexual selection.
scaled_high_allbirds_model <- read_lat_model("scaled_high_allbirds")
scaled_all_highcert_model <- read_lat_model("scaled_all_highcert")


# Figure 2.
raw_all_allbirds_model <- read_lat_model("raw_all_allbirds")
raw_all_certainty_model <- read_lat_model("raw_all_certainty")

# Diet.
raw_all_primary_model <- read_lat_model("raw_all_primary")
raw_all_fruit_model <- read_lat_model("raw_all_fruit")

raw_all_secondary_model <- read_lat_model("raw_all_secondary")
raw_all_invert_model <- read_lat_model("raw_all_invert")

# Migration and territoriality.
raw_all_mig_model <- read_lat_model("raw_all_mig")
raw_all_nomig_model <- read_lat_model("raw_all_nomig")

raw_all_terr_model <- read_lat_model("raw_all_terr")
raw_all_noterr_model <- read_lat_model("raw_all_noterr")

# Proportion of terr models.
raw_all_primterr_model <- read_lat_model("raw_all_primterr")
raw_all_primyearterr_model <- read_lat_model("raw_all_primyearterr")

raw_all_secterr_model <- read_lat_model("raw_all_secterr")
raw_all_secyearterr_model <- read_lat_model("raw_all_secyearterr")

### High certainty maps #####

# Super high certainty map.
raw_high_allbirds_model <- read_lat_model("raw_high_allbirds")
raw_all_highcert_model <- read_lat_model("raw_all_highcert")


##### Models for underlying data #####


# Read in models #

# Directory where models are saved.
first_half <- "Z:/home/sexual_selection/Results/Models/Spatial/500km/"
list.files(first_half)

# Function to read in models.
read_lat_model <- function(model_name){
  model_file <- paste0(first_half, model_name, ".rds")
  readRDS(model_file)
}


# Figure 2.
raw_all_allbirds_model_500 <- read_lat_model("raw_all_allbirds")
raw_all_certainty_model_500 <- read_lat_model("raw_all_certainty")

# Diet.
raw_all_primary_model_500 <- read_lat_model("raw_all_primary")
raw_all_fruit_model_500 <- read_lat_model("raw_all_fruit")

raw_all_secondary_model_500 <- read_lat_model("raw_all_secondary")
raw_all_invert_model_500 <- read_lat_model("raw_all_invert")

# Migration and territoriality.
raw_all_mig_model_500 <- read_lat_model("raw_all_mig")
raw_all_nomig_model_500 <- read_lat_model("raw_all_nomig")

raw_all_terr_model_500 <- read_lat_model("raw_all_terr")
raw_all_noterr_model_500 <- read_lat_model("raw_all_noterr")

# Proportion of terr models.
raw_all_primterr_model_500 <- read_lat_model("raw_all_primterr")
raw_all_primyearterr_model_500 <- read_lat_model("raw_all_primyearterr")

raw_all_secterr_model_500 <- read_lat_model("raw_all_secterr")
raw_all_secyearterr_model_500 <- read_lat_model("raw_all_secyearterr")

### High certainty maps #####

# Super high certainty map.
raw_high_allbirds_model_500 <- read_lat_model("raw_high_allbirds")
raw_all_highcert_model_500 <- read_lat_model("raw_all_highcert")



#################################################################################
                    ####### Try making a fancier plot ######


# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
nbins <- 6

# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
sar_lat_side_plot <- function(ylabel = "", ylimits = c(0,1.2), ybreaks = c(0,0.5,1),
                               lab_x_pos = 60, lab_ypos = 1.2, plot_label = "b",
                               plot_model = sex_model, stats_model = scaled_sex_model,
                              predict_model = raw_all_allbirds_model,
                               sex_score = TRUE, p_year_terr = FALSE, p_include = TRUE){
  
  # Pull out the data from the raw model.
  data_set <-  data.frame(abs_lat = plot_model$X[,2], layer = plot_model$y)

  # Extract model predictions with new data.
  new_data <- data.frame(abs_lat = seq(from = 0, to = round(max(data_set$abs_lat)), by = 1 ))
  new_data$preds <- predict(predict_model, newdata = new_data)[,1]

  # Estimate
  estimate <- last(coef(summary(stats_model))[,1])
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))

  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- last(coef(summary(stats_model))[,1])
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }

  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)

  # # Extract p-values using probability of direction two-tailed test.
  p_value <- last(coef(summary(stats_model))[,4])

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
    scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 70), clip = 'off') +
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
    annotate("text", x = 1, y = ylimits[2], label = plot_label, size = 12, fontface = 2) +
    geom_line(data = new_data, aes(x = abs_lat, y = preds), linetype = "dashed", linewidth = 2)
}


# Figure 2.
sex_lat_plot <- sar_lat_side_plot(stats_model = scaled_all_allbirds_model, 
                                  predict_model = raw_all_allbirds_model,
                                  plot_model = raw_all_allbirds_model_500,
                                  ylabel = "Sexual selection", ylimits = c(0,1.2),
                                  ybreaks =  c(0,0.5,1.0), lab_ypos = 0.12, plot_label = "b")

cert_lat_plot <- sar_lat_side_plot(ylabel = "Data certainty", ylimits = c(1,4),
                   ybreaks =  c(1,2,3,4), lab_ypos = 1.3, plot_label = "d",
                   plot_model = raw_all_certainty_model_500, 
                   stats_model = scaled_all_certainty_model,
                   predict_model = raw_all_certainty_model)

# Primary and secondary consumers.
pri_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(-0.1,2),
                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 2,
                     plot_label = "b", plot_model = raw_all_primary_model_500, 
                     stats_model = scaled_all_primary_model,
                     predict_model = raw_all_primary_model)

sec_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(-0.1,2),
                     ybreaks =  c(0,1.0,2.0), lab_ypos = 2,
                     plot_label = "f", plot_model = raw_all_secondary_model_500, 
                     stats_model = scaled_all_secondary_model,
                     predict_model = raw_all_secondary_model)

# # Trophic niche models.
fruit_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(-0.1,2),
                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 2,
                     plot_label = "d", plot_model = raw_all_fruit_model_500, 
                     stats_model = scaled_all_fruit_model,
                     predict_model = raw_all_fruit_model)

invert_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(-0.1,2),
                     ybreaks =  c(0,1.0, 2.0), lab_ypos = 2,
                     plot_label = "h", plot_model = raw_all_invert_model_500, 
                     stats_model = scaled_all_invert_model,
                     predict_model = raw_all_invert_model)

# Sexual lat gradient for high data certainty.
med_sex_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1.3),
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1.3, plot_label = "b",
                     plot_model = raw_high_allbirds_model_500, 
                     stats_model = scaled_high_allbirds_model,
                     predict_model = raw_high_allbirds_model)

hi_sex_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1.5),
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.5, plot_label = "d",
                     plot_model = raw_all_highcert_model_500, 
                     stats_model = scaled_all_highcert_model,
                     predict_model = raw_all_highcert_model)

# Make the side plots.
mig_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1.5),
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.5,  plot_label = "b",
                     plot_model = raw_all_mig_model_500, 
                     stats_model = scaled_all_mig_model,
                     predict_model = raw_all_mig_model)

no_mig_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1.5),
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.5, plot_label = "d",
                     plot_model = raw_all_nomig_model_500, 
                     stats_model = scaled_all_nomig_model,
                     predict_model = raw_all_nomig_model)

terr_lat_plot <-  sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1.6),
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.6, plot_label = "b",
                     plot_model = raw_all_terr_model_500, 
                     stats_model = scaled_all_terr_model,
                     predict_model = raw_all_terr_model)

no_terr_lat_plot <- sar_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,1.6),
                     ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 1.6, plot_label = "d",
                     plot_model = raw_all_noterr_model_500, 
                     stats_model = scaled_all_noterr_model,
                     predict_model = raw_all_noterr_model)

# Primary terr.
pri_terr_lat_plot <- sar_lat_side_plot(ylabel = "Proportion of species", ylimits = c(0,1),
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1, plot_label = "b",
                     plot_model = raw_all_primterr_model_500, 
                     stats_model = scaled_all_primterr_model,
                     predict_model = raw_all_primterr_model)
# Primary year terr.
pri_yearterr_lat_plot <-  sar_lat_side_plot(ylabel = "Proportion of species", ylimits = c(0,1),
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1, plot_label = "d",
                     plot_model = raw_all_primyearterr_model_500, p_year_terr = TRUE, 
                     stats_model = scaled_all_primyearterr_model,
                     predict_model = raw_all_primyearterr_model)

# Secondary terr.
sec_terr_lat_plot <- sar_lat_side_plot(ylabel = "Proportion of species", ylimits = c(0,1),
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1, plot_label = "f",
                     plot_model = raw_all_secterr_model_500, 
                     stats_model = scaled_all_secterr_model,
                     predict_model = raw_all_secterr_model)

# Secondary year terr.
sec_yearterr_lat_plot <-  sar_lat_side_plot(ylabel = "Proportion of species", ylimits = c(0,1),
                     ybreaks =  c(0,0.5,1.0), lab_ypos = 1, plot_label = "h",
                     plot_model = raw_all_secyearterr_model_500, 
                     stats_model = scaled_all_secyearterr_model,
                     predict_model = raw_all_secyearterr_model)



# # Export the plots.
save(list = ls(pattern =  "lat_plot"), file = "Plots/Maps/behr_500_and_100_latitudinal_sideplots.Rdata")
gc()



################################################################################
                       ##### Extract model coeffs #####


model_names <- c("all_birds",  "primary", "fruit", "secondary", "invert", "migration", 
                 "no_migration", "territorial", "non_territorial")

all_estimates <- data.frame(model = model_names, 
                            estimate = 0,
                            lower = 0,
                            upper = 0)
high_estimates <- all_estimates

extract_coeffs <- function(model){
  estimate <- last(coef(summary(model))[,1])
  se <- last(coef(summary(model))[,2])
  upper <- estimate + 1.96*se
  lower <-  estimate - 1.96*se
  c(estimate, lower, upper)
}


all_estimates[1,2:4] <- extract_coeffs(scaled_all_allbirds_model)

all_estimates[2,2:4] <- extract_coeffs(scaled_all_primary_model)
all_estimates[3,2:4] <- extract_coeffs(scaled_all_fruit_model)
all_estimates[4,2:4] <- extract_coeffs(scaled_all_secondary_model)
all_estimates[5,2:4] <- extract_coeffs(scaled_all_invert_model)

all_estimates[6,2:4] <- extract_coeffs(scaled_all_mig_model)
all_estimates[7,2:4] <- extract_coeffs(scaled_all_nomig_model)
all_estimates[8,2:4] <- extract_coeffs(scaled_all_terr_model)
all_estimates[9,2:4] <- extract_coeffs(scaled_all_noterr_model)


## High models ###


scaled_high_allbirds_model <- read_lat_model("scaled_high_allbirds")
scaled_high_primary_model <- read_lat_model("scaled_high_primary")
scaled_high_fruit_model <- read_lat_model("scaled_high_fruit")
scaled_high_secondary_model <- read_lat_model("scaled_high_secondary")
scaled_high_invert_model <- read_lat_model("scaled_high_invert")
scaled_high_mig_model <- read_lat_model("scaled_high_mig")
scaled_high_nomig_model <- read_lat_model("scaled_high_nomig")
scaled_high_terr_model <- read_lat_model("scaled_high_terr")
scaled_high_noterr_model <- read_lat_model("scaled_high_noterr")

high_estimates[1,2:4] <- extract_coeffs(scaled_high_allbirds_model)
high_estimates[2,2:4] <- extract_coeffs(scaled_high_primary_model)
high_estimates[3,2:4] <- extract_coeffs(scaled_high_fruit_model)
high_estimates[4,2:4] <- extract_coeffs(scaled_high_secondary_model)
high_estimates[5,2:4] <- extract_coeffs(scaled_high_invert_model)
high_estimates[6,2:4] <- extract_coeffs(scaled_high_mig_model)
high_estimates[7,2:4] <- extract_coeffs(scaled_high_nomig_model)
high_estimates[8,2:4] <- extract_coeffs(scaled_high_terr_model)
high_estimates[9,2:4] <- extract_coeffs(scaled_high_noterr_model)


# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(estimate, 2),
  intervals = paste0(round(lower, 2), ", ", 
                     round(upper, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Paste together values for reporting in a table.
high_estimates %<>% mutate(
  round_est = round(estimate, 2),
  intervals = paste0(round(lower, 2), ", ", 
                     round(upper, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

write.csv(all_estimates, "Results/Tables/100k_all_spatial_regression.csv", row.names = FALSE)
write.csv(high_estimates, "Results/Tables/100k_high_spatial_regression.csv", row.names = FALSE)

