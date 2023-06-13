###############################################################################
              # Plotting Latitudinal brms models from HPC  #
###############################################################################

# This script creates the side plots which are used alongside maps produced in
# another script. It uses models that have been run on the HPC.


# Packages to load.
library(magrittr)
library(tictoc)
library(caper)
library(dplyr)
library(effectsize)
library(car)
library(ggplot2)
library(ggpubr)
library(phytools)
library(brms)
library(graph4lg)

# Set the seed.
set.seed(1993)

###############################################################################
                       #### Read in the data #####

# Functions.
source("Code/functions.R")

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)
model_data$abs_lat <- abs(model_data$complete_latitude)


###############################################################################
              #### Prepare predictor variables ######

# Scale continuous predictors to two SD.
model_data %<>% mutate(centroid_z = standardize(centroid_sqrt, two_sd = TRUE))

# Create 5 degree bins for latitude.
bin_range <- seq(0, 80, 5)
bin_labels <- seq(0, 75, 5) 
model_data$binned_lat <- cut(abs(model_data$complete_latitude), breaks=bin_range, labels = bin_labels, include.lowest = TRUE)
model_data$binned_lat %<>% as.character() %>% as.numeric()

# Create grouped data to reuse with different varaibles below.
grouped_lat_data <- model_data %>% filter(binned_lat != 75) %>% group_by(binned_lat)

# Function for grouping sexual selection by trait.
average_lat_bins <- function(grouped_data, predictor = "sexual_score"){
  grouped_data %>% 
    summarise(trait_mean = mean(!!! syms(predictor)),
              trait_sd = sd(!!! syms(predictor)),
              trait_se = sd(!!! syms(predictor))/sqrt(length(!!! syms(predictor))),
              trait_max = trait_mean + trait_se,
              trait_min = trait_mean - trait_se,
              trait_n = length(!!! syms(predictor))) %>% na.omit()
}

# Group sex scores and certainty by lat bins.
lat_data <- grouped_lat_data %>% average_lat_bins()
cert_lat_data <- grouped_lat_data %>% average_lat_bins("cert_reverse")

# Group the data by trophic level and trophic niche.
diet_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% average_lat_bins()
niche_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_niche) %>% average_lat_bins()

# Add in the three monogamous frugivore species that each occupy a single lat bin alone.
fruit_lat_data <- niche_lat_data %>% filter(trophic_niche == "Frugivore")
fruit_lat_data[9,] <- list(40, "Frugivore", 0, 0, 0, 0, 0, 1) 
fruit_lat_data[10,] <- list(45, "Frugivore", 0, 0, 0, 0, 0, 1)
fruit_lat_data[11,] <- list(50, "Frugivore", 0, 0, 0, 0, 0, 1)

# Do same for high and medium certainty.
hi_lat_data <- model_data %>% filter(binned_lat != 75 & sexual_certainty == 1) %>% 
  group_by(binned_lat) %>% average_lat_bins()
med_lat_data <- model_data %>% filter(binned_lat != 75 & sexual_certainty < 3) %>% 
  group_by(binned_lat) %>% average_lat_bins()

# Group the data by trophic level territory and latitude.
terr_diet_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% average_lat_bins("terr_dummy")
terr_diet_lat_data$trait_min[terr_diet_lat_data$trait_min < 0] <- 0   # Get rid of negative error bars.
year_terr_diet_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_binary) %>% average_lat_bins("year_terr_dummy")



###############################################################################
                    #### Read in  brms models ####

# Read in models using raw data.
allbirds_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/allbirds_model.rds")
cert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/cert_model.rds")
primary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/primary_model.rds")
secondary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/secondary_model.rds")
fruit_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/fruit_model.rds")
invert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/invert_model.rds")


# Read in centered models.
centered_allbirds_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_allbirds_model.rds")
centered_cert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_cert_model.rds")
centered_primary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_primary_model.rds")
centered_secondary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_secondary_model.rds")
centered_fruit_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_fruit_model.rds")
centered_invert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_invert_model.rds")


################################################################################
                    #### Export summary tables ####


# Extract relevant coeffcient information.
allbirds_all_estimates <- summary(centered_allbirds_model)$fixed[5,c(1,3,4)]
cert_all_estimates <- summary(centered_cert_model)$fixed[4,c(1,3,4)]
primary_all_estimates <- summary(centered_primary_model)$fixed[5,c(1,3,4)]
fruit_all_estimates <- summary(centered_fruit_model)$fixed[5,c(1,3,4)]
secondary_all_estimates <- summary(centered_secondary_model)$fixed[5,c(1,3,4)]
invert_all_estimates <- summary(centered_invert_model)$fixed[5,c(1,3,4)]
all_estimates <- rbind(allbirds_all_estimates, cert_all_estimates, 
                       primary_all_estimates, fruit_all_estimates,
                       secondary_all_estimates, invert_all_estimates)
row.names(all_estimates) <- c("all_birds", "certainty", "primary",
                              "fruit", "secondary", "invert")
# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0("[", round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2), "]"),
  est_intervals = paste0(round_est, " ", intervals))

# Export the results.
write.csv(all_estimates, "Results/Tables/all_nonphy_lat_regression.csv", row.names = TRUE)




###############################################################################
                    #### side plot function  ######

# Change scipen. (Probs not needed now)
options(scipen = 999)

# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
brms_lat_side_plot_2 <- function(data_set, ylabel = "", ylimits = c(0,1.1), ybreaks = c(0,0.5,1), 
                               lab_x_pos = 60, lab_ypos = 1, plot_label = "b", 
                               plot_model = allbirds_model, stats_model = centered_allbirds_model,
                               sex_score = TRUE, r_include = FALSE, p_include = TRUE){
  
  # Sample size.   # Not sure we need.
  sample_size <- nrow(plot_model$data)
  
  # Extract predictions from brms model.
  predictions <- conditional_effects(plot_model)[[1]] 
  
  # Change back to normal data values that were transformed for ordinal models..
  if (sex_score){
     predictions %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
  }
  
  # Estimate
  estimate <- last(summary(stats_model)$fixed[,1])
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- last(summary(stats_model)$fixed[,1])
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  # Lower. 
  lower <- last(summary(stats_model)$fixed[,3])
  lower <- as.character(format(round(lower, 2), nsmall = 2))
  
  # Upper
  upper <- last(summary(stats_model)$fixed[,4])
  upper <- as.character(format(round(upper, 2), nsmall = 2))
  
  # Change size of CI.
  if (lower == "0.00" | upper == "0.00"){
    lower <- last(summary(stats_model)$fixed[,3])
    lower <- as.character(format(round(lower, 3), nsmall = 3))
    upper <- last(summary(stats_model)$fixed[,4])
    upper <- as.character(format(round(upper, 3), nsmall = 3))
  }
  
  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)
  intervals <- paste0("[", lower, ", ", upper, "]")
  
  # # Extract p-values using probability of direction two-tailed test.
  p_value <- pd_to_p(last(p_direction(stats_model)[,2]))
  
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
  if (p_include){
    stats_label <- paste0(estimate, "\n", p_value)
  } else {
    stats_label <- paste0(estimate, "\n", intervals)
  }

  # ggplot function for sideplots with annotations.
  ggplot(data_set, aes(x = binned_lat, y = trait_mean)) +
    geom_errorbar(aes(ymin = trait_min, ymax = trait_max), 
                  position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
    geom_point(position = position_dodge(width = 1), col = "black") + 
    scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 75), clip = 'off') +
    ylab(ylabel) +
    xlab("Latitude") + theme_classic(base_size = 25) + 
    theme(legend.position = "none", 
          text = element_text(face = "bold"),
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7, fontface = 2) +
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
    geom_line(data = predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1)
}


################################################################################
                     #### Make the side plots  ####

# Main sexual selection plot.
sex_lat_plot <- lat_data %>% 
  brms_lat_side_plot_2(ylabel = "Sexual selection", ylimits = c(0,1.2), 
                       ybreaks =  c(0,0.5,1.0), lab_ypos = 0.2, plot_label = "b")

# Data certainty.
cert_lat_plot <- cert_lat_data %>% 
  brms_lat_side_plot_2(ylabel = "Data certainty", ylimits = c(1,4), 
                     ybreaks =  c(1,2,3,4), lab_ypos = 1.5, plot_label = "d",
                     plot_model = cert_model, stats_model = centered_cert_model,
                     sex_score = FALSE)

# Primary and secondary consumers.
pri_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  brms_lat_side_plot_2(ylabel = "Sexual selection", ylimits = c(0,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_x_pos = 20, lab_ypos = 1.5, 
                     plot_label = "b", plot_model = primary_model,
                     stats_model = centered_primary_model) 

sec_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  brms_lat_side_plot_2(ylabel = "Sexual selection", ylimits = c(0,2), 
                     ybreaks =  c(0,1.0,2.0), lab_x_pos = 20, lab_ypos = 1.5,
                     plot_label = "f", plot_model = secondary_model, 
                     stats_model = centered_secondary_model) 

# Trophic niche models.
fruit_lat_plot <- fruit_lat_data %>% 
  brms_lat_side_plot_2(ylabel = "Sexual selection", ylimits = c(0,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_x_pos = 20, lab_ypos = 1.5, 
                     plot_label = "d", plot_model = fruit_model, 
                     stats_model = centered_fruit_model) 
invert_lat_plot <- niche_lat_data %>% filter(trophic_niche == "Invertivore") %>% 
  brms_lat_side_plot_2(ylabel = "Sexual selection", ylimits = c(0,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_x_pos = 20, lab_ypos = 1.5, 
                     plot_label = "h", plot_model = invert_model, 
                     stats_model = centered_invert_model) 


################################################################################
        #### Add side plots used in supplementary figures ####


# Read in models using raw data.
med_cert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/med_cert_model.rds")
hi_cert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/hi_cert_model.rds")
terr_primary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/terr_primary_model.rds")
yearterr_primary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/yearterr_primary_model.rds")
terr_secondary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/terr_secondary_model.rds")
yearterr_secondary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/yearterr_secondary_model.rds")

# Read in centered models.
centered_med_cert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_med_cert_model.rds")
centered_hi_cert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_hi_cert_model.rds")
centered_terr_primary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_terr_primary_model.rds")
centered_yearterr_primary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_yearterr_primary_model.rds")
centered_terr_secondary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_terr_secondary_model.rds")
centered_yearterr_secondary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Latitude/centered_yearterr_secondary_model.rds")


# Sexual lat gradient for high data certainty.
med_sex_lat_plot <- med_lat_data %>% 
  brms_lat_side_plot_2(ylabel = "Sexual selection", ylimits = c(0,1.2), 
                       ybreaks =  c(0,0.5,1.0), lab_ypos = 0.2, plot_label = "b",
                       plot_model = med_cert_model, stats_model = centered_med_cert_model,
                       sex_score = TRUE)

hi_sex_lat_plot <- hi_lat_data %>% 
  brms_lat_side_plot_2(ylabel = "Sexual selection", ylimits = c(0,1.8), 
                       ybreaks =  c(0,0.5,1.0, 1.5), lab_ypos = 0.3, plot_label = "d",
                       plot_model = hi_cert_model, stats_model = centered_hi_cert_model,
                       sex_score = TRUE)


# Primary terr.
pri_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  brms_lat_side_plot_2(ylabel = "Proportion of species", ylimits = c(0,1), 
                       ybreaks =  c(0,0.5,1.0), lab_ypos = 0.9, plot_label = "b",
                       plot_model = terr_primary_model, stats_model = centered_terr_primary_model,
                       sex_score = FALSE)
# Primary year terr.
pri_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  brms_lat_side_plot_2(ylabel = "Proportion of species", ylimits = c(0,1), 
                       ybreaks =  c(0,0.5,1.0), lab_ypos = 0.9, plot_label = "d",
                       plot_model = yearterr_primary_model, stats_model = centered_yearterr_primary_model,
                       sex_score = FALSE)

# Secondary terr.
sec_terr_lat_plot <- terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  brms_lat_side_plot_2(ylabel = "Proportion of species", ylimits = c(0,1), 
                       ybreaks =  c(0,0.5,1.0), lab_ypos = 0.9, plot_label = "f",
                       plot_model = terr_secondary_model, stats_model = centered_terr_secondary_model,
                       sex_score = FALSE)

# Secondary year terr.
sec_yearterr_lat_plot <- year_terr_diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  brms_lat_side_plot_2(ylabel = "Proportion of species", ylimits = c(0,1), 
                       ybreaks =  c(0,0.5,1.0), lab_ypos = 0.9, plot_label = "h",
                       plot_model = yearterr_secondary_model, stats_model = centered_yearterr_secondary_model,
                       sex_score = FALSE)

# Export the plots.
save(list = ls(pattern =  "lat_plot"), file = "Plots/Maps/latitudinal_sideplots.Rdata")
gc()


