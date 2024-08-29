###############################################################################
                   # Assess the fit of brms models  #
###############################################################################


# Packages to load.
library(magrittr)
library(tictoc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(brms)
library(bayesplot)
library(ggdist)
library(ggbeeswarm)
library(stringr)
library(bayestestR)
library(effectsize)

# Clear the workspace.
rm(list=ls())

###############################################################################
                      #### Read in the data #####

# Functions.
source("Code/functions.R")

# Read in the life history traits.
model_data <- read_ss_data()
  
# Scale continuous predictors to two SD.
model_data %<>% mutate(centroid_z = standardize(sqrt(abs_lat), two_sd = TRUE))

# Create 5 degree bins for latitude.
bin_range <- seq(0, 80, 5)
bin_labels <- seq(0, 75, 5) 
model_data$binned_lat <- cut(model_data$abs_lat, breaks=bin_range, labels = bin_labels, include.lowest = TRUE)
model_data$binned_lat %<>% as.character() %>% as.numeric()

# Create grouped data to reuse with different varaibles below.
grouped_lat_data <- model_data %>% filter(binned_lat != 75) %>% group_by(binned_lat)

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
lat_data <- grouped_lat_data %>% average_lat_bins()
cert_lat_data <- grouped_lat_data %>% average_lat_bins("data_certainty")

# Group the data by trophic level and trophic niche.
diet_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_level_binary) %>% average_lat_bins()
niche_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_niche) %>% average_lat_bins()

# Add in the three monogamous frugivore species that each occupy a single lat bin alone.
fruit_lat_data <- niche_lat_data %>% filter(trophic_niche == "Frugivore")

# Do migration and territoriality for plots.
mig_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, migration_binary) %>% average_lat_bins()
terr_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, territoriality_binary) %>% average_lat_bins()

# Do same for high and medium certainty.
hi_lat_data <- model_data %>% filter(binned_lat != 75 & data_certainty == 4) %>% 
  group_by(binned_lat) %>% average_lat_bins()
med_lat_data <- model_data %>% filter(binned_lat != 75 & data_certainty > 2) %>% 
  group_by(binned_lat) %>% average_lat_bins()

# Group the data by trophic level territory and latitude.
terr_diet_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_level_binary) %>% average_lat_bins("terr_dummy")
year_terr_diet_lat_data <- model_data %>% filter(binned_lat != 75) %>% 
  group_by(binned_lat, trophic_level_binary) %>% average_lat_bins("year_terr_dummy")




###############################################################################
                  #### Read in the models #####

# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Consensus/Latitude/ordinal_"

# Function to read in models. (Haven't set up for this script.)
read_lat_model <- function(model_name, certainty = "all"){
  model_file <- paste0(first_half, model_name, "_uncentered_", certainty,  "_1.rds")
  readRDS(model_file)
}


# Trophic models.
conphy_primary_all_model <- read_lat_model("primary")
conphy_fruit_all_model <- read_lat_model("fruit")
conphy_secondary_all_model <-  read_lat_model("secondary")
conphy_invert_all_model <- read_lat_model("invert")

conphy_mig_all_model <- read_lat_model("mig")
conphy_noterr_all_model <- read_lat_model("no_terr")


# High certainty.
conphy_allbirds_high_model <- read_lat_model("allbirds", "high")

conphy_primary_high_model <- read_lat_model("primary", "high")
conphy_fruit_high_model <- read_lat_model("fruit", "high")
conphy_secondary_high_model <-  read_lat_model("secondary", "high")
conphy_invert_high_model <- read_lat_model("invert", "high")


# Directory where models are saved.
first_half <- "Z:/home/sexual_selection/Results/Models/Old_models/Nonphy_models/Latitude/"
list.files(first_half)

# Function to read in models.
read_nonphy_model <- function(model_name){
  model_file <- paste0(first_half, model_name, "_model.rds")
  readRDS(model_file)
}

nonphy_allbirds_all_model <- read_nonphy_model("allbirds_all")
nonphy_certainty_all_model <- read_nonphy_model("certainty_all")
nonphy_primary_all_model <- read_nonphy_model("primary_all")
nonphy_fruit_all_model <- read_nonphy_model("fruit_all")
nonphy_secondary_all_model <- read_nonphy_model("secondary_all")
nonphy_invert_all_model <- read_nonphy_model("invert_all")

nonphy_allbirds_high_model <- read_nonphy_model("allbirds_high")
nonphy_certainty_high_model <- read_nonphy_model("certainty_high")
nonphy_primary_high_model <- read_nonphy_model("primary_high")
nonphy_fruit_high_model <- read_nonphy_model("fruit_high")
nonphy_secondary_high_model <- read_nonphy_model("secondary_high")
nonphy_invert_high_model <- read_nonphy_model("invert_high")












# Function to read in models. (Haven't set up for this script.)
read_trop_model <- function(model_name, certainty = "all"){
  trop_first_half <- "Z:/home/sexual_selection/Results/Models/Consensus/Tropics/ordinal_"
  model_file <- paste0(trop_first_half, model_name, "_", certainty,  ".rds")
  readRDS(model_file)
}


fruit_all_trop_model <- read_trop_model("fruit")
primary_all_trop_model <- read_trop_model("fruit")



###############################################################################
                         #### side plot function  ######

# Change scipen. (Probs not needed now)
options(scipen = 999)

# Functions to make.

# Change p value to a string, using standard thresholds. 
p_to_string <- function(p_value){
  if (p_value < 0.001 ){
    p_value <- "p < 0.001"
  } else if (p_value < 0.01) {
    p_value <- "p < 0.01"
  } else if (p_value < 0.05) {
    p_value <- "p < 0.05"
  } else {
    p_value <- paste0("p = ", as.character(format(round(p_value, 2), nsmall = 2)))
  }
}



# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
brms_lat_side_plot <- function(data_set, ylabel = "", ylimits = c(0,1.1), ybreaks = c(0,0.5,1), 
                               lab_x_pos = 60, lab_ypos = 1, plot_label = "b", 
                               plot_model = allbirds_all_model,
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
  estimate <- last(summary(plot_model)$fixed[,1])
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- last(summary(plot_model)$fixed[,1])
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  # Lower. 
  lower <- last(summary(plot_model)$fixed[,3])
  lower <- as.character(format(round(lower, 2), nsmall = 2))
  
  # Upper
  upper <- last(summary(plot_model)$fixed[,4])
  upper <- as.character(format(round(upper, 2), nsmall = 2))
  
  # Change size of CI.
  if (lower == "0.00" | upper == "0.00"){
    lower <- last(summary(plot_model)$fixed[,3])
    lower <- as.character(format(round(lower, 3), nsmall = 3))
    upper <- last(summary(plot_model)$fixed[,4])
    upper <- as.character(format(round(upper, 3), nsmall = 3))
  }
  
  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)
  intervals <- paste0("[", lower, ", ", upper, "]")
  
  # # Extract p-values using probability of direction two-tailed test.
  p_value <- p_to_string(pd_to_p(last(p_direction(plot_model)[,2])))
  
  # Create a label 
  if (p_include){
    stats_label <- paste0(estimate, "\n", p_value)
  } else {
    stats_label <- paste0(estimate, "\n", intervals)
  }
  
  # ggplot function for sideplots with annotations.
  ggplot(data_set, aes(x = binned_lat, y = trait_mean)) +
    geom_errorbar(aes(ymin = trait_min, ymax = trait_max),
                  position = position_dodge(width = 1), show.legend = FALSE, col = "darkgrey") + #    col =  "darkgrey") +
    #geom_point(position = position_dodge(width = 1), col = "#442B48") + 
    geom_point(aes(size = trait_n), position = position_dodge(width = 1), alpha = 0.9, col = "#442B48") + 
    scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 75), clip = 'off') +
    scale_size_continuous(range = c(1,6))+
    #scale_alpha_continuous(range = c(0.1,1)) +
    ylab(ylabel) +
    xlab("Latitude") + theme_classic(base_size = 25) + 
    theme(legend.position = "none", 
          #text = element_text(face = "bold"),
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7) + # fontface = 2) +
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
    geom_line(data = predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1)
}


# Primary and secondary consumers.
conphy_pri_lat_plot <- diet_lat_data %>% filter(trophic_level_binary == "Primary") %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(-0.1, 4), 
                     ybreaks =  c(0, 1.0, 2.0, 3.0, 4.0), lab_x_pos = 20, lab_ypos = 1.5, 
                     plot_label = "b", plot_model = conphy_primary_all_model)

conphy_sec_lat_plot <- diet_lat_data %>% filter(trophic_level_binary == "Secondary") %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(-0.1, 4), 
                     ybreaks =  c(0, 1.0, 2.0, 3.0, 4.0), lab_x_pos = 20, lab_ypos = 1.5,
                     plot_label = "f", plot_model = conphy_secondary_all_model)

# Trophic niche models.
conphy_fruit_lat_plot <- fruit_lat_data %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(-0.1, 4), 
                     ybreaks =  c(0, 1.0, 2.0, 3.0, 4.0), lab_x_pos = 20, lab_ypos = 2.5, 
                     plot_label = "d", plot_model = conphy_fruit_all_model)

conphy_invert_lat_plot <- niche_lat_data %>% filter(trophic_niche == "Invertivore") %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(-0.1, 4), 
                     ybreaks =  c(0, 1.0, 2.0, 3.0, 4.0), lab_x_pos = 20, lab_ypos = 1.5, 
                     plot_label = "h", plot_model = conphy_invert_all_model) 












###############################################################################
              #### Actually make the side plots ####



# # Two last models.
# pri_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
#   brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
#                      ybreaks =  c(0,1.0, 2.0), lab_x_pos = 10, lab_ypos = 1.5, 
#                      plot_label = "b", plot_model_data = prim_all_data) 
# 
# sec_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
#   brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
#                      ybreaks =  c(0,1.0,2.0), lab_x_pos = 10, lab_ypos = 1.5,
#                      plot_label = "f", plot_model_data = sec_all_data) 



fruit_lat_data <- niche_lat_data %>% filter(trophic_niche == "Frugivore")
fruit_lat_data[9,] <- list(40, "Frugivore", 0, 0, 0, 0, 0, 1)  # Add in three monogamous frugivore species at higher lat.
fruit_lat_data[10,] <- list(45, "Frugivore", 0, 0, 0, 0, 0, 1)
fruit_lat_data[11,] <- list(50, "Frugivore", 0, 0, 0, 0, 0, 1)

# # Trophic niche models.
# fruit_lat_plot <- fruit_lat_data %>% 
#   brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
#                      ybreaks =  c(0,1.0, 2.0), lab_x_pos = 10, lab_ypos = 1.5, 
#                      plot_label = "d", plot_model_data = frug_all_data) 
# 
# invert_lat_plot <- niche_lat_data %>% filter(trophic_niche == "Invertivore") %>% 
#   brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
#                      ybreaks =  c(0,1.0, 2.0), lab_x_pos = 10, lab_ypos = 1.5, 
#                      plot_label = "h", plot_model_data = invert_all_data) 



################################################################################
                       ##### Plot all together #######

ylabel = "Sexual selection"
ylimits = c(0,4)
ybreaks =  c(0,1.0, 2.0)
lab_x_pos = 10
lab_ypos = 1.5
plot_label = "h"


# Change sexual score back to original values (have to + 1 for ordinal regression)
all_prim_predictions <- prim_all_data[[2]]
all_frug_predictions <- frug_all_data[[2]]
all_sec_predictions <- sec_all_data[[2]]
all_invert_predictions <- invert_all_data[[2]]

all_prim_predictions$estimate__ <- all_prim_predictions$estimate__ - 1
all_frug_predictions$estimate__ <- all_frug_predictions$estimate__ - 1
all_sec_predictions$estimate__ <- all_sec_predictions$estimate__ - 1
all_invert_predictions$estimate__ <- all_invert_predictions$estimate__ - 1

ggplot(lat_data, aes(x = binned_lat, y = trait_mean)) +
  # geom_errorbar(aes(ymin = trait_min, ymax = trait_max), 
  #               position = position_dodge(width = 1), show.legend = FALSE, col =  "darkgrey") + 
  # geom_point(position = position_dodge(width = 1), col = "black") + 
  # scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
  scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
  coord_cartesian(ylim = ylimits, xlim = c(NA, 75), clip = 'off') +
  ylab(ylabel) +
  xlab("Latitude") + theme_classic(base_size = 25) + 
  theme(legend.position = "none", 
        text = element_text(face = "bold"),
        axis.title.y = element_text(size = rel(0.85)),
        axis.title.x = element_text(size = rel(0.85)),
        plot.margin = margin(t = 1, l = 0.2, b = 0.2, unit = "cm")) + 
  #annotate("text", x = lab_x_pos, y =lab_ypos, label = cor_label, size = 7, fontface = 2) +
  annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
  # geom_ribbon(data = plot_predictions, 
  #             aes(x = complete_latitude,  ymin = (lower__ - 1), ymax = (upper__ - 1)), 
  #             fill = "grey70", colour = NA, alpha = 0.2, inherit.aes = FALSE)  + 
  
  geom_line(data = all_prim_predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1, colour = "#05299E") +
  geom_line(data = all_frug_predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1, colour = "#214F4B") +
  geom_line(data = all_sec_predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1, colour = "#8D0801") +
  geom_line(data = all_invert_predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1, colour = "#A88A05")






   
################################################################################
              #### Try family averages ######


# Read in clade function and assign.
source("Code/clade_function.R")
model_data %<>% assign_prum_clades()

model_data$abs_lat <- abs(model_data$complete_latitude)

# Create settings for scatter plots.
point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                      size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                       labs(x = "", y = NULL,  colour = "Clade"),
                       scale_colour_manual(values = prum_clade_colours),
                       scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(0,4.2)),
                       theme_classic(base_size = 20),
                       theme(text = element_text(face = "bold"),
                             legend.position = "none",
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.title = element_text(size = rel(0.8))))


# Filter for each consumer type.
pri_fam_data <- model_data %>% filter(trophic_binary == "Primary") 
sec_fam_data <- model_data %>% filter(trophic_binary == "Secondary") 
fruit_fam_data <- model_data %>% filter(trophic_niche == "Frugivore")
invert_fam_data <- model_data %>% filter(trophic_niche == "Invertivore") 

# Average traits for family.
family_average_data_data <- pri_fam_data %>% group_by(family) %>% 
  summarise(higher_clade = first(higher_clade),
            tree_tip = first(tree_tip),
            sexual_mean = mean(sexual_score),
            clade_sum = length(sexual_score),
            abs_lat = mean(abs_lat))

family_average_data_data <- sec_fam_data %>% group_by(family) %>% 
  summarise(higher_clade = first(higher_clade),
            tree_tip = first(tree_tip),
            sexual_mean = mean(sexual_score),
            clade_sum = length(sexual_score),
            abs_lat = mean(abs_lat))

# Make some nice plots.
family_average_data_data %>% ggplot(aes(x=abs_lat, y=sexual_mean)) + point_settings +
  xlab("Latitude") + 
  geom_ribbon(data = all_sec_predictions, inherit.aes = FALSE,
              aes(x = abs_lat, ymin = lower__-1, ymax = upper__-1), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = all_sec_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1)





