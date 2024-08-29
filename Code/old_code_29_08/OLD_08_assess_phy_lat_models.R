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

# Clear the workspace.
rm(list=ls())

# Functions.
source("Code/functions.R")


###############################################################################
                  #### Read in the models #####

# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Combined_models/Latitude/"

# Read in the models.
allbirds_all_data <- readRDS(paste0(first_half, "all_uncentered_all_data.rds"))
#cert_all_data <- readRDS(paste0(first_half, "certainty_uncentered_all_data.rds"))

# Trophic models.
prim_all_data <- readRDS(paste0(first_half, "primary_uncentered_all_data.rds"))
frug_all_data <- readRDS(paste0(first_half, "frugivore_uncentered_all_data.rds"))
sec_all_data <- readRDS(paste0(first_half, "secondary_uncentered_all_data.rds"))
invert_all_data <- readRDS(paste0(first_half, "invertivore_uncentered_all_data.rds"))



################################################################################
                    #### Export summary tables ####


# Extract relevant coeffcient information.
allbirds_all_estimates <- summary(allbirds_all_data[[1]])$fixed[5,c(1,3,4)]
#cert_all_estimates <- summary(centered_cert_model)$fixed[4,c(1,3,4)]

primary_all_estimates <- summary(prim_all_data[[1]])$fixed[5,c(1,3,4)]
fruit_all_estimates <- summary(frug_all_data[[1]])$fixed[5,c(1,3,4)]
secondary_all_estimates <- summary(sec_all_data[[1]])$fixed[5,c(1,3,4)]
invert_all_estimates <- summary(invert_all_data[[1]])$fixed[5,c(1,3,4)]

all_estimates <- rbind(allbirds_all_estimates, #cert_all_estimates, 
                       primary_all_estimates, fruit_all_estimates,
                       secondary_all_estimates, invert_all_estimates)
row.names(all_estimates) <- c("all_birds", #"certainty", 
                              "primary",
                              "fruit", "secondary", "invert")
# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export the results.
write.csv(all_estimates, "Results/Tables/all_phy_lat_regression.csv", row.names = TRUE)



################################################################################
                        ##### High certainty ######


# Clear the workspace due to memory demands.
rm(list=ls())
gc()

#source("Code/functions.R")

# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Combined_models/Latitude/"

# Read in the models.
allbirds_high_data <- readRDS(paste0(first_half, "all_uncentered_high_data.rds"))
prim_high_data <- readRDS(paste0(first_half, "primary_uncentered_high_data.rds"))
frug_high_data <- readRDS(paste0(first_half, "frugivore_uncentered_high_data.rds"))
sec_high_data <- readRDS(paste0(first_half, "secondary_uncentered_high_data.rds"))
invert_high_data <- readRDS(paste0(first_half, "invertivore_uncentered_high_data.rds"))

# Extract relevant coeffcient information.
allbirds_high_estimates <- summary(allbirds_high_data[[1]])$fixed[5,c(1,3,4)]
primary_high_estimates <- summary(prim_high_data[[1]])$fixed[5,c(1,3,4)]
fruit_high_estimates <- summary(frug_high_data[[1]])$fixed[5,c(1,3,4)]
secondary_high_estimates <- summary(sec_high_data[[1]])$fixed[5,c(1,3,4)]
invert_high_estimates <- summary(invert_high_data[[1]])$fixed[5,c(1,3,4)]

high_estimates <- rbind(allbirds_high_estimates,
                       primary_high_estimates, fruit_high_estimates,
                       secondary_high_estimates, invert_high_estimates)
row.names(high_estimates) <- c("all_birds", 
                              "primary",
                              "fruit", "secondary", "invert")

# Paste together values for reporting in a table.
high_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export the results.
write.csv(high_estimates, "Results/Tables/high_phy_lat_regression.csv", row.names = TRUE)










###############################################################################
                     #### Bin latitude for plots ####

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)

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


###############################################################################
             #### Function to make the side plots ####

brms_lat_side_plot <- function(data_set, ylabel = "", ylimits = c(0,1.1), ybreaks = c(0,0.5,1), 
                               lab_x_pos = 60, lab_ypos = 1, plot_label = "b", 
                               plot_model_data = all_frug_data, sex_score = TRUE,
                               r_include = FALSE){
  
  # Extract predictions from brms model.
  predictions <- plot_model_data[[2]]
  
  # Change sexual score back to original values (have to + 1 for ordinal regression)
  if (sex_score){
    predictions$estimate__ <- predictions$estimate__ - 1
  }
  plot_model <- plot_model_data[[1]]
  # Extract p-values.
  p_values <- brms_pmap(plot_model)[2]
  p_value <- last(p_values[,1])
  
  # Extract r-squared.
  r_squared <- plot_model_data[[3]]
  
  # Sample size.
  sample_size <- nrow(plot_model$data)
  
  # Estimate
  estimate <- summary(plot_model)$fixed["abs_lat",1]
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- summary(plot_model)$fixed["abs_lat",1]
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  
  # To include r-squared values in plot.
  if (r_include){
    # With R squared, which is super low without using phylogeny.
    if(p_value < 0.001){
      cor_label <- paste0("p < 0.001", "\nR\u00b2 = ", as.character(format(round(r_squared, 2), nsmall = 2)))
    } else if (p_value > 0.005) {
      cor_label <- paste0("p = ", as.character(format(round(p_value, 2), nsmall = 2)), "\nR\u00b2 = ", as.character(format(round(r_squared, 2), nsmall = 2)))
    } else {
      cor_label <- paste0("p = ", as.character(format(round(p_value, 3), nsmall = 3)), "\nR\u00b2 = ", as.character(format(round(r_squared, 2), nsmall = 2)))
    }
  } else {
    # Without R squared
    if(p_value < 0.001){
      cor_label <- paste0("p < 0.001", "\n\U03B2 = ", estimate)
      cor_label <- paste0("\U03B2 = ", estimate, "\np < 0.001")
    } else if (p_value > 0.005) {
      cor_label <- paste0("\U03B2 = ", estimate, "\np = ", as.character(format(round(p_value, 2), nsmall = 2)))
    } else {
      cor_label <- paste0("\U03B2 = ", estimate, "\np = ", as.character(format(round(p_value, 3), nsmall = 3)))
    }
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
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = cor_label, size = 7, fontface = 2) +
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
    # geom_ribbon(data = plot_predictions, 
    #             aes(x = complete_latitude,  ymin = (lower__ - 1), ymax = (upper__ - 1)), 
    #             fill = "grey70", colour = NA, alpha = 0.2, inherit.aes = FALSE)  + 
    geom_line(data = predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1)
}

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





