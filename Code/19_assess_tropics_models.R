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
#set.seed(1993)

###############################################################################
                       #### Read in the data #####

# Functions.
source("Code/functions.R")

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")

# Add the tropical non-tropical models.
model_data$trop_non_trop <- NA
model_data$trop_non_trop[abs(model_data$complete_latitude) < 23.43624] <- "trop"
model_data$trop_non_trop[is.na(model_data$trop_non_trop)] <- "non_trop"

# Filter for dietary data.
primary_data <- model_data %>% filter(trophic_binary == "Primary")
secondary_data <- model_data %>% filter(trophic_binary == "Secondary")
fruit_data <- model_data %>% filter(trophic_niche == "Frugivore")
invert_data <- model_data %>% filter(trophic_niche == "Invertivore")

# Filter for eco roles.
mig_data <- model_data %>% filter(migration_binary == "Strong")
non_mig_data <- model_data %>% filter(migration_binary == "Weak")
terr_data <- model_data %>% filter(territory_binary == "Territory")
non_terr_data <- model_data %>% filter(territory_binary == "No territory")


###############################################################################
                    #### Read in  brms models ####

# Read in models using raw data.
primary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/primary_model.rds")
secondary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/secondary_model.rds")
fruit_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/fruit_model.rds")
invert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/invert_model.rds")

mig_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/mig_model.rds")
non_mig_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/non_mig_model.rds")
terr_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/terr_model.rds")
non_terr_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/non_terr_model.rds")

# Read in centered models.
centered_primary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/centered_primary_model.rds")
centered_secondary_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/centered_secondary_model.rds")
centered_fruit_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/centered_fruit_model.rds")
centered_invert_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/centered_invert_model.rds")

centered_mig_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/centered_mig_model.rds")
centered_non_mig_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/centered_non_mig_model.rds")
centered_terr_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/centered_terr_model.rds")
centered_non_terr_model <- readRDS("Z:/home/sexual_selection/Results/Models/Nonphy_models/Tropics/centered_non_terr_model.rds")



################################################################################
                    #### Export summary tables ####

# 
# # Extract relevant coeffcient information.
# primary_all_estimates <- summary(centered_primary_model)$fixed[5,c(1,3,4)]
# fruit_all_estimates <- summary(centered_fruit_model)$fixed[5,c(1,3,4)]
# secondary_all_estimates <- summary(centered_secondary_model)$fixed[5,c(1,3,4)]
# invert_all_estimates <- summary(centered_invert_model)$fixed[5,c(1,3,4)]
# all_estimates <- rbind(allbirds_all_estimates, cert_all_estimates, 
#                        primary_all_estimates, fruit_all_estimates,
#                        secondary_all_estimates, invert_all_estimates)
# row.names(all_estimates) <- c("all_birds", "certainty", "primary",
#                               "fruit", "secondary", "invert")
# # Paste together values for reporting in a table.
# all_estimates %<>% mutate(
#   round_est = round(Estimate, 2),
#   intervals = paste0("[", round(`l-95% CI`, 2), ", ", 
#                      round(`u-95% CI`, 2), "]"),
#   est_intervals = paste0(round_est, " ", intervals))
# 
# # Export the results.
# write.csv(all_estimates, "Results/Tables/all_nonphy_lat_regression.csv", row.names = TRUE)

###############################################################################
                    #### side plot function  ######
options(scipen = 999)

# Function that recreates side plots using both pseudo p-values and credible intervals from centered models.
brms_tropics_side_plot <- function(data_set = primary_data, ylabel = "Sexual selection", 
                                   ylimits = c(0,1.15), ybreaks = c(0,0.5,1.0), 
                                   lab_x_pos = 2, lab_ypos = 1.15, plot_label = "a", 
                                   plot_model = primary_model, stats_model = centered_primary_model,
                                   sex_score = TRUE, r_include = FALSE, x_label = "", 
                                   p_include = TRUE){
  
  # Extract predictions from brms model.
  predictions <- conditional_effects(plot_model)[[1]] 
  predictions %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
  
  # Sample size.
  sample_sizes <- data_set %>% count(trop_non_trop)
  
  # Estimate
  estimate <- summary(plot_model)$fixed[5,1]
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- summary(plot_model)$fixed[5,1]
    estimate <- as.character(format(round(estimate, 3), nsmall = 3))
  }
  
  # Lower. 
  lower <- summary(plot_model)$fixed[5,3]
  lower <- as.character(format(round(lower, 2), nsmall = 2))
  
  # Upper
  upper <- summary(plot_model)$fixed[5,4]
  upper <- as.character(format(round(upper, 2), nsmall = 2))
  
  # Change size of CI.
  if (lower == "0.00" | upper == "0.00"){
    lower <- summary(plot_model)$fixed[5,3]
    lower <- as.character(format(round(lower, 3), nsmall = 3))
    upper <- summary(plot_model)$fixed[5,4]
    upper <- as.character(format(round(upper, 3), nsmall = 3))
  }

  # Paste together estimate and CIs as a string.
  estimate <- paste0("\U03B2 = ", estimate)
  intervals <- paste0("[", lower, ", ", upper, "]")
  
  # # Extract p-values using probability of direction two-tailed test.
  p_value <- pd_to_p(last(p_direction(plot_model)[,2]))
  
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
  ggplot(predictions, aes(x = trop_non_trop, y = estimate__)) +
    geom_errorbar(aes(ymin = lower__, ymax = upper__), 
                  position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) + 
    
    geom_point(position = position_dodge(width = 1), col = "black", size = 3.6) + 
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    scale_x_discrete(limits = c("trop", "non_trop"), labels = c("Tropical", "Non-tropical")) +
    coord_cartesian(ylim = ylimits, clip = 'off') +
    ylab(ylabel) +
    xlab(x_label) + theme_classic(base_size = 25) + 
    theme(legend.position = "none",
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7) +
    annotate("text", x = 0.65, y = ylimits[2], label = plot_label, size = 12, fontface = 2) +
    annotate("text", x = 1, y = predictions$upper__[1] + 0.05, label = sample_sizes[2,2], size = 5) +
    annotate("text", x = 2, y = predictions$upper__[2] + 0.05, label = sample_sizes[1,2], size = 5)
}



################################################################################
              #### Make with centered models for stats ####

# Trophic niche plots.
primary_plot <- brms_tropics_side_plot(
  data_set = primary_data, plot_label = "a", plot_model = primary_model,
  stats_model = centered_primary_model, x_label = expression("1"^ry*" consumers")) 

fruit_plot <- brms_tropics_side_plot(
  data_set = fruit_data, ylabel = NULL, plot_label = "b", plot_model = fruit_model, 
  stats_model = centered_fruit_model, x_label = "Frugivores") + 
  theme(axis.text.y = element_blank())

secondary_plot <- brms_tropics_side_plot(
  data_set = secondary_data, ylabel = NULL, plot_label = "c", plot_model = secondary_model, 
  stats_model = centered_secondary_model, x_label = expression("2"^ry*" consumers")) + 
  theme(axis.text.y = element_blank())

invert_plot <- brms_tropics_side_plot(
  data_set = invert_data, ylabel = NULL, plot_label = "d", plot_model = invert_model, 
  stats_model = centered_invert_model, x_label = "Invertivores") + 
  theme(axis.text.y = element_blank())


# Life history partitions.
mig_plot <- brms_tropics_side_plot(
  data_set = mig_data, plot_label = "e", plot_model = mig_model,
  stats_model = centered_mig_model, x_label = "Migrants") 

non_mig_plot <- brms_tropics_side_plot(
  data_set = non_mig_data, ylabel = NULL, plot_label = "f", plot_model = non_mig_model, 
  stats_model = centered_non_mig_model, x_label = "Non-migrants") + 
  theme(axis.text.y = element_blank())

terr_plot <- brms_tropics_side_plot(
  data_set = terr_data, ylabel = NULL, plot_label = "g", plot_model = terr_model, 
  stats_model = centered_terr_model, x_label = "Territorial") + 
  theme(axis.text.y = element_blank())

non_terr_plot <- brms_tropics_side_plot(
  data_set = non_terr_data, ylabel = NULL, plot_label = "h", plot_model = non_terr_model, 
  stats_model = centered_non_terr_model, x_label = "Non-territorial") + 
  theme(axis.text.y = element_blank())

# Arrange the plots together.
ggarrange(primary_plot, fruit_plot, secondary_plot, invert_plot, 
          mig_plot, non_mig_plot, terr_plot, non_terr_plot,
          nrow = 2, ncol = 4, widths = c(1.2,1,1,1))
ggsave("Plots/Tropics/tropical_comparisons.pdf", height = 10, width = 15, dpi = 600)
ggsave("Plots/Tropics/tropical_comparisons.png", height = 10, width = 15, dpi = 600)



