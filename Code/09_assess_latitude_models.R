###############################################################################
              # Plotting Latitudinal brms models from HPC  #
###############################################################################

# This script extracts model estimates and creates figure S7. 
# It uses models that have been run on the HPC.


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
library(bayestestR)

# Clear the workspace.
rm(list=ls())

# Set the seed.
set.seed(1993)

###############################################################################
                             #### Data #####

# Functions.
source("Code/functions.R")

# Read in the life history traits.
model_data <- read_ss_data()


###############################################################################
              #### Prepare predictor variables ######

# Scale continuous predictors to two SD.
model_data %<>% mutate(centroid_z = standardize(sqrt(abs_lat), two_sd = TRUE))

# Create 5 degree bins for latitude.
bin_range <- seq(0, 80, 5)
bin_labels <- seq(0, 75, 5) 
model_data$binned_lat <- cut(model_data$abs_lat, breaks=bin_range, labels = bin_labels, include.lowest = TRUE)
model_data$binned_lat %<>% as.character() %>% as.numeric()

# Create grouped data to reuse with different varaibles below.
grouped_lat_data <- model_data  %>% group_by(binned_lat)

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
diet_lat_data <- model_data  %>% 
  group_by(binned_lat, trophic_level_binary) %>% average_lat_bins()
niche_lat_data <- model_data  %>% 
  group_by(binned_lat, trophic_niche) %>% average_lat_bins()

# Add in the three monogamous frugivore species that each occupy a single lat bin alone.
fruit_lat_data <- niche_lat_data %>% filter(trophic_niche == "Frugivore")

# Do migration and territoriality for plots.
mig_lat_data <- model_data  %>% 
  group_by(binned_lat, migration_binary) %>% average_lat_bins()
terr_lat_data <- model_data  %>% 
  group_by(binned_lat, territoriality_binary) %>% average_lat_bins()

# Do same for high and medium certainty.
hi_lat_data <- model_data %>% filter(data_certainty == 4) %>% 
  group_by(binned_lat) %>% average_lat_bins()
med_lat_data <- model_data %>% filter(data_certainty > 2) %>% 
  group_by(binned_lat) %>% average_lat_bins()

# Group the data by trophic level territory and latitude.
terr_diet_lat_data <- model_data  %>% 
  group_by(binned_lat, trophic_level_binary) %>% average_lat_bins("terr_dummy")
year_terr_diet_lat_data <- model_data  %>% 
  group_by(binned_lat, trophic_level_binary) %>% average_lat_bins("year_terr_dummy")


###############################################################################
                    #### Read in brms models ####


# Directory where models are saved.
first_half <- "Z:/home/sexual_selection/Results/Models/Latitudinal/Continuous/"
list.files(first_half)

# Function to read in models.
read_lat_model <- function(model_name){
  model_file <- paste0(first_half, model_name, "_model.rds")
  readRDS(model_file)
}

allbirds_all_model <- read_lat_model("allbirds_all")
certainty_all_model <- read_lat_model("certainty_all")
primary_all_model <- read_lat_model("primary_all")
fruit_all_model <- read_lat_model("fruit_all")
secondary_all_model <- read_lat_model("secondary_all")
invert_all_model <- read_lat_model("invert_all")

allbirds_high_model <- read_lat_model("allbirds_high")
certainty_high_model <- read_lat_model("certainty_high")
primary_high_model <- read_lat_model("primary_high")
fruit_high_model <- read_lat_model("fruit_high")
secondary_high_model <- read_lat_model("secondary_high")
invert_high_model <- read_lat_model("invert_high")

# Read in mig and terr models.
mig_all_model <- read_lat_model("mig_all")
no_mig_all_model <- read_lat_model("no_mig_all")
terr_all_model <- read_lat_model("terr_all")
no_terr_all_model <- read_lat_model("no_terr_all")

mig_high_model <- read_lat_model("mig_high")
no_mig_high_model <- read_lat_model("no_mig_high")
terr_high_model <- read_lat_model("terr_high")
no_terr_high_model <- read_lat_model("no_terr_high")

# Read in extra models.
hi_cert_model <- read_lat_model("highcert_all")
prim_allterr_all_model <- read_lat_model("prim_allterr_all")
prim_yearterr_all_model <- read_lat_model("prim_yearterr_all")
sec_allterr_all_model <- read_lat_model("sec_allterr_all")
sec_yearterr_all_model <- read_lat_model("sec_yearterr_all")


################################################################################
                    #### Export summary tables ####


allbirds_all_estimates <- summary(allbirds_all_model)$fixed[5,c(1,3,4)]
cert_all_estimates <- summary(certainty_all_model)$fixed[4,c(1,3,4)]
primary_all_estimates <- summary(primary_all_model)$fixed[5,c(1,3,4)]
fruit_all_estimates <- summary(fruit_all_model)$fixed[5,c(1,3,4)]
secondary_all_estimates <- summary(secondary_all_model)$fixed[5,c(1,3,4)]
invert_all_estimates <- summary(invert_all_model)$fixed[5,c(1,3,4)]

mig_all_estimates <- summary(mig_all_model)$fixed[5,c(1,3,4)]
no_mig_all_estimates <- summary(no_mig_all_model)$fixed[5,c(1,3,4)]
terr_all_estimates <- summary(terr_all_model)$fixed[5,c(1,3,4)]
no_terr_all_estimates <- summary(no_terr_all_model)$fixed[5,c(1,3,4)]

all_estimates <- rbind(allbirds_all_estimates, cert_all_estimates, 
                       primary_all_estimates, fruit_all_estimates,
                       secondary_all_estimates, invert_all_estimates,
                       mig_all_estimates, no_mig_all_estimates,
                       terr_all_estimates, no_terr_all_estimates)
row.names(all_estimates) <- c("all_birds", "certainty", "primary",
                              "fruit", "secondary", "invert", "migration", 
                              "no_migration", "territorial", "non_territorial")
# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export the results.
write.csv(all_estimates, "Results/Tables/all_continuous_lat_regression.csv", row.names = TRUE)


# High certainty estimates
allbirds_high_estimates <- summary(allbirds_high_model)$fixed[5,c(1,3,4)]
cert_high_estimates <- summary(certainty_high_model)$fixed[4,c(1,3,4)]
primary_high_estimates <- summary(primary_high_model)$fixed[5,c(1,3,4)]
fruit_high_estimates <- summary(fruit_high_model)$fixed[5,c(1,3,4)]
secondary_high_estimates <- summary(secondary_high_model)$fixed[5,c(1,3,4)]
invert_high_estimates <- summary(invert_high_model)$fixed[5,c(1,3,4)]

mig_high_estimates <- summary(mig_high_model)$fixed[5,c(1,3,4)]
no_mig_high_estimates <- summary(no_mig_high_model)$fixed[5,c(1,3,4)]
terr_high_estimates <- summary(terr_high_model)$fixed[5,c(1,3,4)]
no_terr_high_estimates <- summary(no_terr_high_model)$fixed[5,c(1,3,4)]

high_estimates <- rbind(allbirds_high_estimates, cert_high_estimates, 
                       primary_high_estimates, fruit_high_estimates,
                       secondary_high_estimates, invert_high_estimates,
                       mig_high_estimates, no_mig_high_estimates,
                       terr_high_estimates, no_terr_high_estimates)
row.names(high_estimates) <- c("all_birds", "certainty", "primary",
                              "fruit", "secondary", "invert", "migration", 
                              "no_migration", "territorial", "non_territorial")
# Paste together values for reporting in a table.
high_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))


# Export the results.
write.csv(high_estimates, "Results/Tables/high_continuous_lat_regression.csv", row.names = TRUE)


###############################################################################
                    #### side plot function  ######

# Change scipen. (Probs not needed now)
options(scipen = 999)

# Function to create panel plots in figure S7
brms_lat_side_plot <- function(data_set, xlabel,  ylabel = "", ylimits = c(-0.1,2.1), ybreaks = c(0,1,2), 
                               lab_x_pos = 40, lab_ypos = 2.1, plot_label = "b", 
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
    stats_label <- paste0(estimate, "  ", p_value)
  } else {
    stats_label <- paste0(estimate, "\n", intervals)
  }

  # ggplot function for sideplots with annotations.
  ggplot(data_set, aes(x = binned_lat, y = trait_mean)) +
  geom_errorbar(aes(ymin = trait_min, ymax = trait_max),
                position = position_dodge(width = 1), show.legend = FALSE, col = "darkgrey") + 
    geom_point(aes(size = trait_n), position = position_dodge(width = 1), alpha = 0.9, col = "#442B48") + 
      scale_x_continuous(breaks = seq(from = 0, to = 70, by = 35)) +
    scale_y_continuous(breaks = ybreaks, labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = ylimits, xlim = c(NA, 75), clip = 'off') +
    scale_size_continuous(range = c(1,5))+
    ylab(ylabel) +
    xlab(xlabel) + theme_classic(base_size = 25) + 
    theme(legend.position = "none", 
          axis.title.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          plot.margin = margin(t = 1, l = 0.2, b = 0.2, r =0.3, unit = "cm")) + 
    annotate("text", x = lab_x_pos, y =lab_ypos, label = stats_label, size = 7) + 
    annotate("text", x = 0, y = ylimits[2], label = plot_label, size = 12, fontface = 2) + 
    geom_line(data = predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1)
}


################################################################################
                     #### Make thhe panels  ####


# Primary and secondary consumers.
pri_lat_plot <- diet_lat_data %>% 
  filter(trophic_level_binary == "Primary") %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", 
                     plot_label = "a", plot_model = primary_all_model,
                     xlabel = expression("1"^ry*" consumers")) 

fruit_lat_plot <- fruit_lat_data %>% 
  brms_lat_side_plot(plot_label = "b", plot_model = fruit_all_model,
                     ylabel = NULL,  xlabel = "Frugivores") + 
  theme(axis.text.y = element_blank())

sec_lat_plot <- diet_lat_data %>% 
  filter(trophic_level_binary == "Secondary") %>% 
  brms_lat_side_plot(plot_label = "c", plot_model = secondary_all_model,
                     ylabel = NULL, xlabel = expression("2"^ry*" consumers")) + 
  theme(axis.text.y = element_blank())


invert_lat_plot <- niche_lat_data %>% 
  filter(trophic_niche == "Invertivore") %>% 
  brms_lat_side_plot(plot_label = "d", plot_model = invert_all_model,
                     ylabel = NULL, xlabel = "Invertivores") + 
  theme(axis.text.y = element_blank())

mig_lat_plot <- mig_lat_data %>% filter(migration_binary == "Strong") %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", plot_label = "e",
                     plot_model = mig_all_model, xlabel = "Migrants")

no_mig_lat_plot <- mig_lat_data %>% filter(migration_binary == "Weak") %>% 
  brms_lat_side_plot(  plot_label = "f", 
                     plot_model = no_mig_all_model,
                     ylabel = NULL, xlabel = "Non-migrants") + 
  theme(axis.text.y = element_blank())

terr_lat_plot <- terr_lat_data %>% 
  filter(territoriality_binary == "Territorial") %>% 
  brms_lat_side_plot(plot_label = "g",
                     plot_model = terr_all_model,
                     ylabel = NULL, xlabel = "Territorial") + 
  theme(axis.text.y = element_blank())

no_terr_lat_plot <- terr_lat_data %>% 
  filter(territoriality_binary == "Non-territorial") %>% 
  brms_lat_side_plot( plot_label = "h",
                     plot_model = no_terr_all_model,
                     ylabel = NULL, xlabel = "Non-territorial") + 
  theme(axis.text.y = element_blank())


ggarrange(pri_lat_plot, fruit_lat_plot, sec_lat_plot, invert_lat_plot, 
          mig_lat_plot, no_mig_lat_plot, terr_lat_plot, no_terr_lat_plot,
          nrow = 2, ncol = 4, widths = c(1.2,1,1,1))


ggsave("Figures/figure_S7.png", height = 10, width = 15, dpi = 600)
ggsave("Figures/figure_S7.pdf", height = 10, width = 15, dpi = 600)


################################################################################
                            #### End ####
################################################################################