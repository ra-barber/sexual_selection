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

# # # Read in the new files.
# lat_files <- list.files(path = c("Z:/home/sexual_selection/Results/Models/Latitude/"),
#                         full.names = T, include.dirs = FALSE, recursive = FALSE)

# Actual models to report in main results section. Centered to aid direct comparison.
# all_frug_models <- combine_brms("frugivore_centered_all") 
# high_frug_models <- combine_brms("frugivore_centered_high") 
# 
# all_invert_models <- combine_brms("invertivore_centered_all")
# 
# primary_all_lat_models <- combine_brms("primary_centered_all")  
# primary_high_lat_models <- combine_brms("primary_centered_high") 
# 
# secondary_all_lat_models <- combine_brms("secondary_centered_all") 
# secondary_high_lat_models <- combine_brms("secondary_centered_high")

# 
# # Models for plotting. Uncentered to see raw relationship with latitude.
# all_frug_uncen_models <- combine_brms("frugivore_uncentered_all") 
# high_frug_uncen_models <- combine_brms("frugivore_uncentered_high")  
# 
# all_invert_uncen_models <- combine_brms("invertivore_uncentered_all") 
# high_invert_uncen_models <- combine_brms("invertivore_uncentered_high") 
# 
# all_primary_uncen_models <- combine_brms("primary_uncentered_all")
# high_primary_uncen_models <- combine_brms("primary_uncentered_high") 
# 
# all_secondary_uncen_models <- combine_brms("secondary_uncentered_all")
# high_secondary_uncen_models <- combine_brms("secondary_uncentered_high") 

# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Combined_models/Latitude/"

# Read in the centered models.
allbirds_uncen_all_data <- readRDS(paste0(first_half, "all_uncentered_all_data.rds"))

#cert_cen_all_data <- readRDS(paste0(first_half, "certainty_centered_all_data.rds"))

prim_cen_all_data <- readRDS(paste0(first_half, "primary_centered_all_data.rds"))
frug_cen_all_data <- readRDS(paste0(first_half, "frugivore_centered_all_data.rds"))
sec_cen_all_data <- readRDS(paste0(first_half, "secondary_centered_all_data.rds"))
invert_cen_all_data <- readRDS(paste0(first_half, "invertivore_centered_all_data.rds"))

# Uncentered models.
prim_uncen_all_data <- readRDS(paste0(first_half, "primary_uncentered_all_data.rds"))
frug_uncen_all_data <- readRDS(paste0(first_half, "frugivore_uncentered_all_data.rds"))
sec_uncen_all_data <- readRDS(paste0(first_half, "secondary_uncentered_all_data.rds"))
invert_uncen_all_data <- readRDS(paste0(first_half, "invertivore_uncentered_all_data.rds"))



################################################################################
                    #### Export summary tables ####


# Extract relevant coeffcient information.
allbirds_all_estimates <- summary(allbirds_cen_all_data[[1]])$fixed[5,c(1,3,4)]
#cert_all_estimates <- summary(centered_cert_model)$fixed[4,c(1,3,4)]

primary_all_estimates <- summary(prim_cen_all_data[[1]])$fixed[5,c(1,3,4)]
fruit_all_estimates <- summary(frug_cen_all_data[[1]])$fixed[5,c(1,3,4)]
secondary_all_estimates <- summary(sec_cen_all_data[[1]])$fixed[5,c(1,3,4)]
invert_all_estimates <- summary(invert_cen_all_data[[1]])$fixed[5,c(1,3,4)]

all_estimates <- rbind(allbirds_all_estimates, #cert_all_estimates, 
                       primary_all_estimates, fruit_all_estimates,
                       secondary_all_estimates, invert_all_estimates)
row.names(all_estimates) <- c("all_birds", #"certainty", 
                              "primary",
                              "fruit", "secondary", "invert")
# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0("[", round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2), "]"),
  est_intervals = paste0(round_est, " ", intervals))

# Export the results.
write.csv(all_estimates, "Results/Tables/all_phy_lat_regression.csv", row.names = TRUE)





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
                               plot_model_data = all_frug_uncen_data, sex_score = TRUE,
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



# Two last models.
pri_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Primary") %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_x_pos = 10, lab_ypos = 1.5, 
                     plot_label = "b", plot_model_data = all_prim_uncen_data) 

sec_lat_plot <- diet_lat_data %>% filter(trophic_binary == "Secondary") %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
                     ybreaks =  c(0,1.0,2.0), lab_x_pos = 10, lab_ypos = 1.5,
                     plot_label = "f", plot_model_data = all_sec_uncen_data) 



fruit_lat_data <- niche_lat_data %>% filter(trophic_niche == "Frugivore")
fruit_lat_data[9,] <- list(40, "Frugivore", 0, 0, 0, 0, 0, 1)  # Add in three monogamous frugivore species at higher lat.
fruit_lat_data[10,] <- list(45, "Frugivore", 0, 0, 0, 0, 0, 1)
fruit_lat_data[11,] <- list(50, "Frugivore", 0, 0, 0, 0, 0, 1)

# Trophic niche models.
fruit_lat_plot <- fruit_lat_data %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_x_pos = 10, lab_ypos = 1.5, 
                     plot_label = "d", plot_model_data = all_frug_uncen_data) 

invert_lat_plot <- niche_lat_data %>% filter(trophic_niche == "Invertivore") %>% 
  brms_lat_side_plot(ylabel = "Sexual selection", ylimits = c(0,2), 
                     ybreaks =  c(0,1.0, 2.0), lab_x_pos = 10, lab_ypos = 1.5, 
                     plot_label = "h", plot_model_data = all_invert_uncen_data) 



################################################################################
                       ##### Plot all together #######

ylabel = "Sexual selection"
ylimits = c(0,2)
ybreaks =  c(0,1.0, 2.0)
lab_x_pos = 10
lab_ypos = 1.5
plot_label = "h"


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
  
  geom_line(data = all_prim_uncen_predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1, colour = "#05299E") +
  geom_line(data = all_frug_uncen_predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1, colour = "#214F4B") +
  geom_line(data = all_sec_uncen_predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1, colour = "#8D0801") +
  geom_line(data = all_invert_uncen_predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1, colour = "#A88A05")





# Change sexual score back to original values (have to + 1 for ordinal regression)
all_prim_uncen_predictions <- all_prim_uncen_data[[2]]
all_frug_uncen_predictions <- all_frug_uncen_data[[2]]
all_sec_uncen_predictions <- all_sec_uncen_data[[2]]
all_invert_uncen_predictions <- all_invert_uncen_data[[2]]

all_prim_uncen_predictions$estimate__ <- all_prim_uncen_predictions$estimate__ - 1
all_frug_uncen_predictions$estimate__ <- all_frug_uncen_predictions$estimate__ - 1
all_sec_uncen_predictions$estimate__ <- all_sec_uncen_predictions$estimate__ - 1
all_invert_uncen_predictions$estimate__ <- all_invert_uncen_predictions$estimate__ - 1
   
################################################################################
              #### Try family averages ######


# Read in clade function and assign.
source("Code/clade_function.R")
model_data %<>% assign_clades()

model_data$abs_lat <- abs(model_data$complete_latitude)

# Create settings for scatter plots.
point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                      size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                       labs(x = "", y = NULL,  colour = "Clade"),
                       scale_colour_manual(values = clade_colours),
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
  geom_ribbon(data = all_sec_uncen_predictions, inherit.aes = FALSE,
              aes(x = abs_lat, ymin = lower__-1, ymax = upper__-1), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = all_sec_uncen_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1)






                         
###############################################################################
                  #### Plot conditonal effects #####


tic()
high_frug_conds <- conditional_effects(high_frug_uncen_models)
toc()

tic()
high_invert_conds <- conditional_effects(high_invert_uncen_models)
toc()




###############################################################################
                      #### Lambda values #####


# Look at lambda values.
hyp <- "sd_tree_tip__Intercept^2 / (sd_tree_tip__Intercept^2 + disc^2) = 0"

(all_lambda <- hypothesis(single_model, hyp, class = NULL))
plot(all_lambda)



###############################################################################
            #### Extract the draws for each model #####



# Simple models first.
simple_all_plotdata <- extract_draws(equi_all_models, simple_names)
terr_simple_all_plotdata <- extract_draws(terr_equi_all_models, terr_simple_names)

# Npp and no npp complex models.
npp_all_plotdata <- extract_draws(npp_cen_all_models, npp_names)
nonpp_all_plotdata <- extract_draws(nonpp_cen_all_models, nonpp_names)

# High certainty data.
npp_high_plotdata <- extract_draws(npp_cen_high_models, npp_names)
terr_simple_high_plotdata <- extract_draws(terr_equi_high_models, terr_simple_names)

terr_simple_all_plotdata %<>% filter(name %in% c("Territoriality", "Migration", "Primary\nconsumer",
                                                  "Seasonality", "Primary\nconsumer\n x \nSeasonality", "Primary\nconsumer\n x \nTerritoriality"))
terr_simple_high_plotdata %<>% filter(name %in% c("Territoriality", "Migration", "Primary\nconsumer",
                       "Seasonality", "Primary\nconsumer\n x \nSeasonality", "Primary\nconsumer\n x \nTerritoriality"))

# Presentation figure.
all_plotdata <- extract_draws(no_inter_cen_all_models, model_colnames)
all_plotdata %<>% filter(name %in%  c("Territoriality", "Migration", "(PC) Primary\nconsumer", "Chick maturity",
                                      "Seasonality", "PC x \nSeasonality"))



###############################################################################
                #### Add axis label and legend information #####


# Vector of predictors in reverse order of plotting.
simple_preds <- simple_names
simple_terr_preds <- terr_simple_names[-4]

simple_terr_preds <- simple_terr_preds[c(4,2,3,1,5,6)]

npp_preds <- npp_names[c(1:7, 9, 8, 10, 12, 11)]
nonpp_preds <- nonpp_names

# Set data as factors for plotting order.
#all_plotdata$name %<>% factor(levels = rev(rev_preds))

# Simple models first.
simple_all_plotdata$name %<>% factor(levels = rev(simple_preds))
terr_simple_all_plotdata$name %<>% factor(levels = rev(simple_terr_preds))

# Npp and no npp complex models.
npp_all_plotdata$name  %<>% factor(levels = rev(npp_preds))
nonpp_all_plotdata$name  %<>% factor(levels = rev(nonpp_preds))

# High certainty data.
npp_high_plotdata$name  %<>% factor(levels = rev(npp_preds))
terr_simple_high_plotdata$name  %<>% factor(levels = rev(simple_terr_preds))



lifehistory <- c("Territoriality",  "Migration",  "Primary\nconsumer", "(Pr) Precocial", "Chick maturity")
resource <-   c("Seasonality", "PC x \nSeasonality")
resource <-   c("Seasonality", "Productivity")
interaction <-  c("Primary\nconsumer\n x \nTerritoriality", "Primary\nconsumer\n x \nSeasonality", "Primary\nconsumer\n x \nProductivity")
#devo_interaction <- c("Pr x \nTerritoriality", "Pr x \nSeasonality", "Pr x \nProductivity")
devo_interaction <- c("DM x PC", "DM x \ngeneration length", "DM x Migration", "DM x \nseasonality",
                      "DM x \nproductivity", "DM x \nterritoriality","CM x \nSeasonality",
                      "CM x \nProductivity", "CM x \nTerritoriality")

# Function for adding colour to legends.
add_legend_info <- function(plotdata){
  plotdata$trait <- NA
  plotdata$trait[plotdata$name %in% lifehistory] <- "Life history"
  plotdata$trait[plotdata$name %in% resource] <- "Environmental"
  plotdata$trait[plotdata$name %in% interaction] <- "Diet interaction"
  plotdata$trait[plotdata$name %in% devo_interaction] <- "Devlopment interaction"
  plotdata
}

# Add the legend info.
simple_all_plotdata %<>% add_legend_info()
terr_simple_all_plotdata %<>% add_legend_info()
npp_all_plotdata %<>% add_legend_info()
nonpp_all_plotdata %<>% add_legend_info()

# High certainty.
npp_high_plotdata %<>% add_legend_info()
terr_simple_high_plotdata %<>% add_legend_info()
# Y axis labels for ordering plot.
#axes_labs <- rev_preds


# Group plot settings.
ggplot_settings <- list(
  # geom_point(alpha = .1, stroke = 0.2, colour = "#586994", size = 0.5, 
  #            position = position_jitter(seed = 1993, height = .35)),
  # stat_pointinterval(size = 9, alpha = rep(c(0.95, 0.6), times = 9), 
  #                    colour = "navy"),
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6), # Add lines.
  #scale_y_discrete(labels = rev(axes_labs)),
  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL),
  theme_classic2(base_size = 12),   # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm')))



###############################################################################
                         #### Colour pal #####

pred_n <- 6

library(ggnewscale)

light_colours <- c( "#77AD78","#7494EA", "#AA4465")
light_colours <- c( "#77AD78","#7494EA", "#C98986")



dark_colours <- c("#136F63", "#05299E", "#8D0801")
dark_colours <- c("#214F4B", "#05299E", "#8D0801")
legend_order <- c("Life history", "Environmental", "Trophic interaction")

# cOLOURS WITH PRECOCIAL INTERACTION.
light_colours <- c( "#77AD78","#7494EA", "#C98986", "#D8C77B")   #F7D460   #F6CE5F
dark_colours <- c("#214F4B", "#05299E", "#8D0801","#A88A05")
light_colours <- c("#7494EA",  "#77AD78", "#C98986", "#D8C77B")   #F7D460   #F6CE5F
dark_colours <- c( "#05299E","#214F4B", "#8D0801","#A88A05")

legend_order <- c( "Environmental", "Life history","Diet interaction", "Devlopment interaction")

###############################################################################
                    #### Quasi random plots #####

gg_quasi_forst <- function(plotdata, pred_n = 6, legend_pos = c(0.15, 0.15), col_pal = light_colours, alpha_val = 0.1){
  # Main publication design for simple models.
  ggplot(plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
    geom_quasirandom(bandwidth = 0.5, size = 0.5,
                     stroke = 0.1, groupOnX = FALSE, alpha = alpha_val, position = position_dodge(2)) +
    scale_colour_manual(values = col_pal, breaks = legend_order) +
    guides(colour = guide_legend(byrow = TRUE, override.aes = list(alpha = 1, shape = 15, size = 4))) +
    new_scale_colour() +
    stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
    guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
    geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
    #scale_y_discrete(labels = rev(axes_labs)) +
    # facet_grid(rows = vars(trait), scales = "free_y") +
    labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
    theme_classic2(base_size = 12) +  # Theme.
    theme(legend.position = legend_pos,
          axis.text.x=element_text(size=rel(1), face = "bold", colour = "black"), 
          axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
          axis.title.x=element_text(size=rel(0.7), face = "bold"),
          legend.text = element_text(size=rel(0.8), face = "bold"), 
          #legend.key.size = unit(0, 'cm'),
          legend.key.height = unit(0.2, 'cm'),
          # legend.key.width = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title = element_blank(),
          strip.text.y = element_blank())
}

gg_quasi_forst(simple_all_plotdata, legend_pos = c(0.85, 0.9), alpha_val = 0.2)
ggsave("Plots/Results/simple_model.tiff", dpi = 600, width = 8, height = 4)

gg_quasi_forst(terr_simple_all_plotdata, pred_n = 6, legend_pos = c(0.125, 0.125), alpha_val = 0.5)
ggsave("Plots/Results/terr_simple_model.tiff", dpi = 600, width = 8, height = 4)  # Changed for presentation size.
ggsave("Plots/Results/figure_5_model.tiff", dpi = 1000, width = 8, height = 5)  # Changed for presentation size.

#gg_quasi_forst(npp_all_plotdata, pred_n = 12, legend_pos = c(0.15, 0.1))
#ggsave("Plots/Results/npp_model.tiff", dpi = 600, width = 8, height = 6)

#gg_quasi_forst(nonpp_all_plotdata, pred_n = 9, legend_pos = c(0.15, 0.125))
#ggsave("Plots/Results/nonpp_model.tiff", dpi = 600, width = 8, height = 5)


# High certainty data.
#gg_quasi_forst(npp_high_plotdata, pred_n = 12, legend_pos = c(0.125, 0.125))
#ggsave("Plots/Results/npp_high_model.tiff", dpi = 600, width = 8, height = 6)

gg_quasi_forst(terr_simple_high_plotdata, pred_n = 6, legend_pos = c(0.125, 0.125), alpha_val = 0.2)
ggsave("Plots/Results/terr_simple_high_model.tiff", dpi = 1000, width = 8, height = 4)


ggsave("Plots/Results/figure_S2_tall.tiff", dpi = 1000, width = 8, height = 8)


# Main publication design for simple models.
ggplot(simple_all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
  scale_colour_manual(values = light_colours[c(2,3, 4)], breaks = legend_order) +
  guides(colour = 
           guide_legend(byrow = TRUE,
                        override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  #scale_y_discrete(labels = rev(axes_labs)) +
  # facet_grid(rows = vars(trait), scales = "free_y") +
  
  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.175, 0.2),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        #legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        # legend.key.width = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.text.y = element_blank())

ggsave("Plots/Results/simple_model.tiff", dpi = 600, width = 8, height = 4)

# Main publication design, but without facet (facet only works if groups are equal size.)
ggplot(simple_all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
  scale_colour_manual(values = light_colours[c(2,3)], breaks = legend_order) +
  guides(colour = 
           guide_legend(byrow = TRUE,
                        override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  # facet_grid(rows = vars(trait), scales = "free_y") +
  
  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.175, 0.2),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        #legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        # legend.key.width = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.text.y = element_blank() 
  )


# Main publication design, but without facet and moving legend when altricial is included.
ggplot(everything_data, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
  scale_colour_manual(values = light_colours, breaks = legend_order) +
  guides(colour = 
           guide_legend(byrow = TRUE,
                        override.aes = list(alpha = 1, shape = 15, size = 4))) +
  new_scale_colour() +
  stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) + #xlim(c(-6,6)) +
  # facet_grid(rows = vars(trait), scales = "free_y") +
  
  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.15, 0.15),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        #legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        # legend.key.width = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.text.y = element_blank() 
  )

ggsave("Plots/Results/everything_model_4.tiff", dpi = 600, width = 8, height = 7)


# With facet grid, but it makes each facet equal so it doesn't work properly.
ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1, position = position_dodge(2)) +
    scale_colour_manual(values = light_colours, breaks = legend_order) +
    guides(colour = 
             guide_legend(byrow = TRUE,
               override.aes = list(alpha = 1, shape = 15, size = 4))) +
    new_scale_colour() +
    stat_pointinterval(size = rep(c(9, 2), times = pred_n)) + 
    guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  #scale_y_discrete(labels = rev(axes_labs), expand = c(5,5)) +
  facet_grid(rows = vars(trait), scales = "free_y") +
  #scale_y_break(c("Diet","Territory")) +
  #scale_x_break(c(2)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.175, 0.9),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        #legend.key.size = unit(0, 'cm'),
         legend.key.height = unit(0.2, 'cm'),
        # legend.key.width = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.text.y = element_blank() 
  )
ggsave("Plots/Results/all_model_5.tiff", dpi = 600, width = 8, height = 4)

unique(all_plotdata$name)




library(ggbreak)


ggplot(all_plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5,
                   stroke = 0.1,# shape=21, 
                   groupOnX = FALSE, alpha = 0.1) +
  scale_colour_manual(values = light_colours, guide="none", breaks = legend_order) +
  new_scale_colour() +
  stat_pointinterval(colour = "black", size = 5, 
                     alpha = rep(c(0.95, 0.8), times = 10)) +
  scale_colour_manual(values = dark_colours, breaks = legend_order) +
  scale_alpha_continuous(range = c(0.8,1)) +
  guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.25, 0.85),
        axis.text.x=element_text(size=rel(0.8), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.8), face = "bold"), 
        legend.key.height = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'))



ggplot(all_plotdata, aes(x = value, y = name)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5, colour = "#586994",
                   stroke = 0.1, shape=21, groupOnX = FALSE, alpha = 0.1) +
  stat_pointinterval(size = 5, 
                     alpha = rep(c(0.95, 0.5), times = 9),
                     colour = "navy") +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm'))


ggplot(all_plotdata, aes(x = value, y = name)) +
  geom_beeswarm(size = 0.5, colour = "#586994", dodge.width = 1, priority = "density",
                stroke = 0.1, shape=21, groupOnX = FALSE, alpha = 1) +
  stat_pointinterval(size = 5, 
                     alpha = rep(c(0.95, 0.5), times = 9),
                     colour = "navy") +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm'))


ggplot(all_plotdata, aes(x = value, y = name)) +
  geom_quasirandom(bandwidth = 0.5, size = 0.5, colour = "#586994",
                   stroke = 0.1, shape=21, groupOnX = FALSE, alpha = 1) +
  stat_pointinterval(size = 5, 
                     alpha = rep(c(0.95, 0.5), times = 9),
                     colour = "navy") +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + # Add lines.
  scale_y_discrete(labels = rev(axes_labs)) +
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL) +
  theme_classic2(base_size = 12) +  # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.8), face = "bold"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold"),
        axis.title.x=element_text(size=rel(0.7), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm'))






##############################################################################
                      #### Leftover code #####


ggplot_settings <- list(
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6), # Add lines.
  #coord_cartesian(xlim = c(-0.25, 0.75)),   # Axes.
  scale_y_discrete(labels = rev(axis_labels)),
  # scale_color_manual(values = rev(plot_palette),  # Colours.
  #                    guide = guide_legend(reverse = TRUE, 
  #                                         override.aes = list(linetype = 0))),
  # scale_fill_manual(values = rev(plot_palette), 
  #                    guide = "none"),
  labs(x = "Standardised Effect Size", y = NULL, fill = NULL, colour = NULL),
  theme_classic2(base_size = 12),   # Theme.
  theme(legend.position = c(0.92, 0.32),
        axis.text.x=element_text(size=rel(0.7), face = "bold"), 
        axis.text.y=element_text(size=rel(0.7), face = "bold"),
        axis.title.x=element_text(size=rel(0.6), face = "bold"),
        legend.text = element_text(size=rel(0.6), face = "bold"), 
        legend.key.height = unit(0.6, 'cm')))



benchmarking <- read.csv("Z:/home/sexual_selection/benchmarking.csv")
benchmarking %<>% filter(elapsed > 0)
plot(elapsed ~ grain_size, data = benchmarking)
plot(elapsed ~ jitter(thread_number), data = benchmarking)
boxplot(log(elapsed) ~ family, data = benchmarking)
boxplot(log(elapsed) ~ thread_number, data = benchmarking)
boxplot(log(elapsed) ~ grain_size, data = benchmarking)

boxplot(elapsed ~ thread_number, data = benchmarking)

hist(benchmarking$elapsed)

linear_test <- lm(elapsed ~ grain_size*thread_number, data = benchmarking)

summary(linear_test)

threads_16 <- benchmarking %>% filter(thread_number == 16)

plot(elapsed ~ grain_size, data = threads_16)

library(interactions)

interact_plot(linear_test, pred = thread_number, modx = grain_size)


# Most basic models.
model_1_all <- readRDS("Z:/home/sexual_selection/Results/Models/Test/sratio_16_100_1.rds")
cum_model <- readRDS("Z:/home/sexual_selection/Results/Models/Test/cumulative_16_100_1.rds")

plot(model_1_all)

mcmc_areas(cum_model, regex_pars = "^b_[a-z]", prob = 0.95)
