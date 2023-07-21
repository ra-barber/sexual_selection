###############################################################################
                # Assess the fit/plot univariate models  #
###############################################################################

# This code won't run properly yet because i change the object name from "deets" to "data"



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
          #### Read in the model data for plotting #####


# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/"
list.files(first_half, pattern = "data")

# Read in the models using all data certainties.
temp_cen_all_data <- readRDS(paste0(first_half, "temp_centered_all_data.rds"))
temp_uncen_all_data <- readRDS(paste0(first_half, "temp_uncentered_all_data.rds"))
mig_cen_all_data <- readRDS(paste0(first_half, "mig_centered_all_data.rds"))
mig_uncen_all_data <- readRDS(paste0(first_half, "mig_uncentered_all_data.rds"))
tro_cen_all_data <- readRDS(paste0(first_half, "tro_centered_all_data.rds"))
tro_uncen_all_data <- readRDS(paste0(first_half, "tro_uncentered_all_data.rds"))
terr_cen_all_data <- readRDS(paste0(first_half, "terr_centered_all_data.rds"))
terr_uncen_all_data <- readRDS(paste0(first_half, "terr_uncentered_all_data.rds"))


temp_cen_all_data[[3]]
mig_cen_all_data[[3]]
tro_cen_all_data[[3]]
terr_cen_all_data[[3]]



################################################################################
                    #### Export summary tables ####


# Extract relevant coeffcient information.
temp_all_estimates <- summary(temp_cen_all_data[[1]])$fixed[5,c(1,3,4)]
mig_all_estimates <- summary(mig_cen_all_data[[1]])$fixed[5,c(1,3,4)]
tro_all_estimates <- summary(tro_cen_all_data[[1]])$fixed[5,c(1,3,4)]
terr_all_estimates <- summary(terr_cen_all_data[[1]])$fixed[5,c(1,3,4)]
all_estimates <- rbind(temp_all_estimates, mig_all_estimates, 
                       tro_all_estimates, terr_all_estimates)

# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0("[", round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2), "]"),
  est_intervals = paste0(round_est, " ", intervals))

# Export the results.
#write.csv(all_estimates, "Results/Tables/all_univariate_regression.csv", row.names = TRUE)



###############################################################################
                #### Prepare family level data #####

# Read in model data.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- model_data$birdtree_name %>% str_replace(" ", "_")
row.names(model_data) <- model_data$tree_tip

# Read in clade function and assign.
source("Code/clade_function.R")
model_data %<>% assign_prum_clades()

# Average traits for family.
family_average_data_data <- model_data %>% group_by(family) %>% 
  summarise(higher_clade = first(higher_clade),
            tree_tip = first(tree_tip),
            sexual_mean = mean(sexual_score),
            clade_sum = length(sexual_score),
            temp_log = mean(temp_log),
            trophic_binary = mean(as.numeric(as.factor(trophic_binary)))-1,
            migration_binary  = mean(as.numeric(as.factor(migration_binary)))-1,
            territory_binary  = mean(as.numeric(as.factor(territory_binary)))-1)

# Group by family and trait for binary traits.
average_family <- function(varible){
  var_data <- model_data %>% dplyr::group_by(family, !!! syms(varible)) %>% 
    summarise(higher_clade = first(higher_clade),
              sexual_mean = mean(sexual_score),
              clade_size = length(sexual_score))
  var_data %>% group_by(family) %>% summarise(clade_sum = sum(clade_size)) %>% right_join(var_data)
}

# Get territory and migration means.
teritory_data <- average_family("territory_binary")
migration_data <- average_family("migration_binary")
trophic_data <- average_family("trophic_binary")

###############################################################################
                     #### Prepare plots #####


# Reorder colours for plotting.
prum_clade_colours <- c("#7494EA", # Palaeognathae
                        "#8D0801", # Galloanserae
                        "#77AD78", # Strisores
                        "#5941A9", # Columbaves
                        "#988F2A", # Gruiformes
                        "#136F63", # Aequorlitornithes
                        "#AA4465", # Opisthocomiformes
                        "#C8B8DB", # Accipitriformes
                        "#05299E", # Strigiformes
                        "#EB6424", # Coraciimorphae
                        "#E0CA3C") # Australaves

names(prum_clade_colours) <- c("Palaeognathae", "Galloanserae", "Strisores",
                               "Columbaves", "Gruiformes", "Aequorlitornithes", 
                               "Opisthocomiformes", "Accipitriformes", "Strigiformes",
                               "Coraciimorphae", "Australaves")







# Change predictions back to same scale as original sex scores.
temp_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
mig_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
tro_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
terr_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)

# Estimate
temp_estimate <- as.character(format(round(summary(temp_cen_all_data[[1]])$fixed[5,1], 2), nsmall = 2))
mig_estimate <- as.character(format(round(summary(mig_cen_all_data[[1]])$fixed[5,1], 2), nsmall = 2))
tro_estimate <- as.character(format(round(summary(tro_cen_all_data[[1]])$fixed[5,1], 2), nsmall = 2))
terr_estimate <- as.character(format(round(summary(terr_cen_all_data[[1]])$fixed[5,1], 2), nsmall = 2))

# Lower CI
temp_lower <- as.character(format(round(summary(temp_cen_all_data[[1]])$fixed[5,3], 2), nsmall = 2))
mig_lower <- as.character(format(round(summary(mig_cen_all_data[[1]])$fixed[5,3], 2), nsmall = 2))
tro_lower <- as.character(format(round(summary(tro_cen_all_data[[1]])$fixed[5,3], 2), nsmall = 2))
terr_lower <- as.character(format(round(summary(terr_cen_all_data[[1]])$fixed[5,3], 2), nsmall = 2))

# Upper CI
temp_upper <- as.character(format(round(summary(temp_cen_all_data[[1]])$fixed[5,4], 2), nsmall = 2))
mig_upper <- as.character(format(round(summary(mig_cen_all_data[[1]])$fixed[5,4], 2), nsmall = 2))
tro_upper <- as.character(format(round(summary(tro_cen_all_data[[1]])$fixed[5,4], 2), nsmall = 2))
terr_upper <- as.character(format(round(summary(terr_cen_all_data[[1]])$fixed[5,4], 2), nsmall = 2))

# Plot labels.
temp_cor_label <- paste0("\U03B2 = ", temp_estimate, "\n[", temp_lower, ", ", temp_upper, "]")
mig_cor_label <- paste0("\U03B2 = ", mig_estimate, "\n[", mig_lower, ", ", mig_upper, "]")
tro_cor_label <- paste0("\U03B2 = ", tro_estimate, "\n[", tro_lower, ", ", tro_upper, "]")
terr_cor_label <- paste0("\U03B2 = ", terr_estimate, "\n[", terr_lower, ", ", terr_upper, "]")

# Use calculated p-values instead.
temp_cor_label <- paste0("\U03B2 = ", temp_estimate, "\np < 0.001")
mig_cor_label <- paste0("\U03B2 = ", mig_estimate, "\np < 0.001")
tro_cor_label <- paste0("\U03B2 = ", tro_estimate, "\np < 0.01")
terr_cor_label <- paste0("\U03B2 = ", terr_estimate, "\np < 0.001")


###############################################################################
                  #### Create continuous side plots ####

# Create settings for scatter plots.
point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                      size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                       labs(x = "", y = NULL,  colour = "Clade"),
                       scale_colour_manual(values = prum_clade_colours),
                       scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(0,4.5)),
                       theme_classic(base_size = 20),
                       theme(text = element_text(face = "bold"),
                             legend.position = "none",
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.title = element_text(size = rel(0.8))))

# Remove low predictions as family averages don't span as far as single species.
temp_plot_predictions <- temp_uncen_all_data[[2]] %>% filter(temp_log > 3.3)
temp_plot_predictions <- temp_plot_predictions %>% filter(temp_log < 7.2)

# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sexual_mean)) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1)

max_temp <- temp_plot_predictions$temp_log %>% max()
min_temp <- temp_plot_predictions$temp_log %>% min()
temp_range <- max_temp-min_temp
temp_x_label <- max_temp - (temp_range*0.2)
temp_x_label_2 <- min_temp + (temp_range*0.05)

###############################################################################
                  #### Create categorical side plots ####

# Plot settings.
categorical_settings <- list(geom_point(aes(group=family, colour = higher_clade, 
               size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
           position = position_jitter(0.2)),
  labs(x = "", y = NULL,  colour = "Clade"),
  scale_colour_manual(values = prum_clade_colours),
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(-0.1,4.7)),
  theme_classic(base_size = 20),
  theme(text = element_text(face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title.x = element_text(size = rel(0.8)),
        axis.title.y = element_text(size = rel(0.8))))



# Migration plots.
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=migration_binary, ymin = estimate__ - se__, 
                    ymax = estimate__ + se__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=migration_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())




# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territory_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=territory_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=territory_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())

# Add labels.
lab_x_pos <- 2.1
lab_y_pos <- 4.4
temp_plot <- temp_plot + annotate("text", x = 6.5, y = lab_y_pos, label = temp_cor_label, size = 7, fontface = 2) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = mig_cor_label, size = 7, fontface = 2) 
tro_plot <- tro_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = tro_cor_label, size = 7, fontface = 2) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = terr_cor_label, size = 7, fontface = 2) 





###############################################################################
                  ### Non-bold text version ####

# Create settings for scatter plots.
point_settings <- list(geom_point(aes(group=higher_clade, colour = higher_clade, 
                                      size = sqrt(clade_sum), alpha = sqrt(clade_sum))),
                       labs(x = "", y = NULL,  colour = "Clade"),
                       scale_colour_manual(values = prum_clade_colours),
                       scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(-0.1,4.7)),
                       theme_classic(base_size = 24),
                       theme(legend.position = "none",
                             axis.line = element_line(size = 1),
                             axis.ticks = element_line(size = 1),
                             axis.title = element_text(size = rel(0.8))))

# Plot settings.
categorical_settings <- list(geom_point(aes(group=family, colour = higher_clade, 
                                            size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
                                        position = position_jitter(0.2)),
                             labs(x = "", y = NULL,  colour = "Clade"),
                             scale_colour_manual(values = prum_clade_colours),
                             scale_y_continuous(breaks = c(0, 1, 2, 3, 4), limits = c(-0.1,4.7)),
                             theme_classic(base_size = 24),
                             theme(legend.position = "none",
                                   axis.line = element_line(size = 1),
                                   axis.ticks = element_line(size = 1),
                                   axis.title.x = element_text(size = rel(0.8)),
                                   axis.title.y = element_text(size = rel(0.8))))

# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sexual_mean)) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = estimate__ - se__, ymax = estimate__ + se__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


# Migration plots.
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=migration_binary, ymin = estimate__ - se__, 
                    ymax = estimate__ + se__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=migration_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())

# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), 
                              expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_binary, ymin = estimate__ - se__, 
                    ymax = estimate__ + se__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territory_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=territory_binary, ymin = estimate__ - se__, 
                    ymax = estimate__ + se__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=territory_binary, y = estimate__), size = 3.6)  + theme(axis.text.y = element_blank())


# Add labels.
lab_x_pos <- 2.1
lab_y_pos <- 4.4

# Non-bold version.
temp_plot <- temp_plot + annotate("text", x = temp_x_label, y = lab_y_pos, label = temp_cor_label, size = 5) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = mig_cor_label, size = 5) 
tro_plot <- tro_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = tro_cor_label, size = 5) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = terr_cor_label, size = 5) 




side_plots <- ggarrange(tro_plot + ylab("Sexual selection"), 
                        mig_plot,
                        terr_plot, 
                        temp_plot,
                        nrow = 1, ncol = 4, align = "h",
                        widths = c(1.2,1,1,1))


side_plots

ggsave("Plots/Presentation/uni_trait_plots.png", height = 5, width = 15, dpi = 600)



                            ######### END #########
###################################################################################


