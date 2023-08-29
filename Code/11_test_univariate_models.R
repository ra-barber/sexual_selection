###############################################################################
               # Assess the fit/plot univariate models  #
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


# Read in the actual models. (This takes up loads of RAM)
all_temp_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/temp_centered_all_models.rds")
all_mig_cen_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/mig_centered_all_models.rds")
all_tro_cen_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/tro_centered_all_models.rds")
all_terr_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/terr_centered_all_models.rds")

all_mig_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/mig_uncentered_all_models.rds")
all_tro_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/tro_uncentered_all_models.rds")

tic()
temp_p_value <- pd_to_p(last(p_direction(all_temp_models)[,2]))
mig_p_value <- pd_to_p(last(p_direction(all_mig_models)[,2]))
tro_p_value <- pd_to_p(last(p_direction(all_tro_models)[,2]))
terr_p_value <- pd_to_p(last(p_direction(all_terr_models)[,2]))
toc()

### Output ###
# > temp_p_value
# [1] 0
# > mig_p_value
# [1] 0
# > tro_p_value
# [1] 0.009375
# > terr_p_value
# [1] 0

rm(all_temp_models, all_mig_models, all_tro_models, all_terr_models)


###############################################################################
#### Read in the model data for plotting #####


# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/"
list.files(first_half, pattern = "data")

# Read in the models using all data certainties.
temp_cen_all_data <- readRDS(paste0(first_half, "temp_centered_all_data.rds"))
#temp_uncen_all_data <- readRDS(paste0(first_half, "temp_uncentered_all_data.rds"))
mig_cen_all_data <- readRDS(paste0(first_half, "mig_centered_all_data.rds"))
#mig_uncen_all_data <- readRDS(paste0(first_half, "mig_uncentered_all_data.rds"))
tro_cen_all_data <- readRDS(paste0(first_half, "tro_centered_all_data.rds"))
#tro_uncen_all_data <- readRDS(paste0(first_half, "tro_uncentered_all_data.rds"))
terr_cen_all_data <- readRDS(paste0(first_half, "terr_centered_all_data.rds"))
#terr_uncen_all_data <- readRDS(paste0(first_half, "terr_uncentered_all_data.rds"))


temp_cen_all_data[[3]]
mig_cen_all_data[[3]]
tro_cen_all_data[[3]]
terr_cen_all_data[[3]]

temp_cen_all_data[[1]]
mig_cen_all_data[[1]]
tro_cen_all_data[[1]]
terr_cen_all_data[[1]]


## Reading in all models takes up too much RAM ## 

# Read in the models using high data certainties.
temp_cen_high_data <- readRDS(paste0(first_half, "temp_centered_high_data.rds"))
# temp_uncen_high_data <- readRDS(paste0(first_half, "temp_uncentered_high_data.rds"))
mig_cen_high_data <- readRDS(paste0(first_half, "mig_centered_high_data.rds"))
# mig_uncen_high_data <- readRDS(paste0(first_half, "mig_uncentered_high_data.rds"))
tro_cen_high_data <- readRDS(paste0(first_half, "tro_centered_high_data.rds"))
# tro_uncen_high_data <- readRDS(paste0(first_half, "tro_uncentered_high_data.rds"))
terr_cen_high_data <- readRDS(paste0(first_half, "terr_centered_high_data.rds"))
# terr_uncen_high_data <- readRDS(paste0(first_half, "terr_uncentered_high_data.rds"))
# 
# 



################################################################################
#### Export summary tables ####


# Extract relevant coeffcient information.
tro_all_estimates <- summary(tro_cen_all_data[[1]])$fixed[5,c(1,3,4)]
mig_all_estimates <- summary(mig_cen_all_data[[1]])$fixed[5,c(1,3,4)]
terr_all_estimates <- summary(terr_cen_all_data[[1]])$fixed[5,c(1,3,4)]
temp_all_estimates <- summary(temp_cen_all_data[[1]])$fixed[5,c(1,3,4)]

all_estimates <- rbind(tro_all_estimates, mig_all_estimates, 
                       terr_all_estimates, temp_all_estimates)

# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export the results.
write.csv(all_estimates, "Results/Tables/all_univariate_regression.csv", row.names = TRUE)


# High certainty,
tro_high_estimates <- summary(tro_cen_high_data[[1]])$fixed[5,c(1,3,4)]
mig_high_estimates <- summary(mig_cen_high_data[[1]])$fixed[5,c(1,3,4)]
terr_high_estimates <- summary(terr_cen_high_data[[1]])$fixed[5,c(1,3,4)]
temp_high_estimates <- summary(temp_cen_high_data[[1]])$fixed[5,c(1,3,4)]

high_estimates <- rbind(tro_high_estimates, mig_high_estimates, 
                        terr_high_estimates, temp_high_estimates)

# Paste together values for reporting in a table.
high_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export the results.
write.csv(high_estimates, "Results/Tables/high_univariate_regression.csv", row.names = TRUE)


###############################################################################
#### Prepare family level data #####

# Read in model data.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- model_data$birdtree_name %>% str_replace(" ", "_")
row.names(model_data) <- model_data$tree_tip


# Set as factor, then re-level for appropriate reference group.
model_data %<>% mutate(
  territory_binary = relevel(as.factor(territory_binary), ref = "No territory"),
  migration_binary = relevel(as.factor(migration_binary), ref = "Weak"),
  trophic_binary = relevel(as.factor(trophic_binary), ref = "Secondary"),
)

# Center categorical predictors.
model_data %<>% mutate(
  terr_bi_c = center_categorical(territory_binary),
  migration_bi_c = center_categorical(migration_binary),
  trophic_level_c = center_categorical(trophic_binary),
)

library(effectsize)
# Scale continuous predictors to two SD.
model_data %<>% mutate(
  temp_seasonality_z = standardize(temp_log, two_sd = TRUE),
)


# Read in clade function and assign.
source("Code/clade_function.R")
model_data %<>% assign_prum_clades()

# Average traits for family.
family_average_data_data <- model_data %>% group_by(family) %>% 
  summarise(higher_clade = first(higher_clade),
            tree_tip = first(tree_tip),
            sexual_mean = mean(sexual_score),
            clade_sum = length(sexual_score),
            temp_seasonality_z = mean(temp_seasonality_z),
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
teritory_data <- average_family("terr_bi_c")
migration_data <- average_family("migration_bi_c")
trophic_data <- average_family("trophic_level_c")

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
temp_cen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
mig_cen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
tro_cen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
terr_cen_all_data[[2]] %<>% 
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
temp_plot_predictions <- temp_cen_all_data[[2]] #%>% filter(temp_seasonality_z > 3.3)
#temp_plot_predictions <- temp_plot_predictions %>% filter(temp_seasonality_z < 7.2)

# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_seasonality_z, y=sexual_mean)) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_seasonality_z, ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1)

max_temp <- temp_plot_predictions$temp_seasonality_z %>% max()
min_temp <- temp_plot_predictions$temp_seasonality_z %>% min()
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
mig_plot <- migration_data %>% ggplot(aes(x=migration_bi_c, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_cen_all_data[[2]], inherit.aes = FALSE,
                aes(x=migration_bi_c, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_cen_all_data[[2]], inherit.aes = FALSE,
             aes(x=migration_bi_c, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())


tro_plot_predictions <- tro_cen_all_data[[2]] %>% filter(trophic_level_c < -0.27 | trophic_level_c > 0.72)

# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_c, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c(-0.2798902, 0.7201098), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_plot_predictions, inherit.aes = FALSE,
                aes(x=trophic_level_c, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_plot_predictions, inherit.aes = FALSE,
             aes(x=trophic_level_c, y = estimate__), size = 3.6)
#

# Trophic plots.


   all_trop_data <- model_data %>% dplyr::group_by(trophic_level_c) %>% 
     summarise(sexual_mean = mean(sexual_score),
               clade_size = length(sexual_score),
               sex_sd = sd(sexual_score),
               sex_se = sd(sexual_score)/sqrt(length(sexual_score)))
   
   
   all_terr_data <- model_data %>% dplyr::group_by(terr_bi_c) %>% 
     summarise(sexual_mean = mean(sexual_score),
               clade_size = length(sexual_score),
               sex_sd = sd(sexual_score),
               sex_se = sd(sexual_score)/sqrt(length(sexual_score)))
   
   

   trophic_data %>% ggplot(aes(x=trophic_level_c, y=sexual_mean)) + 
     categorical_settings +
     scale_x_discrete(limits = c( 0.7201098, -0.2798902), labels = c("Primary", "Secondary")) +
     xlab("Trophic level") +
     geom_errorbar(data = all_trop_data, inherit.aes = FALSE,
                   aes(x=trophic_level_c, ymin = sexual_mean - (sex_se*1.96), 
                       ymax = sexual_mean + (sex_se*1.96)),
                   position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
     geom_point(data = all_trop_data, inherit.aes = FALSE,
                aes(x=trophic_level_c, y = sexual_mean), size = 3.6)
   
   teritory_data %>% ggplot(aes(x=terr_bi_c, y=sexual_mean)) + 
     categorical_settings +
     scale_x_discrete(labels = c("No", "Yes")) +
     xlab("Territoriality") +
     geom_errorbar(data = all_terr_data, inherit.aes = FALSE,
                   aes(x=terr_bi_c, ymin = sexual_mean - (sex_se*1.96), 
                       ymax = sexual_mean + (sex_se*1.96)),
                   position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
     geom_point(data = all_terr_data, inherit.aes = FALSE,
                aes(x=terr_bi_c, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())
   
   
tro_cen_all_data[[1]]

library(brms)
?conditional_effects

conditional_effects(tro_cen_all_data[[1]], re_formula = NULL)

uncen_test <- conditional_effects(tro_uncen_all_data[[1]], re_formula = NULL)
# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territory_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_cen_all_data[[2]], inherit.aes = FALSE,
                aes(x=territory_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_cen_all_data[[2]], inherit.aes = FALSE,
             aes(x=territory_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())

# Add labels.
lab_x_pos <- 2.1
lab_y_pos <- 4.4
temp_plot <- temp_plot + annotate("text", x = 6.5, y = lab_y_pos, label = temp_cor_label, size = 7, fontface = 2) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = mig_cor_label, size = 7, fontface = 2) 
tro_plot <- tro_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = tro_cor_label, size = 7, fontface = 2) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = terr_cor_label, size = 7, fontface = 2) 



###############################################################################
### Group the plots together ####

side_plots <- ggarrange(temp_plot + ylab("Sexual selection"), 
                        mig_plot,# + rremove("y.text"), 
                        tro_plot + ylab("Sexual selection"), 
                        terr_plot,# + rremove("y.text"), 
                        nrow = 2, ncol = 2, align = "h",
                        labels = c("b", "c", "d", "e"), 
                        font.label = list(size = 24),
                        hjust = c(-3.8, -1.7, -3.8, -1.7),
                        widths = c(1.2,1))
side_plots


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
              aes(x = temp_log, ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


# Migration plots.
mig_plot <- migration_data %>% ggplot(aes(x=migration_bi_c, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_cen_all_data[[2]], inherit.aes = FALSE,
                aes(x=migration_bi_c, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_cen_all_data[[2]], inherit.aes = FALSE,
             aes(x=migration_bi_c, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())

# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), 
                              expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_cen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_cen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territory_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_cen_all_data[[2]], inherit.aes = FALSE,
                aes(x=territory_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_cen_all_data[[2]], inherit.aes = FALSE,
             aes(x=territory_binary, y = estimate__), size = 3.6) 


# Non-bold version.
temp_plot <- temp_plot + annotate("text", x = temp_x_label, y = lab_y_pos, label = temp_cor_label, size = 5) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = mig_cor_label, size = 5) 
tro_plot <- tro_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = tro_cor_label, size = 5) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = terr_cor_label, size = 5) 


side_plots <- ggarrange(temp_plot + ylab("Sexual selection"), 
                        mig_plot,# + rremove("y.text"), 
                        tro_plot + ylab("Sexual selection"), 
                        terr_plot,# + rremove("y.text"), 
                        nrow = 2, ncol = 2, align = "h",
                        labels = c("b", "c", "d", "e"), 
                        font.label = list(size = 24),
                        hjust = c(-3.8, -1.7, -3.8, -1.7),
                        widths = c(1.2,1))
side_plots

lab_x_pos_2 <- 0.6
lab_y_pos_2 <- 4.6

# temp_plot <- temp_plot + annotate("text", x = temp_x_label_2, y = lab_y_pos_2, label = "b", size = 10, fontface = 2) 
# mig_plot <- mig_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "c", size = 10, fontface = 2) 
# tro_plot <- tro_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "d", size = 10, fontface = 2) 
# terr_plot <- terr_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "e", size = 10, fontface = 2) 
# 
# side_plots <- ggarrange(temp_plot + ylab("Sexual selection"), 
#                         mig_plot,
#                         tro_plot + ylab("Sexual selection"), 
#                         terr_plot, 
#                         nrow = 2, ncol = 2, align = "h",
#                         widths = c(1.2,1))


## Rearrange for joes order.

tro_plot <- tro_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "b", size = 10, fontface = 2) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "c", size = 10, fontface = 2) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "d", size = 10, fontface = 2) 
temp_plot <- temp_plot + annotate("text", x = temp_x_label_2, y = lab_y_pos_2, label = "e", size = 10, fontface = 2) 

side_plots <- ggarrange(tro_plot + ylab("Sexual selection"), 
                        mig_plot,
                        terr_plot + ylab("Sexual selection"), 
                        temp_plot,
                        nrow = 2, ncol = 2, align = "h",
                        widths = c(1.2,1))
side_plots



side_plots

######### END #########
