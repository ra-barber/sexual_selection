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
library(tidybayes)

# Clear the workspace.
rm(list=ls())

# Functions.
source("Code/functions.R")


###############################################################################
                  #### Read in the models #####


# # Read in the actual models. (This takes up loads of RAM)
# all_temp_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/temp_centered_all_models.rds")
# all_mig_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/mig_uncentered_all_models.rds")
# all_tro_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/tro_uncentered_all_models.rds")
# all_terr_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Univariate/terr_centered_all_models.rds")
# 
# tic()
# temp_p_value <- pd_to_p(last(p_direction(all_temp_models)[,2]))
# mig_p_value <- pd_to_p(last(p_direction(all_temp_models)[,2]))
# tro_p_value <- pd_to_p(last(p_direction(all_tro_models)[,2]))
# terr_p_value <- pd_to_p(last(p_direction(all_terr_models)[,2]))
# toc()
# 
# ### Output ###
# # > temp_p_value
# # [1] 0
# # > mig_p_value
# # [1] 0
# # > tro_p_value
# # [1] 0.00725
# # > terr_p_value
# # [1] 0
# 
# rm(all_temp_models, all_mig_models, all_tro_models, all_terr_models)
# 

###############################################################################
             #### Read in the model data for plotting #####


# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Old_models/Combined_models/Univariate/"

# Read in the models using all data certainties.
temp_uncen_all_data <- readRDS(paste0(first_half, "temp_uncentered__data.rds"))
mig_uncen_all_data <- readRDS(paste0(first_half, "mig_uncentered__data.rds"))
tro_uncen_all_data <- readRDS(paste0(first_half, "tro_uncentered__data.rds"))
terr_uncen_all_data <- readRDS(paste0(first_half, "terr_uncentered__data.rds"))

# Read in centered seasonality to compare.
temp_cen_all_data <- readRDS(paste0(first_half, "temp_centered__data.rds"))

# Read in new consensus models.
consensus_temp_high_model <- readRDS("Z:/home/sexual_selection/Results/Models/Consensus/Univariate/ordinal_temp_high.rds")
consensus_temp_z_high_model <- readRDS("Z:/home/sexual_selection/Results/Models/Consensus/Univariate/ordinal_temp_z_high.rds")
consenus_tro_high_model <-  readRDS("Z:/home/sexual_selection/Results/Models/Consensus/Univariate/ordinal_tro_high.rds")




###############################################################################
                 #### Prepare family level data #####

# Read in model data.
model_data <- read_ss_data("Data/sexual_selection_dataset_12_10.xlsx")
row.names(model_data) <- model_data$tree_tip

# Read in clade function and assign.
source("Code/clade_function.R")
model_data %<>% assign_prum_clades()

# Remove migration and seasonality NA data.
model_data %<>% tidyr::drop_na(migration, seasonality)

# Scale seasonality.
model_data %<>% mutate(temp_seasonality_z = scale(log(seasonality)))

# Average traits for family.
family_seasonality_data <- model_data %>% group_by(family_bird_tree) %>% 
  summarise(higher_clade = first(higher_clade),
            family_ss = mean(sexual_selection),
            clade_sum = length(sexual_selection),
            temp_log = mean(log(seasonality)),
            temp_seasonality_z = mean(temp_seasonality_z))

# Group by family and trait for binary traits.
average_family <- function(varible){
  var_data <- model_data %>% dplyr::group_by(family_bird_tree, !!! syms(varible)) %>% 
    summarise(higher_clade = first(higher_clade),
              family_ss = mean(sexual_selection),
              clade_size = length(sexual_selection))
  var_data %>% group_by(family_bird_tree) %>% summarise(clade_sum = sum(clade_size)) %>% right_join(var_data)
}

# Get territory and migration family summaries.
teritory_data <- average_family("territoriality_binary")
migration_data <- average_family("migration_binary")
trophic_data <- average_family("trophic_level_binary")


##############################################################################
            ######  Extract predicted values ######

# Create new data for extracting predictions.
tro_newdata <- data.frame(trophic_binary = c("Secondary", "Primary"), tree_tip = "typical_species")
mig_newdata <- data.frame(migration_binary = c("Weak", "Strong"), tree_tip = "typical_species")
terr_newdata <- data.frame(territory_binary = c("No territory", "Territory"), tree_tip = "typical_species")
temp_newdata <- data.frame(temp_log = c(3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7), tree_tip = "typical_species")

# Function to convert ordinal expected predictions back to response variable.
post_epred_ord_means <- function(fit, newdata = fit$data, draw_number = ndraws(fit), new_level_method = "uncertainty"){
  dta <- fit$data
  y <- all.vars(as.formula(fit$formula))[1]
  lvls <- sort(unique(as.numeric(as.character(dta[[y]]))))
  pep <- posterior_epred(fit, newdata, seed = 1234, 
                         allow_new_levels = TRUE, 
                         sample_new_levels = new_level_method, 
                         ndraws = draw_number)  #incl_autocor = FALSE
  nsamp <- dim(pep)[1]
  npts <- dim(pep)[2]
  m <- matrix(nrow = nsamp, ncol = npts)
  for (i in 1:nsamp){
    for (j in 1:npts){
      m[i, j] <- sum(pep[i, j, ] * lvls)
    }
  }
  m
}


# Model predictions using a sample of 1000 draws to increase speed. (Still takes ~60mins)
tro_preds <- post_epred_ord_means(tro_uncen_all_data[[1]], newdata = tro_newdata, draw_number = 1000)
mig_preds <- post_epred_ord_means(mig_uncen_all_data[[1]], newdata = mig_newdata, draw_number = 1000)
terr_preds <- post_epred_ord_means(terr_uncen_all_data[[1]], newdata = terr_newdata, draw_number = 1000)
temp_preds <- post_epred_ord_means(temp_uncen_all_data[[1]], newdata = temp_newdata, draw_number = 1000)

all_predictions <- list(tro_preds, mig_preds, terr_preds, temp_preds)
saveRDS(all_predictions, "Results/Models/Univariate/model_predictions.rds")

# Extract temperature summaries.
temp_est <- apply(temp_preds, 2, mean)
temp_err <- apply(temp_preds, 2, stanard_error)
temp_summaries <- data.frame(temp_log = c(3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7),
                             estimate = temp_est,
                             error = temp_err,
                             lower = temp_est - (1.96*temp_err),
                             upper = temp_est + (1.96*temp_err))

# Extract summaries of categorical variables.
extract_summaries <- function(matrix, levels =  c("Reference", "Main")){
  ref_mean <- mean(matrix[,1])
  main_mean <- mean(matrix[,2])
  
  ref_se <- stanard_error(matrix[,1])
  main_se <- stanard_error(matrix[,2])
  
  main_lower <- main_mean - (1.96*main_se)
  main_upper <- main_mean + (1.96*main_se)
  
  ref_lower <- ref_mean - (1.96*ref_se)
  ref_upper <- ref_mean + (1.96*ref_se)
  
  data.frame(effect_level = levels,
             estimate = c(ref_mean, main_mean),
             error = c(ref_se, main_se),
             lower = c(ref_lower, main_lower),
             upper = c(ref_upper, main_upper))
  
}

# Extract categorical variables.
tro_summaries <- extract_summaries(tro_preds, levels =  c("Secondary", "Primary"))
mig_summaries <- extract_summaries(mig_preds, levels = c("Weak", "Strong"))
terr_summaries <- extract_summaries(terr_preds, levels = c("Non-territorial", "Territorial"))


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

# Extract estimates.
extract_estimate <- function(model){
  as.character(format(round(summary(model)$fixed[5,1], 2), nsmall = 2))
}

temp_estimate <- extract_estimate(temp_cen_all_data[[1]])
mig_estimate <- extract_estimate(mig_uncen_all_data[[1]])
tro_estimate <- extract_estimate(tro_uncen_all_data[[1]])
terr_estimate <- extract_estimate(terr_uncen_all_data[[1]])

# Use calculated p-values from beginning of script.
temp_cor_label <- paste0("\U03B2 = ", temp_estimate, "\np < 0.001")
mig_cor_label <- paste0("\U03B2 = ", mig_estimate, "\np < 0.001")
tro_cor_label <- paste0("\U03B2 = ", tro_estimate, "\np < 0.01")
terr_cor_label <- paste0("\U03B2 = ", terr_estimate, "\np < 0.001")

# Plot settings.
categorical_settings <- list(
  geom_point(aes(group=family_bird_tree, colour = higher_clade, 
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

# Label positions.
lab_x_pos <- 2.1
lab_y_pos <- 1.8
lab_x_pos_2 <- 0.6
lab_y_pos_2 <- 1.8

# Temp labels.
max_temp <- family_seasonality_data$temp_log %>% max()
min_temp <- family_seasonality_data$temp_log %>% min()
temp_range <- max_temp-min_temp
temp_x_label <- max_temp - (temp_range*0.2)
temp_x_label_2 <- min_temp + (temp_range*0.05)


################################################################################
                      ##### Create plots ######


add_errorbars <- function(plot, summary_data){
  plot +
    geom_errorbar(data = summary_data, inherit.aes = FALSE,
                  aes(x=effect_level, ymin = log(lower), 
                      ymax = log(upper)),
                  position = position_dodge(width = 1), show.legend = FALSE, width = 0.4, linewidth = 1) +
    geom_point(data = summary_data, inherit.aes = FALSE,
               aes(x=effect_level, y = log(estimate)), size = 2) 
}

# Trophic level.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=log(family_ss + 1))) + 
  categorical_settings + ylim(c(-0.01, 1.9)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") 
tro_plot <- add_errorbars(tro_plot, tro_summaries)

# Migration.
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=log(family_ss+1))) + 
  categorical_settings + ylim(c(-0.01, 1.9)) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") + theme(axis.text.y = element_blank())
mig_plot <- add_errorbars(mig_plot, mig_summaries)

# Territory.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=log(family_ss+1))) + 
  categorical_settings +ylim(c(-0.01, 1.9)) +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality")
terr_plot <- add_errorbars(terr_plot, terr_summaries)

# Seasonality.
temp_plot <- family_seasonality_data %>% ggplot(aes(x=temp_log, y=log(family_ss+1))) + 
  point_settings + xlab("Seasonality") +ylim(c(-0.01, 1.9)) + xlim(c(3.3,7.2)) +
  geom_ribbon(data = temp_summaries, inherit.aes = FALSE,
              aes(x = temp_log, ymin = log(lower), ymax = log(upper)), 
              fill = "grey70", colour = NA, size = 2, alpha = 0.5) +
  geom_line(data = temp_summaries, aes(x = temp_log, y = log(estimate)), 
            linetype = "dashed", linewidth = 1.5, colour = "black") + 
  theme(axis.text.y = element_blank())


# Add estimates.
temp_plot <- temp_plot + annotate("text", x = temp_x_label, y = lab_y_pos, label = temp_cor_label, size = 5) 
mig_plot <-  mig_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = mig_cor_label, size = 5) 
tro_plot <- tro_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = tro_cor_label, size = 5) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = terr_cor_label, size = 5) 

# Add plot labels.
tro_plot <- tro_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "b", size = 10, fontface = 2) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "c", size = 10, fontface = 2) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "d", size = 10, fontface = 2) 
temp_plot <- temp_plot + annotate("text", x = temp_x_label_2, y = lab_y_pos_2, label = "e", size = 10, fontface = 2) 

# Plot together.
side_plots_preds <- ggarrange(tro_plot + ylab("Sexual selection (log)"), 
                              mig_plot,
                              terr_plot + ylab("Sexual selection (log)"), 
                              temp_plot,
                              nrow = 2, ncol = 2, align = "h",
                              widths = c(1.2,1))

# Export.
saveRDS(side_plots_preds, "Plots/Trees/side_plots.rds")


################################################################################
                      ######## END #########
################################################################################



# Get overall averages.
tropic_se <- model_data %>% group_by(trophic_level_binary) %>%
  summarise(sexual_mean = mean(sexual_selection),
            sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))

terr_se <- model_data %>% group_by(territoriality_binary) %>%
  summarise(sexual_mean = mean(sexual_selection),
            sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))

mig_se <- model_data %>% group_by(migration_binary) %>%
  summarise(sexual_mean = mean(sexual_selection),
            sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))





################################################################################
            ###### Create categorical plots ###### 



# Migration plot.
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=family_ss)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration")  +
  geom_errorbar(data = mig_se, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = sexual_mean - (1.96*sexual_se), 
                    ymax = sexual_mean + (1.96*sexual_se)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_se, inherit.aes = FALSE,
             aes(x=migration_binary, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())

# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=family_ss)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sexual_mean - (1.96*sexual_se), 
                    ymax = sexual_mean + (1.96*sexual_se)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=family_ss)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_se, inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = sexual_mean - (1.96*sexual_se), 
                    ymax = sexual_mean + (1.96*sexual_se)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_se, inherit.aes = FALSE,
             aes(x=territoriality_binary, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())

################################################################################
                 #### Create temperature side plot ####


# Change predictions back to same scale as original sex scores.
temp_plot_predictions <- temp_uncen_all_data[[2]] %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1) %>% 
  filter(temp_log > 3.3 & temp_log < 7.2)

# Make some nice plots.
temp_plot <- family_seasonality_data %>% ggplot(aes(x=temp_log, y=family_ss)) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


################################################################################
                   ##### Add the labels ######

# Temp labels.
max_temp <- temp_plot_predictions$temp_log %>% max()
min_temp <- temp_plot_predictions$temp_log %>% min()
temp_range <- max_temp-min_temp
temp_x_label <- max_temp - (temp_range*0.2)
temp_x_label_2 <- min_temp + (temp_range*0.05)

lab_x_pos <- 2.1
lab_y_pos <- 4.4

lab_x_pos_2 <- 0.6
lab_y_pos_2 <- 4.6

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


###############################################################################
                     #####  Log scale ########


mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=log(family_ss+1))) + 
  categorical_settings + ylim(c(-0.01, 1.9)) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_se, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = log(sexual_mean + 1 - (1.96*sexual_se)), 
                    ymax = log(sexual_mean + 1 + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_se, inherit.aes = FALSE,
             aes(x=migration_binary, y = log(sexual_mean + 1)), size = 3.6) + theme(axis.text.y = element_blank())

tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=log(family_ss+1))) + 
  categorical_settings +ylim(c(-0.01, 1.9)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = log(sexual_mean + 1 - (1.96*sexual_se)), 
                    ymax = log(sexual_mean + 1 + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = log(sexual_mean + 1)), size = 3.6) + theme(axis.text.y = element_blank())

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=log(family_ss+1))) + 
  categorical_settings +ylim(c(-0.01, 1.9)) +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_se, inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = log(sexual_mean + 1 - (1.96*sexual_se)), 
                    ymax = log(sexual_mean + 1 + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_se, inherit.aes = FALSE,
             aes(x=territoriality_binary, y = log(sexual_mean + 1)), size = 3.6) + theme(axis.text.y = element_blank())


# Label positions.
lab_x_pos <- 2.1
lab_y_pos <- 1.8
lab_x_pos_2 <- 0.6
lab_y_pos_2 <- 1.8

# Add estimates.
temp_plot <- temp_plot + annotate("text", x = temp_x_label, y = lab_y_pos, label = temp_cor_label, size = 5) 
mig_plot <-  mig_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = mig_cor_label, size = 5) 
tro_plot <- tro_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = tro_cor_label, size = 5) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos, y = lab_y_pos, label = terr_cor_label, size = 5) 

# Add plot labels.
tro_plot <- tro_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "b", size = 10, fontface = 2) 
mig_plot <- mig_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "c", size = 10, fontface = 2) 
terr_plot <- terr_plot + annotate("text", x = lab_x_pos_2, y = lab_y_pos_2, label = "d", size = 10, fontface = 2) 
temp_plot <- temp_plot + annotate("text", x = temp_x_label_2, y = lab_y_pos_2, label = "e", size = 10, fontface = 2) 

side_plots_log <- ggarrange(tro_plot + ylab("Sexual selection"), 
                        mig_plot,
                        terr_plot + ylab("Sexual selection"), 
                        temp_plot,
                        nrow = 2, ncol = 2, align = "h",
                        widths = c(1.2,1))



################################################################################
                  ###### Normal scale #######


add_errorbars <- function(plot, summary_data){
  plot +
    geom_errorbar(data = summary_data, inherit.aes = FALSE,
                  aes(x=effect_level, ymin = lower, 
                      ymax = upper),
                  position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 2) +
    geom_point(data = summary_data, inherit.aes = FALSE,
               aes(x=effect_level, y = estimate), size = 5.6) 
}

# Trophic level.
tro_plot_2 <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=family_ss)) + 
  categorical_settings + ylim(c(-0.01, 4.7)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") 
tro_plot_2 <- add_errorbars(tro_plot_2, tro_summaries)

# Migration.
mig_plot_2 <- migration_data %>% ggplot(aes(x=migration_binary, y=family_ss)) + 
  categorical_settings + ylim(c(-0.01, 4.7)) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") + theme(axis.text.y = element_blank())
mig_plot_2 <- add_errorbars(mig_plot_2, mig_summaries)

# Territory.
terr_plot_2 <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=family_ss)) + 
  categorical_settings +ylim(c(-0.01, 4.7)) +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality")
terr_plot_2 <- add_errorbars(terr_plot_2, terr_summaries)

# Seasonality.
temp_plot_2 <- family_seasonality_data %>% ggplot(aes(x=temp_log, y=family_ss)) + point_settings +
  xlab("Seasonality") +ylim(c(-0.01, 4.7)) + xlim(c(3.3,7.2)) +
  geom_ribbon(data = temp_summaries, inherit.aes = FALSE,
              aes(x = temp_log, ymin = lower, ymax = upper), fill = "grey70", colour = NA, size = 1, alpha = 0.5) +
  geom_line(data = temp_summaries, aes(x = temp_log, y = estimate), linetype = "dashed", linewidth = 2, colour = "black") + 
  theme(axis.text.y = element_blank())


# Plot together.
side_plots_preds_2 <- ggarrange(tro_plot_2 + ylab("Sexual selection (log)"), 
                              mig_plot_2,
                              terr_plot_2 + ylab("Sexual selection (log)"), 
                              temp_plot_2,
                              nrow = 2, ncol = 2, align = "h",
                              widths = c(1.2,1))


################################################################################
                    ##### Try log scale ####


trophic_data %>% ggplot(aes(x=trophic_level_binary, y=log(family_ss+1))) + 
  categorical_settings +ylim(c(0, 1.7)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = log(sexual_mean + 1 - (2*sexual_se)), 
                    ymax = log(sexual_mean + 1 + (2*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = log(sexual_mean + 1)), size = 2.6) + theme(axis.text.y = element_blank())

tro_plot

###############################################################################

consensus_temp_preds <- conditional_effects(consensus_temp_high_model)[[1]]
consensus_tempz_preds <- conditional_effects(consensus_tempz_high_model)[[1]]

consensus_temp_preds

temp_plot_predictions_2 <- consensus_temp_preds %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1) %>% 
  filter(temp_log > 3.3 & temp_log < 7.2)

# Make some nice plots.
temp_plot <- family_seasonality_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions_2, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions_2, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())

new_temp <- data.frame(temp_log = seq(from = 3, to = 6, by = 0.5))
tic()
test_new_preds <- posterior_epred(consensus_temp_high_model, newdata = newdata, re_formula = NA)
toc()

str(test_new_preds)



library(emmeans)
tic()
testing_emmeans <- consensus_temp_high_model %>% emmeans( ~ temp_log, epred = TRUE, newdata = new_temp)
toc()

tic()
tro_emmeans <- tro_uncen_all_data[[1]] %>% emmeans( ~ trophic_binary, epred = TRUE)
toc()

epred_draws()

new_tro_data <- data.frame(trophic_level_binary = c("Primary", "Secondary"), tree_tip = NA)
consensus_tro <- conditional_effects(consenus_tro_high_model)

consensus_tro_NA <- posterior_epred(consenus_tro_high_model, new_data = new_tro_data, allow_new_levels = TRUE)

epred_predicted <- consenus_tro_high_model %>%
  epred_draws(newdata = new_tro_data)

epred_predicted %>% mean_hdi()

epred_predicted %>% group_by(trophic_level_binary) %>% summarise(pred_mean = mean(.epred*as.numeric(as.character(.category))))

consensus_tro_NA[[1]]

mean(consensus_tro_NA[,3,1])

(consensus_tro_NA)
# Expected values with tidybayes
tic()
tidy_epred <- tro_uncen_all_data[[1]] %>% 
  epred_draws(newdata = new_tro_data, re_formula = NA)
toc()

# Summarise.
tidy_epred %>% group_by(trophic_binary) %>% median_hdi(.epred)
tidy_epred %>% group_by(trophic_binary) %>% mean_hdi(.epred)

tidy_epred %>% filter(trophic_binary == "Primary") %>% pull(.epred) %>% mean()
tidy_epred %>% filter(trophic_binary == "Secondary") %>% pull(.epred) %>% mean()

# Try predicting just the one factor.
tic()
tidy_epred_2 <- tro_uncen_all_data[[1]] %>% 
  epred_draws(newdata = data.frame(trophic_binary = "Primary"), re_formula = NA)
toc()

tidy_epred_2 %>% group_by(trophic_binary) %>% mean_hdi(.epred)

mean_hdi()
plot_preds <- bind_rows(
  "Expectation of predicted draws" = rename(tidy_epred, .prediction = .epred),
  .id = "draw_type") 

plot_preds %>% 
  group_by(trophic_binary) %>% 
  median_hdi(.prediction)
library(emmeans)
tic()
temp_trends <- temp_uncen_all_data[[1]] %>% 
  emtrends(~ temp_log,
           var = "temp_log",
            at = list(temp_log = seq(from = 3.5, to = 7, by = 0.5)),
           epred = TRUE, re_formula = NA)
toc()
temp_trends %>% gather_emmeans_draws()

# Posterior predictions across civil_liberties
grand_mean_temp_uncen <- temp_uncen_all_data[[1]] %>% 
  epred_draws(newdata = new_temp, 
              re_formula = NA)

grand_mean_temp_uncen %<>% mutate(
  both = as.numeric(as.character(.category))*.epred
)

plot_grand_mean_civlib <- ggplot(grand_mean_temp_uncen, 
                                 aes(x = temp_log, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Civil liberties index", y = "Predicted media freedom index",
       fill = "Credible interval") +
  ggthemes::theme_clean() +
  theme(legend.position = "bottom")

library(tidyr)
grand_mean_temp_uncen %>% select(.category, .epred) %>%  pivot_wider(names_from = .category, values_from = .epred)

grand_mean_temp_uncen %>% mean_hdi()




consensus_temp_high_model

grand_mean_temp_uncen <- temp_uncen_all_data[[1]] %>% 
  epred_draws(newdata = new_temp, 
              re_formula = NA)

test_new_level <-  tro_uncen_all_data[[1]] %>% 
  emmeans(~ trophic_binary,
          at = list(trophic_binary = c("Secondary", "Primary"), 
                    tree_tip = "Generic"),
          epred = TRUE, re_formula = NULL, 
          allow_new_levels = TRUE, sample_new_levels = "gaussian") 

epred_draws(newdata = data.frame(trophic_binary = "Primary"), re_formula = NA)

test_new_level %>%  gather_emmeans_draws() %>% median_hdi()



post_epred_ord_means <- function(fit, newdata = fit$data, new_level_method = "uncertainty"){
  dta <- fit$data
  y <- all.vars(as.formula(fit$formula))[1]
  lvls <- sort(unique(as.numeric(as.character(dta[[y]]))))
  pep <- posterior_epred(fit, newdata, allow_new_levels = TRUE, sample_new_levels = new_level_method)
  nsamp <- dim(pep)[1]
  npts <- dim(pep)[2]
  m <- matrix(nrow = nsamp, ncol = npts)
  for (i in 1:nsamp){
    for (j in 1:npts){
      m[i, j] <- sum(pep[i, j, ] * lvls)
    }
  }
  m
}


post_epred_ord_means_10 <- function(fit, newdata = fit$data, draw_number = 10){
  dta <- fit$data
  y <- all.vars(as.formula(fit$formula))[1]
  lvls <- sort(unique(as.numeric(as.character(dta[[y]]))))
  pep <- posterior_epred(fit, newdata, ndraws = draw_number, allow_new_levels = TRUE)
  nsamp <- dim(pep)[1]
  npts <- dim(pep)[2]
  m <- matrix(nrow = nsamp, ncol = npts)
  for (i in 1:nsamp){
    for (j in 1:npts){
      m[i, j] <- sum(pep[i, j, ] * lvls)
    }
  }
  m
}







# Testing the new custom function with different predictions and inputs.

test_epred_function <- post_epred_ord_means(consenus_tro_high_model)
test_epred_function_10 <- post_epred_ord_means_10(consenus_tro_high_model)
test_epred_function_100 <- post_epred_ord_means_10(consenus_tro_high_model, draw_number = 100)


new_tro_data_2 <- data.frame(trophic_level_binary = c("Secondary", "Primary"), tree_tip = NA)
new_tro_data_3 <- data.frame(trophic_level_binary = c("Secondary", "Primary"), tree_tip = "new_species")

test_epred_function_new <- post_epred_ord_means(consenus_tro_high_model, newdata = new_tro_data_2)
test_epred_function_10_new <- post_epred_ord_means_10(consenus_tro_high_model, newdata = new_tro_data_2)
test_epred_function_100_new <- post_epred_ord_means_10(consenus_tro_high_model, draw_number = 100, newdata = new_tro_data_2)

test_epred_function_new_new <- post_epred_ord_means(consenus_tro_high_model, newdata = new_tro_data_3)
test_epred_function_10_new_sp <- post_epred_ord_means_10(consenus_tro_high_model, newdata = new_tro_data_3)
test_epred_function_100_new_sp <- post_epred_ord_means_10(consenus_tro_high_model, draw_number = 100, newdata = new_tro_data_3)

# Test different methods for estimating typical random effect.
test_epred_new_sp_uncen <- post_epred_ord_means(consenus_tro_high_model, newdata = new_tro_data_3)
test_epred_new_sp_gaus <- post_epred_ord_means(consenus_tro_high_model, newdata = new_tro_data_3, new_level_method = "gaussian")
test_epred_function_sp_old <- post_epred_ord_means(consenus_tro_high_model, newdata = new_tro_data_3, new_level_method = "old_levels")

test_epred_function_new

quantile(test_epred_function_100_new[, 1], probs = c(0.5, 0.025, 0.975))
quantile(test_epred_function_100_new[, 2], probs = c(0.5, 0.025, 0.975))



extract_summaries <- function(matrix){
  sec_mean <- mean(matrix[,1])
  prim_mean <- mean(matrix[,2])
  
  sec_se <- stanard_error(matrix[,1])
  prim_se <- stanard_error(matrix[,2])
  
  prim_lower <- prim_mean - (1.96*prim_se)
  prim_upper <- prim_mean + (1.96*prim_se)
  
  sec_lower <- sec_mean - (1.96*sec_se)
  sec_upper <- sec_mean + (1.96*sec_se)
  
  data.frame(trophic_level_binary = c("Secondary", "Primary"),
             epred_mean = c(sec_mean, prim_mean),
             epred_se = c(sec_se, prim_se),
             epred_lower = c(sec_lower, prim_lower),
             epred_upper = c(sec_upper, prim_upper))
  
}

extract_summaries(test_epred_function_new)
extract_summaries(test_epred_function_new_new)

extract_summaries(test_epred_function_new)
extract_summaries(test_epred_function_10_new)
extract_summaries(test_epred_function_100_new)


mean(test_epred_function_new[, 1])
mean(test_epred_function_new[, 2])


stanard_error(test_epred_function_new[, 1])
stanard_error(test_epred_function_new[, 2])

test_tro_model <- data.frame(trophic_level_binary = c("Secondary", "Primary"),
  epred_mean = c(mean(test_epred_function_new[, 1]), mean(test_epred_function_new[, 2])),
                             epred_se = stanard_error(test_epred_function_new[, 1]), stanard_error(test_epred_function_new[, 2]))

test_error_bars <- function(error_bar_data){
  error_bar_data %>% ggplot(aes(x = trophic_level_binary, y = sqrt(epred_mean-1), ymin = sqrt(epred_lower-1), ymax = sqrt(epred_upper-1))) +
    geom_errorbar() + geom_point() + ylim(c(0,2.3)) + theme_classic()
}

test_error_bars(extract_summaries(test_epred_function_new))
#test_error_bars(extract_summaries(test_epred_function_100_new))
test_error_bars(extract_summaries(test_epred_function_new_new))

test_error_bars <- function(error_bar_data){
  error_bar_data %>% ggplot(aes(x = trophic_level_binary, y = sqrt(epred_mean-1), ymin = sqrt(epred_lower-1), ymax = sqrt(epred_upper-1))) +
    geom_errorbar() + geom_point() + ylim(c(0,2.3)) + theme_classic()
}

extract_summaries_sqrt <- function(matrix){
  
  matrix[,1] <- sqrt(matrix[,1])
  matrix[,2] <- sqrt(matrix[,2])
  
  sec_mean <- mean(matrix[,1])
  prim_mean <- mean(matrix[,2])
  
  sec_se <- stanard_error(matrix[,1])
  prim_se <- stanard_error(matrix[,2])
  
  prim_lower <- prim_mean - (1.96*prim_se)
  prim_upper <- prim_mean + (1.96*prim_se)
  
  sec_lower <- sec_mean - (1.96*sec_se)
  sec_upper <- sec_mean + (1.96*sec_se)
  
  data.frame(trophic_level_binary = c("Secondary", "Primary"),
             epred_mean = c(sec_mean, prim_mean),
             epred_se = c(sec_se, prim_se),
             epred_lower = c(sec_lower, prim_lower),
             epred_upper = c(sec_upper, prim_upper))
  
}


extract_summaries(test_epred_function_new)
extract_summaries_sqrt(test_epred_function_new)
test_error_bars_sqrt(extract_summaries_sqrt(test_epred_function_new))
test_error_bars(extract_summaries(test_epred_function_new))


test_error_bars_sqrt <- function(error_bar_data){
  error_bar_data %>% ggplot(aes(x = trophic_level_binary, y = epred_mean-1, ymin = epred_lower-1, ymax = epred_upper-1)) +
    geom_errorbar() + geom_point() + ylim(c(0,2.3)) + theme_classic()
}


hist(test_epred_function_new[,1])
hist(test_epred_function_new[,2])
test_epred_function_new_new

mean(test_epred_function_new_new[, 1])
mean(test_epred_function_new_new[, 2])
stanard_error(test_epred_function_new_new[, 1])
stanard_error(test_epred_function_new_new[, 2])
qi(test_epred_function_new_new[, 1])


trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sqrt(family_ss))) + 
  categorical_settings +ylim(c(0, 2.3)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = test_tro_model, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sqrt(epred_mean-1 - (1.96*epred_se)), 
                    ymax = sqrt(epred_mean-1 + (1.96*epred_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = test_tro_model, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sqrt(epred_mean-1)), size = 3.6) + theme(axis.text.y = element_blank())












test_lin_pred <- posterior_linpred(consenus_tro_high_model)
test_lin_pred_2 <- posterior_linpred(consenus_tro_high_model, transform = TRUE)
str(test_lin_pred)
test_custom_function <- post_epred_ord_means(consenus_tro_high_model)

fit_data$new_function <- apply(test_lin_pred_2, 2, mean)



str(test_custom_function)
fit_data$new_function <- apply(test_custom_function, 2, mean)
fit_data <- consenus_tro_high_model$data


new_model_preds <- fit_data %>% group_by(trophic_level_binary) %>% summarise(mean_nf = mean(new_function-1),
                                                          se_nf = sd(new_function-1)/sqrt(length(new_function)))
fit_data %>% group_by(trophic_level_binary) %>% qi()
fit_data %>% filter(trophic_level_binary == "Primary") %>% pull(new_function) %>% qi()
trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sqrt(family_ss))) + 
  categorical_settings +ylim(c(0, 2.7)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = new_model_preds, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sqrt(mean_nf - (1.96*se_nf)), 
                    ymax = sqrt(mean_nf + (1.96*se_nf))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = new_model_preds, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sqrt(mean_nf)), size = 3.6) + theme(axis.text.y = element_blank())



# 1*pnorm(0.337) + 
#   2*(pnorm(0.359) - pnorm(0.337)) +
#   3*(pnorm(0.186) - pnorm(0.359)) +
#   4*(pnorm(0.117) - pnorm(0.186)) +
#   5*(pnorm(0.0000967) - pnorm(0.117)) +
#   6*(1 - pnorm(0.0000967))
# 
# test <- grand_mean_temp_uncen %>% filter(temp_log == 6)  %>% mean_hdi()
# 
# 1*pnorm(test$.epred[1]) + 
#   2*(pnorm(test$.epred[2]) - pnorm(test$.epred[1])) +
#   3*(pnorm(test$.epred[3]) - pnorm(test$.epred[2])) +
#   4*(pnorm(test$.epred[4]) - pnorm(test$.epred[3])) +
#   5*(pnorm(test$.epred[5]) - pnorm(test$.epred[4])) +
#   6*(1 - pnorm(test$.epred[5]))




################################################################################
                  ######                #######


mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sqrt(family_ss))) + 
  categorical_settings + ylim(c(0, 2.7)) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration")  +
  geom_errorbar(data = mig_se, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = sqrt(sexual_mean - (1.96*sexual_se)), 
                    ymax = sqrt(sexual_mean + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_se, inherit.aes = FALSE,
             aes(x=migration_binary, y = sqrt(sexual_mean)), size = 3.6) + theme(axis.text.y = element_blank())

tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sqrt(family_ss))) + 
  categorical_settings +ylim(c(0, 2.7)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sqrt(sexual_mean - (1.96*sexual_se)), 
                    ymax = sqrt(sexual_mean + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sqrt(sexual_mean)), size = 3.6) + theme(axis.text.y = element_blank())

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=sqrt(family_ss))) + 
  categorical_settings +ylim(c(0, 2.7)) +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_se, inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = sqrt(sexual_mean - (1.96*sexual_se)), 
                    ymax = sqrt(sexual_mean + (1.96*sexual_se))),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_se, inherit.aes = FALSE,
             aes(x=territoriality_binary, y = sqrt(sexual_mean)), size = 3.6) + theme(axis.text.y = element_blank())

# Make some nice plots.
temp_plot <- family_seasonality_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + ylim(c(0,2.7)) +
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())











# temp_uncen_all_data[[2]]
# temp_plot_predictions$lower__[temp_plot_predictions$lower__ < 0] <- 0


# Temp_predictions
centered_linear_preds <- conditional_effects(temp_cen_all_data[[1]], method = "posterior_linpred")
?conditional_effects



temp_test <- temp_uncen_all_data[[2]] %>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)

temp_test <- temp_cen_all_data[[2]] %>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)






temp_model <- temp_uncen_all_data[[1]]

temp_plot_predictions <- temp_test %>% filter(temp_log > 3.3 & temp_log < 7.2)
temp_plot_predictions$lower__[temp_plot_predictions$lower__ < 0] <- 0


family_seasonality_data %>% ggplot(aes(x=temp_seasonality_z, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + ylim(c(0,2)) +
  geom_ribbon(data = temp_test, inherit.aes = FALSE,
              aes(x = temp_seasonality_z, ymin = sqrt(lower__), ymax = sqrt(upper__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_test, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())




# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(estimate__-(1.96*se__)), ymax = sqrt(estimate__+(1.96*se__))), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=family_ss)) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


temp_plot <- family_seasonality_data %>% ggplot(aes(x=temp_log, y=family_ss)) + point_settings +
  xlab("Seasonality") + 
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())











family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + geom_smooth(method = "lm", data = temp_plot_predictions, aes(y = sqrt(estimate__))) +
  geom_ribbon(data = temp_plot_predictions, mapping = aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), inherit.aes = FALSE, fill = "grey70", colour = NA, alpha = 0.2)





# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + ylim(c(0,2)) +
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


# family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
#   xlab("Seasonality") + 
#   geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
#               aes(x = temp_log, ymin = sqrt(estimate__)-sqrt((1.96*se__)), ymax = sqrt(estimate__)+sqrt(1.96*se__)), fill = "grey70", colour = NA, alpha = 0.2) +
#   geom_line(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


family_average_data_data %>% ggplot(aes(x=temp_log, y=family_ss)) + point_settings +
  xlab("Seasonality") +
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = lower__, ymax = upper__), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())


family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + geom_smooth(method = "lm", data = temp_plot_predictions, aes(y = sqrt(estimate__))) +
  geom_ribbon(data = temp_plot_predictions, mapping = aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), inherit.aes = FALSE, fill = "grey70", colour = NA, alpha = 0.2)

family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + geom_smooth(method = "lm", data = temp_lin_preds, aes(y = sqrt(estimate__))) +
  geom_ribbon(data = temp_lin_preds, mapping = aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), inherit.aes = FALSE, fill = "grey70", colour = NA, alpha = 0.2)




?geom_smooth
temp_lin_preds

################################################################################
                 #### Using family averages ######


mig_plot +
  geom_errorbar(data = migration_average, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = mean_ss - (1.96*error_ss), 
                    ymax = mean_ss + (1.96*error_ss)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = migration_average, inherit.aes = FALSE,
             aes(x=migration_binary, y = mean_ss), size = 3.6) + theme(axis.text.y = element_blank())

tro_plot +
  geom_errorbar(data = trophic_average, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = mean_ss - (1.96*error_ss), 
                    ymax = mean_ss + (1.96*error_ss)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = trophic_average, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = mean_ss), size = 3.6)

terr_plot  +
  geom_errorbar(data = territory_average, inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = mean_ss - (1.96*error_ss), 
                    ymax = mean_ss + (1.96*error_ss)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = territory_average, inherit.aes = FALSE,
             aes(x=territoriality_binary, y = mean_ss), size = 3.6) 


################################################################################
                    #### Overall averages ######


mig_plot <- mig_plot + ylim(c(0,2)) +
  geom_errorbar(data = mig_se, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = sexual_mean - (1.96*sexual_se), 
                    ymax = sexual_mean + (1.96*sexual_se)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_se, inherit.aes = FALSE,
             aes(x=migration_binary, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())

tro_plot <- tro_plot + ylim(c(0,2)) +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sexual_mean - (1.96*sexual_se), 
                    ymax = sexual_mean + (1.96*sexual_se)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())

terr_plot <- terr_plot + ylim(c(0,2)) + 
  geom_errorbar(data = terr_se, inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = sexual_mean - (1.96*sexual_se), 
                    ymax = sexual_mean + (1.96*sexual_se)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_se, inherit.aes = FALSE,
             aes(x=territoriality_binary, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())



################################################################################
                #### sqrt Overall averages ######



# Migration plot.
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=family_ss)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = migration_average, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = mean_ss - (1.96*error_ss), 
                    ymax = mean_ss + (1.96*error_ss)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = migration_average, inherit.aes = FALSE,
             aes(x=migration_binary, y = mean_ss), size = 3.6) + theme(axis.text.y = element_blank())

# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=family_ss)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = trophic_average, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = mean_ss - (1.96*error_ss), 
                    ymax = mean_ss + (1.96*error_ss)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = trophic_average, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = mean_ss), size = 3.6)

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=family_ss)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = territory_average, inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = mean_ss - (1.96*error_ss), 
                    ymax = mean_ss + (1.96*error_ss)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = territory_average, inherit.aes = FALSE,
             aes(x=territoriality_binary, y = mean_ss), size = 3.6) 










# tropic_se <- model_data %>% group_by(trophic_level_binary) %>% 
#   summarise(sexual_mean = mean(sexual_selection),
#             sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))
# 
# terr_se <- model_data %>% group_by(territoriality_binary) %>% 
#   summarise(sexual_mean = mean(sexual_selection),
#             sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))
# 
# mig_se <- model_data %>% group_by(migration_binary) %>% 
#   summarise(sexual_mean = mean(sexual_selection),
#             sexual_se = sd(sexual_selection)/sqrt(length(sexual_selection)))




################################################################################
                 #### Create temperature side plot ####


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

temp_uncen_all_data[[2]]
temp_linear_predictions[[1]]
temp_linear_predictions <- conditional_effects(temp_uncen_all_data[[1]], method = "posterior_linpred")

temp_posterior_predictions <- conditional_effects(temp_uncen_all_data[[1]], method = "posterior_predict")
?conditional_effects

temp_linear_predictions[[1]]
?posterior_epred()

pnorm(temp_lin_preds$estimate__)

post

# started at 1.12 pm.

temp_lin_preds <- temp_linear_predictions[[1]]



# Change predictions back to same scale as original sex scores.
temp_lin_preds %<>% 
  mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1) %>% 
  filter(temp_log > 3.3 & temp_log < 7.2)

# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + ylim(c(0,2)) +
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(estimate__-(1.96*se__)), ymax = sqrt(estimate__+(1.96*se__))), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())



# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + ylim(c(0,2)) +
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(estimate__-(1.96*se__)), ymax = sqrt(estimate__+(1.96*se__))), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1) + theme(axis.text.y = element_blank())

# Make some nice plots.
temp_plot <- family_average_data_data %>% ggplot(aes(x=temp_log, y=sqrt(family_ss))) + point_settings +
  xlab("Seasonality") + ylim(c(0,2)) +
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = sqrt(lower__), ymax = sqrt(upper__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_smooth(data = temp_plot_predictions, aes(y = sqrt(estimate__)), linetype = "dashed", linewidth = 1, method = "lm", stat = "summary") + theme(axis.text.y = element_blank())

temp_plot_predictions$lower__[temp_plot_predictions$lower__ < 0] <- 0

?geom_smooth

family_average_data_data %>% ggplot(aes(x=temp_log, y=family_ss)) + point_settings +
  xlab("Seasonality") + ylim(c(0,2)) +
  geom_ribbon(data = temp_plot_predictions, inherit.aes = FALSE,
              aes(x = temp_log, ymin = estimate__-(1.96*se__), ymax = estimate__+(1.96*se__)), fill = "grey70", colour = NA, alpha = 0.2) +
  geom_line(data = temp_plot_predictions, aes(y = estimate__), linetype = "dashed", linewidth = 1) #+ theme(axis.text.y = element_blank())






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
temp_plot_predictions <- temp_uncen_all_data[[2]] %>% filter(temp_log > 3.3 & temp_log < 7.2)
# temp_plot_predictions <- temp_plot_predictions %>% filter(temp_log < 7.2)

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
categorical_settings <- list(geom_point(aes(group=family_bird_tree, colour = higher_clade, 
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
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=family_ss)) + 
  categorical_settings + #ylim(c(0,2)) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") + 
  geom_errorbar(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=migration_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=migration_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())


migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_pointrange(
    stat = "summary",
    fun.min = standard_error,
    fun.max = standard_error,
    fun = median
  )


mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_se, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = sexual_mean - sexual_se*1.96, 
                    ymax = sexual_mean + sexual_se*1.96),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data =mig_se, inherit.aes = FALSE,
             aes(x=migration_binary, y = sexual_mean), size = 3.6) + theme(axis.text.y = element_blank())



# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)


tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_binary, ymin = estimate__-(se__*1.96), 
                    ymax = estimate__+(se__*1.96)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)


trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sqrt(sexual_mean))) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_pred_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sqrt(pred_mean - se_mean*1.96), 
                    ymax = sqrt(pred_mean + se_mean*1.96)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_pred_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sqrt(pred_mean)), size = 3.6)


trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sexual_mean_2 - sexual_se*2, 
                    ymax = sexual_mean_2 + sexual_se*2),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sexual_mean_2), size = 3.6)




trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sqrt(sexual_mean))) + 
  categorical_settings + ylim(c(0,2)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = sqrt(sexual_mean_2 - sexual_se*2), 
                    ymax = sqrt(sexual_mean_2 + sexual_se*2)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = sqrt(sexual_mean_2)), size = 3.6)



trophic_data %>% ggplot(aes(x=trophic_level_binary, y=log(sexual_mean+1))) + 
  categorical_settings + ylim(c(0,1.5)) +
  scale_x_discrete(limits = c("Primary", "Secondary"), labels = c("Primary", "Secondary")) +
  xlab("Trophic level") +
  geom_errorbar(data = tropic_se, inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = log((sexual_mean_2 - sexual_se)+1), 
                    ymax = log((sexual_mean_2 + sexual_se)+1)),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tropic_se, inherit.aes = FALSE,
             aes(x=trophic_level_binary, y = log(sexual_mean_2+1)), size = 3.6)


# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=territoriality_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())

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
categorical_settings <- list(geom_point(aes(group=family_bird_tree, colour = higher_clade, 
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
mig_plot <- migration_data %>% ggplot(aes(x=migration_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") +
  geom_errorbar(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=migration_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = mig_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=migration_binary, y = estimate__), size = 3.6) + theme(axis.text.y = element_blank())

# Trophic plots.
tro_plot <- trophic_data %>% ggplot(aes(x=trophic_level_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(limits = c("Primary", "Secondary"), 
                   labels = c(expression("1"^ry*""), 
                              expression("2"^ry*""))) +
  xlab("Trophic level") +
  geom_errorbar(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=trophic_level_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = tro_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=trophic_binary, y = estimate__), size = 3.6)

# Territory plots.
terr_plot <- teritory_data %>% ggplot(aes(x=territoriality_binary, y=sexual_mean)) + 
  categorical_settings +
  scale_x_discrete(labels = c("No", "Yes")) +
  xlab("Territoriality") +
  geom_errorbar(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
                aes(x=territoriality_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = terr_uncen_all_data[[2]], inherit.aes = FALSE,
             aes(x=territoriality_binary, y = estimate__), size = 3.6) 


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
###################################################################################






brms_univar_plot <- function(data_set, ylabel = "", ylimits = c(0,1.1), ybreaks = c(0,0.5,1), 
                             lab_x_pos = 2.2, lab_ypos = 2.2, plot_label = "b", 
                             plot_model = allbirds_model, sex_score = TRUE,
                             r_include = TRUE){
  
  # Extract predictions from brms model.
  predictions <- conditional_effects(plot_model)[[1]]
  
  # Change sexual score back to original values (have to + 1 for ordinal regression)
  if (sex_score){
    predictions$estimate__ <- predictions$estimate__ - 1
    predictions$lower__ <- predictions$lower__ - 1
    predictions$upper__ <- predictions$upper__ - 1
  }
  
  # Extract p-values.
  p_values <- brms_pmap(plot_model)[2]
  p_value <- last(p_values[,1])
  
  # Extract r-squared.
  r_squared <- Bayes_R2_MZ(plot_model)[[1]]
  
  # Sample size.
  sample_size <- nrow(plot_model$data)
  
  # Estimate
  estimate <- summary(plot_model)$fixed[5,1]
  estimate <- as.character(format(round(estimate, 2), nsmall = 2))
  
  # Redo estimate if it's too small.
  if (estimate == "0.00"){
    estimate <- summary(plot_model)$fixed[5,1]
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


migration_data %>% ggplot(aes(x = migration_binary, y = trait_mean)) +
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

# Settings for binary traits.
plot_settings <- list(geom_point(aes(group=family, colour = higher_clade, 
                                     size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
                                 position = position_jitter(0.2)),
                      labs(x = "", y = NULL,  colour = "Clade"),
                      scale_colour_manual(values = clade_colours),
                      scale_y_continuous(breaks = c(0, 1, 2, 3)),
                      theme_classic(base_size = 20),
                      theme(text = element_text(face = "bold"),
                            legend.position = "none",
                            axis.line = element_line(size = 1),
                            axis.ticks = element_line(size = 1),
                            axis.title.x = element_text(size = rel(0.8)),
                            axis.title.y = element_text(size = rel(0.8))))
migration_data %>% ggplot(aes(x=migration_binary, y=log(sexual_mean + 1))) + 
  geom_point(aes(group=family, colour = higher_clade, 
                 size = sqrt(clade_size), alpha = sqrt(clade_sum)), 
             position = position_jitter(0.2)) +
  labs(x = "", y = NULL,  colour = "Clade") +
  scale_colour_manual(values = clade_colours) +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  theme_classic(base_size = 20) +
  theme(text = element_text(face = "bold") +
          legend.position = "none",
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title.x = element_text(size = rel(0.8)),
        axis.title.y = element_text(size = rel(0.8))) +
  scale_x_discrete(limits = c("Weak", "Strong"), labels = c("No", "Yes")) +
  xlab("Migration") + theme(axis.text.y = element_blank()) +
  geom_errorbar(data = predictions, inherit.aes = FALSE,
                aes(x=migration_binary, ymin = lower__, 
                    ymax = upper__),
                position = position_dodge(width = 1), show.legend = FALSE, width = 0.2, linewidth = 1) +
  geom_point(data = predictions, inherit.aes = FALSE,
             aes(x=migration_binary, y = estimate__), size = 3.6) +
  annotate("text", x = lab_x_pos, y =lab_ypos, label = cor_label, size = 7, fontface = 2) 

data_set <- migration_data
