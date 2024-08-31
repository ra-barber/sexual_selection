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
             #### Read in the model data for plotting #####


# First half of the pathway.
first_half <- "Z:/home/sexual_selection/Results/Combined_models/Univariate/"

# Read unscaled models to plot against untransformed data.
temp_uncen_all_data <- readRDS(paste0(first_half, "temp_data.rds"))
mig_uncen_all_data <- readRDS(paste0(first_half, "mig_data.rds"))
tro_uncen_all_data <- readRDS(paste0(first_half, "tro_data.rds"))
terr_uncen_all_data <- readRDS(paste0(first_half, "terr_data.rds"))

# Read in scaled seasonality effect size to compare against multivariate models.
temp_cen_all_data <- readRDS(paste0(first_half, "temp_scaled_data.rds"))


###############################################################################
                 #### Prepare family level data #####


# Read in model data.
model_data <- read_ss_data()
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


# Create dataframe for extracting predictions.
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
                         ndraws = draw_number)
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
temp_err <- apply(temp_preds, 2, standard_error)
temp_summaries <- data.frame(temp_log = c(3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7),
                             estimate = temp_est,
                             error = temp_err,
                             lower = temp_est - (1.96*temp_err),
                             upper = temp_est + (1.96*temp_err))

# Extract summaries of categorical variables.
extract_summaries <- function(matrix, levels =  c("Reference", "Main")){
  ref_mean <- mean(matrix[,1])
  main_mean <- mean(matrix[,2])
  
  ref_se <- standard_error(matrix[,1])
  main_se <- standard_error(matrix[,2])
  
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


# Function for plotting models.
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