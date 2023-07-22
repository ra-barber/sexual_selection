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

all_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Multivariate/multivariate_all_models.rds")
high_models <- readRDS("Z:/home/sexual_selection/Results/Models/Combined_models/Multivariate/multivariate_high_models.rds")

###############################################################################
                    #### Model Diagnostics #####


# Look at the summaries.
summary(all_models)
summary(high_models)

# # Plot the model.
plot(all_models)
plot(high_models)

# Check the posterior distribution.
pp_check(all_models)
ggsave("Plots/Results/all_multivariate_ppcheck.pdf", width = 8, height = 5.5, device = cairo_pdf)
ggsave("Plots/Results/all_multivariate_ppcheck.png", width = 8, height = 5.5)
pp_check(high_models)
ggsave("Plots/Results/high_multivariate_ppcheck.pdf", width = 8, height = 5.5, device = cairo_pdf)
ggsave("Plots/Results/high_multivariate_ppcheck.png", width = 8, height = 5.5)

###############################################################################
                       #### Forest Plots #####


# Quick function to avoid repeating lines.
brms_forest <- function(model){
  mcmc_areas(model, regex_pars = "^b_[a-z]", prob = 0.95, prob_outer = 0.99)
}

# Feb models.
brms_forest(all_models)
brms_forest(high_models)

#brms_forest(log_high_models)
#brms_forest(test_model)
###############################################################################
                          #### P mapping #####


# Check p-values.
brms_pmap <- function(model){
  p_map(model, regex_pars = "^b_[a-z]")
}

# Feb models.
brms_pmap(all_models)
brms_pmap(high_models)

###############################################################################
           #### Lambda values (Shouldn't work for this family) #####
# 
# 
# # Look at lambda values.
# hyp <- "sd_tree_tip__Intercept^2 / (sd_tree_tip__Intercept^2 + disc^2) = 0"
# 
# (all_lambda <- hypothesis(all_models, hyp, class = NULL))
# (high_models <- hypothesis(high_models, hyp, class = NULL))
# plot(all_lambda)


###############################################################################
                         #### R2 values #####
gc()
# library(tictoc)
# 
# # This takes forever and smashes the ram.
# tic()
# all_r_squared <- bayes_R2(all_models, categorical = TRUE)
# toc()
# 
# tic()
# high_r_squared <- bayes_R2(high_models, categorical = TRUE)
# toc()


################################################################################
                ##### Export summary tables ####

# Extract relevant coeffcient information.
all_estimates <- summary(all_models)$fixed[5:10,c(1,3,4)]
high_estimates <- summary(high_models)$fixed[5:10,c(1,3,4)]

# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, " ", intervals))

high_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ",
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, " ", intervals))

# Reorder.
order <- c("trophic_level_c", "migration_bi_c","terr_bi_c","temp_seasonality_z","terr_bi_c:trophic_level_c","trophic_level_c:temp_seasonality_z") 

all_estimates <- all_estimates[order,]
high_estimates <- high_estimates[order,]

# Export the results.
write.csv(all_estimates, "Results/Tables/all_multivariate_regression.csv", row.names = TRUE)
write.csv(high_estimates, "Results/Tables/high_multivariate_regression.csv", row.names = TRUE)


###############################################################################
                  #### Making the plot names #####

# Predictor names.
predictor_names <- c("Territoriality", "Migration", "Primary\nconsumer",
                  "Seasonality", "Primary\nconsumer\n x \nSeasonality",
                  "Primary\nconsumer\n x \nTerritoriality")
predictor_names <- c("Territoriality", "Migration", "1ry consumer",
                     "Seasonality", "1ry consumer\nx seasonality",
                     "1ry consumer\nx territoriality")


###############################################################################
            #### Extract the draws for each model #####

# Extract draws function.
extract_draws <- function(model, column_names){
  # Pull out the draws.
  model_draws <- as_draws_df(model, c("^b"), regex = TRUE)[,-c(1:4)] %>% 
    as.data.frame() %>% dplyr::select(-c(.chain, .iteration, .draw))
  # Change column names.
  colnames(model_draws) <- column_names
  # Pivot longer.
  model_draws %>% tidyr::pivot_longer(cols = column_names)
}

# # Extract draws function.
# log_extract_draws <- function(model, column_names){
#   # Pull out the draws.
#   model_draws <- as_draws_df(model, c("^b"), regex = TRUE)[,-c(1)] %>% 
#     as.data.frame() %>% dplyr::select(-c(.chain, .iteration, .draw))
#   # Change column names.
#   colnames(model_draws) <- column_names
#   # Pivot longer.
#   model_draws %>% tidyr::pivot_longer(cols = column_names)
# }

# Simple models first.
all_plotdata <- extract_draws(all_models, predictor_names)
high_plotdata <- extract_draws(high_models, predictor_names)
#log_high_plotdata <- log_extract_draws(log_high_models, predictor_names)



###############################################################################
                #### Add axis label and legend information #####

# Change order for plotting.
predictor_order <- c("Seasonality", "Migration", "Primary\nconsumer", "Territoriality",
                     "Primary\nconsumer\n x \nSeasonality",
                     "Primary\nconsumer\n x \nTerritoriality")
predictor_order <- c("Seasonality", "Migration", "1ry consumer", "Territoriality",
                     "1ry consumer\nx seasonality",
                     "1ry consumer\nx territoriality")

predictor_order <- c("1ry consumer", "Migration", "Territoriality","Seasonality",
                     "1ry consumer\nx territoriality", "1ry consumer\nx seasonality")

# predictor_order <- c("Seasonality", "Migration", "expression(bold("1"^ry*"")) consumer", "Territoriality",
#                      "expression(bold("1"^ry*"")) consumer\nx seasonality",
#                      "expression(bold("1"^ry*"")) consumer\nx territoriality")


# Vector of predictors in reverse order of plotting.
all_plotdata$name %<>% factor(levels = rev(predictor_order))
high_plotdata$name %<>% factor(levels = rev(predictor_order))
#log_high_plotdata$name %<>% factor(levels = rev(predictor_order))

# Add info on models.
lifehistory <- c("Territoriality",  "Migration",  "Primary\nconsumer", "1ry consumer",
                 "(Pr) Precocial", "Chick maturity")
resource <-   c("Seasonality", "Productivity")
interaction <-  c("Primary\nconsumer\n x \nTerritoriality", 
                  "Primary\nconsumer\n x \nSeasonality",
                  "Primary\nconsumer\n x \nProductivity",
                  "1ry consumer\nx seasonality",
                  "1ry consumer\nx territoriality")

devo_interaction <- c("DM x PC", "DM x \ngeneration length", "DM x Migration", "DM x \nseasonality",
                      "DM x \nproductivity", "DM x \nterritoriality","CM x \nSeasonality",
                      "CM x \nProductivity", "CM x \nTerritoriality")

# Function for adding colour to legends.
add_legend_info <- function(plotdata){
  plotdata$trait <- NA
  plotdata$trait[plotdata$name %in% lifehistory] <- "Life history"
  plotdata$trait[plotdata$name %in% resource] <- "Environment"
  plotdata$trait[plotdata$name %in% interaction] <- "Diet interaction"
  plotdata$trait[plotdata$name %in% devo_interaction] <- "Devlopment interaction"
  plotdata
}

# Add the legend info.
all_plotdata %<>% add_legend_info()
high_plotdata %<>% add_legend_info()
#log_high_plotdata %<>% add_legend_info()

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

light_colours <- c( "#C98986",  "#77AD78", "#7494EA", "#D8C77B")   #F7D460   #F6CE5F
light_colours <- c("#77AD78", "#C98986", "#7494EA", "#D8C77B") 
dark_colours <- c( "#05299E","#214F4B", "#8D0801","#A88A05")

legend_order <- c(  "Life history","Environment","Diet interaction", "Devlopment interaction")

y_axis_labs <- rev(expression("Territoriality", "Migration", "1"^ry*" consumer",
                              "Seasonality",
                              "1"^ry*" consumer\nx seasonality",
                              "1"^ry*" consumer\nx territoriality"))

y_axis_labs <- rev(expression("Seasonality", "Migration", "1"^ry*" consumer",
                     "Territoriality",
                     atop(NA,atop(textstyle("1"^ry*" consumer"),textstyle("× seasonality"))),
                     atop(NA,atop(textstyle("1"^ry*" consumer"), textstyle("× territoriality")))))

y_axis_labs <- rev(expression("1"^ry*" consumer", "Migration", 
                              "Territoriality", "Seasonality",
                              atop(NA,atop(textstyle("1"^ry*" consumer"),textstyle("× territoriality"))),
                              atop(NA,atop(textstyle("1"^ry*" consumer"), textstyle("× seasonality")))))

###############################################################################
                    #### Quasi random plots #####

gg_quasi_forst <- function(plotdata, pred_n = 6, legend_pos = c(0.15, 0.15), col_pal = light_colours, alpha_val = 0.1){
  # Main publication design for simple models.
  ggplot(plotdata, aes(x = value, y = name, colour = trait, fill = trait)) +
    geom_quasirandom(bandwidth = 0.5, 
                     size = 0.5, 
                     #shape=19,  # Can't figure out how to get filled in circles that look good.
                     stroke = 0.05,  # Probs gotta play around with stroke too.
                     groupOnX = FALSE, 
                     alpha = alpha_val, position = position_dodge(2)) +
    scale_colour_manual(values = col_pal, breaks = legend_order) +
    scale_fill_manual(values = col_pal, breaks = legend_order) +
    guides(fill = "none",
      colour = guide_legend(byrow = TRUE, 
                                 override.aes = list(alpha = 1, 
                                                     shape = 15, size = 4))) +
    scale_y_discrete(labels = y_axis_labs) +
    #new_scale_colour() +
    stat_pointinterval(size = rep(c(9, 2), times = pred_n), colour = "black") + 
    #guides(fill = "none", colour = guide_legend(byrow = TRUE)) +
    geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + 
    labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
    theme_classic2(base_size = 15) +  # Theme.
    theme(legend.position = legend_pos,
          axis.text.x=element_text(size=rel(1), face = "bold", colour = "black"), 
          axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
          axis.title.x=element_text(size=rel(0.7), face = "bold"),
          legend.text = element_text(size=rel(0.8), face = "bold"), 
          legend.key.height = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title = element_blank(),
          strip.text.y = element_blank())
}

# Export model figures.
gg_quasi_forst(all_plotdata, legend_pos = c(0.14, 0.9), alpha_val = 0.1, col_pal = light_colours) + xlim(c(-3.75,3))
ggsave("Plots/Results/all_model.tiff", dpi = 1200, width = 8, height = 4, compression = "lzw")

gg_quasi_forst(high_plotdata, legend_pos = c(0.15, 0.125), alpha_val = 0.2)
ggsave("Plots/Results/high_model.tiff", dpi = 600, width = 8, height = 4)


###############################################################################
                           #### Eye plots #####


library(ggdist)
gg_stat_eye <- function(plotdata, pred_n = 6, legend_pos = c(0.15, 0.15), col_pal = light_colours, alpha_val = 0.1){
ggplot(plotdata, aes(x = value, y = name, fill = trait, height =  0.6)) +
    stat_slab(alpha = alpha_val, side = "both", normalize = "none", fill_type = "gradient") +
  stat_pointinterval(size = rep(c(12, 4), times = pred_n), colour = "black", show.legend = FALSE) + 
  scale_colour_manual(values = col_pal, breaks = legend_order) +
  scale_fill_manual(values = col_pal, breaks = legend_order) +
    scale_y_discrete(labels = y_axis_labs) +
  geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + 
  labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
    guides(colour = "none", size = "none", fill = guide_legend(byrow = TRUE, override.aes = list(alpha = 1))) +
  theme_classic2(base_size = 20) +  # Theme.
  theme(legend.position = legend_pos,
        axis.text.x=element_text(size=rel(1), face = "bold", colour = "black"), 
        axis.text.y=element_text(size=rel(0.9), face = "bold", colour = "black"),
        axis.title.x=element_text(size=rel(0.8), face = "bold"),
        legend.text = element_text(size=rel(0.7), face = "bold"), 
        #legend.key.height = unit(0.2, 'cm'),
        #legend.key.size = unit(0.2, 'cm'),
        #legend.key.width = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.title = element_blank(),
        strip.text.y = element_blank())
}


all_plot <- gg_stat_eye(all_plotdata, legend_pos = c(0.25, 0.9), alpha_val = 0.5, col_pal = light_colours) + xlim(c(-3.5,3))
ggsave("Plots/Results/all_model.tiff", dpi = 1200, width = 8, height = 5, compression = "lzw")

high_plot <- gg_stat_eye(high_plotdata, legend_pos = c(0.25, 0.9), alpha_val = 0.5, col_pal = light_colours) + xlim(c(-4.5,3))
ggsave("Plots/Results/high_model.tiff", dpi = 1200, width = 8, height = 5, compression = "lzw")

#gg_stat_eye(log_high_plotdata, legend_pos = c(0.25, 0.9), alpha_val = 0.5, col_pal = light_colours)
#ggsave("Plots/Results/loghigh_model.tiff", dpi = 1200, width = 8, height = 4, compression = "lzw")


###############################################################################
                          #### Add arrows #####


# Create arrow data for gradient.
arrow_data <- data.frame(x = seq(from = -3.5, to = 2.5, by = 0.01), xend = 3, y = 1, yend = 1)

# Create arrow plot with sexual selection label.
arrow_plot <- ggplot(arrow_data) + geom_segment(mapping = aes(x = x, xend = 3, y = 1, yend = 1, colour = after_stat(x)),
                        lineend = "butt", linejoin = "mitre",
                        size = 3, arrow = arrow(length = unit(0.1, "inches"), type = "closed")) + 
  annotate("text", label = "Sexual selection", x = 0, y = 1, vjust = -1, size = 7.5, fontface = 2) +
  scale_colour_gradient(high = "black", low = "white") + theme_void() + 
  theme(legend.position = "none", plot.margin = margin(b = -40),
        plot.background=element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = 'white', colour = NA))


ggarrange(arrow_plot, all_plot, nrow = 2, heights = c(0.2,1), align = "v")
ggsave("Plots/Results/all_model_arrow.tiff", dpi = 1200, width = 8, height = 5.5, compression = "lzw")
ggsave("Plots/Results/all_model.pdf", width = 8, height = 5.5, device = cairo_pdf)

ggarrange(arrow_plot, high_plot, nrow = 2, heights = c(0.2,1), align = "v")
ggsave("Plots/Results/high_model.pdf", width = 8, height = 5.5, device = cairo_pdf)
ggsave("Plots/Results/high_model_arrow.tiff", dpi = 1200, width = 8, height = 5.5, compression = "lzw")



################################################################################
                     ##### Non-bold version ####



gg_stat_eye <- function(plotdata, pred_n = 6, legend_pos = c(0.15, 0.15), col_pal = light_colours, alpha_val = 0.1){
  ggplot(plotdata, aes(x = value, y = name, fill = trait, height =  0.6)) +
    stat_slab(alpha = alpha_val, side = "both", normalize = "none", fill_type = "gradient") +
    stat_pointinterval(size = rep(c(12, 4), times = pred_n), colour = "black", show.legend = FALSE) + 
    scale_colour_manual(values = col_pal, breaks = legend_order) +
    scale_fill_manual(values = col_pal, breaks = legend_order) +
    scale_y_discrete(labels = y_axis_labs) +
    geom_vline(xintercept=0, lty=2, size = 0.5, alpha = 0.6) + 
    labs(x = "Standardised effect size", y = NULL, fill = NULL, colour = NULL) +
    guides(colour = "none", size = "none", fill = guide_legend(byrow = TRUE, override.aes = list(alpha = 1))) +
    theme_classic2(base_size = 20) +  # Theme.
    theme(legend.position = legend_pos,
          axis.text.x=element_text(size=rel(1), colour = "black"), 
          axis.text.y=element_text(size=rel(0.9), colour = "black"),
          axis.title.x=element_text(size=rel(0.8)),
          #axis.line = element_line(linewidth = 0.5),
          line = element_line(linewidth = 0.5),
          legend.text = element_text(size=rel(0.7)), 
          #legend.key.height = unit(0.2, 'cm'),
          #legend.key.size = unit(0.2, 'cm'),
          #legend.key.width = unit(0.5, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title = element_blank(),
          strip.text.y = element_blank())
}

# Create arrow plot with sexual selection label.
arrow_plot <- ggplot(arrow_data) + geom_segment(mapping = aes(x = x, xend = 3, y = 1, yend = 1, colour = after_stat(x)),
                                                lineend = "butt", linejoin = "mitre",
                                                size = 3, arrow = arrow(length = unit(0.1, "inches"), type = "closed")) + 
  annotate("text", label = "Increasing sexual selection", x = 0, y = 1, vjust = -1, size = 6) +
  scale_colour_gradient(high = "black", low = "white") + theme_void() + 
  theme(legend.position = "none", plot.margin = margin(b = -40),
        plot.background=element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = 'white', colour = NA))


all_plot <- gg_stat_eye(all_plotdata, legend_pos = c(0.2, 0.9), alpha_val = 0.5, col_pal = light_colours) + xlim(c(-3.5,3))
ggarrange(arrow_plot, all_plot, nrow = 2, heights = c(0.2,1), align = "v")
ggsave("Plots/Results/figure_5.pdf", width = 8, height = 5.5, device = cairo_pdf)
ggsave("Plots/Results/figure_5.png", width = 8, height = 5.5)

high_plot <- gg_stat_eye(high_plotdata, legend_pos = c(0.2, 0.9), alpha_val = 0.5, col_pal = light_colours) + xlim(c(-4.5,3))
ggarrange(arrow_plot, high_plot, nrow = 2, heights = c(0.2,1), align = "v")
ggsave("Plots/Results/figure_EDF2.pdf", width = 8, height = 5.5, device = cairo_pdf)
ggsave("Plots/Results/figure_EDF2.png", width = 8, height = 5.5)


                           #### END #####
###############################################################################

