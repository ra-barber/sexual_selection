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

all_models <- readRDS("Z:/home/sexual_selection/Results/Models/Old_models/Combined_models/Multivariate/multivariate_all_models.rds")
high_models <- readRDS("Z:/home/sexual_selection/Results/Models/Old_models/Combined_models/Multivariate/multivariate_high_models.rds")

# Look at the summaries.
summary(all_models)
summary(high_models)

gc()

###############################################################################
                    #### Model Diagnostics #####


# Look at the summaries.
summary(all_models)
summary(high_models)

# # # Plot the model.
# plot(all_models)
# plot(high_models)


###############################################################################
                       #### Forest Plots #####

# 
# # Quick function to avoid repeating lines.
# brms_forest <- function(model){
#   mcmc_areas(model, regex_pars = "^b_[a-z]", prob = 0.95, prob_outer = 0.95)
# }

# # Feb models.
# brms_forest(all_models)
# brms_forest(high_models)

# ###############################################################################
#                           #### P mapping #####
# 
# 
# # Check p-values.
# brms_pmap <- function(model){
#   p_map(model, regex_pars = "^b_[a-z]")
# }
# 
# # Feb models.
# brms_pmap(all_models)
# brms_pmap(high_models)


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
            #### Extract the draws for each model #####

# Predictor names.
predictor_names <- c("Territoriality", "Migration", "1ry consumer",
                     "Seasonality", "1ry consumer\nx seasonality",
                     "1ry consumer\nx territoriality")

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

# Simple models first.
all_plotdata <- extract_draws(all_models, predictor_names)
high_plotdata <- extract_draws(high_models, predictor_names)


###############################################################################
                #### Add axis label and legend information #####

# Change order for plotting.
predictor_order <- c("1ry consumer", "Migration", "Territoriality","Seasonality",
                     "1ry consumer\nx territoriality", "1ry consumer\nx seasonality")

# Vector of predictors in reverse order of plotting.
all_plotdata$name %<>% factor(levels = rev(predictor_order))
high_plotdata$name %<>% factor(levels = rev(predictor_order))

# Add info on models.
lifehistory <- c("Territoriality",  "Migration", "1ry consumer")
resource <-   c("Seasonality")
interaction <-  c("1ry consumer\nx seasonality",
                  "1ry consumer\nx territoriality")


# Function for adding colour to legends.
add_legend_info <- function(plotdata){
  plotdata$trait <- NA
  plotdata$trait[plotdata$name %in% lifehistory] <- "Life history"
  plotdata$trait[plotdata$name %in% resource] <- "Environment"
  plotdata$trait[plotdata$name %in% interaction] <- "Diet interaction"
  plotdata
}

# Add the legend info.
all_plotdata %<>% add_legend_info()
high_plotdata %<>% add_legend_info()



###############################################################################
                         #### Colour pal #####


library(ggnewscale)

# Dark and light colours (not all used anymore)
#light_colours <- c( "#77AD78","#7494EA", "#C98986", "#D8C77B")   #F7D460   #F6CE5F
#dark_colours <- c("#214F4B", "#05299E", "#8D0801","#A88A05")
#light_colours <- c("#7494EA",  "#77AD78", "#C98986", "#D8C77B")   #F7D460   #F6CE5F

#light_colours <- c( "#C98986",  "#77AD78", "#7494EA", "#D8C77B")   #F7D460   #F6CE5F
light_colours <- c("#77AD78", "#C98986", "#7494EA", "#D8C77B") 
#dark_colours <- c( "#05299E","#214F4B", "#8D0801","#A88A05")

# Manual order for legend.
legend_order <- c("Life history","Environment","Diet interaction", "Devlopment interaction")

# Fancy y axis labels to have superscript
y_axis_labs <- rev(expression("1"^ry*" consumer", "Migration", 
                              "Territoriality", "Seasonality",
                              atop(NA,atop(textstyle("1"^ry*" consumer"),textstyle("× territoriality"))),
                              atop(NA,atop(textstyle("1"^ry*" consumer"), textstyle("× seasonality")))))


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
    labs(x = "Standardized effect size", y = NULL, fill = NULL, colour = NULL) +
    guides(colour = "none", size = "none", fill = guide_legend(byrow = TRUE, override.aes = list(alpha = 1))) +
    theme_classic2(base_size = 20) +  # Theme.
    theme(legend.position = legend_pos,
          axis.text.x=element_text(size=rel(1), colour = "black"), 
          axis.text.y=element_text(size=rel(0.9), colour = "black"),
          axis.title.x=element_text(size=rel(0.8)),
          line = element_line(linewidth = 0.5),
          legend.text = element_text(size=rel(0.7)), 
          legend.spacing.y = unit(0.1, 'cm'),
          legend.title = element_blank(),
          strip.text.y = element_blank())
}

# Create arrow data for gradient.
arrow_data <- data.frame(x = seq(from = -3.5, to = 2.5, by = 0.01), xend = 3, y = 1, yend = 1)

# Create arrow plot with sexual selection label.
arrow_plot <- ggplot(arrow_data) + geom_segment(mapping = aes(x = x, xend = 3, y = 1, yend = 1, colour = after_stat(x)),
                                                lineend = "butt", linejoin = "mitre",
                                                size = 3, arrow = arrow(length = unit(0.1, "inches"), type = "closed")) + 
  annotate("text", label = "Increasing sexual selection", x = 0, y = 1, vjust = -1, size = 6) +
  scale_colour_gradient(high = "black", low = "white") + theme_void() + 
  theme(legend.position = "none", plot.margin = margin(b = -40),
        plot.background=element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = 'white', colour = NA))

# Export the plots.
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

