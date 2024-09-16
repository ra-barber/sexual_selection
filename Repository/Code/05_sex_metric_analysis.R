###############################################################################
                    #### Sex metric model assessment ####
###############################################################################

# This script summarises models comparing sexual selection scores against
# alternative metrics commonly used to estimate sexual selection.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(caper)
library(dplyr)
library(effectsize)
library(janitor)
library(ggplot2)
library(ggpubr)
library(phytools)
library(brms)
library(readxl)

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                          #### Prepare data ####


# Read in the sex metric data.
pathway <- "Data/supplementary_dataset_1.xlsx"
metric_data <- read_excel(pathway, sheet = 4, na = "NA") %>% clean_names()
birdtree_data <-  read_excel(pathway, sheet = 2, na = "NA") %>% clean_names()
model_data <- left_join(metric_data, birdtree_data)
model_data$tree_tip <- gsub(" ", "_", model_data$scientific_name_bird_tree)

# Select testes mass.
testes_data <- model_data %>% 
  dplyr::select(scientific_name_bird_tree, tree_tip, sexual_selection, data_certainty, 
                sex_role_reversal, residual_testes_mass) %>% na.omit()
testes_data %<>% filter(data_certainty > 2)

# Select Bateman gradient metrics.
bateman_data <- model_data %>% 
  dplyr::select(scientific_name_bird_tree, tree_tip, sexual_selection, data_certainty, 
                bateman_gradient, bateman_individuals) %>% na.omit()

# Select OSS metrics.
oss_data <- model_data %>% 
  dplyr::select(scientific_name_bird_tree, tree_tip, sexual_selection, data_certainty, 
                sex_role_reversal, opportunity_for_sexual_selection, 
                oss_individuals, oss_estimated) %>% na.omit()


###############################################################################
                      #### Predictors ####


# Scale continuous predictors to two SD.
testes_data %<>% mutate(
  testes_z = standardize(residual_testes_mass, two_sd = TRUE))

bateman_data %<>% mutate(
  bateman_z = standardize(bateman_gradient, two_sd = TRUE))

oss_data %<>% mutate(
  oss_z = standardize(opportunity_for_sexual_selection, two_sd = TRUE),
  oss_log = log(opportunity_for_sexual_selection + 1),
  oss_log_z = standardize(oss_log, two_sd = TRUE))

# Make a sexual selection factor variable for plotting colours.
testes_data$sexual_selection_plot <- testes_data$sexual_selection %>% as.factor()
bateman_data$sexual_selection_plot <- bateman_data$sexual_selection %>% as.factor()
oss_data$sexual_selection_plot <- oss_data$sexual_selection %>% as.factor()


###############################################################################
                          #### brms models ####


# Pathway.
first_half <- "Z:/home/sexual_selection/Results/Models/Metrics/"

# Read models. (Centered for comparing stats. Uncentered for plotting)
cen_bateman_model <- readRDS(paste0(first_half, "bateman_centered_phylo_models.rds"))
bateman_model <- readRDS(paste0(first_half, "bateman_uncentered_phylo_models.rds"))

cen_oss_model <- readRDS(paste0(first_half, "oss_centered_phylo_models.rds"))
oss_model <- readRDS(paste0(first_half, "oss_uncentered_phylo_models.rds"))
sens_oss_model  <- readRDS(paste0(first_half, "sensoss_centered_phylo_models.rds"))

cen_teste_model <- readRDS(paste0(first_half, "testes_centered_phylo_models.rds"))
teste_model <- readRDS(paste0(first_half, "testes_uncentered_phylo_models.rds"))


# Summaries.
summary(cen_bateman_model)
summary(cen_oss_model)
summary(cen_teste_model)

# Extract and round the estimate from each model.
extract_estimate <- function(model){
  estimate <- last(summary(model)$fixed[,1])
  as.character(format(round(estimate, 2), nsmall = 2))
}

# Extract estimate.
bateman_estimate <- extract_estimate(cen_bateman_model)
oss_estimate <- extract_estimate(cen_oss_model)
teste_estimate <- extract_estimate(cen_teste_model)

# Extract p value.
library(bayestestR)

phy_bateman_p_value <- pd_to_p(last(p_direction(cen_bateman_model)[,2])) %>% round(2)
phy_oss_p_value <- pd_to_p(last(p_direction(cen_oss_model)[,2])) 
phy_teste_p_value <- pd_to_p(last(p_direction(cen_teste_model)[,2]))

# Extract r squared.
phy_bateman_r2 <- round(Bayes_R2_MZ(cen_bateman_model)[[1]], 2)
phy_oss_r2 <- round(Bayes_R2_MZ(cen_oss_model)[[1]], 2)
phy_teste_r2 <- round(Bayes_R2_MZ(cen_teste_model)[[1]], 2)

phy_bateman_r2 <- as.character(format(round(marginal_R2_MZ(cen_bateman_model)[[1]], 2), nsmall = 2))
phy_oss_r2 <- as.character(format(round(marginal_R2_MZ(cen_oss_model)[[1]], 2), nsmall = 2))
phy_teste_r2 <- as.character(format(round(marginal_R2_MZ(cen_teste_model)[[1]], 2), nsmall = 2))

# Label without line breaks.
phy_bateman_label <-  paste0("\U03B2 = ", bateman_estimate, "   p < 0.001   R\u00b2 = ", phy_bateman_r2)
phy_oss_label <-  paste0("\U03B2 = ", oss_estimate, "   p < 0.001   R\u00b2 = ", phy_oss_r2)
phy_teste_label <-  paste0("\U03B2 = ", teste_estimate, "   p < 0.001   R\u00b2 = ", phy_teste_r2)


################################################################################
                 #### Export summary tables ####


# Extract relevant coefficient information.
bateman_estimates <- summary(cen_bateman_model)$fixed[5,c(1,3,4)]
oss_estimates <- summary(cen_oss_model)$fixed[5,c(1,3,4)]
teste_estimates <- summary(cen_teste_model)$fixed[5,c(1,3,4)]
sens_oss_estimates <- summary(sens_oss_model)$fixed[5,c(1,3,4)]

all_estimates <- rbind(teste_estimates, bateman_estimates, 
                       oss_estimates, sens_oss_estimates)

# Paste together values for reporting in a table.
all_estimates %<>% mutate(
  round_est = round(Estimate, 2),
  intervals = paste0(round(`l-95% CI`, 2), ", ", 
                     round(`u-95% CI`, 2)),
  est_intervals = paste0(round_est, ", ", intervals))

# Export the results.
write.csv(all_estimates, "Results/Tables/sex_metric_comparison.csv", row.names = TRUE)


################################################################################
               ##### Make the EDF figure version #####


# Get the position for the stats.
get_x_pos <- function(response){
  range <- (max(response) - min(response))*0.45
  max(response) - range
}

bateman_x <- get_x_pos(bateman_data$bateman_gradient)
oss_x <- get_x_pos(oss_data$oss_log)
teste_x <- get_x_pos(testes_data$residual_testes_mass)

# Get the position for the sample size label. 
get_x_pos <- function(response){
  range <- (max(response) - min(response))*0.2
  max(response) - range
}

bateman_sample_x <- get_x_pos(bateman_data$bateman_gradient)
oss_sample_x <- get_x_pos(oss_data$oss_log)
teste_sample_x <- get_x_pos(testes_data$residual_testes_mass)

# Create the labels with italic "n"
bateman_sample_n <-  expression(~italic(n)~'= 14 species')
oss_sample_n <-  expression(~italic(n)~'= 79 species')
teste_sample_n <-  expression(~italic(n)~'= 977 species')

# Get predictions.
bateman_preds <- conditional_effects(bateman_model)[[1]] 
oss_preds <- conditional_effects(oss_model)[[1]] 
teste_preds <- conditional_effects(teste_model)[[1]]


# Function for OSS and BG plots.
pred_plot <- function(response = "bateman_gradient", dataset = bateman_data, 
                             plot_label = phy_bateman_label, 
                             x_label = "Bateman gradient",
                             y_lab_pos = 6.1, point_size = "bateman_individuals", 
                             x_lab_pos = bateman_x, sample_x = bateman_sample_x, sample_label = bateman_sample_n,
                             predict_data = bateman_preds){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", y = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", x = response)) + 
    geom_jitter(aes_string(size = point_size), 
                alpha = 0.75, width = 0, height = 0.2) + 
    scale_size_continuous(range = c(5,8)) +
    theme_classic(base_size = 18) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    ylab("Sexual selection") + 
    xlab(x_label) +
    scale_y_discrete(limits = as.factor(0:4), expand = expansion(add = c(1,1.5))) +
    scale_x_continuous(expand = expansion(mult = c(0.1))) +
    annotate("text", y = y_lab_pos, x = x_lab_pos, label = plot_label, size = 5.25) +
    annotate("text", y = 0.3, x = sample_x, label = sample_label, size = 5.25) + 
    geom_smooth(data = predict_data, inherit.aes = FALSE, 
                aes_string(x = response, y = "estimate__"),
                linetype = "dashed", linewidth = 1, col = "black", se = FALSE)
}

# Function for teste plot.
pred_plot_teste <- function(response = "residual_testes_mass", dataset = testes_data, 
                                   plot_label = phy_teste_label, 
                                   x_label = "Residual testes mass",
                                   y_lab_pos = 6.1, 
                                   x_lab_pos = teste_x,
                                   predict_data = teste_preds){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", y = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", x = response)) + 
    geom_jitter(size = 3,
                alpha = 0.5, width = 0, height = 0.2) + 
    scale_size_continuous(range = c(5,8)) +
    theme_classic(base_size = 18) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    ylab("Sexual selection") + 
    xlab(x_label) +
    scale_y_discrete(limits = as.factor(0:4), expand = expansion(add = c(1,1.5))) +
    scale_x_continuous(expand = expansion(mult = c(0.1))) +
    annotate("text", y = y_lab_pos, x = x_lab_pos, label = plot_label, size = 5.25) +
    annotate("text", y = 0.3, x = teste_sample_x, label = teste_sample_n, size = 5.25) + 
    geom_smooth(data = predict_data, inherit.aes = FALSE, 
                aes_string(x = response, y = "estimate__"),
                linetype = "dashed", linewidth = 1, col = "black", se = FALSE)
}

# Colour pal.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

# Testes mass plot.
teste_plot <- pred_plot_teste()

# Opportunity for sexual selection plot.
oss_plot <- pred_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                      x_label = "Opportunity for sexual selection", 
                      point_size = "oss_individuals", predict_data = oss_preds, 
                      x_lab_pos = oss_x, sample_x = oss_sample_x,
                      sample_label = oss_sample_n) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

# Bateman plot.
pal <- c('#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
bateman_plot <- pred_plot() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

# Combine plots together.
ggarrange(teste_plot,  bateman_plot, oss_plot,
          nrow = 1, ncol = 3, widths = c(1.1,1,1), labels = c("a", "b", "c"),
          hjust = c(-3.5, -1.4, -1.4), font.label = list(size = 28))

# Export.
ggsave("Figures/Fig_2_base.tiff", width = 12, height = 4.5, compression = "lzw", )


################################################################################
                               #### End ####
################################################################################