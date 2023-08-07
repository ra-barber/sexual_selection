###############################################################################
                    ##### Bateman gradient analysis #####
###############################################################################

# This script performs the supplementary analysis for testes mass with sexual 
# selection scores.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(tictoc)
library(caper)
library(dplyr)
library(effectsize)
library(janitor)
library(car)
library(ggplot2)
library(ggpubr)
library(phytools)
library(brms)

# Read in the functions. 
source("Code/functions.R")




###############################################################################
                  #### Prepare bateman data ####


# Read in some data.
avian_data <- read.csv("Data/bateman_gradients/tim_avian_data.csv") %>% clean_names()
avian_data$tree_tip <- gsub(" ", "_", avian_data$birdtree_name)

oss_data <- read.csv("Data/bateman_gradients/Tim_OSS_database.csv") %>% clean_names()
oss_data$tree_tip <- gsub(" ", "_", oss_data$birdtree_name)
oss_data %<>% filter(exclude != "YES")

# Select batemman graident metrics.
bateman_data <- avian_data %>% 
  dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, 
                n_f, beta_f, r_f, n_m, beta_m, r_m, beta_t, beta_g) %>% na.omit()
bateman_data[bateman_data$beta_f > bateman_data$beta_m, "birdtree_name"]

oss_data <- oss_data %>% 
  dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, sex_role_reversal,
                n_f, is_f, n_m, is_m)

# Create a single column of OSS, swapping over Jacana jacana as a role reversed species.
oss_data$oss <- oss_data$is_m
oss_data$n_oss <- oss_data$n_m

# Species to take female OSS from.
oss_data[oss_data$is_f > oss_data$is_m,]
female_species <- c("Clamator glandarius",  "Jacana jacana", "Zonotrichia albicollis", "Zonotrichia leucophrys")

# Get data from females for these species.
oss_data$oss[oss_data$birdtree_name %in% female_species] <- oss_data$is_f[oss_data$birdtree_name %in% female_species]
oss_data$n_oss[oss_data$birdtree_name %in% female_species] <- oss_data$n_f[oss_data$birdtree_name %in% female_species]

# Change n to 572 as stated in the paper. 
oss_data$n_oss[oss_data$birdtree_name == "Dendroica caerulescens"] <- 572 # Number of males included in the Germain et al. 2021 study.

# Try removing Malurus for now (looks like OSS was estimated incorrectly (altho could be correct))
oss_data %<>% filter(birdtree_name != "Malurus cyaneus")

# Remove extra columns and NA species.
oss_data %<>% 
  dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, sex_role_reversal, oss, n_oss) %>% na.omit()
 

zero_mating_variance <- oss_data %>% filter(oss == 0)

bateman_data
zero_mating_variance %>% select(birdtree_name, tree_tip, sexual_selection, sex_role_reversal, data_certainty, oss, n_oss)

# 
# # Select opportunity metrics.
# oss_data <- avian_data %>% 
#   dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, 
#                 n_f, is_f, i_f, n_m, is_m, i_m) %>% na.omit()

# Weight metrics by sample size for combining studies.
bateman_data$weighted_beta_f <- bateman_data$beta_f * bateman_data$n_f   # Bateman gradient.
bateman_data$weighted_beta_m <- bateman_data$beta_m * bateman_data$n_m

bateman_data$weighted_corr_f <- bateman_data$r_f * bateman_data$n_f   # Bateman gradient.
bateman_data$weighted_corr_m <- bateman_data$r_m * bateman_data$n_m

oss_data$weighted_oss <- oss_data$oss * oss_data$n_oss  

# oss_data$weighted_oss_f <- oss_data$is_f * oss_data$n_f  # Opportunity for sexual selection.
# oss_data$weighted_oss_m <- oss_data$is_m * oss_data$n_m
# 
# oss_data$weighted_selection_f <- oss_data$i_f * oss_data$n_f # Opportunity for selection.
# oss_data$weighted_selection_m <- oss_data$i_m * oss_data$n_m

# Combine multiple studies.
bateman_data %<>% group_by(birdtree_name) %>% 
  summarise(tree_tip = first(tree_tip),
            sexual_selection = first(sexual_selection),
            data_certainty = first(data_certainty),
            
            beta_f = mean(beta_f),
            weight_beta_f = sum(weighted_beta_f)/ sum(n_f),
            beta_sqrt_f = sqrt(weight_beta_f + 0.06),    # sqrt females and add minimum size.
            
            beta_m = mean(beta_m),
            weight_beta_m = sum(weighted_beta_m)/ sum(n_m),
            beta_log_m = log(weight_beta_m),
            
            corr_f = mean(r_f),
            weight_corr_f = sum(weighted_corr_f)/ sum(n_f),
            
            corr_m = mean(r_m),
            weight_corr_m = sum(weighted_corr_m)/ sum(n_m),
            
            corr_mean = mean(weight_corr_f, weight_corr_m),
            
            beta_t = mean(beta_t),
            beta_g = mean(beta_g),
            
            n_f = sum(n_f),
            n_m = sum(n_m))

# # Combine multiple studies.
# oss_data %<>% group_by(birdtree_name) %>% 
#   summarise(tree_tip = first(tree_tip),
#             sexual_selection = first(sexual_selection),
#             data_certainty = first(data_certainty),
#             
#             oss_f = mean(is_f),
#             weight_oss_f = sum(weighted_oss_f)/ sum(n_f),
#             oss_log_f = log(weight_oss_f),
#             
#             oss_m = mean(is_m),
#             weight_oss_m = sum(weighted_oss_m)/ sum(n_m),
#             oss_log_m = log(weight_oss_m),
#             
#             selection_f = mean(i_f),
#             weight_selection_f = sum(weighted_selection_f)/ sum(n_f),
#             selection_log_f = log(weight_selection_f),
#             
#             selection_m = mean(i_m),
#             weight_selection_m = sum(weighted_selection_m)/ sum(n_m),
#             selection_log_m = log(weight_selection_m),
#             
#             n_f = sum(n_f),
#             n_m = sum(n_m))


oss_data %<>% group_by(birdtree_name) %>% 
  summarise(tree_tip = first(tree_tip),
            sexual_selection = first(sexual_selection),
            data_certainty = first(data_certainty),
            
            oss = mean(oss),
            weight_oss = sum(weighted_oss)/ sum(n_oss),
            oss_log = log(weight_oss + 1),
            n_oss = sum(n_oss))

# Make negative reciprocal if we need a straighter line.
oss_data$oss_recip <- (1/(oss_data$weight_oss + 1))*-1

hist(oss_data$oss_recip) # Still skewed so not worth the extra complexity of a weird transformation.
hist(oss_data$oss_log)

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]

# Drop tips on the tree.
bateman_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, bateman_data$tree_tip))
oss_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, oss_data$tree_tip))

# Make a covariance matrix, and order data the same.
bateman_covar <- ape::vcv.phylo(bateman_tree)
oss_covar <- ape::vcv.phylo(oss_tree)


###############################################################################
                      #### Predictors ####


# Scale continuous predictors to two SD.
bateman_data %<>% mutate(
  beta_f_z = standardize(weight_beta_f, two_sd = TRUE),
  beta_sqrt_f_z = standardize(beta_sqrt_f, two_sd = TRUE),
  corr_f_z = standardize(weight_corr_f, two_sd = TRUE),
  beta_m_z = standardize(weight_beta_m, two_sd = TRUE),
  beta_log_m_z = standardize(beta_log_m, two_sd = TRUE),
  corr_m_z = standardize(weight_corr_m, two_sd = TRUE),
  beta_t_z = standardize(beta_t, two_sd = TRUE),
  beta_g_z = standardize(beta_g, two_sd = TRUE),
  corr_mean_z = standardise(corr_mean, two_sd = TRUE)
)

oss_data %<>% mutate(
  oss_z = standardize(weight_oss, two_sd = TRUE),
  oss_log_z = standardize(oss_log, two_sd = TRUE),
)


# Make sexual selection a factor for plotting colours.
bateman_data$sexual_selection_plot <- bateman_data$sexual_selection %>% as.factor()
oss_data$sexual_selection_plot <- oss_data$sexual_selection %>% as.factor()

# Prepare response variables.
bateman_data$sexual_selection <- bateman_data$sexual_selection + 1
oss_data$sexual_selection <- oss_data$sexual_selection + 1


################################################################################
                  ##### And in teste data #####

# Read in some data.
teste_data <- read.csv("Data/teste_mass/both_testes_datasets_01_08.csv") %>% clean_names()
teste_data %<>% filter(data_certainty > 2)

# Scale continuous predictors to two SD.
teste_data %<>% mutate(
  teste_z = standardize(residual_testes_mass, two_sd = TRUE)
)

# Prepare response variables.
teste_data$sexual_selection_plot <- teste_data$sexual_selection %>% as.factor()
teste_data$sexual_selection <- teste_data$sexual_selection + 1

# Read in phy model.
phy_teste_model <- readRDS("Z:/home/sexual_selection/Results/Models/Testes/teste_model_1.rds")

# Brms formula.
teste_formula <- brmsformula("sexual_selection ~ residual_testes_mass", family = cumulative())

# Run models.
teste_model <- brm(
  teste_formula, data = teste_data, 
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")




###############################################################################
                    #### Run brms models ######

# Brms formula.
bateman_male_formula <- brmsformula("sexual_selection ~ corr_m_z + (1|gr(tree_tip, cov=A))", family = cumulative())
oss_formula <- brmsformula("sexual_selection ~ oss_log_z + (1|gr(tree_tip, cov=A))", family = cumulative())

# Phy models.
phy_bateman_male_model <- brm(
  bateman_male_formula, data = bateman_data, data2 = list(A=bateman_covar),
  iter = 25000, warmup = 20000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")

phy_oss_model <- brm(
  oss_formula, data = oss_data, data2 = list(A=oss_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
  init = 0, normalize = FALSE, backend = "cmdstanr")

###############################################################################
                  #### Extract summary stats ######


# Summaries.
summary(phy_bateman_male_model)
summary(phy_oss_model)
summary(phy_teste_model)

# Extract estimate.
phy_bateman_male_estimate <- last(summary(phy_bateman_male_model)$fixed[,1])
phy_oss_estimate <- last(summary(phy_oss_model)$fixed[,1])
phy_teste_estimate <- last(summary(phy_teste_model)$fixed[,1])

phy_bateman_male_estimate <- as.character(format(round(phy_bateman_male_estimate, 2), nsmall = 2))
phy_oss_estimate <- as.character(format(round(phy_oss_estimate, 2), nsmall = 2))
phy_teste_estimate <- as.character(format(round(phy_teste_estimate, 2), nsmall = 2))


# Extract p value.
library(bayestestR)

phy_bateman_male_p_value <- pd_to_p(last(p_direction(phy_bateman_male_model)[,2])) %>% round(2)
phy_oss_p_value <- pd_to_p(last(p_direction(phy_oss_model)[,2])) 
phy_teste_p_value <- pd_to_p(last(p_direction(phy_teste_model)[,2]))

# Extract r squared.
phy_bateman_male_r2 <- round(Bayes_R2_MZ(phy_bateman_male_model)[[1]], 2)
phy_oss_r2 <- round(Bayes_R2_MZ(phy_oss_model)[[1]], 2)
phy_teste_r2 <- round(Bayes_R2_MZ(phy_teste_model)[[1]], 2)

# Label with line breaks.
phy_bateman_male_label <-  paste0("\U03B2 = ", phy_bateman_male_estimate, "\np = ", phy_bateman_male_p_value, "\nr\u00b2 = ", phy_bateman_male_r2)
phy_oss_label <-  paste0("\U03B2 = ", phy_oss_estimate, "\np < 0.01\nr\u00b2 = ", phy_oss_r2)

# Label without line breaks.
phy_bateman_male_label <-  paste0("\U03B2 = ", phy_bateman_male_estimate, "   p = ", phy_bateman_male_p_value, "   r\u00b2 = ", phy_bateman_male_r2)
phy_oss_label <-  paste0("\U03B2 = ", phy_oss_estimate, "   p < 0.001   r\u00b2 = ", phy_oss_r2)
phy_teste_label <-  paste0("\U03B2 = ", phy_teste_estimate, "   p < 0.001   r\u00b2 = ", phy_teste_r2)




################################################################################
                      ##### Non-phylo models ####


# Brms formula.
bateman_male_formula <- brmsformula("sexual_selection ~ weight_corr_m", family = cumulative())
oss_formula <- brmsformula("sexual_selection ~ oss_log", family = cumulative())

# Run models.
bateman_male_model <- brm(
  bateman_male_formula, data = bateman_data, 
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")

oss_model <- brm(
  oss_formula, data = oss_data,
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")


# Extract predictions from brms model.
bateman_male_preds <- conditional_effects(bateman_male_model)[[1]] 
oss_preds <- conditional_effects(oss_model)[[1]] 
teste_preds <- conditional_effects(teste_model)[[1]]

# Change back to normal data values that were transformed for ordinal models..
bateman_male_preds %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
oss_preds %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
teste_preds %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)


################################################################################
                        ##### Plotting #####


# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
pal <- c('#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')


# Get the position for the stats. 90% of max value. Reponse is a vector of values for y axis.
get_x_pos <- function(response){
  range <- (max(response) - min(response))*0.45
  max(response) - range
}

bateman_m_x <- get_x_pos(bateman_data$corr_m)
oss_x <- get_x_pos(oss_data$oss_log)
teste_x <- get_x_pos(teste_data$residual_testes_mass)

# Create a function to plot the data.
prediction_plot <- function(response = "beta_f", dataset = bateman_data, 
                        plot_label = bateman_female_label, 
                        x_label = "Female bateman gradient",
                        y_lab_pos = 6.1, point_size = "n_f", 
                        x_lab_pos = 1.5,
                        predict_data = bateman_female_preds){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", y = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", x = response)) + 

    geom_jitter(aes_string(size = point_size), 
                alpha = 0.75, width = 0, height = 0.2) + 
    scale_size_continuous(range = c(5,10)) +
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
    geom_line(data = predict_data, inherit.aes = FALSE, 
              aes_string(x = response, y = "estimate__"), linetype = "dashed", linewidth = 1)
}

# Bateman gradient plot.
pal <- c('#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
bateman_male_preds <- conditional_effects(bateman_male_model)[[1]] 
male_bateman_plot <- prediction_plot("weight_corr_m", plot_label = phy_bateman_male_label, 
                                     x_label = "Bateman Gradient", x_lab_pos = bateman_m_x,
                                     point_size = "n_m", predict_data = bateman_male_preds) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

# Opportunity for SS plot
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
oss_preds <- conditional_effects(oss_model)[[1]] 
oss_plot <- prediction_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                                 x_label = "Opportunity for sexual selection", point_size = "n_oss", 
                                 predict_data = oss_preds, x_lab_pos = oss_x) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

# Testes mass plot.
teste_data$teste_n <- 1
teste_preds <- conditional_effects(teste_model)[[1]]
teste_plot <- prediction_plot("residual_testes_mass", teste_data, plot_label = phy_teste_label, 
                x_label = "Residual testes mass", point_size = "teste_n",
                predict_data = teste_preds, x_lab_pos = teste_x)




ggarrange(teste_plot, oss_plot, male_bateman_plot, 
          nrow = 1, ncol = 3, widths = c(1.1,1,1), labels = c("a", "b", "c"),
          hjust = c(-3.5, -1.4, -1.4), font.label = list(size = 28))

ggsave("Plots/Diagnostics/ss_metric_analysis.png", width = 12, height = 4)



ssd_data <- read.csv("../Size_dimorphism/Data/trait_averages.csv") %>% clean_names()
colnames(ssd_data)[1] <- "birdtree_name"

ssd_data %<>% select(birdtree_name, hand_wing_index)

teste_wing_data <- teste_data %>% left_join(ssd_data)

wing_teste_model <- lm(residual_testes_mass ~ hand_wing_index, data =teste_wing_data)

summary(wing_teste_model)

plot(jitter(teste_wing_data$sexual_selection) ~ wing_teste_model$residuals)

teste_wing_data$wing_resids <- wing_teste_model$residuals

# Log relative size.
teste_wing_data %>% ggplot(aes(group = as.factor(sexual_selection), x = sexual_selection,
                          colour = as.factor(sexual_selection),
                          fill = as.factor(sexual_selection), y = wing_resids)) + 
  geom_jitter(size = 3, alpha = 0.5) + 
  geom_boxplot(alpha = 1, colour = "black", fill = NA, outlier.shape = NA, size = 0.75) + 
  theme_classic(base_size = 20) +  
  scale_fill_manual(values = pal) + 
  scale_colour_manual(values = pal) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5)) + 
  xlab("Sexual selection") + ylab("Residual testes mass")

hist(teste_wing_data$hand_wing_index)


teste_data$resid_teste <- teste_data$teste/teste_data$body_mass

plot(teste_data$resid_teste ~ teste_data$teste)

plot(teste_data$sexual_selection ~ teste_data$teste)
plot(teste_data$sexual_selection ~ teste_data$resid_teste)

hist(teste_data$resid_teste)
hist(teste_data$resid_teste_2)
teste_data$resid_teste_2 <- exp(teste_data$teste)/exp(teste_data$body_mass)

mean_teste_data <- teste_data %>% group_by(sexual_selection) %>% summarise(mean_teste = mean(exp(resid_teste)))


prediction_plot("resid_teste", teste_data, plot_label = phy_teste_label, 
                x_label = "Opportunity for sexual selection", sex = "male", 
                point_size = "teste_n", symbol_y = 1.05, symbol_scale = 0.85, symbol_x = 3.9,
                predict_data = teste_preds, x_lab_pos = teste_x,  y_lab_pos = 5.4) + xlim(c(-2, 0.5)) + 
  geom_point(data = mean_teste_data, aes(x = mean_teste, y = sexual_selection), size = 10, inherit.aes = FALSE)

#################################################################################
           ##### Create the plots using smoothed predictions #####


# Create a function to plot the data.
smooth_pred_plot <- function(response = "beta_f", dataset = bateman_data, 
                            plot_label = bateman_female_label, 
                            x_label = "Female bateman gradient",
                            y_lab_pos = 4.4, point_size = "n_f", 
                            x_lab_pos = 1.5,
                            sex = "female", 
                            symbol_x = 4, symbol_y = 1.1, symbol_scale = 0.65,
                            predict_data = bateman_female_preds){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", y = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", x = response)) + 
    
    geom_jitter(aes_string(size = point_size), 
                alpha = 0.75, width = 0, height = 0.2) + 
    scale_size_continuous(range = c(5,10)) +
    theme_classic(base_size = 20) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    
    ylab("Sexual selection") + 
    xlab(x_label) +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    annotate("text", y = y_lab_pos, x = x_lab_pos, label = plot_label, size = 7) +
    geom_smooth(data = predict_data, inherit.aes = FALSE, 
              aes_string(x = response, y = "estimate__"),
              linetype = "dashed", linewidth = 1, col = "black", se = FALSE)
}


oss_plot_2 <- smooth_pred_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                              x_label = "Opportunity for sexual selection", sex = "male", 
                              point_size = "n_oss", symbol_y = 1.05, symbol_scale = 0.85, symbol_x = 3.9,
                              predict_data = oss_preds, x_lab_pos = oss_x,  y_lab_pos = 5.4) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())
























## Creating it on the same scale.

female_bateman_plot <- female_bateman_plot + scale_x_continuous(limits = c(-0.1, 0.9))
male_bateman_plot <- male_bateman_plot + scale_x_continuous(limits = c(-0.1, 0.9))

female_oss_plot <- female_oss_plot + scale_x_continuous(limits = c(-3, 2))
male_oss_plot <- male_oss_plot + scale_x_continuous(limits = c(-3, 2))

ggarrange(female_bateman_plot, female_oss_plot, male_bateman_plot, male_oss_plot,
          nrow = 2, ncol = 2, heights = c(1,1.1), labels = c("a", "b", "c", "d"),
          hjust = c(-4.9, -2.5, -4.9, -2.5), font.label = list(size = 28))

ggsave("Plots/Diagnostics/ss_correlations_same_scale.png", width = 12, height = 12)





library(ggrepel)



# Get the position for the stats. 90% of max value. Reponse is a vector of values for y axis.
get_y_pos <- function(response){
  range <- (max(response) - min(response))*0.1
  max(response) - range
}

# Extract the y value for the labels.
bateman_f_y <- get_y_pos(bateman_data$corr_f)
bateman_m_y <- get_y_pos(bateman_data$corr_m)
oss_f_y <-  get_y_pos(oss_data$oss_log_f)
oss_m_y <- get_y_pos(oss_data$oss_log_m)


################################################################################
##### Label repel plots #######


# Function to make label plot.
label_repel_plot <- function(response = "beta_f", dataset = bateman_data, 
                             plot_label = bateman_female_label, 
                             y_label = "Female bateman gradient",
                             y_lab_pos = 1.5){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", x = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", y = response)) + 
    theme_classic(base_size = 20) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    xlab("Sexual selection") + ylab(y_label) +
    geom_label_repel(aes(label = birdtree_name), colour = "black") +
    annotate("text", x = 1, y = y_lab_pos, label = plot_label, size = 7)
}


# Create label repel plots.
female_bateman_plot <- label_repel_plot("weight_corr_f", plot_label = bateman_female_label, y_lab_pos = bateman_f_y)
male_bateman_plot <- label_repel_plot("weight_corr_m", plot_label = phy_bateman_male_label, 
                                      y_label = "Male Bateman Gradient", y_lab_pos = bateman_m_y)
mean_bateman_plot <- label_repel_plot("corr_mean", plot_label = phy_bateman_male_label, 
                                      y_label = "Male Bateman Gradient", y_lab_pos = bateman_m_y)


female_oss_plot <- label_repel_plot("oss_log_f", oss_data, plot_label = oss_female_label, 
                                    y_label = "Female OSS", y_lab_pos = oss_f_y) 
male_oss_plot <- label_repel_plot("oss_log_m", oss_data, plot_label = oss_label, 
                                  y_label = "Male OSS", y_lab_pos = oss_m_y)   #ylab("Male opportunity for sexual selection")

oss_plot <- label_repel_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                             y_label = "Male OSS")

# Arrange them together.
ggarrange(female_bateman_plot + rremove("xlab"), male_bateman_plot + rremove("xlab"), female_oss_plot, male_oss_plot)

# Export them.
ggsave("Plots/Diagnostics/ss_metrics.png", width = 12, height = 12)
ggsave("Plots/Diagnostics/ss_metrics.pdf", width = 8, height = 8)






##################################################################################

library(jpeg)
library(cowplot)

# Size of images 2480 x 3508
male_symbol <- readJPEG("Data/Images/male_symbol.jpg")
female_symbol <- readJPEG("Data/Images/female_symbol.jpg")

# # Sex specific plot settings.
# male_image <- draw_image(male_symbol, x = 0.25, y= 0.625, scale = 0.7)
# female_image <- list(draw_image(female_symbol, x = 0.2, y= 0.60, scale = 0.65))


jitter_plot <- function(response = "beta_f", dataset = bateman_data, 
                        plot_label = bateman_female_label, 
                        y_label = "Female bateman gradient",
                        y_lab_pos = 1.5, point_size = "n_f", 
                        sex = "female", 
                        symbol_x = 3.9, symbol_y = 1, symbol_scale = 0.65){
  
  if (sex == "female"){
    symbol <- female_symbol
  } else {
    symbol <- male_symbol
  }
  
  # Expand the range.
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", x = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", y = response)) + 
    draw_image(symbol, x = symbol_x, y = symbol_y, scale = symbol_scale, clip = "off") +
    geom_jitter(aes_string(size = point_size), #size = 5, 
                alpha = 0.75, width = 0.2, height = 0) + 
    scale_size_continuous(range = c(5,10)) +
    theme_classic(base_size = 20) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    #xlab("Sexual selection \U2640") + 
    xlab("Sexual selection") + 
    ylab(y_label) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    annotate("text", x = 1, y = y_lab_pos, label = plot_label, size = 7)
}



female_bateman_plot <- jitter_plot("weight_corr_f", plot_label = bateman_female_label, 
                                   y_label = "Bateman gradient",  y_lab_pos = bateman_f_y, symbol_y = 0.2) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_smooth(method = "lm") 

male_bateman_plot <- jitter_plot("weight_corr_m", plot_label = bateman_male_label, 
                                 y_label = "Male Bateman Gradient", y_lab_pos = bateman_m_y, sex = "male", 
                                 point_size = "n_m",
                                 symbol_y = 0.4, symbol_x = 3.9) +
  theme(axis.title = element_blank(), axis.text.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = 0.2), breaks = c(0,0.5, 1.0))

female_oss_plot <- jitter_plot("oss_log_f", oss_data, plot_label = oss_female_label, 
                               y_label = "Opportunity for sexual selection", y_lab_pos = oss_f_y, symbol_y = -1.1)
#draw_image(female_symbol, x = 4, y = -1.1, scale = 0.65, clip = "off")

male_oss_plot <- jitter_plot("oss_log_m", oss_data, plot_label = oss_label, 
                             y_label = "Male OSS", y_lab_pos = oss_m_y, sex = "male", 
                             point_size = "n_m", symbol_y = 1.05, symbol_scale = 0.85, symbol_x = 3.9) + 
  theme(axis.title.y = element_blank())



ggarrange(female_bateman_plot, male_bateman_plot, female_oss_plot, male_oss_plot,
          nrow = 2, ncol = 2, heights = c(1,1.1), labels = c("a", "b", "c", "d"),
          hjust = c(-4.9, -2.5, -4.9, -2.5), font.label = list(size = 28))
ggsave("Plots/Diagnostics/ss_jitter_metrics.png", width = 12, height = 12)



###############################################################################
##### Points with labels #######

# Function to make label plot.
point_label_plot <- function(response = "beta_f", dataset = bateman_data, 
                             plot_label = bateman_female_label, 
                             y_label = "Female bateman gradient",
                             y_lab_pos = 1.5){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", x = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", y = response)) + 
    geom_point(size = 5, alpha = 0.99) + 
    theme_classic(base_size = 20) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    xlab("Sexual selection") + ylab(y_label) +
    geom_text_repel(aes(label = birdtree_name), colour = "black",
                    force = 10, max.time = 10, max.iter = 50000, 
                    min.segment.length = 0, seed = 1993,
                    segment.linetype = 5,
                    point.padding = 5, point.size = 5) +
    annotate("text", x = 1, y = y_lab_pos, label = plot_label, size = 7)
}

# Make the plots.
female_bateman_plot <- point_label_plot("weight_corr_f", plot_label = bateman_female_label, y_lab_pos = bateman_f_y)
male_bateman_plot <- point_label_plot("weight_corr_m", plot_label = bateman_male_label, y_label = "Male Bateman Gradient", y_lab_pos = bateman_m_y)
female_oss_plot <- point_label_plot("oss_log_f", oss_data, plot_label = oss_female_label, y_label = "Female OSS", y_lab_pos = oss_f_y) 
male_oss_plot <- point_label_plot("oss_log_m", oss_data, plot_label = oss_label, y_label = "Male OSS", y_lab_pos = oss_m_y)

# Export the plots.
ggarrange(female_bateman_plot + rremove("xlab"), male_bateman_plot + rremove("xlab"), female_oss_plot, male_oss_plot)
ggsave("Plots/Diagnostics/ss_point_label_metrics.png", width = 12, height = 12)



###############################################################################
###### Jitter #######


jitter_plot <- function(response = "beta_f", dataset = bateman_data, 
                        plot_label = bateman_female_label, 
                        y_label = "Female bateman gradient",
                        y_lab_pos = 1.5, point_size = "n_f"){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", x = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", y = response)) + 
    geom_jitter(aes_string(size = point_size), #size = 5, 
                alpha = 0.75, width = 0.2) + 
    scale_size_continuous(range = c(5,10)) +
    theme_classic(base_size = 20) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    xlab("Sexual selection") + ylab(y_label) +
    annotate("text", x = 1, y = y_lab_pos, label = plot_label, size = 7)
}


female_bateman_plot <- jitter_plot("weight_corr_f", plot_label = bateman_female_label, y_lab_pos = bateman_f_y) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

male_bateman_plot <- jitter_plot("weight_corr_m", plot_label = bateman_male_label, 
                                 y_label = "Male Bateman Gradient", y_lab_pos = bateman_m_y) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

female_oss_plot <- jitter_plot("oss_log_f", oss_data, plot_label = oss_female_label, 
                               y_label = "Female OSS", y_lab_pos = oss_f_y) 
male_oss_plot <- jitter_plot("oss_log_m", oss_data, plot_label = oss_label, 
                             y_label = "Male OSS", y_lab_pos = oss_m_y)

jitter_plot("oss_log", oss_data, plot_label = "", 
            y_label = "Male OSS", point_size = "n_oss")

oss_plot <- label_repel_plot("oss_log", oss_data, plot_label = "", 
                             y_label = "Male OSS")



ggarrange(female_bateman_plot, male_bateman_plot, female_oss_plot, male_oss_plot,
          nrow = 2, ncol = 2, heights = c(1,1.1))
ggsave("Plots/Diagnostics/ss_jitterl_metrics.png", width = 12, height = 12)

























##################################################################################
               #######m      END                   ########

predictions <- conditional_effects(bateman_female_model)[[1]] 


conditional_effects(bateman_female_model)


bateman_data %>% ggplot(aes_string(group = "sexual_selection_plot", y = "sexual_selection",
                              colour = "sexual_selection_plot",
                              fill = "sexual_selection_plot", x = "weight_corr_f")) +
  geom_jitter(aes_string(size = "n_f"), 
              alpha = 0.75, width = 0, height = 0.2) +
  scale_size_continuous(range = c(5,10)) +
  theme_classic(base_size = 20) +  
  scale_fill_manual(values = pal) + 
  scale_colour_manual(values = pal) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5)) + 
  #xlab("Sexual selection \U2640") + 
  ylab("Sexual selection") + 
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  geom_line(data = predictions, inherit.aes = FALSE, 
            aes(x = weight_corr_f, y = (estimate__)), linetype = "dashed", linewidth = 1)

colnames(bateman_data)
# Extract predictions from brms model.
predictions <- conditional_effects(bateman_female_model)[[1]] 



bateman_female_p_value <- pd_to_p(last(p_direction(bateman_female_model)[,2]))
bateman_male_p_value <- pd_to_p(last(p_direction(bateman_male_model)[,2])) %>% round(2)

oss_female_p_value <- pd_to_p(last(p_direction(oss_female_model)[,2])) %>% round(2)
oss_male_p_value <- pd_to_p(last(p_direction(oss_male_model)[,2])) 

# Extract r squared.
bateman_female_r2 <- round(Bayes_R2_MZ(bateman_female_model)[[1]], 2)
bateman_male_r2 <- round(Bayes_R2_MZ(bateman_male_model)[[1]], 2)

oss_female_r2 <- round(Bayes_R2_MZ(oss_female_model)[[1]], 2)
oss_male_r2 <- round(Bayes_R2_MZ(oss_male_model)[[1]], 2)














################################################################################
                        #### Opportunity for selection ######




label_repel_plot("selection_m", oss_data)
label_repel_plot("selection_f", oss_data)

label_repel_plot("weight_selection_m", oss_data)
label_repel_plot("weight_selection_f", oss_data)

label_repel_plot("selection_log_f", oss_data)
label_repel_plot("selection_log_m", oss_data)




