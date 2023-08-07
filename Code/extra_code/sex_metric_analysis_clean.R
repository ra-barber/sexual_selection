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
                  #### Read in data ####


# Read in some data.
avian_data <- read.csv("Data/bateman_gradients/tim_avian_data.csv") %>% clean_names()
avian_data$tree_tip <- gsub(" ", "_", avian_data$birdtree_name)

oss_data <- read.csv("Data/bateman_gradients/Tim_OSS_database.csv") %>% clean_names()
oss_data$tree_tip <- gsub(" ", "_", oss_data$birdtree_name)
oss_data %<>% filter(exclude != "YES")

# Select batemman graident metrics.
bateman_data <- avian_data %>% 
  dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, 
                n_m, r_m) %>% na.omit()
#bateman_data[bateman_data$beta_f > bateman_data$beta_m, "birdtree_name"]

# Select OSS metrics.
oss_data <- oss_data %>% 
  dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, sex_role_reversal,
                n_f, is_f, n_m, is_m)


###############################################################################
                  ##### Prepare single OSS metric #####


# Create a single column of OSS, swapping over Jacana jacana as a role reversed species.
oss_data$n_oss <- oss_data$n_m
oss_data$oss <- oss_data$is_m

# Species to take female OSS from.
oss_data[oss_data$is_f > oss_data$is_m,]
female_species <- c("Clamator glandarius",  "Jacana jacana", "Zonotrichia albicollis", "Zonotrichia leucophrys")

# Get data from females for these species.
oss_data$oss[oss_data$birdtree_name %in% female_species] <- oss_data$is_f[oss_data$birdtree_name %in% female_species]
oss_data$n_oss[oss_data$birdtree_name %in% female_species] <- oss_data$n_f[oss_data$birdtree_name %in% female_species]

# Try removing Malurus for now (looks like OSS was estimated incorrectly (altho could be correct))
oss_data %<>% filter(birdtree_name != "Malurus cyaneus")

# Remove extra columns and NA species.
oss_data %<>% 
  dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, sex_role_reversal, n_oss, oss) %>% na.omit()
 

###############################################################################
                  ##### Add zero bateman species #####

# Filter for species with 0 oss, and therefore 0 bateman (in theory).
zero_mating_variance <- oss_data %>% filter(oss == 0)

# Select columns and merge.
extra_bateman_data <- zero_mating_variance %>% select(birdtree_name, tree_tip, sexual_selection, data_certainty, n_oss, oss)
colnames(extra_bateman_data)[5:6] <- colnames(bateman_data)[5:6]
bateman_data <- rbind(bateman_data, extra_bateman_data)

###############################################################################
                   ##### Average studies  #####

# Weight metrics by sample size for combining studies.
bateman_data$weighted_corr_m <- bateman_data$r_m * bateman_data$n_m
oss_data$weighted_oss <- oss_data$oss * oss_data$n_oss  

# Combine multiple studies.
bateman_data %<>% group_by(birdtree_name) %>% 
  summarise(tree_tip = first(tree_tip),
            sexual_selection = first(sexual_selection),
            data_certainty = first(data_certainty),
            
            corr_m = mean(r_m),
            weight_corr_m = sum(weighted_corr_m)/ sum(n_m),
            n_m = sum(n_m))

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
  corr_m_z = standardize(weight_corr_m, two_sd = TRUE)
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

phy_bateman_male_r2 <- as.character(format(round(marginal_R2_MZ(phy_bateman_male_model)[[1]], 2), nsmall = 2))
phy_oss_r2 <- as.character(format(round(marginal_R2_MZ(phy_oss_model)[[1]], 2), nsmall = 2))
phy_teste_r2 <- as.character(format(round(marginal_R2_MZ(phy_teste_model)[[1]], 2), nsmall = 2))

# Label without line breaks.
phy_bateman_male_label <-  paste0("\U03B2 = ", phy_bateman_male_estimate, "   p < 0.001   R\u00b2 = ", phy_bateman_male_r2)
phy_oss_label <-  paste0("\U03B2 = ", phy_oss_estimate, "   p < 0.001   R\u00b2 = ", phy_oss_r2)
phy_teste_label <-  paste0("\U03B2 = ", phy_teste_estimate, "   p < 0.001   R\u00b2 = ", phy_teste_r2)




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





###############################################################################
                    #### remove extra zeroes. ######

# Brms formula.
bateman_male_formula <- brmsformula("sexual_selection ~ corr_m_z + (1|gr(tree_tip, cov=A))", family = cumulative())
oss_formula <- brmsformula("sexual_selection ~ oss_log_z + (1|gr(tree_tip, cov=A))", family = cumulative())

sens_bateman_data <- bateman_data %>% filter(corr_m != 0)
sens_oss_data <- oss_data %>% filter(weight_oss != 0)

# Scale continuous predictors to two SD.
sens_bateman_data %<>% mutate(corr_m_z = standardize(weight_corr_m, two_sd = TRUE))
sens_oss_data %<>% mutate(oss_log_z = standardize(oss_log, two_sd = TRUE))

# Phy models.
sens_bateman_model <- brm(
  bateman_male_formula, data = sens_bateman_data, data2 = list(A=bateman_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")

sens_oss_model <- brm(
  oss_formula, data = oss_data, data2 = list(A=oss_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
  init = 0, normalize = FALSE, backend = "cmdstanr")


Bayes_R2_MZ(sens_bateman_model)
Bayes_R2_MZ(sens_oss_model)

marginal_R2_MZ(sens_bateman_model)
marginal_R2_MZ(sens_oss_model)


################################################################################
                 #### Export summary tables ####

# Extract relevant coefficient information.
bateman_estimates <- summary(phy_bateman_male_model)$fixed[5,c(1,3,4)]
oss_estimates <- summary(phy_oss_model)$fixed[5,c(1,3,4)]
teste_estimates <- summary(phy_teste_model)$fixed[5,c(1,3,4)]
sens_bateman_estimates <- summary(sens_bateman_model)$fixed[5,c(1,3,4)]
sens_oss_estimates <- summary(sens_oss_model)$fixed[5,c(1,3,4)]


all_estimates <- rbind(teste_estimates, bateman_estimates, sens_bateman_estimates, 
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
                       ##### Raw correlations #####


cor(teste_data$sexual_selection, teste_data$residual_testes_mass, method = "spearman")
cor(oss_data$sexual_selection, oss_data$oss_log, method = "spearman")
cor(bateman_data$sexual_selection, bateman_data$weight_corr_m, method = "spearman")

cor(teste_data$sexual_selection, teste_data$residual_testes_mass, method = "pearson")
cor(oss_data$sexual_selection, oss_data$oss_log, method = "pearson")
cor(bateman_data$sexual_selection, bateman_data$weight_corr_m, method = "pearson")






################################################################################
                        ##### Plotting #####

# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

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
    geom_line(data = predict_data, inherit.aes = FALSE, 
              aes_string(x = response, y = "estimate__"), linetype = "dashed", linewidth = 1)
}

# Bateman gradient plot.
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


################################################################################
               ##### Make the EDF figure version #####


# Get the position for the sample size label. 80% of max value. Reponse is a vector of values for y axis.
get_x_pos <- function(response){
  range <- (max(response) - min(response))*0.2   # 0.05 with just n = 65
  max(response) - range
}

bateman_sample_x <- get_x_pos(bateman_data$corr_m)
oss_sample_x <- get_x_pos(oss_data$oss_log)
teste_sample_x <- get_x_pos(teste_data$residual_testes_mass)

# Create the labels with italic "n"
bateman_sample_n <-  expression(~italic(n)~'= 65 species')
oss_sample_n <-  expression(~italic(n)~'= 79 species')
teste_sample_n <-  expression(~italic(n)~'= 977 species')

# Function for OSS and BG plots.
smooth_pred_plot <- function(response = "beta_f", dataset = bateman_data, 
                             plot_label = bateman_female_label, 
                             x_label = "Female bateman gradient",
                             y_lab_pos = 6.1, point_size = "n_f", 
                             x_lab_pos = 1.5, sample_x = bateman_sample_x, sample_label = bateman_sample_n,
                             predict_data = bateman_female_preds){
  
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
smooth_pred_plot_teste <- function(response = "beta_f", dataset = bateman_data, 
                                   plot_label = bateman_female_label, 
                                   x_label = "Female bateman gradient",
                                   y_lab_pos = 6.1, point_size = "n_f", 
                                   x_lab_pos = 1.5,
                                   predict_data = bateman_female_preds){
  
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


# Testes mass plot.
teste_plot_4 <- smooth_pred_plot_teste("residual_testes_mass", teste_data, plot_label = phy_teste_label, 
                                       x_label = "Residual testes mass", point_size = "teste_n",
                                       predict_data = teste_preds, x_lab_pos = teste_x)

# Bateman plot.
male_bateman_plot_3 <- smooth_pred_plot("weight_corr_m", plot_label = phy_bateman_male_label, 
                                        x_label = "Bateman gradient", x_lab_pos = bateman_m_x,
                                        point_size = "n_m", predict_data = bateman_male_preds) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

# Opportunity for sexual selection plot.
oss_plot_3 <- smooth_pred_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                               x_label = "Opportunity for sexual selection", point_size = "n_oss", 
                               predict_data = oss_preds, x_lab_pos = oss_x, sample_x = oss_sample_x,
                               sample_label = oss_sample_n) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

# Combine plots together.
ggarrange(teste_plot_4,  male_bateman_plot_3, oss_plot_3,
          nrow = 1, ncol = 3, widths = c(1.1,1,1), labels = c("a", "b", "c"),
          hjust = c(-3.5, -1.4, -1.4), font.label = list(size = 28))

# Export.
ggsave("Plots/Diagnostics/ss_metric_analysis_edf_2.png", width = 12, height = 4.5)









#################################################################################
           ##### Create the plots using smoothed predictions #####


# Create a function to plot the data.
smooth_pred_plot <- function(response = "beta_f", dataset = bateman_data, 
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
    geom_smooth(data = predict_data, inherit.aes = FALSE, 
                aes_string(x = response, y = "estimate__"),
                linetype = "dashed", linewidth = 1, col = "black", se = FALSE)
}




# Testes mass plot.
teste_plot_2 <- smooth_pred_plot("residual_testes_mass", teste_data, plot_label = phy_teste_label, 
                              x_label = "Residual testes mass", point_size = "teste_n",
                              predict_data = teste_preds, x_lab_pos = teste_x)

male_bateman_plot_2 <- smooth_pred_plot("weight_corr_m", plot_label = phy_bateman_male_label, 
                                     x_label = "Bateman gradient", x_lab_pos = bateman_m_x,
                                     point_size = "n_m", predict_data = bateman_male_preds) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

oss_plot_2 <- smooth_pred_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                            x_label = "Opportunity for sexual selection", point_size = "n_oss", 
                            predict_data = oss_preds, x_lab_pos = oss_x) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

ggarrange(teste_plot_2,  male_bateman_plot_2, oss_plot_2,
          nrow = 1, ncol = 3, widths = c(1.1,1,1), labels = c("a", "b", "c"),
          hjust = c(-3.5, -1.4, -1.4), font.label = list(size = 28))

ggsave("Plots/Diagnostics/ss_metric_analysis_smoothed_all.png", width = 12, height = 4)



####################################################################################
                        ##### Add sample size ########




# Testes mass plot.
teste_plot_2 <- smooth_pred_plot("residual_testes_mass", teste_data, plot_label = phy_teste_label, 
                                 x_label = "Residual testes mass\nn = 977", point_size = "teste_n",
                                 predict_data = teste_preds, x_lab_pos = teste_x)

male_bateman_plot_2 <- smooth_pred_plot("weight_corr_m", plot_label = phy_bateman_male_label, 
                                        x_label = "Bateman gradient\nn = 65", x_lab_pos = bateman_m_x,
                                        point_size = "n_m", predict_data = bateman_male_preds) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

oss_plot_2 <- smooth_pred_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                               x_label = "Opportunity for sexual selection\nn = 79", point_size = "n_oss", 
                               predict_data = oss_preds, x_lab_pos = oss_x) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

ggarrange(teste_plot_2,  male_bateman_plot_2, oss_plot_2,
          nrow = 1, ncol = 3, widths = c(1.1,1,1), labels = c("a", "b", "c"),
          hjust = c(-3.5, -1.4, -1.4), font.label = list(size = 28))

#ggsave("Plots/Diagnostics/ss_metric_analysis_edf.png", width = 12, height = 4.5)




# Create a function to plot the data.
smooth_pred_plot_teste <- function(response = "beta_f", dataset = bateman_data, 
                             plot_label = bateman_female_label, 
                             x_label = "Female bateman gradient",
                             y_lab_pos = 6.1, point_size = "n_f", 
                             x_lab_pos = 1.5,
                             predict_data = bateman_female_preds){
  
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
    geom_smooth(data = predict_data, inherit.aes = FALSE, 
                aes_string(x = response, y = "estimate__"),
                linetype = "dashed", linewidth = 1, col = "black", se = FALSE)
}


# Testes mass plot.
teste_plot_3 <- smooth_pred_plot_teste("residual_testes_mass", teste_data, plot_label = phy_teste_label, 
                 x_label = "Residual testes mass\nn = 977", point_size = "teste_n",
                 predict_data = teste_preds, x_lab_pos = teste_x)



ggarrange(teste_plot_3,  male_bateman_plot_2, oss_plot_2,
          nrow = 1, ncol = 3, widths = c(1.1,1,1), labels = c("a", "b", "c"),
          hjust = c(-3.5, -1.4, -1.4), font.label = list(size = 28))

ggsave("Plots/Diagnostics/ss_metric_analysis_edf.png", width = 12, height = 4.5)




###############################################################################
                   ##### add sample size inside plot #####

# Get the position for the stats. 90% of max value. Reponse is a vector of values for y axis.
get_x_pos <- function(response){
  range <- (max(response) - min(response))*0.2   # 0.05 with just n = 65
  max(response) - range
}

bateman_sample_x <- get_x_pos(bateman_data$corr_m)
oss_sample_x <- get_x_pos(oss_data$oss_log)
teste_sample_x <- get_x_pos(teste_data$residual_testes_mass)

# bateman_sample_n <- paste0("*n* = ", nrow(bateman_data))
# oss_sample_n <- paste0("n = ", nrow(oss_data))
# teste_sample_n <- paste0("n = ", nrow(teste_data))

bateman_sample_n <-  expression(~italic(n)~'= 65 species')
oss_sample_n <-  expression(~italic(n)~'= 79 species')
teste_sample_n <-  expression(~italic(n)~'= 977 species')

# 
# 
# test_label <-   paste0("~italic(n)~ = ", nrow(oss_data))
# expression(test_label)
# expression('No. of'~italic(bacteria X)~'isolates with corresponding types')
# test_label <-  expression(~italic(n)~'= 69')
# 
# install.packages("mdthemes")
# 
# expression(paste("X comes from ",italic("normal distribution")))
# test_label <- expression(paste0(italic("n", " = ", nrow(bateman_data))))
# test_label <- expression(italic("n", " = ", nrow(bateman_data)))
# Create a function to plot the data.
smooth_pred_plot <- function(response = "beta_f", dataset = bateman_data, 
                             plot_label = bateman_female_label, 
                             x_label = "Female bateman gradient",
                             y_lab_pos = 6.1, point_size = "n_f", 
                             x_lab_pos = 1.5, sample_x = bateman_sample_x, sample_label = bateman_sample_n,
                             predict_data = bateman_female_preds){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", y = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", x = response)) + 
    
    geom_jitter(aes_string(size = point_size), 
                alpha = 0.75, width = 0, height = 0.2) + 
    scale_size_continuous(range = c(5,8)) +
    theme_classic(base_size = 18) +  
    #mdthemes::md_theme_classic(base_size = 18) +
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

# Create a function to plot the data.
smooth_pred_plot_teste <- function(response = "beta_f", dataset = bateman_data, 
                                   plot_label = bateman_female_label, 
                                   x_label = "Female bateman gradient",
                                   y_lab_pos = 6.1, point_size = "n_f", 
                                   x_lab_pos = 1.5,
                                   predict_data = bateman_female_preds){
  
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


# Testes mass plot.
teste_plot_4 <- smooth_pred_plot_teste("residual_testes_mass", teste_data, plot_label = phy_teste_label, 
                                       x_label = "Residual testes mass", point_size = "teste_n",
                                       predict_data = teste_preds, x_lab_pos = teste_x)


male_bateman_plot_3 <- smooth_pred_plot("weight_corr_m", plot_label = phy_bateman_male_label, 
                                        x_label = "Bateman gradient", x_lab_pos = bateman_m_x,
                                        point_size = "n_m", predict_data = bateman_male_preds) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

oss_plot_3 <- smooth_pred_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                               x_label = "Opportunity for sexual selection", point_size = "n_oss", 
                               predict_data = oss_preds, x_lab_pos = oss_x, sample_x = oss_sample_x,
                               sample_label = oss_sample_n) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

ggarrange(teste_plot_4,  male_bateman_plot_3, oss_plot_3,
          nrow = 1, ncol = 3, widths = c(1.1,1,1), labels = c("a", "b", "c"),
          hjust = c(-3.5, -1.4, -1.4), font.label = list(size = 28))

ggsave("Plots/Diagnostics/ss_metric_analysis_edf_test.png", width = 12, height = 4.5)


###############################################################################
            ##### With phy line #####

# Brms formula.
bateman_male_formula <- brmsformula("sexual_selection ~ weight_corr_m + (1|gr(tree_tip, cov=A))", family = cumulative())
oss_formula <- brmsformula("sexual_selection ~ oss_log + (1|gr(tree_tip, cov=A))", family = cumulative())

# Phy models.
phy_raw_bateman_male_model <- brm(
  bateman_male_formula, data = bateman_data, data2 = list(A=bateman_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")

phy_raw_oss_model <- brm(
  oss_formula, data = oss_data, data2 = list(A=oss_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
  init = 0, normalize = FALSE, backend = "cmdstanr")

teste_model

bateman_male_preds <- conditional_effects(phy_raw_bateman_male_model, re_formula = NULL)[[1]] 
oss_preds <- conditional_effects(phy_raw_oss_model)[[1]] 
#teste_preds <- conditional_effects(phy_teste_model)[[1]]




male_bateman_plot_3 <- smooth_pred_plot("weight_corr_m", plot_label = phy_bateman_male_label, 
                                        x_label = "Bateman gradient\nn = 65", x_lab_pos = bateman_m_x,
                                        point_size = "n_m", predict_data = bateman_male_preds) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

oss_plot_3 <- smooth_pred_plot("oss_log", oss_data, plot_label = phy_oss_label, 
                               x_label = "Opportunity for sexual selection\nn = 79", point_size = "n_oss", 
                               predict_data = oss_preds, x_lab_pos = oss_x) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 

# Testes mass plot.
teste_plot_4 <- smooth_pred_plot_teste("residual_testes_mass", teste_data, plot_label = phy_teste_label, 
                                       x_label = "Residual testes mass\nn = 977", point_size = "teste_n",
                                       predict_data = teste_preds, x_lab_pos = teste_x)



ggarrange(teste_plot_4,  male_bateman_plot_3, oss_plot_3,
          nrow = 1, ncol = 3, widths = c(1.1,1,1), labels = c("a", "b", "c"),
          hjust = c(-3.5, -1.4, -1.4), font.label = list(size = 28))


################################################################################
                    ###### Try all three #####

multivariate_data <- bateman_data %>% left_join(oss_data) %>% left_join(teste_data) %>% na.omit()


# Brms formula.
multivariate_formula <- brmsformula("sexual_selection ~ corr_m_z + oss_log_z*teste_z + (1|gr(tree_tip, cov=A))", family = cumulative(), decomp = "QR")

# Phy models.
phy_multivariate_model <- brm(
  multivariate_formula, data = multivariate_data, data2 = list(A=bateman_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, 
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  normalize = FALSE, backend = "cmdstanr")


# Summaries.
summary(phy_multivariate_model)
plot(phy_multivariate_model)
# Extract estimate.
phy_multivariate_estimate <- last(summary(phy_multivariate_model)$fixed[,1])
phy_multivariate_estimate <- as.character(format(round(phy_multivariate_estimate, 2), nsmall = 2))

# Extract p value.
phy_multivariate_p_value <- pd_to_p(last(p_direction(phy_multivariate_model)[,2])) %>% round(2)

# Extract r squared.
phy_multivariate_r2 <- round(Bayes_R2_MZ(phy_multivariate_model)[[1]], 2)



# and Mitchell, 1992; Veall and Zimmermann, 1994).
Bayes_R2_MZ <- function(fit, ...) {
  y_pred <- fitted(fit, scale = "linear", summary = FALSE, ...)
  var_fit <- apply(y_pred, 1, var)
  if (fit$formula$family$family == "cumulative" ||
      fit$formula$family$family == "bernoulli") {
    if (fit$formula$family$link == "probit" || 
        fit$formula$family$link == "probit_approx") {
      var_res <- 1
    }
    else if (fit$formula$family$link == "logit") {
      var_res <- pi^2 / 3 
    }
  } 
  else {
    sum_fit <- summary(fit)
    sig_res <- sum_fit$spec_pars["sigma", "Estimate"]
    var_res <- sig_res^2
  } 
  R2_MZ <- var_fit / (var_fit + var_res)
  print(
    data.frame(
      Estimate = mean(R2_MZ), 
      Est.Error = sd(R2_MZ), 
      "l-95% CI" = quantile(R2_MZ, 0.025),
      "u-95% CI" = quantile(R2_MZ, 0.975),
      row.names = "Bayes_R2_MZ", 
      check.names = FALSE), 
    digits = 6)
}


marginal_R2_MZ(phy_multivariate_model)
Bayes_R2_MZ(phy_oss_model)









###################################################################################
              ###### Adding in hand wing index scores for teste mass #####


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




