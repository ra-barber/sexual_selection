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
#library(graph4lg)

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

oss_data <- oss_data %>% 
  dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, sex_role_reversal,
                n_f, is_f, n_m, is_m)

# Create a single column of OSS, swapping over Jacana jacana as a role reversed species.
oss_data$oss <- oss_data$is_m
oss_data$oss_n <- oss_data$n_m

# Species to take female OSS from.
oss_data[oss_data$is_f > oss_data$is_m,]
female_species <- c("Clamator glandarius",  "Jacana jacana", "Zonotrichia albicollis", "Zonotrichia leucophrys")

# Get data from females for these species.
oss_data$oss[oss_data$birdtree_name %in% female_species] <- oss_data$is_f[oss_data$birdtree_name %in% female_species]
oss_data$oss_n[oss_data$birdtree_name %in% female_species] <- oss_data$n_f[oss_data$birdtree_name %in% female_species]

# Change n to 1 for 
oss_data$oss_n[oss_data$birdtree_name == "Dendroica caerulescens"] <- 572 # Number of males included in the Germain et al. 2021 study.

# Remove extra columns and NA species.
oss_data %<>% 
  dplyr::select(birdtree_name, tree_tip, sexual_selection, data_certainty, sex_role_reversal, oss, oss_n) %>% na.omit()
 

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

oss_data$weighted_oss <- oss_data$oss * oss_data$oss_n  

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
            
            
            beta_t = mean(beta_t),
            beta_g = mean(beta_g),
            
            n_f = sum(n_f),
            n_m = sum(n_m))

# Combine multiple studies.
oss_data %<>% group_by(birdtree_name) %>% 
  summarise(tree_tip = first(tree_tip),
            sexual_selection = first(sexual_selection),
            data_certainty = first(data_certainty),
            
            oss_f = mean(is_f),
            weight_oss_f = sum(weighted_oss_f)/ sum(n_f),
            oss_log_f = log(weight_oss_f),
            
            oss_m = mean(is_m),
            weight_oss_m = sum(weighted_oss_m)/ sum(n_m),
            oss_log_m = log(weight_oss_m),
            
            selection_f = mean(i_f),
            weight_selection_f = sum(weighted_selection_f)/ sum(n_f),
            selection_log_f = log(weight_selection_f),
            
            selection_m = mean(i_m),
            weight_selection_m = sum(weighted_selection_m)/ sum(n_m),
            selection_log_m = log(weight_selection_m),
            
            n_f = sum(n_f),
            n_m = sum(n_m))

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
  beta_g_z = standardize(beta_g, two_sd = TRUE)
)

oss_data %<>% mutate(
  oss_f_z = standardize(weight_oss_f, two_sd = TRUE),
  oss_log_f_z = standardize(oss_log_f, two_sd = TRUE),
  oss_m_z = standardize(weight_oss_m, two_sd = TRUE),
  oss_log_m_z = standardize(oss_log_m, two_sd = TRUE),
)

# Make sexual selection a factor for plotting colours.
bateman_data$sexual_selection_plot <- bateman_data$sexual_selection %>% as.factor()
oss_data$sexual_selection_plot <- oss_data$sexual_selection %>% as.factor()

# Prepare response variables.
bateman_data$sexual_selection <- bateman_data$sexual_selection + 1
oss_data$sexual_selection <- oss_data$sexual_selection + 1


###############################################################################
                    #### Run brms models ######

# Brms formula.
bateman_female_formula <- brmsformula("sexual_selection ~ beta_sqrt_f_z + (1|gr(tree_tip, cov=A))", family = cumulative())
bateman_male_formula <- brmsformula("sexual_selection ~ beta_log_m_z + (1|gr(tree_tip, cov=A))", family = cumulative())

bateman_female_formula <- brmsformula("sexual_selection ~ corr_f_z + (1|gr(tree_tip, cov=A))", family = cumulative())
bateman_male_formula <- brmsformula("sexual_selection ~ corr_m_z + (1|gr(tree_tip, cov=A))", family = cumulative())

oss_female_formula <- brmsformula("sexual_selection ~ oss_log_f_z + (1|gr(tree_tip, cov=A))", family = cumulative())
oss_male_formula <- brmsformula("sexual_selection ~ oss_log_m_z + (1|gr(tree_tip, cov=A))", family = cumulative())


# Run phy brms models.
bateman_female_model <- brm(
  bateman_female_formula, data = bateman_data, data2 = list(A=bateman_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
  init = 0, normalize = FALSE, backend = "cmdstanr")

bateman_male_model <- brm(
  bateman_male_formula, data = bateman_data, data2 = list(A=bateman_covar),
  iter = 25000, warmup = 15000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")

oss_female_model <- brm(
  oss_female_formula, data = oss_data, data2 = list(A=oss_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
  init = 0, normalize = FALSE, backend = "cmdstanr")

oss_male_model <- brm(
  oss_male_formula, data = oss_data, data2 = list(A=oss_covar),
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")


###############################################################################
                  #### Extract summary stats ######


# Summaries.
summary(bateman_female_model)
summary(bateman_male_model)

summary(oss_female_model)
summary(oss_male_model)

# Extract estimate.
# Estimate
bateman_female_estimate <- last(summary(bateman_female_model)$fixed[,1])
bateman_male_estimate <- last(summary(bateman_male_model)$fixed[,1])

oss_female_estimate <- last(summary(oss_female_model)$fixed[,1])
oss_male_estimate <- last(summary(oss_male_model)$fixed[,1])

bateman_female_estimate <- as.character(format(round(bateman_female_estimate, 2), nsmall = 2))
bateman_male_estimate <- as.character(format(round(bateman_male_estimate, 2), nsmall = 2))
oss_female_estimate <- as.character(format(round(oss_female_estimate, 2), nsmall = 2))
oss_male_estimate <- as.character(format(round(oss_male_estimate, 2), nsmall = 2))


# Extract p value.
library(bayestestR)

bateman_female_p_value <- pd_to_p(last(p_direction(bateman_female_model)[,2]))
bateman_male_p_value <- pd_to_p(last(p_direction(bateman_male_model)[,2])) %>% round(2)

oss_female_p_value <- pd_to_p(last(p_direction(oss_female_model)[,2])) %>% round(2)
oss_male_p_value <- pd_to_p(last(p_direction(oss_male_model)[,2])) 

# Extract r squared.
bateman_female_r2 <- round(Bayes_R2_MZ(bateman_female_model)[[1]], 2)
bateman_male_r2 <- round(Bayes_R2_MZ(bateman_male_model)[[1]], 2)

oss_female_r2 <- round(Bayes_R2_MZ(oss_female_model)[[1]], 2)
oss_male_r2 <- round(Bayes_R2_MZ(oss_male_model)[[1]], 2)


bateman_female_label <- paste0("\U03B2 = ", bateman_female_estimate, "\np < 0.05\nr\u00b2 = ", bateman_female_r2)
bateman_male_label <-  paste0("\U03B2 = ", bateman_male_estimate, "\np = ", bateman_male_p_value, "\nr\u00b2 = ", bateman_male_r2)

oss_female_label <-  paste0("\U03B2 = ", oss_female_estimate, "\np = ", oss_female_p_value, "\nr\u00b2 = ", oss_female_r2)
oss_male_label <-  paste0("\U03B2 = ", oss_male_estimate, "\np < 0.01\nr\u00b2 = ", oss_male_r2)


bateman_female_label <- paste0("\U03B2 = ", bateman_female_estimate, "   p < 0.01   r\u00b2 = ", bateman_female_r2)
bateman_male_label <-  paste0("\U03B2 = ", bateman_male_estimate, "   p = ", bateman_male_p_value, "   r\u00b2 = ", bateman_male_r2)

oss_female_label <-  paste0("\U03B2 = ", oss_female_estimate, "   p = ", oss_female_p_value, "   r\u00b2 = ", oss_female_r2)
oss_male_label <-  paste0("\U03B2 = ", oss_male_estimate, "   p < 0.01   r\u00b2 = ", oss_male_r2)


################################################################################
                      ##### Plotting #####


# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')
pal <- c('#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')


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
male_bateman_plot <- label_repel_plot("weight_corr_m", plot_label = bateman_male_label, 
                                      y_label = "Male Bateman Gradient", y_lab_pos = bateman_m_y)

female_oss_plot <- label_repel_plot("oss_log_f", oss_data, plot_label = oss_female_label, 
                                    y_label = "Female OSS", y_lab_pos = oss_f_y) 
male_oss_plot <- label_repel_plot("oss_log_m", oss_data, plot_label = oss_male_label, 
                                  y_label = "Male OSS", y_lab_pos = oss_m_y)   #ylab("Male opportunity for sexual selection")

# Arrange them together.
ggarrange(female_bateman_plot + rremove("xlab"), male_bateman_plot + rremove("xlab"), female_oss_plot, male_oss_plot)

# Export them.
ggsave("Plots/Diagnostics/ss_metrics.png", width = 12, height = 12)
ggsave("Plots/Diagnostics/ss_metrics.pdf", width = 8, height = 8)



################################################################################
                  ##### Standard point plots  #######

# Function to make label plot.
point_plot <- function(response = "beta_f", dataset = bateman_data, 
                             plot_label = bateman_female_label, 
                             y_label = "Female bateman gradient",
                             y_lab_pos = 1.5){
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", x = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", y = response)) + 
    geom_point(size = 10, alpha = 0.99) + 
    theme_classic(base_size = 20) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    xlab("Sexual selection") + ylab(y_label) +
    annotate("text", x = 1, y = y_lab_pos, label = plot_label, size = 7)
}


female_bateman_plot <- point_plot("weight_corr_f", plot_label = bateman_female_label, y_lab_pos = bateman_f_y)

male_bateman_plot <- point_plot("weight_corr_m", plot_label = bateman_male_label, 
                                      y_label = "Male Bateman Gradient", y_lab_pos = bateman_m_y)

female_oss_plot <- point_plot("oss_log_f", oss_data, plot_label = oss_female_label, 
                              y_label = "Female OSS", y_lab_pos = oss_f_y) 
male_oss_plot <- point_plot("oss_log_m", oss_data, plot_label = oss_male_label, 
                            y_label = "Male OSS", y_lab_pos = oss_m_y)


ggarrange(female_bateman_plot + rremove("xlab"), male_bateman_plot + rremove("xlab"), female_oss_plot, male_oss_plot)
ggsave("Plots/Diagnostics/ss_point_metrics.png", width = 12, height = 12)



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
                    #min.segment.length = Inf, 
                    min.segment.length = 0, seed = 1993,
                    segment.linetype = 5,
                    point.padding = 5, point.size = 5) +
    annotate("text", x = 1, y = y_lab_pos, label = plot_label, size = 7)
}


female_bateman_plot <- point_label_plot("weight_corr_f", plot_label = bateman_female_label, y_lab_pos = bateman_f_y)

male_bateman_plot <- point_label_plot("weight_corr_m", plot_label = bateman_male_label, 
                                y_label = "Male Bateman Gradient", y_lab_pos = bateman_m_y)

female_oss_plot <- point_label_plot("oss_log_f", oss_data, plot_label = oss_female_label, 
                              y_label = "Female OSS", y_lab_pos = oss_f_y) 
male_oss_plot <- point_label_plot("oss_log_m", oss_data, plot_label = oss_male_label, 
                            y_label = "Male OSS", y_lab_pos = oss_m_y)


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
male_oss_plot <- jitter_plot("oss_log_m", oss_data, plot_label = oss_male_label, 
                                  y_label = "Male OSS", y_lab_pos = oss_m_y)


ggarrange(female_bateman_plot, male_bateman_plot, female_oss_plot, male_oss_plot,
          nrow = 2, ncol = 2, heights = c(1,1.1))
ggsave("Plots/Diagnostics/ss_jitterl_metrics.png", width = 12, height = 12)


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

male_oss_plot <- jitter_plot("oss_log_m", oss_data, plot_label = oss_male_label, 
                             y_label = "Male OSS", y_lab_pos = oss_m_y, sex = "male", 
                             point_size = "n_m", symbol_y = 1.05, symbol_scale = 0.85, symbol_x = 3.9) + 
  theme(axis.title.y = element_blank())



ggarrange(female_bateman_plot, male_bateman_plot, female_oss_plot, male_oss_plot,
          nrow = 2, ncol = 2, heights = c(1,1.1), labels = c("a", "b", "c", "d"),
          hjust = c(-4.9, -2.5, -4.9, -2.5), font.label = list(size = 28))
ggsave("Plots/Diagnostics/ss_jitter_metrics.png", width = 12, height = 12)


# ###############################################################################
#                ###### Add correlations #####
# 
# # Brms formula.
# bateman_female_formula <- brmsformula("sexual_selection ~ weight_corr_f + (1|gr(tree_tip, cov=A))", family = cumulative())
# bateman_male_formula <- brmsformula("sexual_selection ~ weight_corr_m + (1|gr(tree_tip, cov=A))", family = cumulative())
# 
# oss_female_formula <- brmsformula("sexual_selection ~ oss_log_f + (1|gr(tree_tip, cov=A))", family = cumulative())
# oss_male_formula <- brmsformula("sexual_selection ~ oss_log_m + (1|gr(tree_tip, cov=A))", family = cumulative())
# 
# 
# # Run phy brms models.
# bateman_female_model <- brm(
#   bateman_female_formula, data = bateman_data, data2 = list(A=bateman_covar),
#   iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
#   init = 0, normalize = FALSE, backend = "cmdstanr")
# 
# bateman_male_model <- brm(
#   bateman_male_formula, data = bateman_data, data2 = list(A=bateman_covar),
#   iter = 25000, warmup = 15000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
#   normalize = FALSE, backend = "cmdstanr")
# 
# oss_female_model <- brm(
#   oss_female_formula, data = oss_data, data2 = list(A=oss_covar),
#   iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
#   init = 0, normalize = FALSE, backend = "cmdstanr")
# 
# oss_male_model <- brm(
#   oss_male_formula, data = oss_data, data2 = list(A=oss_covar),
#   iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
#   normalize = FALSE, backend = "cmdstanr")
# 
# 
# 
# # Extract predictions from brms model.
# predictions <- conditional_effects(bateman_female_model)[[1]] 
# 
# # Change back to normal data values that were transformed for ordinal models..
# predictions %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
# 
# 
# geom_line(data = predictions, aes(x = abs_lat, y = (estimate__)), linetype = "dashed", linewidth = 1)
# 
# 
# 
# # Plot with x and y axis reversed.
# jitter_plot <- function(response = "beta_f", dataset = bateman_data, 
#                         plot_label = bateman_female_label, 
#                         y_label = "Female bateman gradient",
#                         y_lab_pos = 1.5, point_size = "n_f", 
#                         sex = "female", 
#                         symbol_x = 4, symbol_y = 1.1, symbol_scale = 0.65){
#   
#   if (sex == "female"){
#     symbol <- female_symbol
#   } else {
#     symbol <- male_symbol
#   }
#   
#   # Expand the range.
#   
#   
#   
#   dataset %>% ggplot(aes_string(group = "sexual_selection_plot", y = "sexual_selection_plot",
#                                 colour = "sexual_selection_plot",
#                                 fill = "sexual_selection_plot", x = response)) + 
#     draw_image(symbol, x = symbol_x, y = symbol_y, scale = symbol_scale, clip = "off") +
#     geom_jitter(aes_string(size = point_size), #size = 5, 
#                 alpha = 0.75, width = 0, height = 0.2) + 
#     scale_size_continuous(range = c(5,10)) +
#     theme_classic(base_size = 20) +  
#     scale_fill_manual(values = pal) + 
#     scale_colour_manual(values = pal) + 
#     theme(legend.position = "none",
#           line = element_line(linewidth = 0.5)) + 
#     #xlab("Sexual selection \U2640") + 
#     ylab("Sexual selection") + 
#     xlab(y_label) +
#     scale_x_continuous(expand = expansion(mult = 0.2)) +
#     annotate("text", x = 1, y = y_lab_pos, label = plot_label, size = 7)
# }
# 
# 
# jitter_plot() + geom_line(data = predictions, inherit.aes = FALSE, aes(x = weight_corr_f, y = (estimate__)), linetype = "dashed", linewidth = 1)
# bateman_data
# 


################################################################################
                      ##### Non-phylo models ####


# Brms formula.
bateman_female_formula <- brmsformula("sexual_selection ~ weight_corr_f", family = cumulative())
bateman_male_formula <- brmsformula("sexual_selection ~ weight_corr_m", family = cumulative())

oss_female_formula <- brmsformula("sexual_selection ~ oss_log_f", family = cumulative())
oss_male_formula <- brmsformula("sexual_selection ~ oss_log_m", family = cumulative())


# Run phy brms models.
bateman_female_model <- brm(
  bateman_female_formula, data = bateman_data, 
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
  init = 0, normalize = FALSE, backend = "cmdstanr")

bateman_male_model <- brm(
  bateman_male_formula, data = bateman_data, 
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")

oss_female_model <- brm(
  oss_female_formula, data = oss_data, 
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, control = list(adapt_delta = 0.99),
  init = 0, normalize = FALSE, backend = "cmdstanr")

oss_male_model <- brm(
  oss_male_formula, data = oss_data,
  iter = 10000, warmup = 5000, chains = 8, thin = 20, cores = 8, init = 0, control = list(adapt_delta = 0.99),
  normalize = FALSE, backend = "cmdstanr")


# Extract predictions from brms model.
bateman_female_preds <- conditional_effects(bateman_female_model)[[1]]
bateman_male_preds <- conditional_effects(bateman_male_model)[[1]] 
oss_female_preds <- conditional_effects(oss_female_model)[[1]]
oss_male_preds <- conditional_effects(oss_male_model)[[1]] 

# Change back to normal data values that were transformed for ordinal models..
bateman_female_preds %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
bateman_male_preds %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
oss_female_preds %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)
oss_male_preds %<>% mutate(estimate__ = estimate__ - 1, lower__ = lower__ - 1, upper__ = upper__ - 1)



bateman_female_label <- paste0("\U03B2 = ", bateman_female_estimate, "   p < 0.01   r\u00b2 = ", bateman_female_r2)
bateman_male_label <-  paste0("\U03B2 = ", bateman_male_estimate, "   p = ", bateman_male_p_value, "   r\u00b2 = ", bateman_male_r2)

oss_female_label <-  paste0("\U03B2 = ", oss_female_estimate, "   p = ", oss_female_p_value, "   r\u00b2 = ", oss_female_r2)
oss_male_label <-  paste0("\U03B2 = ", oss_male_estimate, "   p < 0.01   r\u00b2 = ", oss_male_r2)


# Get the position for the stats. 90% of max value. Reponse is a vector of values for y axis.
get_x_pos <- function(response){
  range <- (max(response) - min(response))*0.5
  max(response) - range
}

# Extract the y value for the labels.
bateman_f_x <- get_x_pos(bateman_data$corr_f)
bateman_m_x <- get_x_pos(bateman_data$corr_m)
oss_f_x <-  get_x_pos(oss_data$oss_log_f)
oss_m_x <- get_x_pos(oss_data$oss_log_m)





prediction_plot <- function(response = "beta_f", dataset = bateman_data, 
                        plot_label = bateman_female_label, 
                        y_label = "Female bateman gradient",
                        y_lab_pos = 1.5, point_size = "n_f", 
                        sex = "female", 
                        symbol_x = 4, symbol_y = 1.1, symbol_scale = 0.65,
                        predict_data = bateman_female_preds){
  
  if (sex == "female"){
    symbol <- female_symbol
    x_label <- "Female sexual selection"
  } else {
    symbol <- male_symbol
    x_label <- "Male sexual selection"
  }
  
  dataset %>% ggplot(aes_string(group = "sexual_selection_plot", y = "sexual_selection_plot",
                                colour = "sexual_selection_plot",
                                fill = "sexual_selection_plot", x = response)) + 
    draw_image(symbol, x = symbol_x, y = symbol_y, scale = symbol_scale, clip = "off") +
    geom_jitter(aes_string(size = point_size), #size = 5, 
                alpha = 0.75, width = 0, height = 0.2) + 
    scale_size_continuous(range = c(5,10)) +
    theme_classic(base_size = 20) +  
    scale_fill_manual(values = pal) + 
    scale_colour_manual(values = pal) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5)) + 
    #xlab("Sexual selection \U2640") + 
    ylab(x_label) + 
    xlab(y_label) +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    annotate("text", y = 4.4, x = y_lab_pos, label = plot_label, size = 7) +
    geom_line(data = predict_data, inherit.aes = FALSE, 
              aes_string(x = response, y = "estimate__"), linetype = "dashed", linewidth = 1)
}


female_bateman_plot <- prediction_plot("weight_corr_f", plot_label = bateman_female_label, 
                                   y_label = "Bateman gradient",  y_lab_pos = bateman_f_x, symbol_y = 4) +
  theme(axis.title.x = element_blank())

female_oss_plot <- prediction_plot("oss_log_f", oss_data, plot_label = oss_female_label, 
                               y_label = "Opportunity for sexual selection", 
                               y_lab_pos = oss_f_x, symbol_y = -1.1, predict_data = oss_female_preds) +
  theme(axis.title = element_blank(), axis.text.y = element_blank())

male_bateman_plot <- prediction_plot("weight_corr_m", plot_label = bateman_male_label, 
                                     y_label = "Bateman Gradient", y_lab_pos = bateman_m_x, sex = "male", 
                                     point_size = "n_m",
                                     symbol_y = 0.4, symbol_x = 3.9, predict_data = bateman_male_preds)
  #scale_x_continuous(expand = expansion(mult = 0.2), breaks = c(0,0.5, 1.0))



male_oss_plot <- prediction_plot("oss_log_m", oss_data, plot_label = oss_male_label, 
                             y_label = "Opportunity for sexual selection", y_lab_pos = oss_m_x, sex = "male", 
                             point_size = "n_m", symbol_y = 1.05, symbol_scale = 0.85, symbol_x = 3.9,
                             predict_data = oss_male_preds) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())



ggarrange(female_bateman_plot, female_oss_plot, male_bateman_plot, male_oss_plot,
          nrow = 2, ncol = 2, heights = c(1,1.1), labels = c("a", "b", "c", "d"),
          hjust = c(-4.9, -2.5, -4.9, -2.5), font.label = list(size = 28))


ggsave("Plots/Diagnostics/ss_correlations.png", width = 12, height = 12)
ggsave("Plots/Diagnostics/ss_correlations.pdf", width = 12, height = 12)


## Creating it on the same scale.

female_bateman_plot <- female_bateman_plot + scale_x_continuous(limits = c(-0.1, 0.9))
male_bateman_plot <- male_bateman_plot + scale_x_continuous(limits = c(-0.1, 0.9))

female_oss_plot <- female_oss_plot + scale_x_continuous(limits = c(-3, 2))
male_oss_plot <- male_oss_plot + scale_x_continuous(limits = c(-3, 2))

ggarrange(female_bateman_plot, female_oss_plot, male_bateman_plot, male_oss_plot,
          nrow = 2, ncol = 2, heights = c(1,1.1), labels = c("a", "b", "c", "d"),
          hjust = c(-4.9, -2.5, -4.9, -2.5), font.label = list(size = 28))

ggsave("Plots/Diagnostics/ss_correlations_same_scale.png", width = 12, height = 12)




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




