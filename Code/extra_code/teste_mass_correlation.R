###############################################################################
                        ##### Testes mass analysis #####
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
library(car)
library(ggplot2)
library(ggpubr)
library(phytools)
library(brms)
library(graph4lg)

# Read in the functions. 
source("Code/functions.R")




###############################################################################
                       #### Prepare Baker data ####


# Read in some data.
model_data <- read.csv("Data/teste_mass/baker_teste_mass.csv") %>% clean_names()
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)

# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]

# Drop tips on the tree.
model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Make a covariance matrix, and order data the same.
model_covar <- ape::vcv.phylo(model_tree)

# Reorder the matrix so it's random, to maximise parallel processing speed.
mat_order <- sample(1:nrow(model_covar), size = nrow(model_covar), replace = FALSE)
model_covar <- reorder_mat(model_covar, rownames(model_covar)[mat_order])
row.names(model_data) <- model_data$tree_tip
model_data <- model_data[row.names(model_covar),]

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  teste_z = standardize(teste, two_sd = TRUE)
)

# Prepare response variables.
model_data$sexual_selection <- model_data$sexual_selection + 1


###############################################################################
                    #### Set model formula ######


# Brms formula.
model_formula <- "sexual_selection ~ teste_z + (1|gr(tree_tip, cov=A))"
brms_formula <- brmsformula(model_formula, family = cumulative())

# # Add un-informative priors.
normal_priors <- c(prior(normal(0,1), class="Intercept"),
                   prior(normal(0,1), class="b"),
                   prior(gamma(2,1), "sd")) # Gamma 2,1 probs seems to work well given all previous models end up with values around 1

# Run brms models.
testes_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  prior = normal_priors,
  iter = 10000,
  warmup = 5000,
  chains = 2,
  thin = 2,
  cores = 4,
  init = 0,
  #file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  #control = list(adapt_delta = 0.6),
  threads = threading(2),
)


ordinal_bayes_r2 <- Bayes_R2_MZ(brms_model)

plot(conditional_effects(brms_model), points = TRUE, point_args = list(height = 0.3))


# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')


model_data %>% ggplot(aes(group = as.factor(sexual_selection), x = sexual_selection,
                          colour = as.factor(sexual_selection),
                          fill = as.factor(sexual_selection), y = teste)) + 
  geom_jitter(size = 3, alpha = 0.5) + 
  geom_boxplot(alpha = 1, colour = "black", fill = NA, outlier.shape = NA, size = 0.75) + 
  theme_classic(base_size = 20) +  
  scale_fill_manual(values = pal) + 
  scale_colour_manual(values = pal) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5)) + 
  xlab("Sexual selection") + ylab("Residual testes mass")










# Read in some data.
clements_data <- read.csv("Data/bateman_gradients/test_mass.csv") %>% clean_names()
gbd_teste_mass <- read.csv("Data/teste_mass/gbd_teste_mass.csv") %>% clean_names()
gbd_teste_mass %<>% na.omit()
colnames(gbd_teste_mass)[1] <- "birdtree_name"


# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)

gbd_teste_mass %<>% left_join(model_data)

# Read in a tree.
tree <- read.tree("Data/Trees/ebird_24may.tre")
jetz_tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]



data$tree_tip <-  data$scientific_name_clements %>% str_replace(" ", "_")
gbd_teste_mass$tree_tip <-  gbd_teste_mass$new_jetz_name_tree_tip_name %>% str_replace(" ", "_")


tree <- drop.tip(tree, setdiff(tree$tip.label, data$tree_tip))
jetz_tree <- drop.tip(jetz_tree, setdiff(jetz_tree$tip.label, gbd_teste_mass$tree_tip))


cor(data$sexual_selection, data$teste_mass, method = "spearman")

cor(gbd_teste_mass$sexual_selection, data$teste_mass, method = "spearman")


data$tree_tip <-  data$scientific_name_clements %>% str_replace(" ", "_")

row.names(data) <- data$tree_tip


library(phylolm)

library(ape)
library(phytools)

tree <- force.ultrametric(tree)

sex_model <- phylolm(sexual_selection ~ teste_mass, phy = tree, data = data)

hist(sex_model$residuals)

teste_model <- phylolm(teste_mass ~ sexual_selection, phy = tree, data = data)

hist(teste_model$residuals)

summary(teste_model)




###############################################################################
                           #### BRMS models ####

# Packages to load.
library(magrittr)
library(tictoc)
library(caper)
library(dplyr)
library(effectsize)
library(car)
library(ggplot2)
library(ggpubr)
library(phytools)
library(brms)
library(graph4lg)


# Read in the tree.
model_tree <- read.tree("Data/Trees/prum_trees.tre")[[1]]

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)

# Filter for high cert.
if (data_type == "high"){
  model_data %<>% filter(sexual_certainty < 3)
}

# Read in gbd teste mass data.
gbd_teste_mass <- read.csv("Data/teste_mass/gbd_teste_mass.csv") %>% clean_names()
gbd_teste_mass %<>% na.omit()
colnames(gbd_teste_mass)[1] <- "birdtree_name"

gbd_teste_mass %<>% distinct()

model_data %<>% left_join(gbd_teste_mass)

model_data <- tidyr::drop_na(model_data, residual_testes_mass)


# Drop tips on the tree.
model_tree <- drop.tip(model_tree, setdiff(model_tree$tip.label, model_data$tree_tip))

# Make a covariance matrix, and order data the same.
model_covar <- ape::vcv.phylo(model_tree)

# Reorder the matrix so it's random, to maximise parallel processing speed.
mat_order <- sample(1:nrow(model_covar), size = nrow(model_covar), replace = FALSE)
model_covar <- reorder_mat(model_covar, rownames(model_covar)[mat_order])
row.names(model_data) <- model_data$tree_tip
model_data <- model_data[row.names(model_covar),]

###############################################################################
                           #### Section 2 ####

# Scale continuous predictors to two SD.
model_data %<>% mutate(
  teste_z = standardize(residual_testes_mass, two_sd = TRUE)
)

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1


###############################################################################
#### Set model formula ######


model_formula <- "sexual_score ~ teste_z + (1|gr(tree_tip, cov=A))"


# brms formula.
brms_formula <- brmsformula(model_formula, family = cumulative())


# # Add un-informative priors.
normal_priors <- c(prior(normal(0,1), class="Intercept"),
                   prior(normal(0,1), class="b"),
                   prior(gamma(2,1), "sd")) # Gamma 2,1 probs seems to work well given all previous models end up with values around 1

# Run brms models.
brms_model <- brm(
  brms_formula,
  data = model_data,
  data2 = list(A=model_covar),
  prior = normal_priors,
  iter = 1000,
  warmup = 500,
  chains = 2,
  thin = 2,
  cores = 4,
  init = 0,
  #file = model_pathway,
  normalize = FALSE,
  backend = "cmdstanr",
  #control = list(adapt_delta = 0.6),
  threads = threading(2),
)




ordinal_bayes_r2 <- Bayes_R2_MZ(brms_model)

linear_bayes_r2 <- bayes_R2(brms_model)

plot(conditional_effects(brms_model), points = TRUE, jitter_width = 0.1)
conditional_effects(brms_model, surface = TRUE)


model_data %>% ggplot(aes(group = as.factor(sexual_score), fill = as.factor(sexual_score), y = teste_z)) + geom_boxplot(alpha = 0.5)

###############################################################################
                           ####  Clements data ####

# Read in some data.
clements_data <- read.csv("Data/bateman_gradients/test_mass.csv") %>% clean_names()

# Read in a tree.
clements_tree <- read.tree("Data/Trees/ebird_24may.tre")


clements_data %>% ggplot(aes(group = as.factor(sexual_selection), fill = as.factor(sexual_selection), y = teste_mass)) + geom_boxplot(alpha = 0.5)


###############################################################################
                         ####  Baker data ####


# Read in some data.
baker_data <- read.csv("Data/teste_mass/baker_teste_mass.csv") %>% clean_names()


baker_data %>% ggplot(aes(group = as.factor(sexual_selection), x = sexual_selection,
                          colour = as.factor(sexual_selection),
                          fill = as.factor(sexual_selection), y = teste)) + 
  geom_jitter(size = 3, alpha = 0.5) + 
  geom_boxplot(alpha = 1, colour = "black", fill = NA, outlier.shape = NA, size = 0.75) + 
  theme_classic(base_size = 20) +  
  scale_fill_manual(values = pal) + 
  scale_colour_manual(values = pal) + 
  theme(legend.position = "none",
        line = element_line(linewidth = 0.5)) + 
  xlab("Sexual selection") + ylab("Residual teste size")


# Create a palette to match bin length.
pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

ss_barplot <- function(dataset = full_data){
  dataset %>% ggplot(aes(x = sexual_score, fill = as.factor(sexual_score))) +
    geom_bar() +  scale_fill_manual(values = pal) + 
    theme_classic(base_size = 20) + 
    theme(legend.position = "none",
          line = element_line(linewidth = 0.5))
}



###############################################################################
                           #### Section 4 ####

###############################################################################
                           #### Section 5 ####

###############################################################################
                           #### Section 6 ####


###############################################################################
                           #### Section 7 ####


###############################################################################
                           #### Section 8 ####

# Look Rob, you've had your fun with the sectioning. 
# They'll be no more sectioning today.


###############################################################################
                             #### END ####
###############################################################################


###############################################################################
              #### All the stuff I'm afraid to delete ####


