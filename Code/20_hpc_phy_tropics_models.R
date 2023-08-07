###############################################################################
                  # Tropical brms models on cluster  #
###############################################################################

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

# Set the seed.
set.seed(1993)


################################################################################
                     #### Set up array iteration ####


# Get the array number from the job script.
array_number <- as.numeric(Sys.getenv("ARRAY_NUMBER"))
array_number

# Create types for hpc jobs.
tree_number <- 1:25

# Model (In order of size and therefore speed)
model_type <- c("frugivore", "primary", "invertivore", "secondary",
                "mig", "no_mig", "terr", "no_terr")

# Set the data types.
data_type <- c("all", "high")


# Expand the grid.
all_combos <- expand.grid(tree_number, data_type, model_type)

# Tree type.
tree_number <- all_combos[array_number, 1] %>% as.numeric()

# Data type.
data_type <- all_combos[array_number, 2] %>% as.character()

# Model type.
model_type <- all_combos[array_number, 3] %>% as.character()


###############################################################################
                      #### Read in the data #####

# Functions.
source("Code/functions.R")

# Read in the life history traits.
model_data <- read.csv("Data/sexual_traits.csv")
model_data$tree_tip <- gsub(" ", "_", model_data$birdtree_name)


###############################################################################
                #### Prepare predictor variables ######

# Prepare response variables.
model_data$sexual_score <- model_data$sexual_score + 1

# Add the tropical non-tropical models.
model_data$trop_non_trop <- NA
model_data$trop_non_trop[abs(model_data$complete_latitude) < 23.43624] <- "trop"
model_data$trop_non_trop[is.na(model_data$trop_non_trop)] <- "non_trop"

# Filter for high cert.
if (data_type == "high"){
  model_data %<>% filter(sexual_certainty < 3)
}

# Filter for trophic level.
if (model_type == "primary"){
  model_data %<>% filter(trophic_binary == "Primary")
}
if (model_type == "secondary"){
  model_data %<>% filter(trophic_binary == "Secondary")
}
if (model_type == "frugivore"){
  model_data %<>% filter(trophic_niche == "Frugivore")
}
if (model_type == "invertivore"){
  model_data %<>% filter(trophic_niche == "Invertivore")
}

# Filter for ecological partition.
if (model_type == "mig"){
  model_data %<>% filter(migration_binary == "Strong")
}
if (model_type == "no_mig"){
  model_data %<>% filter(migration_binary == "Weak")
}
if (model_type == "terr"){
  model_data %<>% filter(territory_binary == "Territory")
}
if (model_type == "no_terr"){
  model_data %<>% filter(territory_binary == "No territory")
}

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
                    #### Run brms models ######

# Tropical brms function.
trop_brms_model <- function(data_set = model_data, response = "sexual_score",  
                            predictor = "trop_non_trop", family = "cumulative"){
  
  # Center categorical predictors.
  data_set %<>%  mutate(
    trop_non_trop = relevel(as.factor(trop_non_trop), ref = "trop"),
    trop_non_trop_c = center_categorical(trop_non_trop))
  
  # Create model formula.
  model_formula <- formula(paste0(response, " ~ ", predictor))
  
  # Add the model family.
  if (family == "cumulative"){
    brms_formula <- brmsformula(model_formula, family = cumulative())
  } else {
    brms_formula <- brmsformula(model_formula, family = bernoulli())
  }
  
  # Priors
  linear_priors <- c(prior(normal(0,1), class="Intercept"),
                     prior(normal(0,1), class="b"))
  # Run brms models.
  brm(
    brms_formula,
    data = data_set,
    prior = linear_priors,
    iter = 10000,
    warmup = 5000,
    chains = 2,
    thin = 20,
    cores = 32,
    init = 0,
    file = model_pathway,
    normalize = FALSE,
    backend = "cmdstanr",
    threads = threading(16)
  )
}


# Create model pathway.
first_half <- "Results/Models/Tropics/"
model_pathway <- paste0(first_half, model_type, "_", data_type, "_", tree_number, ".rds")

