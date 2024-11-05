###############################################################################
                ##### Converting BirdTree to eBird #####
###############################################################################

# Script to facilitate converting traits from outdated BirdTree taxonomy to more
# recent eBird (Clements taxonomy).

# This script will be updated over time when a consistent eBird phylogeny is 
# published. Please refer to the github url for the most recent version:

# https://github.com/ra-barber/sexual_selection


# In this example script, I use trophic level as a trait to convert. This is 
# because it fairly stable between sub-species, so should be the simplest 
# example using species-level data from S1 Data. 

# The process to convert traits is partitioned into different steps depending if
# species are easy 1:1 matches, or instances when birdtree species are split or 
# lumped. The conversion assumptions and methods are noted in more detail in the
# methods section of S2 Text. Fundamentally, the conversion is based on:

## Splits ##
# For splits, traits are simply duplicated. As more datasets are published
# for eBird species, this step will become increasingly inaccurate. It will work
# best with conserved traits (such as trophic level), but will be inaccurate for 
# plastic traits that vary geographically (i.e. body mass or dichromatism).

## Lumps ##
# For lumps, traits are simply averaged. This works well for categorical traits,
# but will be dodgy for continuous traits. This is because range sizes and 
# populations will be vastly different between merged species. One option is to 
# to take means, weighted by range size. << This step will be added on github.

## Many to many ##
# Because the bird tree phylogeny is quite old, some species have been lumped
# and subsequently split again. As the ebird taxonomy is updated, 
# this step will become longer. The best option is to manually assess these 
# species for each trait to make sure it makes sense. 
  

# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(stringr)
library(dplyr)
library(janitor)
library(readxl)


###############################################################################
                   #### Traits and cross walk ####

# Read in the most recent version of the cross walk.
clements_crosswalk <- read_excel("Data/S1 Data.xlsx", sheet = 5) %>% clean_names()

# Read in Bird 
birdtree_data <- read_excel("Data/S1 Data.xlsx", sheet = 2) %>% clean_names()

# Select the traits you want, plus species name.
trait_to_merge <- birdtree_data %>% select(scientific_name_bird_tree, trophic_level)

# Merge together.
crosswalk_traits <- clements_crosswalk %>% left_join(trait_to_merge)


###############################################################################
                       #### 1:1 matches + New ####
 
# These examples are the easiest. 

# Filter for the 1:1 matches. This is the simplest option.
one_to_one_clean <- crosswalk_traits %>% filter(match_type == "1BT to 1CL")

# The new species to science won't have traits to convert currently. They should
# be done manually. When an eBird phylogeny is published, one 
# option will be to use the nearest neighbour for new species until original 
# data becomes available.
new_clean <- crosswalk_traits %>% filter(match_type == "New")


################################################################################
                        #### Birdtree Splits ####

# List of match types to filter.
unique(crosswalk_traits$match_type)

# Filter for splits. Data will be simply duplicated for these species.
one_to_many <- crosswalk_traits %>% filter(match_type == "1BT to many CL")

# Clean dataset.
one_to_many_clean <- one_to_many


################################################################################
                       #### Birdtree lumps ####

# Lumps are harder. For categorical data, we need to check when two Birdtree 
# species have conflicting life history traits.

# Filter for lumps.
many_to_one <- crosswalk_traits %>% filter(match_type == "Many BT to 1CL")

# Function for finding conflicts when lumping birdtree species.
get_lump_conflicts <- function(dataset = many_to_one, predictor = "trophic_level"){
  
  # Empty vector to save problem species.
  problem_lumps <- c()
  
  # Loop through all ebird species, print out conflicts and save problem species.
  for(ebird_spec in unique(dataset$scientific_name_clements)){
    
    # Filter for ebird species.
    loop_spec <- dataset %>% filter(scientific_name_clements == ebird_spec)
    
    # Pull out trait column.
    trait_column <- loop_spec %>% pull(predictor)
    
    # If there are conflicts in the trait, 
    if (n_distinct(trait_column) > 1){
      print(loop_spec)
      problem_lumps <- c(problem_lumps, ebird_spec)
    }
    
  }
  return(problem_lumps)
}


# Default values are for sexual selection scores.
conflict_species <- get_lump_conflicts()

# Filter traits for conflict species.
conflict_data <- crosswalk_traits %>% 
  filter(scientific_name_clements %in% conflict_species)

# In this example, we will take the nominate species for simplicity. Normally, 
# daughter birdtree species are from the smaller range, and tend to display 
# behaviours which are not representative of the larger population. 

# However!  

# This is not always the case, and any conflict species should be manually checked 
# when converting new traits from BirdTree to Clements. In reality, few species
# should appear in this section. But once traits have been converted, they will
# probably be duplicated by other researchers without checking, so this step is 
# important! 

# Just pull out nominate species in this example.
many_to_one_clean <- many_to_one %>% filter(daughter == 0)


###############################################################################
                  #### Many to many matches  ####

# This are instances when multiple lumps and splits have occurred. For categorical
# traits, we can remove the complex birdtree species. This turns them into 
# regular splits.

# Filter for instances matching complex species complexes. 
many_to_many <- crosswalk_traits %>% filter(match_type == "Many BT to many CL")

# Removing these rows makes the many to many dataset into normal splits.
bird_tree_species_to_remove <- c("Aratinga rubritorquis",
                                 "Copsychus stricklandii", 
                                 "Epinecrophylla fjeldsaai", 
                                  "Zosterops montanus", 
                                 "Zosterops salvadorii")

# Many to many clean.
many_to_many_clean <- many_to_many %>% 
  filter(!scientific_name_bird_tree %in% bird_tree_species_to_remove)


###############################################################################
               #### Combine sections and export ####

# Combine all the different steps.
all_ebird_clean <- rbind(new_clean, one_to_one_clean, one_to_many_clean, 
                         many_to_one_clean, many_to_many_clean)

# Custom pathway for exporting.
pathway <- "Data/new_clements_trait.csv"

# Export.
write.csv(all_ebird_clean, pathway, row.names = FALSE)


###############################################################################
                             #### END ####
###############################################################################