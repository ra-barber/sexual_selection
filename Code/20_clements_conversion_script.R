###############################################################################
        ##### Converting jetz traits to clements phylogeny #####
###############################################################################

# Script to set up to automate trait matching across taxonomies.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(skimr)
library(tictoc)
library(stringr)
library(caper)
library(dplyr)
library(janitor)
library(ggpubr)
library(readxl)


###############################################################################
                   #### Traits and crosswalk ####

# Read in crosswalk.

## Need to use most up to date version of crosswalk.
ebird_tree <- read.tree("Data/Trees/ebird_24may.tre")
ebird_species <- ebird_tree$tip.label %>% str_replace("_", " ") %>% as.data.frame()

# Most up to date version of the crosswalk.
clements_crosswalk <- read_excel("Data/supplementary_dataset_1.xlsx", sheet = 5) %>% clean_names()

clements_crosswalk <- read_excel("Data/Clements_conversion/jetz_to_clements_crosswalk_07_06.xlsx") %>% clean_names()

# Read in ecotraits.
eco_traits <- read.csv("Data/ecotraits_09_05_without_invalids.csv") %>% clean_names()
eco_traits <- read.csv("Data/sexual_selection_cleaned_01_08.csv") %>% clean_names()  # We have since cleaned this data even more.
eco_traits <- read_excel("Data/birdtree_ecotraits_06_09_2023.xlsx") %>% clean_names()

colnames(eco_traits)[1] <- "scientific_name_bird_tree"

# # # This section you can change for different traits # # #

# Select the traits you want, plus species name.
trait_to_merge <- eco_traits %>% select(scientific_name_bird_tree, head_dichro, 
                                        upper_dichro, under_dichro, tail_dichro, 
                                        wing_dichro, dichro_total)

# Merge together.
crosswalk_traits <- clements_crosswalk %>% left_join(trait_to_merge)


###############################################################################
                    #### One to one matches + new ####


# Filter for lumps.
one_to_one_clean <- crosswalk_traits %>% filter(match_type == "1BT to 1CL")

new_clean <- crosswalk_traits %>% filter(match_type == "New")

################################################################################
                       #### Do splits ####

# List of match types to filter.
unique(crosswalk_traits$match_type)

# Filter for lumps.
one_to_many <- crosswalk_traits %>% filter(match_type == "1BT to many CL")



## Write a loop to lower certainty scores for high cert species.

# 
# # Loop through all ebird species, print out conflicts and save problem species.
# for(birdtree_spec in unique(one_to_many$bird_tree_species)){
#   
#   # Filter for ebird species.
#   loop_spec <- one_to_many %>% filter(bird_tree_species == birdtree_spec)
#   
#   # Nominate cert score.
#   nom_cert <- loop_spec$data_certainty[loop_spec$nominate == "Yes"]
#   
#   # If the nominate sexual certainty score is 1, change it 2 for now due to inferrence.
#   if (nom_cert == 1){
#     one_to_many$data_certainty[one_to_many$bird_tree_species == birdtree_spec & one_to_many$nominate == "No"] <- 2
#   }
#   
# }

# Clean dataset.
one_to_many_clean <- one_to_many



###############################################################################
           #### Create function to check conflicting lumps ####

# List of match types to filter.
unique(crosswalk_traits$match_type)

# Filter for lumps.
many_to_one <- crosswalk_traits %>% filter(match_type == "Many BT to 1CL")

# Function for finding conflicts when lumping jetz species.
get_lump_conflicts <- function(dataset = many_to_one, predictor = "dichro_total"){
  
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

###############################################################################
            #### Check sexual selection lumps ####

# Default values are for sexual selection scores.
conflict_species <- get_lump_conflicts()

# Check certainty scores (undoubtedly will be many more conflicts.)
cert_conflict_species <- get_lump_conflicts(predictor = "data_certainty")

# Combine lists.
conflict_species <- c(score_conflict_species, cert_conflict_species)

# Filter traits for conflict species.
conflict_data <- crosswalk_traits %>% filter(scientific_name_clements %in% conflict_species)


### Green pheasant should be scored as 3 for SS due to conflict.
# Note to remmber for later.

# The jetz species where the nominate has worse data certainty: 
# Catharacta antarctica 
# Catharacta lonnbergi   - better known
#
# Sporophila melanops   - better known 
# Sporophila ruficollis
#
# Laterallus jamaicensis 
# Laterallus tuerosi - better known 
# 
# Lamprotornis chloropterus
# Lamprotornis elisabeth  - better known 
#
#
#
#
#
#
#
#

# Export species to manually change.
write.csv(conflict_data, "Data/Clements_conversion/sexual_selection_lump_conflicts.csv", row.names = FALSE)

# Species to manually change.
manual_changes <- c("Ramphocelus passerinii", "Stercorarius antarcticus")


###############################################################################
                 #### Combine non-problem lumps ####

# Remove manual changes.
many_to_one_clean <- many_to_one %>% filter(!ebird_species %in% manual_changes)

# Just pull out nominate species.
many_to_one_clean <- many_to_one %>% filter(daughter == 0)

###############################################################################
#                ### Combine problem lumps ####
# 
# 
# # Remove manual changes.
# many_to_one_manual <- many_to_one %>% filter(ebird_species %in% manual_changes)
# 
# # Pick this species because it had EPP data.
# manual_1 <- many_to_one_manual %>% filter(ebird_species == "Ramphocelus passerinii" & nominate == "No")
# 
# # Pick this species because it had correct sources, but update sex score to match nominate.
# manual_2 <- many_to_one_manual %>% filter(ebird_species == "Stercorarius antarcticus" & nominate == "No")
# manual_2$sexual_score <- 2
# manual_2$data_certainty <- 1


###############################################################################
       #### Many to many matches (probably will be manual) ####

# Filter for lumps.
many_to_many <- crosswalk_traits %>% filter(match_type == "Many BT to many CL")

# Removing these rows makes the many to many dataset clean.
bird_tree_species_to_remove <- c("Epinecrophylla fjeldsaai", "Copsychus stricklandii", 
                                 "Aratinga rubritorquis", "Zosterops montanus", 
                                 "Zosterops salvadorii")

# Many to many clean.
many_to_many_clean <- many_to_many %>% filter(!scientific_name_bird_tree %in% bird_tree_species_to_remove)

# Check that we haven't lost any ebird species.
n_distinct(many_to_many$scientific_name_clements)
n_distinct(many_to_many_clean$scientific_name_clements)




###############################################################################
                           #### Combine scores. ####

all_ebird_clean <- rbind(new_clean, one_to_one_clean, one_to_many_clean, many_to_one_clean, many_to_many_clean)



write.csv(all_ebird_clean, "Data/Clements_conversion/clements_dichro_scores.csv", row.names = FALSE)



###############################################################################
                             #### END ####
###############################################################################