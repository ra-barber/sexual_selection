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
                             #### Data ####

# Read in crosswalk.
clements_crosswalk <- read_excel("Data/jetz_to_clements_crosswalk_07_06.xlsx") %>% clean_names()

# Read in ecotraits.
eco_traits <- read.csv("Data/ecotraits_09_05_without_invalids.csv") %>% clean_names()
colnames(eco_traits)[1] <- "bird_tree_species"




###############################################################################
                 #### Pull out specific trait ####

trait_to_merge <- eco_traits %>% select(bird_tree_species, sexual_score)

crosswalk_traits <- clements_crosswalk %>% left_join(trait_to_merge)


###############################################################################
                  #### Check conflicting lumps ####
unique(crosswalk_traits$match_type)
many_to_one <- crosswalk_traits %>% filter(match_type == "Many BT to 1eBird")

for(ebird_spec in unique(many_to_one$ebird_species)){
  loop_spec <- many_to_one %>% filter(ebird_species == ebird_spec)
  
  if (n_distinct(loop_spec$sexual_score) > 1){
    print(loop_spec)
  }
  
}


###############################################################################
                           #### Section 3 ####


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


