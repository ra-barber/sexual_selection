###############################################################################
               ##### Testing the ebird crosswalk #####
###############################################################################

# This script does something new. Pretty sick eh.


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


###############################################################################
                       #### Read in crosswalk ####

# Read in the ebird crosswalk.
ebird_crosswalk <- read.csv("Data/crosswalk_BT_EB.csv") %>% clean_names()

# Read in the ebird tree.
ebird_tree <- read.tree("Data/Trees/ebird_24may.tre")
ebird_tree_species <- ebird_tree$tip.label %>% str_replace("_", " ")

# Print unique list of match types.
unique(ebird_crosswalk$match_notes)


###############################################################################
                     #### Extinct ebird species ####

# Pull out extinct species.
extinct_ebird <- ebird_crosswalk %>% filter(match_notes == "Extinct")

# Check to make sure none are in the tree. (all NA so obvs not in the tree.)
intersect(extinct_ebird$ebird_species, ebird_tree_species)

# Check possib extinct.
poss_extinct_ebird <- ebird_crosswalk %>% filter(match_notes == "Poss. Extinct")

# Remove extinct species.
ebird_crosswalk_2 <- ebird_crosswalk %>% filter(match_notes != "Extinct")
ebird_crosswalk_2 <- ebird_crosswalk_2 %>% filter(match_notes != "Poss. Extinct")

# There are extra extinct ebird species included in the tree.
setdiff(ebird_tree_species, ebird_crosswalk$ebird_species)
intersect(poss_extinct_ebird$ebird_species, ebird_tree_species)


###############################################################################
                #### Look at newly described species ####

# Print unique list of match types.
unique(ebird_crosswalk_2$match_notes)

# Pull out newly described species.
new_ebird <- ebird_crosswalk_2 %>% filter(match_notes == "Newly described species")

# Remove newly described birdlife species.
new_ebird %<>% tidyr::drop_na(ebird_species)

# Remove newly described species.
ebird_crosswalk_3 <- ebird_crosswalk_2 %>% filter(!ebird_species %in% new_ebird$ebird_species)



###############################################################################
                  #### Remove 1:1 matches and check ####

# Print unique list of match types.
unique(ebird_crosswalk_3$match_notes)

# Filter for 1 to 1 matches.
one_bt_to_one_eb <-  ebird_crosswalk_3 %>% filter(match_notes == "1BT to 1eBird")

# Pull out non-matches.
non_one_to_one_matches <-  ebird_crosswalk_3 %>% filter(match_notes != "1BT to 1eBird")

# Check overlap.
intersect(one_bt_to_one_eb$bird_tree_species, non_one_to_one_matches$bird_tree_species)
intersect(one_bt_to_one_eb$ebird_species, non_one_to_one_matches$ebird_species)

# Create new ebird crosswalk list.
ebird_crosswalk_4 <-  ebird_crosswalk_3 %>% filter(match_notes != "1BT to 1eBird")


###############################################################################
              #### Pull out 1 1BT to many eBird splits ####

# Print unique list of match types.
unique(ebird_crosswalk_4$match_notes)

# Filter for 1 to 1 matches.
one_bt_to_many_eb <-  ebird_crosswalk_4 %>% filter(match_notes == "1BT to many eBird")

# Pull out non-matches.
none_matches <-  ebird_crosswalk_4 %>% filter(match_notes != "1BT to many eBird")

# Check overlap.
intersect(one_bt_to_many_eb$bird_tree_species, none_matches$bird_tree_species)
intersect(one_bt_to_many_eb$ebird_species, none_matches$ebird_species)



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


