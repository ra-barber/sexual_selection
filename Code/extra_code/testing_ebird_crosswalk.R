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

# Read in some data.
ebird_crosswalk <- read.csv("Data/crosswalk_BT_EB.csv") %>% clean_names()

# Read in a tree.
ebird_tree <- read.tree("Data/Trees/ebird_24may.tre")
ebird_tree_species <- ebird_tree$tip.label %>% str_replace("_", " ")

###############################################################################
                #### Remove 1:1 matches and check ####

ebird_crosswalk %>% count(match_notes)


bt1_eb1 <-  ebird_crosswalk %>% filter(match_notes == "1BT to 1eBird")
nbt1_eb1 <-  ebird_crosswalk %>% filter(match_notes != "1BT to 1eBird")

intersect(bt1_eb1$bird_tree_species, nbt1_eb1$bird_tree_species)

ebird_dupes <- ebird_crosswalk %>% get_dupes(ebird_species)


ebird_dupes %>% count(match_notes)

###############################################################################
                           #### Extinct ebird species ####


extinct_ebird <-  ebird_crosswalk %>% filter(match_notes == "Extinct")

intersect(extinct_ebird$ebird_species, ebird_tree_species)
intersect("Acrocephalus luscinius", ebird_tree_species)


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


