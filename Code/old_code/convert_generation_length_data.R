###############################################################################
             ##### Convert generation length to jetz #####
###############################################################################

# This script converts generation data to jetz phylogeny.


# Clean the environment.
rm(list=ls())

# Load packages.
library(magrittr)
library(tictoc)
library(stringr)
library(caper)
library(dplyr)
library(janitor)
library(ggpubr)
library(skimr)

# Read in the functions.
source("Code/functions.R")

###############################################################################
                          #### Data #####

# Read in the life history traits.
bird_traits <- read.csv("Data/GBD_BirdTreeTaxo_ecotraits_2022_November_7th.csv") %>% clean_names()

# Read in the generation data.
generation_data <- read.csv("Data/generation_length_birdlife.csv") %>% clean_names()
colnames(generation_data)[7] <- "birdlife_name"

# Read in the crosswalk.
crosswalk <- read.csv("Data/birdlife_jetz_crosswalk.csv") %>% clean_names()
colnames(crosswalk)[1:2] <- c("birdlife_name", "bird_tree_name")

# Now look at any NA species.
crosswalk %>% filter(is.na(birdlife_name))
crosswalk %>% filter(is.na(bird_tree_name))
crosswalk %<>% na.omit()
crosswalk %<>% filter(match_type != "")



###############################################################################
                    #### Pull out 1:1 matches #####


# Pull out all the species that are 1:1 match.
one_to_one <- crosswalk %>% filter(match_type == "1BL to 1BT") %>% 
  select(birdlife_name, bird_tree_name)

# Potentially make sure none of them are in the other categories.

# Filter the generation data.
one_to_one_gen <- generation_data %>% filter(birdlife_name %in% one_to_one$birdlife_name)

# Join the datasets.
one_to_one_both <- left_join(one_to_one, one_to_one_gen)

# Select just the columns we need.
one_to_one_both %<>% dplyr::select(bird_tree_name, birdlife_name, adult_survival, age_at_first_breeding, 
                                   maximum_longevity, z, gen_fls, gen_fl, gen_fs, gen_length)



###############################################################################
                 #### Pull out many to 1 BT #####

# Update the crosswalk.
crosswalk_2 <- crosswalk %>% filter(!birdlife_name %in% one_to_one_both$birdlife_name)

# Pull out all the species that are many:1 match.
many_to_one <- crosswalk_2 %>% filter(match_type == "Many BL to 1BT") %>% 
  select(birdlife_name, bird_tree_name)

# Filter the generation data.
many_to_one_gen <- generation_data %>% filter(birdlife_name %in% many_to_one$birdlife_name)

# Join the datasets.
many_to_one_both <- left_join(many_to_one, many_to_one_gen)

# Get the averages.
many_to_one_both %<>% group_by(bird_tree_name) %>% 
  summarise(birdlife_name = first(birdlife_name),
            adult_survival = mean(adult_survival),
            age_at_first_breeding = mean(age_at_first_breeding),
            maximum_longevity = mean(maximum_longevity),
            z = mean(z),
            gen_fls = mean(gen_fls),
            gen_fl = mean(gen_fl),
            gen_fs = mean(gen_fs),
            gen_length = mean(gen_length))

# Looking good so far.
rbind(one_to_one_both, many_to_one_both) %>% pull(bird_tree_name) %>% n_distinct()
rbind(one_to_one_both, many_to_one_both) %>% get_dupes(bird_tree_name)


###############################################################################
                  #### Pull out 1 to many BT #####

# Update the crosswalk.
crosswalk_3 <- crosswalk_2 %>% filter(!birdlife_name %in% many_to_one_both$birdlife_name)

# Pull out all the species that are many:1 match.
one_to_many <- crosswalk_3 %>% filter(match_type == "1BL to many BT") %>% 
  select(birdlife_name, bird_tree_name)

# Filter the generation data.
one_to_many_gen <- generation_data %>% filter(birdlife_name %in% one_to_many$birdlife_name)

# Join the datasets.
one_to_many_both <- left_join(one_to_many, one_to_many_gen)

# Select just the traits we want.
one_to_many_both %<>% dplyr::select(bird_tree_name, birdlife_name, adult_survival, age_at_first_breeding, 
                            maximum_longevity, z, gen_fls, gen_fl, gen_fs, gen_length)

# Join all together.
joined_dataset <- rbind(one_to_one_both, many_to_one_both, one_to_many_both)

# Check number of distinct.
n_distinct(joined_dataset$bird_tree_name)

# Export the new data.
write.csv(joined_dataset, "Data/jetz_generation_length.csv", row.names = FALSE)

# Look at the duplicate species we need to manually edit in excel.
dupe_species <- get_dupes(joined_dataset, bird_tree_name)


###############################################################################
        #### Fill in the gaps where there was mismatch #####

# Read in the clean data.
clean_gen_data <- read.csv("Data/jetz_generation_length_clean.csv")

# Check we cleaned it.
get_dupes(clean_gen_data, bird_tree_name)

# See the NAs left.
clean_gen_data %>% filter(is.na(gen_length))

# Check if we're missing any species.
sexual_traits <- read.csv("Data/sexual_traits.csv")
eco_traits <- read.csv("Data/GBD_BirdTreeTaxo_ecotraits_2022_November_7th.csv") %>% clean_names()


clean_gen_data %>% filter(bird_tree_name %in% sexual_traits$bird_tree_name)
setdiff(clean_gen_data$bird_tree_name, sexual_traits$bird_tree_name)
setdiff(sexual_traits$bird_tree_name, clean_gen_data$bird_tree_name)
setdiff(eco_traits$bird_tree_name, clean_gen_data$bird_tree_name)

# In the generation data, but using the old jetz name under birdlife taxo.
# Excellens has since been merged so use same data for both.
"Campylopterus curvipennis" 
"Campylopterus excellens"  

# Invalid taxon so ditch.
"Phyllastrephus leucolepis"

# Subspecies of Polioptila guianensis according to birdlife.
"Polioptila clementsi" 

# Now check which missing birdtree species are in the gen data under birdlife.
NA_data <- clean_gen_data %>% filter(is.na(gen_length))

# Replace those species.
NA_data_spec_wrong <- NA_data %>% filter(bird_tree_name %in% generation_data$birdlife_name)
match_fixed_gen_data <- left_join(NA_data_spec_wrong, generation_data, by = c("bird_tree_name" = "birdlife_name"))
match_fixed_gen_data %<>% select(bird_tree_name, birdlife_name, adult_survival.y,
                                age_at_first_breeding.y, maximum_longevity.y, z.y,
                                gen_fls.y, gen_fl.y, gen_fs.y, gen_length.y)
colnames(match_fixed_gen_data) <- colnames(clean_gen_data)

# Replace the NA species with the correct ones.
clean_gen_data_2 <- clean_gen_data %>% filter(!bird_tree_name %in% match_fixed_gen_data$bird_tree_name)

# Merge the new species back in.
clean_gen_data_2 <- rbind(match_fixed_gen_data, clean_gen_data_2)

###############################################################################
               #### Fill in the gaps leftover  #####


# See the NAs left.
clean_gen_data_2 %>% filter(is.na(gen_length)) %>% pull(bird_tree_name)

setdiff(sexual_traits$bird_tree_name, clean_gen_data_2$bird_tree_name)
setdiff(eco_traits$bird_tree_name, clean_gen_data_2$bird_tree_name)

subset_names <- colnames(clean_gen_data_2)[3:10]

# Alophoixus affinis is under Thapsinillas affinis
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Alophoixus affinis", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Thapsinillas affinis", subset_names]

#  "Aramides cajanea" is under = Aramides cajaneus   
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Aramides cajanea", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Aramides cajaneus", subset_names]

# "Claravis godefrida" is under = Claravis geoffroyi
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Claravis godefrida", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Claravis geoffroyi", subset_names]

# "Damophila julie" is under = Juliamyia julie       
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Damophila julie", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Juliamyia julie", subset_names]

# "Dryocopus galeatus"  is under = Hylatomus galeatus     
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Dryocopus galeatus", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Hylatomus galeatus", subset_names]

# "Hirundo fuliginosa"  is under = Petrochelidon fuliginosa     
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Hirundo fuliginosa", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Petrochelidon fuliginosa", subset_names]

# "Hypogramma hypogrammicum" is under = Arachnothera hypogrammica
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Hypogramma hypogrammicum", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Arachnothera hypogrammica", subset_names]

# "Hypsipetes thompsoni" is under = Cerasophila thompsoni
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Hypsipetes thompsoni", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Cerasophila thompsoni", subset_names]

# "Meliphaga gracilis" is under = Microptilotis gracilis
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Meliphaga gracilis", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Microptilotis gracilis", subset_names]

# "Myrmeciza ferruginea" is under = Myrmoderus ferrugineus
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Myrmeciza ferruginea", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Myrmoderus ferrugineus", subset_names]

# "Parus montanus" is under = Poecile montanus     
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Parus montanus", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Poecile montanus", subset_names]

# "Porzana flaviventer" is under = Hapalocrex flaviventer
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Porzana flaviventer", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Hapalocrex flaviventer", subset_names]

# "Psittacula calthropae"   is under = Psittacula calthrapae     
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Psittacula calthropae", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Psittacula calthrapae", subset_names]

# "Strophocincla cachinnans" is under = Trochalopteron cachinnans
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Strophocincla cachinnans", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Trochalopteron cachinnans", subset_names]

# "Strophocincla fairbanki" is under = Trochalopteron fairbanki
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Strophocincla fairbanki", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Trochalopteron fairbanki", subset_names]

# "Suiriri islerorum" is under = Suiriri affinis
clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Suiriri islerorum", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Suiriri affinis", subset_names]


# "Haematopus finschi" appears to be missing but can take gen average.
haem_genus <- grep("Haematopus", clean_gen_data_2$bird_tree_name, value = TRUE)

gen_averages <- clean_gen_data_2 %>% filter(bird_tree_name %in% haem_genus) %>% 
  select(-c(bird_tree_name, birdlife_name)) %>% colMeans(na.rm = TRUE)

clean_gen_data_2[clean_gen_data_2$bird_tree_name == "Haematopus finschi", subset_names] <- gen_averages

###############################################################################
                    #### Final checks  #####

clean_gen_data_2 %>% filter(is.na(gen_length))
clean_gen_data_2 %>% get_dupes(bird_tree_name)

setdiff(sexual_traits$bird_tree_name, clean_gen_data_2$bird_tree_name)
setdiff(eco_traits$bird_tree_name, clean_gen_data_2$bird_tree_name)

# In the generation data, but using the old jetz name under birdlife taxo.
# Excellens has since been merged so use same data for both.
"Campylopterus curvipennis" 
"Campylopterus excellens"

# Subspecies of Polioptila guianensis according to birdlife.
"Polioptila clementsi"

new_data <- data.frame(bird_tree_name = c("Campylopterus curvipennis", "Campylopterus excellens", "Polioptila clementsi"), 
                       birdlife_name = c("Campylopterus curvipennis", "Campylopterus curvipennis", "Polioptila guianensis"), 
                       adult_survival = NA, 
                       age_at_first_breeding = NA, 
                       maximum_longevity = NA, 
                       z = NA, 
                       gen_fls = NA, 
                       gen_fl = NA, 
                       gen_fs = NA, 
                       gen_length = NA)

new_data[new_data$birdlife_name == "Campylopterus curvipennis", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Campylopterus curvipennis", subset_names]

new_data[new_data$birdlife_name == "Polioptila guianensis", subset_names] <- 
  generation_data[generation_data$birdlife_name == "Polioptila guianensis", subset_names]

 

clean_gen_data_3 <- rbind(clean_gen_data_2, new_data)

###############################################################################
                      #### Final final checks  #####

clean_gen_data_3 %>% filter(is.na(gen_length))
clean_gen_data_3 %>% get_dupes(bird_tree_name)

setdiff(sexual_traits$bird_tree_name, clean_gen_data_3$bird_tree_name)
setdiff(eco_traits$bird_tree_name, clean_gen_data_3$bird_tree_name)

# Fucking nailed it.
write.csv(clean_gen_data_3, "Data/jetz_generation_length_clean_05_01.csv", row.names = FALSE)













###############################################################################
                            ##### END ######