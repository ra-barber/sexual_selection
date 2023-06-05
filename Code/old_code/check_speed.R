###############################################################################
                  ##### Check BRMS output speeds #####
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

# Read in the functions. 
source("Code/functions.R")


###############################################################################
                             #### Data ####

eight_threads_outputs <- list.files(path = c("Z:/home/sexual_selection/Jobs/feb/eight_threads/"), full.names = T, include.dirs = FALSE, recursive = FALSE)
sixteen_threads_outputs <- list.files(path = c("Z:/home/sexual_selection/Jobs/feb/sixteen_threads/"), full.names = T, include.dirs = FALSE, recursive = FALSE)

benchmarking <- read.csv("Z:/home/sexual_selection/thread_benchmarking.csv")
adapt_benchmarking <- read.csv("Z:/home/sexual_selection/testing_adapt_delta.csv")
treedepth_benchmarking <- read.csv("Z:/home/sexual_selection/testing_treedepth.csv")
stepsize_benchmarking <- read.csv("Z:/home/sexual_selection/testing_stepsize.csv")

###############################################################################
                           #### Set up array info ####

# Create types for hpc jobs.
tree_number <- 1:15

# Model
model_type <- c("npp", "no_npp", "no_inter")

# Centered or uncentered.
center <- c("centered", "uncentered")

# Set the data types.
#data_type <- c("all", "high")
data_type <- c("all", "high")

# Expand the grid.
all_combos <- expand.grid(tree_number, model_type, center, data_type)

###############################################################################
               #### Add time finished for eight threads ####

all_combos$eight_threads <- NA
#x <- 104
# Loop through array numbers.
for (x in 1:length(eight_threads_outputs)){
  # Pull out array number.
  array_number <- str_extract(eight_threads_outputs[x], "(\\d+)\\.Rout$") %>% str_extract("\\d+") %>% as.numeric()
  # Read  in lines.
  output_file <- readLines(eight_threads_outputs[x])
  
  # Check it finished and if so add time to dataframe.
if(!str_detect(last(output_file), "Chain")){
  time_elapsed <- last(output_file) %>% str_split(pattern = " ", simplify = TRUE)
  time_elapsed <- time_elapsed[-c(length(time_elapsed))] %>% last() %>% as.numeric()
  #time_elapsed <- time_elapsed[[2]] %>% as.numeric()
  all_combos$eight_threads[array_number]  <- time_elapsed/60/60
} else {
  # Read  in the last line.
  chain_1_match <- output_file %>% str_detect("Chain 1 Iter")
  chain_1_last <- last(output_file[chain_1_match])
  
  chain_2_match <- output_file %>% str_detect("Chain 2 Iter")
  chain_2_last <- last(output_file[chain_2_match])
  
  chain_1_percent <- chain_1_last %>% str_extract("[0-9]+%") %>% str_extract("[0-9]+") %>% as.numeric()/100
  chain_2_percent <- chain_2_last %>% str_extract("[0-9]+%") %>% str_extract("[0-9]+") %>% as.numeric()/100
  time_elapsed <- 71.75/min(c(chain_1_percent, chain_2_percent))
  
  
  # percent_done <- last(output_file) %>% str_extract("[1-9]+%") %>% str_extract("[1-9]+") %>% as.numeric()/100
  # time_elapsed <- 72/percent_done
  all_combos$eight_threads[array_number]  <- time_elapsed
}
}


###############################################################################
            #### Add time finished for sixteen threads ####

all_combos$sixteen_threads <- NA
x <- 33
# Loop through array numbers.
for (x in 1:length(sixteen_threads_outputs)){
  # Pull out array number.
  array_number <- str_extract(sixteen_threads_outputs[x], "(\\d+)\\.Rout$") %>% str_extract("\\d+") %>% as.numeric()
  # Read  in lines.
  output_file <- readLines(sixteen_threads_outputs[x])
  
  # Check it finished and if so add time to dataframe.
  if(!str_detect(last(output_file), "Chain")){
    time_elapsed <- last(output_file) %>% str_split(pattern = " ", simplify = TRUE)
    time_elapsed <- time_elapsed[-c(length(time_elapsed))] %>% last() %>% as.numeric()
    #time_elapsed <- time_elapsed[[2]] %>% as.numeric()
    all_combos$sixteen_threads[array_number]  <- time_elapsed/60/60
  } else {
    # Read  in the last line.
    chain_1_match <- output_file %>% str_detect("Chain 1 Iter")
    chain_1_last <- last(output_file[chain_1_match])
    
    chain_2_match <- output_file %>% str_detect("Chain 2 Iter")
    chain_2_last <- last(output_file[chain_2_match])
    
    chain_1_percent <- chain_1_last %>% str_extract("[0-9]+%") %>% str_extract("[0-9]+") %>% as.numeric()/100
    chain_2_percent <- chain_2_last %>% str_extract("[0-9]+%") %>% str_extract("[0-9]+") %>% as.numeric()/100
    time_elapsed <- 71.75/min(c(chain_1_percent, chain_2_percent))
    
    
    # percent_done <- last(output_file) %>% str_extract("[1-9]+%") %>% str_extract("[1-9]+") %>% as.numeric()/100
    # time_elapsed <- 72/percent_done
    all_combos$sixteen_threads[array_number]  <- time_elapsed
  }
}

both_types <- all_combos[all_combos$sixteen_threads > 1,] %>% na.omit()
library(tidyr)
both_types %<>%  pivot_longer(cols = c("eight_threads", "sixteen_threads"), names_to = "threads", values_to = "time_elapsed")


###############################################################################
                           #### Plot the results ####

colnames(both_types)[1:4] <- c("tree_number", "model", "center", "data_type")
colnames(all_combos)[1:4] <- c("tree_number", "model", "center", "data_type")


ggplot(both_types, aes(x = threads, y = log(time_elapsed), col = data_type)) + geom_boxplot()
ggplot(all_combos, aes(x = as.factor(tree_number),
                       y = sixteen_threads, col = data_type)) + 
  geom_boxplot() + geom_hline(yintercept = 72)


###############################################################################
                           #### Section 4 ####

# Try for models with new gamma prior.
ten_threads_outputs <- list.files(path = c("Z:/home/sexual_selection/Jobs/feb/ten_threads/"), full.names = T, include.dirs = FALSE, recursive = FALSE)

# Create types for hpc jobs.
tree_number <- 1:15

# Model
model_type <- c("chicksqrt", "logistic", "logistic_chicksqrt", "equidistant",
                "no_terr", "no_terr_chicksqrt" ) #"npp", "no_npp", "no_inter")

# Centered or uncentered.
center <- c("centered", "uncentered")

# Set the data types.
#data_type <- c("all", "high")
data_type <- c("all", "high")

# Expand the grid.
ten_combos <- expand.grid(model_type, center, data_type, tree_number)


ten_combos$ten_threads <- NA

# Loop through array numbers.
for (x in 1:length(ten_threads_outputs)){
  # Pull out array number.
  array_number <- str_extract(ten_threads_outputs[x], "(\\d+)\\.Rout$") %>% str_extract("\\d+") %>% as.numeric()
  # Read  in lines.
  output_file <- readLines(ten_threads_outputs[x])
  
  # Check it finished and if so add time to dataframe.
  if(!str_detect(last(output_file), "Chain")){
    time_elapsed <- last(output_file) %>% str_split(pattern = " ", simplify = TRUE)
    time_elapsed <- time_elapsed[-c(length(time_elapsed))] %>% last() %>% as.numeric()
    #time_elapsed <- time_elapsed[[2]] %>% as.numeric()
    ten_combos$ten_threads[array_number]  <- time_elapsed/60/60
  } else {
    # Read  in the last line.
    chain_1_match <- output_file %>% str_detect("Chain 1 Iter")
    chain_1_last <- last(output_file[chain_1_match])
    
    chain_2_match <- output_file %>% str_detect("Chain 2 Iter")
    chain_2_last <- last(output_file[chain_2_match])
    
    chain_1_percent <- chain_1_last %>% str_extract("[0-9]+%") %>% str_extract("[0-9]+") %>% as.numeric()/100
    chain_2_percent <- chain_2_last %>% str_extract("[0-9]+%") %>% str_extract("[0-9]+") %>% as.numeric()/100
    time_elapsed <- 21/min(c(chain_1_percent, chain_2_percent))
    
    
    # percent_done <- last(output_file) %>% str_extract("[1-9]+%") %>% str_extract("[1-9]+") %>% as.numeric()/100
    # time_elapsed <- 72/percent_done
    ten_combos$ten_threads[array_number]  <- time_elapsed
  }
}

###############################################################################
                           #### Benchmarking experiment ####

benchmarking %>% count(thread_number)

benchmarking %>% filter(thread_number < 100) %>% 
  ggplot(aes(x = thread_number, y =elapsed)) + geom_jitter() + geom_smooth()



benchmarking %>% 
  ggplot(aes(x = as.factor(thread_number), y =elapsed/60/60)) + 
  geom_boxplot() + geom_smooth(method = "lm") + geom_hline(yintercept = 1.44)


###############################################################################
                           #### Section 6 ####

adapt_benchmarking %>% count(adapt_delta)

adapt_benchmarking %>% distinct() %>% na.omit() %>%  count(adapt_delta)

adapt_benchmarking %>% distinct() %>% na.omit()  %>% 
  ggplot(aes(x = as.factor(adapt_delta), y =elapsed)) + geom_boxplot() + geom_jitter(height = 0)

adapt_benchmarking %>% distinct() %>% na.omit()  %>% 
  ggplot(aes(x = as.factor(adapt_delta), y =log(elapsed))) + geom_boxplot() + geom_jitter(height = 0)

adapt_benchmarking %>% distinct() %>% na.omit()  %>% 
  ggplot(aes(fill = as.factor(adapt_delta), col = as.factor(adapt_delta), x =elapsed)) + geom_density(alpha = 0.3)

adapt_benchmarking %>% distinct() %>% na.omit()  %>% 
  ggplot(aes(fill = as.factor(adapt_delta), col = as.factor(adapt_delta), x =log(elapsed))) + geom_density(alpha = 0.3)

adapt_benchmarking %>% distinct() %>% na.omit()  %>% 
  ggplot(aes(x = as.factor(adapt_delta), y =elapsed)) + geom_boxplot() + geom_jitter(height = 0) + expand_limits(y = 0)

stepsize_benchmarking %>% distinct() %>% na.omit()  %>% 
  ggplot(aes(x = as.factor(step_size), y =elapsed)) + geom_boxplot() + geom_jitter(height = 0)+ expand_limits(y = 0)

treedepth_benchmarking %>% distinct() %>% na.omit()  %>% 
  ggplot(aes(x = as.factor(max_treedepth), y =elapsed)) + geom_boxplot() + geom_jitter(height = 0)+ expand_limits(y = 0)

treedepth_benchmarking %>% filter(max_treedepth == 5) %>% pull(elapsed) %>% mean()
treedepth_benchmarking %>% filter(max_treedepth == 10) %>% pull(elapsed) %>% mean()

adapt_benchmarking %>% 
  ggplot(aes(x = as.factor(thread_number), y =elapsed/60/60)) + 
  geom_boxplot() + geom_smooth(method = "lm") + geom_hline(yintercept = 1.44)
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


