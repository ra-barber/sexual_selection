dyn.load("C:/Windows/System32/OpenCL.dll", FALSE, TRUE)

install.packages("OpenCL", INSTALL_opts = "--no-clean-on-error --no-test-load")

library(OpenCL)


ggplot(phylo_data, aes(x=body_mass_log, y=jitter(sexual_score), group = sexual_score)) + 
  geom_density(alpha=0.1) + ylab("Body Mass (Logged)") + theme_classic()


breaks <- quantile(phylo_data$body_mass_log, seq(0,1,length.out=11), na.rm=TRUE)

# Cut the cell values into bins using the breaks.
phylo_data$cuts <- cut(phylo_data$body_mass_log, breaks = breaks, include.lowest = T)

phylo_data %>% group_by(cuts) %>% summarise(sex_mean = mean(sexu))



time_taken <- read.csv("Data/all_2000_time_taken.csv")

str(time_taken)


plot(log(elapsed) ~ jitter(thread_number), data = time_taken)
library(magrittr)
lm(elapsed ~ thread_number, data = time_taken) %>% summary()

library(ggplot2)
library(tidyr)
library(dplyr)
ggplot(time_taken %>% filter(elapsed < 125000), aes(x = thread_number, y = elapsed)) + geom_point() + geom_smooth(method = "lm")

boxplot(elapsed ~ normalize, data = time_taken)
boxplot(elapsed ~ randomise, data = time_taken)


lm(elapsed ~ (normalize + thread_number + randomise + priors)^2 , data = time_taken) %>% summary()
library(brms)

brm(pri)

density <- read.csv("Data/density_estimates.csv")


bird_tree_dens <- density %>% group_by(birdtree_sciname) %>% summarise(mean_dens = mean(density))
colnames(bird_tree_dens)[1] <- "bird_tree_name"

model_data <- read.csv("Data/sexual_traits.csv")

dens_test <- left_join(model_data, bird_tree_dens)

boxplot(mean_dens ~ sexual_score, data = dens_test)

mod1 <- lm(sexual_score ~ mean_dens*body_mass_log, data = dens_test)

library(interactions)

interact_plot(mod1, mean_dens, body_mass_log)

flock_size <- read.csv("Data/flock_size.csv")

head(flock_size)

flock_species <- flock_size %>% group_by(SCIENTIFIC_NAME) %>% 
  summarise(max_flock = max(median_abund))

colnames(flock_species)[1] <- "bird_tree_name"

flock_test <-  left_join(model_data, flock_species)

boxplot(log(max_flock) ~ sexual_score, data = flock_test)

flock_test$log_flock <- log(flock_test$max_flock)

mod2 <- lm(sexual_score ~ log_flock*body_mass_log, data = flock_test)

interact_plot(mod2, log_flock, body_mass_log)

summary(mod1)


summary(mod2)

mod3 <- lm(sexual_score ~ body_mass_log, data = model_data)

summary(mod3)

library(phylolm)

mod4 <- phylolm(sexual_score ~ body_mass_log, data = model_data, phy = model_tree)
summary(mod4)

rownames(flock_test) <- flock_test$tree_tip

mod6 <- phylolm(sexual_score ~ log_flock, data = flock_test, phy = model_tree)
summary(mod6)

colnames(model_data)
head(gen_data)
colnames(gen_data)[7] <- "bird_tree_name"

gen_test <- left_join(model_data, gen_data)

rownames(gen_test) <- gen_test$tree_tip

mod7 <- phylolm(sexual_score ~ log(GenLength), data = gen_test, phy = model_tree)
mod8 <- phylolm(sexual_score ~ Adult.survival, data = gen_test, phy = model_tree)
mod9 <- phylolm(sexual_score ~ log(Maximum.longevity), data = gen_test, phy = model_tree)

summary(mod7)
summary(mod8)
summary(mod9)

hist(gen_test$Adult.survival)

plot(Adult.survival ~ log(Maximum.longevity), data = gen_test)
colnames(gen_test)

cor(gen_data$Adult.survival, log(gen_data$Maximum.longevity))

mod10 <- lm(sexual_score ~ log(GenLength), data = gen_test)
mod11 <- lm(sexual_score ~ Adult.survival, data = gen_test)
mod12 <- lm(sexual_score ~ log(Maximum.longevity), data = gen_test)


summary(mod10)
summary(mod11)
summary(mod12)


mod13 <- phylolm(sexual_score ~ log(Maximum.longevity)*Adult.survival, data = gen_test, phy = model_tree)

summary(mod13)




plot(jitter(sexual_score, 2) ~ log(GenLength), data = gen_test)
plot(sexual_score ~ log(GenLength), data = gen_test)
max(na.omit(gen_test$GenLength))
lines(predict(mod10, newdata = data.frame(GenLength = seq(from = 1, to = 28, by = 0.1))))
plot(jitter(sexual_score, 2) ~ Adult.survival, data = gen_test)
plot(sexual_score ~ Adult.survival, data = gen_test)

lines(predict(mod11))

plot(jitter(sexual_score, 2) ~ log(Maximum.longevity), data = gen_test)