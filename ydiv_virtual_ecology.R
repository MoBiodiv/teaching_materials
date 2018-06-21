# load packages
rm(list=ls())
library(mobr)
library(mobsim)
library(tidyverse)

# set regional S, N and the 'evenness' of the regional SAD
S_regional <- 50
N_regional <- 2000

# shape of lognormal SAD can be parameterised using the cv = sd(N)/mean(N)
cv_abund <- list(list(cv_abund = 1), # cv is negatively correlated with evenness
                 list(cv_abund = 4))

# generate regional communities
maps = NULL
maps$agg = lapply(cv_abund, function(x) 
  sim_thomas_community(S_regional, N_regional, 'lnorm', x, fix_s_sim = F))
maps$poi = lapply(cv_abund, function(x) 
  sim_poisson_community(S_regional, N_regional, 'lnorm', x, fix_s_sim = F))

# plot the regional communities
par(mfrow=c(2,2))
for(i in seq_along(cv_abund)) {
  plot(maps$agg[[i]], axes=F, xlab='', ylab='',
       main=paste('Aggregated (CV = ', cv_abund[[i]]$cv_abund,')'))
  plot(maps$poi[[i]], axes=F, xlab='', ylab='', 
       main=paste('Random (CV = ', cv_abund[[i]]$cv_abund,')'))
}

# set the number of quadrats to sample
n_quadrats_low <- 5
n_quadrats_high <- 20

comms_low = lapply(maps, function(x) 
lapply(x, function(y) 
sample_quadrats(y, n_quadrats_low, plot = F, quadrat_area = 0.01,
method = 'random', avoid_overlap = T)))

comms_high = lapply(maps, function(x) 
lapply(x, function(y) 
sample_quadrats(y, n_quadrats_high, plot = F, quadrat_area = 0.01,
method = 'random', avoid_overlap = T)))

head(comms_low)
# illustrate sampling process
par(mfrow=c(2,2))
sample_quadrats(maps$agg[[1]], n_quadrats_low, plot = T, quadrat_area = 0.01,
                method = 'random', avoid_overlap = T)
title(paste('Aggregated, n_quadrats = ', n_quadrats_low))
sample_quadrats(maps$agg[[1]], n_quadrats_high, plot = T, quadrat_area = 0.01,
                method = 'grid', avoid_overlap = T)
title(paste('Aggregated, n_quadrats = ', n_quadrats_high))
sample_quadrats(maps$poi[[2]], n_quadrats_low, plot = T, quadrat_area = 0.01,
                method = 'transect', avoid_overlap = T)
title(paste('Random, n_quadrats = ', n_quadrats_low))
sample_quadrats(maps$poi[[2]], n_quadrats_high, plot = T, quadrat_area = 0.01,
                method = 'random', avoid_overlap = T)
title(paste('Random, n_quadrats = ', n_quadrats_high))


# aggregate comms data into a community and attributes dataframes ready for mobr
spdat = bind_rows(bind_rows(lapply(comms_low$agg, function(x) x$spec_dat)),
  bind_rows(lapply(comms_high$agg, function(x) x$spec_dat)),
  bind_rows(lapply(comms_low$poi, function(x) x$spec_dat)),
  bind_rows(lapply(comms_high$poi, function(x) x$spec_dat)))

# replace the NAs with zeroes
spdat[is.na(spdat)] <- 0

coords = rbind(bind_rows(lapply(comms_low$agg, function(x) x$xy_dat)),
  bind_rows(lapply(comms_high$agg, function(x) x$xy_dat)),
  bind_rows(lapply(comms_low$poi, function(x) x$xy_dat)),
  bind_rows(lapply(comms_high$poi, function(x) x$xy_dat)))

plot_attr = data.frame(coords, 
                       sample_effort = rep(rep(c('low', 'high'), times = c(n_quadrats_low*2, n_quadrats_high*2)), times = 2),
                       spatial = rep(c('agg', 'poi'), each = (n_quadrats_low+n_quadrats_high)*2),
                       SAD_CV = c(rep(c(1,4), each = n_quadrats_low), rep(c(1,4), each = n_quadrats_high)),
                       Replicate = rep(c(1:n_quadrats_low, 1:n_quadrats_low, 1:n_quadrats_high, 1:n_quadrats_high), times = 2))

# create a group variable (for use with mobr)
plot_attr$group = paste(plot_attr$sample_effort, plot_attr$spatial, plot_attr$SAD_CV, sep='_')


# join the community data frame with the plot attributes
comm_dat <- bind_cols(spdat, plot_attr)
head(comm_dat)

# create some mob-in objects for use in mobr
mob_in_agg1 <- make_mob_in(spdat[c(1:10, 21:40), ], plot_attr[c(1:10, 21:40), ]) 
mob_in_agg4 <- make_mob_in(spdat[c(11:20, 41:60), ], plot_attr[c(11:20, 41:60), ])
mob_in_poi1 <- make_mob_in(spdat[c(61:70, 81:100), ], plot_attr[c(61:70, 81:100), ]) 
mob_in_poi4 <- make_mob_in(spdat[c(71:80, 101:120), ], plot_attr[c(71:80, 101:120), ])


# plot some different types of rarefaction curves
par(mfrow=c(2,2))
plot_rarefaction(mob_in_agg1, env_var = 'group', method = 'indiv', pooled = F)
plot_rarefaction(mob_in_agg4, env_var = 'group', method = 'indiv', pooled = F)
plot_rarefaction(mob_in_poi1, env_var = 'group', method = 'indiv', pooled = F)
plot_rarefaction(mob_in_poi4, env_var = 'group', method = 'indiv', pooled = F)

par(mfrow=c(2,2))
plot_rarefaction(mob_in_agg1, env_var = 'group', method = 'indiv', pooled = T)
plot_rarefaction(mob_in_agg4, env_var = 'group', method = 'indiv', pooled = T)
plot_rarefaction(mob_in_poi1, env_var = 'group', method = 'indiv', pooled = T)
plot_rarefaction(mob_in_poi4, env_var = 'group', method = 'indiv', pooled = T)

par(mfrow=c(2,2))
plot_rarefaction(mob_in_agg1, env_var = 'group', method = 'spat', pooled = T)
plot_rarefaction(mob_in_agg4, env_var = 'group', method = 'spat', pooled = T)
plot_rarefaction(mob_in_poi1, env_var = 'group', method = 'spat', pooled = T)
plot_rarefaction(mob_in_poi4, env_var = 'group', method = 'spat', pooled = T)

# calculate discrete biodiversity metrics
mob_out_agg1 <- get_mob_stats(mob_in_agg1, group_var = 'group')
plot(mob_out_agg1, multi_panel = T)
mob_out_agg4 <- get_mob_stats(mob_in_agg4, group_var = 'group')
plot(mob_out_agg4, multi_panel = T)

mob_out_poi1 <- get_mob_stats(mob_in_poi1, group_var = 'group')
plot(mob_out_poi1, multi_panel = T)

mob_out_poi4 <- get_mob_stats(mob_in_poi4, group_var = 'group')
plot(mob_out_poi4, multi_panel = T)


# wide to long (for tidyverse)
comm_long <- comm_dat %>%
gather(species, abundance, species1:species9) %>%
as_tibble()


# calculate some biodiversity metrics
comm_alpha_summary <- comm_long %>%
  # remove zeroes
  filter(abundance > 0) %>%
  group_by(sample_effort, spatial, SAD_CV, Replicate) %>%
  summarise(
  N = sum(abundance),
  S = sum(abundance>0),
  PIE = calc_PIE(abundance),
  S_PIE = calc_PIE(abundance, ENS = TRUE),
  S_chao = calc_chao1(abundance)) %>%
  ungroup()

comm_alpha_summary$sample_effort <- factor(comm_alpha_summary$sample_effort, levels = c('low', 'high'))


ggplot() +
  facet_wrap(spatial ~ SAD_CV, scales = 'free') +
  geom_boxplot(data = comm_alpha_summary,
  aes(x = sample_effort, y = S))


ggplot() +
  facet_wrap(spatial ~ SAD_CV, scales = 'free') +
  geom_boxplot(data = comm_alpha_summary,
               aes(x = sample_effort, y = S_chao))

ggplot() +
  facet_wrap(spatial ~ SAD_CV, scales = 'free') +
  geom_boxplot(data = comm_alpha_summary,
               aes(x = sample_effort, y = S_PIE))

