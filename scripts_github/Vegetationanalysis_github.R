#Scripts for analysis and plots for main manuscript (H3.2 compare vegetation and surface water)

library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(tidybayes)
library(readr)

####loading, cleaning, and joining data####

ponds_landscape <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/LandscapeFeatures_ponds.csv")%>%  
  #makes a key column for the different pond/date
  mutate(key = paste(date, pond_id, sep = "_")) %>%  
  select("key", "size_ac", "mean_depth_ft", "age", "percent_duckweed")


#load Pb SW for each pond
SW <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/SurfaceWater_lead_ponds.csv")%>%
  #makes a key column for the different pond/date
  mutate(key = paste(date, pond_id, sep = "_"),
         SW_CONCENTRATION_Pb = as.numeric(SW_CONCENTRATION_Pb))%>%
  #and then select for only the columns of interest
  select("key", "pond_id", "SW_CONCENTRATION_Pb")

ponds_allfeatures <- left_join(SW, ponds_landscape, by="key")

#load Pb data for plants
veg <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/Vegetation_lead_ponds.csv")%>%
  #makes a key column for the different pond/date
  mutate(key = paste(date, pond_id, sep = "_"))%>%
  #and then select for only the columns of interest
  select("key", "vegetation_type", "veg_Pb_ppm")

#Add SW to Plants
sw_plants <- left_join(veg, SW, by="key")


####analysis vegetation lead content####
sw_plants %>% 
  ggplot(aes(x = veg_Pb_ppm))+
  theme_bw() +
  stat_bin(bins=30) +
  geom_histogram()

#note: very small sample size

m_veg <- brm(data = sw_plants, family = "gamma",
               veg_Pb_ppm ~ SW_CONCENTRATION_Pb + vegetation_type,
               iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
               control = list(adapt_delta = 0.99999, max_treedepth = 18))

summary(m_veg)
