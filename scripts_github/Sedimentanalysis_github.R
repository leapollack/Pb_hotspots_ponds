#Scripts for analysis and plots for main manuscript (H1 and H3.3)

library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(readr)

####loading, cleaning, and joining data####
#load ponds info
ponds_landscape <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/LandscapeFeatures_ponds.csv")%>%
    mutate(key_pond = paste(date, pond_id, sep = "_"))
  #filter on columns of interest
  select("key_pond", "date", "pond_id", "size_ac", "mean_depth_ft", "age", "land_use", "mean_imp", "mean_builtyear", "total_traffic_500m")


#load Pb data for each sediment core depth
sediment <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/Sediment_lead_ponds.csv")%>%
  mutate(key_pond = paste(date, pond_id, sep = "_"))  %>%
  #filter on columns of interest
  select( "key_pond", "location", "sampling_run","sediment_depth", "Pb_ppm")
 

#add ponds info to sediments
sediment_details <- left_join(sediment, ponds_landscape, by="key_pond")

####model selection -- hypotheses comparison for surface sediment -- excluding lakes ####

#filter out lakes and focusing on top 10 cm of sediment
sediment_ponds <- sediment_details%>% 
  filter(sediment_depth == "0") %>%
  filter(name != "Powderhorn Lake") %>%
  filter(name != "Loeb Lake") %>%
  filter(name != "Spring Lake")

#compare the impact of imp, parcel age, traffic, while controlling for pond size and the age of the basin

#imp
m_imp <- brm(data = sediment_ponds, family = gaussian(),
             Pb_ppm ~ size_ac + age + mean_imp + (1 | pond_id) ,
             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
             control = list(adapt_delta = 0.999, max_treedepth = 17))

#parcel age
m_year <- brm(data = sediment_ponds, family = gaussian(),
              Pb_ppm ~ size_ac + age + mean_builtyear  + (1 | pond_id ) ,
              iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
              control = list(adapt_delta = 0.999, max_treedepth = 17))

#traffic
#zscore traffic volume
sediment_ponds <- sediment_ponds %>%
  mutate(total_traffic_500m_z = scale(total_traffic_500m, center=TRUE, scale = TRUE)) 

m_traffic <- brm(data = sediment_ponds, family = gaussian(),
                 Pb_ppm ~ size_ac + age + total_traffic_500m_z + (1 | pond_id) ,
                 iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                 control = list(adapt_delta = 0.999, max_treedepth = 17))

#basin age alone
m_basin <- brm(data = sediment_ponds, family = gaussian(),
               Pb_ppm ~ size_ac + age + (1 | pond_id ) ,
               iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
               control = list(adapt_delta = 0.999, max_treedepth = 17))


m_imp <- add_criterion(m_imp,"loo", moment_match = TRUE)
m_year <- add_criterion(m_year,"loo", moment_match = TRUE)
m_traffic <- add_criterion(m_traffic,"loo", moment_match = TRUE)
m_basin <- add_criterion(m_basin,"loo", moment_match = TRUE)

loo(m_imp, m_year, m_traffic, m_basin)

#try taking out size

m_year_2 <- brm(data = sediment_ponds, family = gaussian(),
              Pb_ppm ~ age + mean_builtyear  + (1 | pond_id ) ,
              iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
              control = list(adapt_delta = 0.999, max_treedepth = 17))


m_year_2 <- add_criterion(m_year_2,"loo")
m_year_2 <- add_criterion(m_year_2,"loo", moment_match=TRUE)

loo(m_imp, m_year, m_traffic, m_basin, m_year_2)

m_yeartraffic <- brm(data = sediment_ponds, family = gaussian(),
              Pb_ppm ~ age + mean_builtyear + total_traffic_500m_z + (1 | pond_id ) ,
              iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
              control = list(adapt_delta = 0.999, max_treedepth = 17))

m_yeartraffic <- add_criterion(m_yeartraffic,"loo")
m_yeartraffic <- add_criterion(m_yeartraffic,"loo", moment_match = TRUE)

loo(m_imp, m_year, m_traffic, m_basin, m_year_2, m_yeartraffic)

m_yearimp<- brm(data = sediment_ponds, family = gaussian(),
                     Pb_ppm ~ age + mean_builtyear + mean_imp + (1 | pond_id ) ,
                     iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                     control = list(adapt_delta = 0.999, max_treedepth = 17))

m_yearimp<- add_criterion(m_yearimp,"loo")
m_yearimp<- add_criterion(m_yearimp,"loo", moment_match = TRUE)

loo(m_imp, m_year, m_traffic, m_basin, m_year_2, m_yeartraffic, m_yearimp)

summary(m_year_2)

#plot figure 2

sediment_ponds %>% 
  ggplot(aes(x = mean_builtyear, y = Pb_ppm))+
  theme_bw() +
  geom_point()+
  xlab("Mean Parcel Construction Age")+
  ylab("Surface Pb (ppm)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_smooth(method=lm)

p <- conditional_effects(m_year_2) %>% plot(print = FALSE)
p$mean_builtyear$data %>%
  ggplot(aes(x = mean_builtyear)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = sediment_ponds, aes(y = Pb_ppm))+
  xlab("Mean Parcel Age")+
  ylab("Surface Sediment Pb (ppm)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave("~/Desktop/Landscape-level Pond Analysis/plots/Fig2.jpg", width= 6, height = 5)

####comparing soil as a predictor -- subset of the data set####


#add mean soil lead content
#load soil
soil <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/Soil_lead_ponds.csv")%>% mutate(key_pond = paste(date, pond_id, sep = "_"))

#calculate mean soil values for each pond
soil <- soil %>%
  group_by(pond_id) %>%
  mutate(mean_soil_pb = mean(soil_Pb_ppm))%>%
  select("pond_id", "mean_soil_pb")%>%
  distinct(pond_id, mean_soil_pb, .keep_all = TRUE)


#add soil to landscape details
landscape_soil <- left_join(ponds_landscape, soil, by="pond_id")

#add landscape_soil info to sediments
soil_sediment <- left_join(sediment, landscape_soil, by="key_pond")%>%
#filter to only ponds with known sediment lead values
  filter (mean_soil_pb != "NA") %>%
  filter(sediment_depth == "0") %>%
  filter(name != "Powderhorn Lake") %>%
  filter(name != "Loeb Lake") %>%
  filter(name != "Spring Lake")


#rerun best fit model with subset of data that only have soil values
m_year_3 <- brm(data = soil_sediment, family = gaussian(),
                               Pb_ppm ~ age + mean_builtyear+ (1 | pond_id ) ,
                               iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                               control = list(adapt_delta = 0.999, max_treedepth = 17), save_pars = save_pars(all = TRUE))

m_year_3 <- add_criterion(m_year_3, criterion="loo")
m_year_3 <- add_criterion(m_year_3, criterion="loo", moment_match = TRUE)

m_year_soil <- brm(data = soil_sediment, family = gaussian(),
                Pb_ppm ~ age + mean_builtyear+ mean_soil_pb + (1 | pond_id ) ,
                iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                control = list(adapt_delta = 0.999, max_treedepth = 17), save_pars = save_pars(all = TRUE))

m_year_soil <- add_criterion(m_year_soil, criterion="loo")
m_year_soil <- add_criterion(m_year_soil, criterion="loo", moment_match = TRUE)


loo(m_year_3, m_year_soil)

m_soil <- brm(data = soil_sediment, family = gaussian(),
                   Pb_ppm ~ age + mean_soil_pb + (1 | pond_id ) ,
                   iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
                   control = list(adapt_delta = 0.999, max_treedepth = 17), save_pars = save_pars(all = TRUE))

m_soil <- add_criterion(m_soil, criterion="loo")
m_soil <- add_criterion(m_soil, criterion="loo", moment_match = TRUE)


loo(m_year_3, m_year_soil, m_soil)

#does not improve fit, year is still a better predictor, and they are highly correlated 

####analysis impact of depth on Pb####

#negative binomial, includes lakes

m_depth <- brm(data = sediment_details, family = negbinomial(),
               Pb_ppm ~ sediment_depth + (1 | pond_id / location) ,
               iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(),
               control = list(adapt_delta = 0.999, max_treedepth = 17))

summary(m_depth)
conditional_effects(m_depth) #impact of depth on Pb

#plot figure 4
p <- conditional_effects(m_depth) %>% plot(print = FALSE)
p$sediment_depth$data %>%
  ggplot(aes(x = sediment_depth)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = sediment_details, aes(y = Pb_ppm))+
  xlab("Sediment depth from surface (cm)")+
  ylab("Sediment Pb (ppm)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_cartesian(ylim = c(0, 500))

ggsave("~/Desktop/Landscape-level Pond Analysis/plots/Fig4_zoomedin.jpg", width= 6, height = 5)
