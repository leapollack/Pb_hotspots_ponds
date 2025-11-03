#Scripts for analysis and plots for main manuscript (H2 bottom water)

library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(readr)

####functions for unscaling data when plotting####
unscale <- function(scaled_vector, unscaled_vector = NULL) {
  # Extract scaling parameters from the scaled vector attributes if available
  center <- attr(scaled_vector, "scaled:center")
  scale_param <- attr(scaled_vector, "scaled:scale")
  
  # If attributes are not available and unscaled_vector is provided,
  # calculate parameters from the unscaled vector
  if (is.null(center) && !is.null(unscaled_vector)) {
    center <- mean(unscaled_vector, na.rm = TRUE)
  }
  if (is.null(scale_param) && !is.null(unscaled_vector)) {
    scale_param <- sd(unscaled_vector, na.rm = TRUE)
  }
  
  # Check if we have the necessary parameters
  if (is.null(center) || is.null(scale_param)) {
    stop("Unable to determine scaling parameters. Either the scaled vector must have attributes or unscaled_vector must be provided.")
  }
  
  # Reverse the scaling: multiply by scale and add center
  unscaled <- as.vector(scaled_vector) * scale_param + center
  
  return(unscaled)
}

unscale_dataframe <- function(scaled, original) {
  # Make a copy of the scaled dataframe to avoid modifying the original
  result <- scaled
  
  # Find columns ending with '_scaled'
  scaled_cols <- names(scaled)[grepl("_scaled$", names(scaled))]
  
  if (length(scaled_cols) == 0) {
    message("No columns ending with '_scaled' found in the scaled dataframe.")
    return(result)
  }
  
  # Process each scaled column
  for (col in scaled_cols) {
    # Extract the base column name (remove '_scaled' suffix)
    base_name <- sub("_scaled$", "", col)
    
    # Check if the base column exists in the original dataframe
    if (base_name %in% names(original)) {
      # Unscale the column using the original column as reference
      unscaled_values <- unscale(scaled[[col]], original[[base_name]])
      
      # Add the unscaled column to the result dataframe
      # Use the base name (without '_scaled') for the new column
      result[[base_name]] <- unscaled_values
      
      message(sprintf("Unscaled column '%s' added as '%s'", col, base_name))
    } else {
      warning(sprintf("Column '%s' not found in original dataframe. Skipping '%s'.", 
                      base_name, col))
    }
  }
  
  return(result)
}

####loading, cleaning, and joining data####
#load ponds info
ponds_landscape <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/LandscapeFeatures_ponds.csv")%>%
  mutate(key_pond = paste(date, pond_id, sep = "_")) %>%  
  select("key_pond", "size_ac", "mean_depth_ft", "age", "percent_duckweed")

#load Pb data for each sediment core depth
sediment <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/Sediment_lead_ponds.csv")  %>%
  #makes a key column for the different locations
  mutate(key = paste(date, pond_id, location, sep = "_"))

#calculate mean sediment surface values for each location
sediment <- sediment %>%
  filter(sediment_depth == "0") %>%
  group_by(key) %>%
  summarise(sediment_Pb = mean(Pb_ppm))%>%
  ungroup()%>%
#and then select for only the columns of interest
  select("key", "sediment_Pb")

#load Pb bottom water for each core location
BW <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/BottomWater_lead_ponds.csv") %>%
  #makes a key column for the different locations
  mutate(key = paste(date, pond_id, location, sep = "_"),
         BW_CONCENTRATION_Pb = as.numeric(BW_CONCENTRATION_Pb))%>%
  #and then select for only the columns of interest
  select("key", "BW_CONCENTRATION_Pb", "PH")

#load internal features for each pond
ponds_internal <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/InternalFeatures_ponds.csv") %>%
  #makes a key column for the different locations
  mutate(key = paste(date, pond_id, location, sep = "_"),
         key_pond = paste(date, pond_id, sep = "_"))%>%
  #and then select for only the columns of interest
  select("key", "key_pond", "pond_id", "location", "mean_do", "rtrm_topbtm")

#add ponds info to internal features
ponds_allfeatues <- left_join(ponds_landscape, ponds_internal, by="key_pond")

#grab bottom DO, conductivity, and temp values from profiles
profile <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/WaterProfiles_ponds.csv") %>%
  #makes a key column for the different locations
  mutate(key = paste(date, pond_id, location, sep = "_"))

#make new dataframe with max_depth_cm only and associated values
profile <- profile %>%
  group_by(key)%>%
  mutate(bottom_do = ifelse(max_depth == depth_cm, do_mgL, ifelse(NA)),
         bottom_cond = ifelse(max_depth == depth_cm, rev_cond_uscm, ifelse(NA)),
         bottom_temp = ifelse(max_depth == depth_cm, temp_c, ifelse(NA)))%>%
  filter(!is.na(bottom_do),
         !is.na(bottom_cond),
         !is.na(bottom_temp))%>%
  ungroup()%>%
  select("key","bottom_do", "bottom_cond", "bottom_temp")

#add bottom profile values to the pond features
ponds_allfeatues <- left_join(ponds_allfeatues, profile, by="key")


#add sediment to internal features
ponds_features_sediment <- left_join(ponds_allfeatues, sediment, by="key")

#add BW to everything
ponds_BW <- left_join(BW, ponds_features_sediment, by="key")%>%
  filter(key != "NA_NA_NA")

#add AF
ponds_AF <- read_csv("~/Desktop/Landscape-level Pond Analysis/cleaned data/2024_anoxic_factor_May302025.csv")%>%
  select("pond_id", "af_mean_seas")

ponds_BW_test <- left_join(ponds_BW, ponds_AF, by="pond_id")

####model selection -- hypotheses comparison -- unfrozen ponds only ####

#filter out frozen ponds and scale
ponds_BW3 <- ponds_BW %>%
  filter(percent_duckweed != "NA") %>%
  mutate(BW_CONCENTRATION_Pb = ifelse(is.na(BW_CONCENTRATION_Pb), 
                                      0.01, 
                                      BW_CONCENTRATION_Pb),
         censored = ifelse(BW_CONCENTRATION_Pb == 0.01, -1, 0),
         size_ac_scaled = scale(size_ac)[,1],
         mean_depth_ft_scaled = scale(mean_depth_ft)[,1],
         mean_do_scaled = scale(mean_do)[,1],
         rtrm_topbtm_scaled = scale(rtrm_topbtm)[,1],
         bottom_do_scaled = scale(bottom_do)[,1],
         bottom_cond_scaled = scale(bottom_cond)[,1],
         bottom_temp_scaled = scale(bottom_temp)[,1],
         sediment_Pb_scaled = scale(sediment_Pb)[,1],
         percent_duckweed_scaled = scale(percent_duckweed)[,1])

#filter for initial model comparison
ponds_BW4 <- ponds_BW3 %>%
  filter(rtrm_topbtm != "NA")

#compare the impact of DO, RTRM (Density Impacted Mixing), Duckweed Cover, Conductivity (Salt), and Temp while controlling for pond depth and the local sediment lead content

#mean do
model_cens_gamma_do <- brm(data = ponds_BW4, family = "gamma",
                           bf(BW_CONCENTRATION_Pb | cens(censored) ~ mean_depth_ft_scaled + mean_do_scaled + sediment_Pb_scaled + (1 | pond_id)),
                           iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))


#bottom do
model_cens_gamma_bottomdo <- brm(data = ponds_BW4, family = "gamma",
                                  bf(BW_CONCENTRATION_Pb | cens(censored) ~ mean_depth_ft_scaled + bottom_do_scaled + sediment_Pb_scaled + (1 | pond_id)),
                                  iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

#rtrm
model_cens_gamma_rtrm <- brm(data = ponds_BW4, family = "gamma",
                             bf(BW_CONCENTRATION_Pb | cens(censored) ~  mean_depth_ft_scaled + rtrm_topbtm_scaled + sediment_Pb_scaled + (1 | pond_id)),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

#cond
model_cens_gamma_cond <- brm(data = ponds_BW4, family = "gamma",
                             bf(BW_CONCENTRATION_Pb | cens(censored) ~  mean_depth_ft_scaled + bottom_cond_scaled + sediment_Pb_scaled + (1 | pond_id)),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))


#cover
model_cens_gamma_cover <- brm(data = ponds_BW4, family = "gamma",
                              bf(BW_CONCENTRATION_Pb | cens(censored) ~  mean_depth_ft_scaled + percent_duckweed_scaled + sediment_Pb_scaled + (1 | pond_id)),
                              iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

#temp
model_cens_gamma_temp <- brm(data = ponds_BW4, family = "gamma",
                                 bf(BW_CONCENTRATION_Pb | cens(censored) ~  mean_depth_ft_scaled + bottom_temp_scaled + sediment_Pb_scaled + (1 | pond_id)),
                                 iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores())

model_cens_gamma_do <- add_criterion(model_cens_gamma_do, "loo")
model_cens_gamma_do <- add_criterion(model_cens_gamma_do, "loo", moment_match = TRUE)
model_cens_gamma_bottomdo <- add_criterion(model_cens_gamma_bottomdo, "loo")
model_cens_gamma_bottomdo <- add_criterion(model_cens_gamma_bottomdo, "loo", moment_match = TRUE)
model_cens_gamma_rtrm <- add_criterion(model_cens_gamma_rtrm, "loo")
model_cens_gamma_rtrm <- add_criterion(model_cens_gamma_rtrm, "loo", moment_match = TRUE)
model_cens_gamma_cond <- add_criterion(model_cens_gamma_cond, "loo")
model_cens_gamma_cond <- add_criterion(model_cens_gamma_cond, "loo", moment_match = TRUE)
model_cens_gamma_cover <- add_criterion(model_cens_gamma_cover, "loo")
model_cens_gamma_cover <- add_criterion(model_cens_gamma_cover, "loo", moment_match = TRUE)
model_cens_gamma_temp <- add_criterion(model_cens_gamma_temp, "loo")
model_cens_gamma_temp <- add_criterion(model_cens_gamma_temp, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_bottomdo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover,model_cens_gamma_temp)

#cond, temp, rtrm fits best 

#include both cond and temp

model_cens_gamma_condtemp <- brm(data = ponds_BW4, family = "gamma",
                                 bf(BW_CONCENTRATION_Pb | cens(censored) ~ mean_depth_ft_scaled + bottom_temp_scaled + bottom_cond_scaled + sediment_Pb_scaled + (1 | pond_id)),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))

model_cens_gamma_condtemp <- add_criterion(model_cens_gamma_condtemp, "loo")
model_cens_gamma_condtemp <- add_criterion(model_cens_gamma_condtemp, "loo", moment_match = TRUE)


loo(model_cens_gamma_do, model_cens_gamma_bottomdo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover,model_cens_gamma_temp, model_cens_gamma_condtemp)

#including both doesn't really improve fit

#what about adding rtrm instead?

model_cens_gamma_condrtrm <- brm(data = ponds_BW4, family = "gamma",
                                 bf(BW_CONCENTRATION_Pb | cens(censored) ~ mean_depth_ft_scaled + bottom_cond_scaled + sediment_Pb_scaled + rtrm_topbtm_scaled +(1 | pond_id)),
                                 iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))

model_cens_gamma_condrtrm <- add_criterion(model_cens_gamma_condrtrm, "loo")
model_cens_gamma_condrtrm <- add_criterion(model_cens_gamma_condrtrm, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_bottomdo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover,model_cens_gamma_temp, model_cens_gamma_condtemp, model_cens_gamma_condrtrm)

#what about adding cover + conductivity

model_cens_gamma_condcover <- brm(data = ponds_BW4, family = "gamma",
                                 bf(BW_CONCENTRATION_Pb | cens(censored) ~ mean_depth_ft_scaled + percent_duckweed_scaled + bottom_cond_scaled + sediment_Pb_scaled +(1 | pond_id)),
                                 iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))

model_cens_gamma_condcover <- add_criterion(model_cens_gamma_condcover, "loo")
model_cens_gamma_condcover <- add_criterion(model_cens_gamma_condcover, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_bottomdo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover,model_cens_gamma_temp, model_cens_gamma_condtemp, model_cens_gamma_condrtrm, model_cens_gamma_condcover)

#does size matter here as well?

model_cens_gamma_size <- brm(data = ponds_BW4, family = "gamma",
                                  bf(BW_CONCENTRATION_Pb | cens(censored) ~ size_ac_scaled + mean_depth_ft_scaled + sediment_Pb_scaled +(1 | pond_id)),
                                  iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))

model_cens_gamma_size <- add_criterion(model_cens_gamma_size, "loo")
model_cens_gamma_size <- add_criterion(model_cens_gamma_size, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_bottomdo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover,model_cens_gamma_temp, model_cens_gamma_condtemp, model_cens_gamma_condrtrm, model_cens_gamma_condcover, model_cens_gamma_size)


loo(model_cens_gamma_do, model_cens_gamma_bottomdo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover,model_cens_gamma_temp, model_cens_gamma_condtemp, model_cens_gamma_condrtrm, model_cens_gamma_condcover, model_cens_gamma_condcover)

#try cutting out depth on best fit model

model_cens_gamma_condcover_2 <- brm(data = ponds_BW4, family = "gamma",
                                  bf(BW_CONCENTRATION_Pb | cens(censored) ~  percent_duckweed_scaled + bottom_cond_scaled + sediment_Pb_scaled + (1 | pond_id)),
                                  iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))

model_cens_gamma_condcover_2 <- add_criterion(model_cens_gamma_condcover_2, "loo")
model_cens_gamma_condcover_2 <- add_criterion(model_cens_gamma_condcover_2, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_bottomdo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover,model_cens_gamma_temp, model_cens_gamma_condtemp, model_cens_gamma_condrtrm, model_cens_gamma_condcover, model_cens_gamma_condcover_2, model_cens_gamma_size)

#run best fit model with full data set

model_cens_gamma_full <- brm(data = ponds_BW3, family = "gamma",
                             bf(BW_CONCENTRATION_Pb | cens(censored) ~  percent_duckweed_scaled + bottom_cond_scaled + sediment_Pb_scaled + (1 | pond_id)),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))

summary(model_cens_gamma_full)

##plots and figures 3b

#conditional effects of duckweed

ponds_BW3 %>% 
  ggplot(aes(x = percent_duckweed, y = BW_CONCENTRATION_Pb))+
  theme_bw() +
  geom_point()+
  xlab("Percent Duckweed Cover")+
  ylab("Bottom Water Pb")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_smooth(method=lm)

p <- conditional_effects(model_cens_gamma_full) %>% plot(print = FALSE)

model_cens_gamma_full$data$bottom_cond_scaled %>% range()

p$percent_duckweed_scaled$data %>%
  unscale_dataframe(original = ponds_BW3) %>%
  ggplot(aes(x = percent_duckweed)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = ponds_BW3, aes(y = BW_CONCENTRATION_Pb))+
  xlab("Percent Duckweed Cover")+
  ylab("Bottom Water Pb (ppb)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#conditional effects of conductivity

ponds_BW3 %>% 
  ggplot(aes(x = bottom_cond, y = BW_CONCENTRATION_Pb))+
  theme_bw() +
  geom_point()+
  xlab("Conductivity of Bottom Water")+
  ylab("Bottom Water Pb")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_smooth(method=lm)

p$bottom_cond_scaled$data %>%
  unscale_dataframe(original = ponds_BW3) %>%
  distinct(bottom_cond, bottom_cond_scaled)

p$bottom_cond_scaled$data %>%
  unscale_dataframe(original = ponds_BW3) %>%
  ggplot(aes(x = bottom_cond)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = ponds_BW3, aes(y = BW_CONCENTRATION_Pb))+
  xlab("Conductivity of Bottom Water (µS/cm)")+
  ylab("Bottom Water Pb (ppb)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave("~/Desktop/Landscape-level Pond Analysis/plots/Fig3b.jpg", width= 4, height = 5)


####model selection -- hypotheses comparison -- frozen ponds only for appendix S1 ####

#filter to only frozen ponds
ponds_BW1 <- ponds_BW %>%
  filter(is.na(percent_duckweed))%>%
  filter(BW_CONCENTRATION_Pb != "NA")


ponds_BW2 <- ponds_BW1 %>%
  mutate(BW_CONCENTRATION_Pb = ifelse(is.na(BW_CONCENTRATION_Pb), 
                                      0.01, 
                                      BW_CONCENTRATION_Pb),
         censored = ifelse(BW_CONCENTRATION_Pb == 0.01, -1, 0),
         size_ac_scaled = scale(size_ac)[,1],
         mean_depth_ft_scaled = scale(mean_depth_ft)[,1],
         mean_do_scaled = scale(mean_do)[,1],
         rtrm_topbtm_scaled = scale(rtrm_topbtm)[,1],
         bottom_do_scaled = scale (bottom_do)[,1],
         bottom_cond_scaled = scale (bottom_cond)[,1],
         bottom_temp_scaled = scale (bottom_temp)[,1],
         sediment_Pb_scaled = scale (sediment_Pb)[,1],
         percent_duckweed_scaled = scale(percent_duckweed)[,1],
         PH_scaled = scale(PH)[,1])

#filter for initial model comparison
ponds_BW5 <- ponds_BW2 %>%
  filter(sediment_Pb != "NA")

#compare the impact of DO, RTRM (Density Impacted Mixing), Conductivity (Salt), and Temp while controlling for pond depth and the local sediment lead content

#mean do
model_gamma_do <- brm(data = ponds_BW5, family = "gamma",
                           bf(BW_CONCENTRATION_Pb ~ mean_depth_ft_scaled + mean_do_scaled + sediment_Pb_scaled + (1 | pond_id)),
                           iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))


#bottom do
model_gamma_bottomdo <- brm(data = ponds_BW5, family = "gamma",
                                 bf(BW_CONCENTRATION_Pb  ~ mean_depth_ft_scaled + bottom_do_scaled + sediment_Pb_scaled + (1 | pond_id)),
                                 iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

#rtrm
model_gamma_rtrm <- brm(data = ponds_BW5, family = "gamma",
                             bf(BW_CONCENTRATION_Pb ~  mean_depth_ft_scaled + rtrm_topbtm_scaled + sediment_Pb_scaled + (1 | pond_id)),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 20, adapt_delta = 0.9999))

#cond
model_gamma_cond <- brm(data = ponds_BW5, family = "gamma",
                             bf(BW_CONCENTRATION_Pb  ~  mean_depth_ft_scaled + bottom_cond_scaled + sediment_Pb_scaled + (1 | pond_id)),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 20))

#temp
model_gamma_temp <- brm(data = ponds_BW5, family = "gamma",
                             bf(BW_CONCENTRATION_Pb  ~  mean_depth_ft_scaled + bottom_temp_scaled + sediment_Pb_scaled + (1 | pond_id)),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 20))


model_gamma_do <- add_criterion(model_gamma_do, "loo")
model_gamma_do <- add_criterion(model_gamma_do, "loo", moment_match = TRUE)
model_gamma_bottomdo <- add_criterion(model_gamma_bottomdo, "loo")
model_gamma_bottomdo <- add_criterion(model_gamma_bottomdo, "loo", moment_match = TRUE)
model_cens_gamma_rtrm <- add_criterion(model_gamma_rtrm, "loo")
model_cens_gamma_rtrm <- add_criterion(model_gamma_rtrm, "loo", moment_match = TRUE)
model_cens_gamma_cond <- add_criterion(model_gamma_cond, "loo")
model_cens_gamma_cond <- add_criterion(model_gamma_cond, "loo", moment_match = TRUE)
model_cens_gamma_temp <- add_criterion(model_gamma_temp, "loo")
model_cens_gamma_temp <- add_criterion(model_gamma_temp, "loo", moment_match = TRUE)

loo(model_gamma_do, model_gamma_bottomdo, model_gamma_rtrm, model_gamma_cond,model_gamma_temp)


#Double check best fit going down from here
#temp, rtrm, cond best fit

model_gamma_temprtrm <- brm(data = ponds_BW5, family = "gamma",
                           bf(BW_CONCENTRATION_Pb ~ mean_depth_ft_scaled + bottom_temp_scaled + rtrm_topbtm_scaled + sediment_Pb_scaled + (1 | pond_id)),
                           iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

model_gamma_temprtrm <- add_criterion(model_gamma_temprtrm, "loo")
model_gamma_temprtrm <- add_criterion(model_gamma_temprtrm, "loo", moment_match = TRUE)

#cond and temp

model_gamma_condtemp <- brm(data = ponds_BW5, family = "gamma",
                                bf(BW_CONCENTRATION_Pb ~ mean_depth_ft_scaled + bottom_temp_scaled + bottom_cond_scaled + sediment_Pb_scaled + (1 | pond_id)),
                                iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

model_gamma_condtemp <- add_criterion(model_gamma_condtemp , "loo")
model_gamma_condtemp  <- add_criterion(model_gamma_condtemp , "loo", moment_match = TRUE)


loo(model_gamma_do, model_gamma_bottomdo, model_gamma_rtrm, model_gamma_cond,model_gamma_temp, model_gamma_condtemp, model_gamma_temprtrm)

#rtrm and cond

model_gamma_rtrmcond <- brm(data = ponds_BW5, family = "gamma",
                               bf(BW_CONCENTRATION_Pb ~ mean_depth_ft_scaled + rtrm_topbtm_scaled + bottom_cond_scaled + sediment_Pb_scaled + (1 | pond_id)),
                               iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

model_gamma_rtrmcond <- add_criterion(model_gamma_rtrmcond, "loo")
model_gamma_rtrmcond <- add_criterion(model_gamma_rtrmcond, "loo", moment_match = TRUE)


loo(model_gamma_do, model_gamma_bottomdo, model_gamma_rtrm, model_gamma_cond,model_gamma_temp, model_gamma_condtemp, model_gamma_temprtrm, model_gamma_rtrmcond)

#Model comparisons:
#elpd_diff se_diff
#model_gamma_temprtrm  0.0       0.0   
#model_gamma_temp     -3.4       1.6   
#model_gamma_condtemp -3.6       1.4   
#model_gamma_rtrm     -4.4       1.1   
#model_gamma_cond     -4.8       1.0   
#model_gamma_do       -5.2       1.5   
#model_gamma_rtrmcond -5.2       1.2   
#model_gamma_bottomdo -5.5       1.3

summary(model_gamma_temprtrm)