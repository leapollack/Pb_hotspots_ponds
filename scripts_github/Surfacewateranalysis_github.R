#Scripts for analysis and plots for main manuscript (H2 surface water)

library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(readr)

####functions for unscaling data when plotting####
#unscale functions
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
  #makes a key column for the different pond/date
  mutate(key = paste(date, pond_id, sep = "_")) %>%  
  select("key", "size_ac", "mean_depth_ft", "age", "percent_duckweed")

#load Pb data for each pond 
sediment <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/Sediment_lead_ponds.csv")%>%  
  #makes a key column for the different pond/date
  mutate(key = paste(date, pond_id, sep = "_"))

#calculate mean sediment surface values for each pond
sediment <- sediment %>%
  filter(sediment_depth == "0") %>%
  group_by(key) %>%
  summarize(sediment_Pb_pond = mean(Pb_ppm))%>%
  ungroup()%>%
  #and then select for only the columns of interest
  select("key", "sediment_Pb_pond")

#load Pb surface water for each pond
SW <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/SurfaceWater_lead_ponds.csv")%>%
  #makes a key column for the different pond/date
  mutate(key = paste(date, pond_id, sep = "_"),
       SW_CONCENTRATION_Pb = as.numeric(SW_CONCENTRATION_Pb))%>%
  #and then select for only the columns of interest
  select("key", "pond_id", "SW_CONCENTRATION_Pb", "PH")


#load internal features for each pond
ponds_internal <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/InternalFeatures_ponds.csv") %>%
  #makes a key column for the different pond/date
  mutate(key = paste(date, pond_id, sep = "_"))%>%
  #and then select for only the columns of interest         
  select("key", "location", "mean_do", "rtrm_topbtm")

#take mean of landscape features per pond_date
internal_perpond <- ponds_internal %>%
  group_by(key) %>%
  summarize(pond_mean_do = mean(mean_do),
            pond_rtrm_topbtm = mean(rtrm_topbtm)) %>%
  ungroup() %>%
  select("key", "pond_mean_do", "pond_rtrm_topbtm")

##grab surface DO, Conductivity, and Temp value from profiles
profile <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/WaterProfiles_ponds.csv")%>%
  #makes a key column for the different pond/date
  mutate(key = paste(date, pond_id, sep = "_"))

#make new dataframe with min_depth_cm only and associated values
profile_update <- profile %>%
  filter(depth_cm == "0")%>%
  group_by(key)%>%
  summarise(surface_do = mean(do_mgL),
            surface_cond = mean (rev_cond_uscm),
            surface_temp = mean (temp_c))%>%
  ungroup()%>%
  #and then select for only the columns of interest
  select("key", "surface_do", "surface_cond", "surface_temp")

#add SW features to landscape features ##already an issue with multiples
ponds_allfeatures <- left_join(SW, ponds_landscape, by="key")

#add profile
ponds_allfeatures <- left_join(ponds_allfeatures, profile_update, by="key") 

#add sediment
ponds_allfeatures <- left_join(ponds_allfeatures, sediment, by="key") 

#add internal features
ponds_SW <- left_join(ponds_allfeatures, internal_perpond, by="key") 

####model selection -- hypotheses comparison -- excluding frozen ponds#### 

#filter out frozen ponds
ponds_SW2 <- ponds_SW %>%
  filter(percent_duckweed != "NA")

#set up data for censored model and scale continuous variables
ponds_SW3 <- ponds_SW2 %>%
  mutate(SW_CONCENTRATION_Pb = ifelse(is.na(SW_CONCENTRATION_Pb), 
                                      0.01, 
                                      SW_CONCENTRATION_Pb),
         censored = ifelse(SW_CONCENTRATION_Pb == 0.01, -1, 0),
         size_ac_scaled = scale(size_ac)[,1],
         mean_depth_ft_scaled = scale(mean_depth_ft)[,1],
         pond_mean_do_scaled = scale(pond_mean_do)[,1],
         pond_rtrm_topbtm_scaled = scale(pond_rtrm_topbtm)[,1],
         surface_do_scaled = scale(surface_do)[,1],
         surface_cond_scaled = scale(surface_cond)[,1],
         surface_temp_scaled = scale(surface_temp)[,1],
         sediment_Pb_pond_scaled = scale(sediment_Pb_pond)[,1],
         percent_duckweed_scaled = scale(percent_duckweed)[,1])


#filter for initial model comparison
ponds_SW4 <- ponds_SW3 %>%
  filter(pond_rtrm_topbtm != "NA",
         surface_do != "NA")

#compare the impact of DO, RTRM (Density Impacted Mixing), Duckweed Cover, Conductivity (Salt), and Temp while controlling for pond depth and the local sediment lead content

#mean do
model_cens_gamma_do <- brm(data = ponds_SW4, family = "gamma",
                           bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + pond_mean_do_scaled),
                           iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores())


#surface do
model_cens_gamma_surfacedo <- brm(data = ponds_SW4, family = "gamma",
                                  bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + surface_do_scaled),
                                  iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores())

#rtrm
model_cens_gamma_rtrm <- brm(data = ponds_SW4, family = "gamma",
                                  bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + pond_rtrm_topbtm_scaled),
                                  iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

#cond
model_cens_gamma_cond <- brm(data = ponds_SW4, family = "gamma",
                                  bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + surface_cond_scaled),
                                  iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))


#cover
model_cens_gamma_cover <- brm(data = ponds_SW4, family = "gamma",
                                  bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + percent_duckweed_scaled),
                                  iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores())

#temp
model_cens_gamma_temp <- brm(data = ponds_SW4, family = "gamma",
                              bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + surface_temp_scaled),
                              iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores())


model_cens_gamma_do <- add_criterion(model_cens_gamma_do,"loo", moment_match = TRUE)
model_cens_gamma_surfacedo <- add_criterion(model_cens_gamma_surfacedo,"loo", moment_match = TRUE)
model_cens_gamma_rtrm <- add_criterion(model_cens_gamma_rtrm, "loo", moment_match = TRUE)
model_cens_gamma_cond <- add_criterion(model_cens_gamma_cond, "loo", moment_match = TRUE)
model_cens_gamma_cover <- add_criterion(model_cens_gamma_cover, "loo", moment_match = TRUE)
model_cens_gamma_temp <- add_criterion(model_cens_gamma_temp, "loo", moment_match = TRUE)


loo(model_cens_gamma_do, model_cens_gamma_surfacedo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover, model_cens_gamma_temp)

##rtrm is a better predictor than conductivity, followed by DO

#size
model_cens_gamma_size <- brm(data = ponds_SW4, family = "gamma",
                             bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + size_ac_scaled + surface_temp_scaled),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores())


model_cens_gamma_size <- add_criterion(model_cens_gamma_size, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_surfacedo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover, model_cens_gamma_temp, model_cens_gamma_size)

#rtrm + temp
model_cens_gamma_rtrmtemp <- brm(data = ponds_SW4, family = "gamma",
                             bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + pond_rtrm_topbtm_scaled + surface_temp_scaled),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))


model_cens_gamma_rtrmtemp <- add_criterion(model_cens_gamma_rtrmtemp, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_surfacedo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover, model_cens_gamma_temp, model_cens_gamma_size, model_cens_gamma_rtrmtemp)

#rtrm + do
model_cens_gamma_rtrmdo <- brm(data = ponds_SW4, family = "gamma",
                                 bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + pond_rtrm_topbtm_scaled + pond_mean_do_scaled),
                                 iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

model_cens_gamma_rtrmdo <- add_criterion(model_cens_gamma_rtrmdo, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_surfacedo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover, model_cens_gamma_temp, model_cens_gamma_size, model_cens_gamma_rtrmtemp, model_cens_gamma_rtrmdo)

#try removing depth
model_cens_gamma_rtrm_2 <- brm(data = ponds_SW4, family = "gamma",
                               bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + pond_rtrm_topbtm_scaled ),
                               iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

model_cens_gamma_rtrm_2 <- add_criterion(model_cens_gamma_rtrm_2, "loo", moment_match = TRUE)

loo(model_cens_gamma_do, model_cens_gamma_surfacedo, model_cens_gamma_rtrm, model_cens_gamma_cond, model_cens_gamma_cover, model_cens_gamma_temp, model_cens_gamma_size, model_cens_gamma_rtrmtemp, model_cens_gamma_rtrmdo, model_cens_gamma_rtrm_2 ) #best without depth

#run best fit model with full data set

model_cens_gamma_full <- brm(data = ponds_SW3, family = "gamma",
                             bf(SW_CONCENTRATION_Pb | cens(censored) ~ sediment_Pb_pond_scaled + pond_rtrm_topbtm_scaled ),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))


summary(model_cens_gamma_full)

##plots and figures 3a

ponds_SW3 %>% 
  ggplot(aes(x = pond_rtrm_topbtm, y = SW_CONCENTRATION_Pb))+
  theme_bw() +
  geom_point()+
  xlab("RTRM")+
  ylab("Surface Water Pb")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_smooth(method=lm)

p <- conditional_effects(model_cens_gamma_full) %>% plot(print = FALSE)

p$pond_rtrm_topbtm_scaled$data %>%
  unscale_dataframe(original = ponds_SW3) %>%
  ggplot(aes(x = pond_rtrm_topbtm)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = ponds_SW4, aes(y = SW_CONCENTRATION_Pb))+
  xlab("Relative Thermal Resistance to Mixing (RTRM)")+
  ylab("Surface Water Pb (ppb)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave("~/Desktop/Landscape-level Pond Analysis/plots/Fig3a.jpg", width= 4, height = 5)

####model selection -- hypotheses comparison -- frozen ponds only for appendix S1 ####

#filter to only frozen ponds
ponds_SW1 <- ponds_SW %>%
  filter(is.na(percent_duckweed))%>%
  filter(SW_CONCENTRATION_Pb != "NA")

ponds_SW2 <- ponds_SW1 %>%
  mutate(SW_CONCENTRATION_Pb = ifelse(is.na(SW_CONCENTRATION_Pb), 
                                      0.01, 
                                      SW_CONCENTRATION_Pb),
         censored = ifelse(SW_CONCENTRATION_Pb == 0.01, -1, 0),
         size_ac_scaled = scale(size_ac)[,1],
         mean_depth_ft_scaled = scale(mean_depth_ft)[,1],
         pond_mean_do_scaled = scale(pond_mean_do)[,1],
         pond_rtrm_topbtm_scaled = scale(pond_rtrm_topbtm)[,1],
         surface_do_scaled = scale(surface_do)[,1],
         surface_cond_scaled = scale(surface_cond)[,1],
         surface_temp_scaled = scale (surface_temp)[,1],
         sediment_Pb_pond_scaled = scale (sediment_Pb_pond)[,1],
         percent_duckweed_scaled = scale(percent_duckweed)[,1],
         PH_scaled = scale(PH)[,1])


#filter for initial model comparison
ponds_SW5 <- ponds_SW2 %>%
  filter(sediment_Pb_pond != "NA")


#compare the impact of DO, RTRM (Density Impacted Mixing), Conductivity (Salt), and Temp while controlling for pond depth and the local sediment lead content

#mean do
model_gamma_do <- brm(data = ponds_SW5, family = "gamma",
                           bf(SW_CONCENTRATION_Pb ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + pond_mean_do_scaled),
                           iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.99999))


#rtrm
model_gamma_rtrm <- brm(data = ponds_SW5, family = "gamma",
                             bf(SW_CONCENTRATION_Pb  ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + pond_rtrm_topbtm_scaled),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))

#cond
model_gamma_cond <- brm(data = ponds_SW5, family = "gamma",
                             bf(SW_CONCENTRATION_Pb ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + surface_cond_scaled),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))


#temp
model_gamma_temp <- brm(data = ponds_SW5, family = "gamma",
                             bf(SW_CONCENTRATION_Pb ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + surface_temp_scaled),
                             iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.9999))


model_gamma_do <- add_criterion(model_gamma_do, "loo")
model_gamma_do <- add_criterion(model_gamma_do, "loo", moment_match = TRUE)
model_gamma_rtrm <- add_criterion(model_gamma_rtrm, "loo")
model_gamma_rtrm <- add_criterion(model_gamma_rtrm, "loo", moment_match = TRUE)
model_gamma_cond <- add_criterion(model_gamma_cond, "loo")
model_gamma_cond <- add_criterion(model_gamma_cond, "loo", moment_match = TRUE)
model_gamma_temp <- add_criterion(model_gamma_temp, "loo")
model_gamma_temp <- add_criterion(model_gamma_temp, "loo", moment_match = TRUE)

loo(model_gamma_do, model_gamma_rtrm, model_gamma_cond, model_gamma_temp)

#temp, do, rtrm are better fits than the other variables

#temp and do
model_gamma_tempdo <- brm(data = ponds_SW5, family = "gamma",
                      bf(SW_CONCENTRATION_Pb ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + pond_mean_do_scaled + surface_temp_scaled),
                      iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.99999))

#temp and rtrm
model_gamma_temprtrm <- brm(data = ponds_SW5, family = "gamma",
                      bf(SW_CONCENTRATION_Pb ~ sediment_Pb_pond_scaled + mean_depth_ft_scaled + surface_temp_scaled + pond_rtrm_topbtm_scaled),
                      iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(max_treedepth = 19, adapt_delta = 0.99999))

model_gamma_tempdo <- add_criterion(model_gamma_tempdo, "loo")
model_gamma_tempdo <- add_criterion(model_gamma_tempdo, "loo", moment_match = TRUE)
model_gamma_temprtrm <- add_criterion(model_gamma_temprtrm, "loo")
model_gamma_temprtrm <- add_criterion(model_gamma_temprtrm, "loo", moment_match = TRUE)

loo(model_gamma_do, model_gamma_rtrm, model_gamma_cond, model_gamma_temp, model_gamma_tempdo, model_gamma_temprtrm)

summary(model_gamma_temp )
