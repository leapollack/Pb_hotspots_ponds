#Scripts for analysis and plots for main manuscript (H3.1 compare across water column depths)

library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(tidybayes)
library(readr)

####loading, cleaning, and joining data####

SW <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/SurfaceWater_lead_ponds.csv")

SW <- rename(SW, water_pb = SW_CONCENTRATION_Pb)

BW <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/BottomWater_lead_ponds.csv") %>%
  #makes a key column for the different locations
  mutate(key = paste(date, pond_id, location, sep = "_"),
         BW_CONCENTRATION_Pb = as.numeric(BW_CONCENTRATION_Pb))

BW <- rename(BW, water_pb = BW_CONCENTRATION_Pb)

PW <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/PoreWater_lead_ponds.csv") %>%
  #makes a key column for the different locations
  mutate(key = paste(date, pond_id, location, sep = "_"),
         PW_CONCENTRATION_Pb = as.numeric(PW_CONCENTRATION_Pb))

PW <- PW %>%
  group_by(date, pond_id, location, sediment_depth) %>%
  mutate(mean_PW_pb = mean(PW_CONCENTRATION_Pb))%>%
  ungroup() %>%
  select("key", "sediment_depth", "mean_PW_pb", "water_depth", "pond_id", "location")%>%
  distinct(pond_id, mean_PW_pb, .keep_all = TRUE)

PW_surface <- PW %>%
  filter(sediment_depth == "0") %>%
  select("key", "mean_PW_pb", "water_depth", "pond_id", "location")

PW_surface <- rename(PW_surface, water_pb = mean_PW_pb)

#load Pb data for each sediment core depth
sediment <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/Sediment_lead_ponds.csv")

#calculate mean sediment surface values for each pond
sediment <- sediment %>%
  filter(sediment_depth == "0") %>%
  group_by(pond_id) %>%
  summarise(sediment_Pb = mean(Pb_ppm))%>%
  ungroup()%>%
  #and then select for only the columns of interest
  select("pond_id", "sediment_Pb")

#grab pond depth as well to control for height of pond
ponds_landscape <- read_csv("~/Desktop/Landscape-level Pond Analysis/data_github/LandscapeFeatures_ponds.csv")%>%
  select("pond_id", "mean_depth_ft")

#join all together

all_pb <- PW_surface %>%
  bind_rows(BW %>%
              select(pond_id, location, water_depth, water_pb)) %>%
  bind_rows(SW %>%
              select(pond_id, water_depth, water_pb))

all_pb <- left_join(all_pb, sediment)

all_pb <- left_join(all_pb, ponds_landscape)

all_pb %>%
  print(n=Inf)


####analysis across water column depth types####

#add censoring
all_pb2 <- all_pb %>%
  mutate(water_pb = ifelse(is.na(water_pb), 
                                  0.01, 
                                  water_pb),
         censored = ifelse(water_pb == 0.01, -1, 0))

model_depth_cens_gamma <- brm(data = all_pb2, family = "gamma",
                           bf(water_pb | cens(censored) ~ water_depth + sediment_Pb + mean_depth_ft+ (1|pond_id)),
                           iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores())

#median odds ratios
#estimating the differences between water column depths
model_depth_cens_gamma_emm <- model_depth_cens_gamma %>% 
  emmeans::emmeans(~ water_depth) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
#values
model_depth_cens_gamma_emm  %>% 
  median_hdci()

####figure 1####

#final form for fig 1
all_pb %>% 
  filter(!is.na(water_depth),
         !is.na(water_pb)) %>% 
  mutate(water_depth_label = case_when(
    water_depth == "SW" ~ "Surface water",
    water_depth == "BW" ~ "Bottom water",
    water_depth == "PW" ~ "Pore water",
    TRUE ~ as.character(water_depth)
  )) %>% 
  ggplot(aes(x = factor(pond_id), y = water_pb, 
             color = water_depth_label, 
             shape = water_depth_label)) +
  geom_jitter(width = 0.1, height = 0, size = 2.5) +
  geom_hline(yintercept = 2, color = "red", linetype = "dashed", linewidth = 1) +  # red line
  annotate("text", x = "81", y = 2, label = "Chronic Aquatic Hazard Level", 
           vjust = -0.5, hjust = 1.02, color = "red", size = 4) +
  scale_y_log10() +
  scale_color_manual(values = c(
    "Surface water" = "lightblue",
    "Bottom water" = "blue",
    "Pore water" = "black"),
  limits = c("Surface water", "Bottom water", "Pore water")  # sets order
  ) +
  scale_shape_manual(values = c(
    "Surface water" = 16,  # circle
    "Bottom water" = 17,   # triangle
    "Pore water" = 15     # square
  ),
  limits = c("Surface water", "Bottom water", "Pore water")  # sets order
  ) +
  labs(
    x = "Pond ID",  # swapped labels
    y = "Water Pb concentration",
    color = NULL,
    shape = NULL
  ) +
  guides(color = guide_legend(ncol = 1),
         shape = guide_legend(ncol = 1)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_blank(),  # removes major horizontal gridlines
        panel.grid.minor.y = element_blank() ) 

ggsave("~/Desktop/Landscape-level Pond Analysis/plots/Fig1.jpg", width= 9, height = 5)
