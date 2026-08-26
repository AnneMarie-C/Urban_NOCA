# Script based on Kristen Lalla, Vanessa Poirier and Barbara Frei

#.libPaths("D:/R Packages")

library(lme4)
library(ggplot2)
library(car)
library(tidyverse)
library(corrplot)
library(MuMIn)
library(cowplot)
library(ResourceSelection)
library(corrplot)
library(dplyr)
library(performance) # for Nakagawa conditional/marginal R2
library(partR2)
library(splines)
library(pROC)
library(logistf)

# =============================================================================================== #
# 1. LOADING DATA & CHECKING FOR COLINEARITY FOR ALL YEARS COMBINED                               #
# =============================================================================================== #
#### Load data and check for colinearity ####

# Anne-Marie's filepaths.
obs.2223 <- read.csv("03_processed_data/noca_obs_df_mea_sd_2223.csv")
random.2223 <- read.csv("03_processed_data/noca_ran_df_mea_sd_2223.csv")

obs.2324 <- read.csv("03_processed_data/noca_obs_df_mea_sd_2324.csv")
random.2324 <- read.csv("03_processed_data/noca_ran_df_mea_sd_2324.csv")

obs.2425 <- read.csv("03_processed_data/noca_obs_df_mea_sd_2425.csv")
random.2425 <- read.csv("03_processed_data/noca_ran_df_mea_sd_2425.csv")


random.2223 <- random.2223 %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)
random.2324 <- random.2324 %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)
random.2425 <- random.2425 %>% 
  dplyr::rename(Lat = randLat,
                Lon = randLon)

# Combine used and available data
obs.2223 <- obs.2223 %>% 
  mutate(Used = 1)
random.2223 <- random.2223 %>% 
  mutate(Used = 0)
# Combine used and available data
obs.2324 <- obs.2324 %>% 
  mutate(Used = 1)
random.2324 <- random.2324 %>% 
  mutate(Used = 0)
# Combine used and available data
obs.2425 <- obs.2425 %>% 
  mutate(Used = 1)
random.2425 <- random.2425 %>% 
  mutate(Used = 0)

obs.2223$ID <- as.factor(obs.2223$ID)
obs.2324$ID <- as.factor(obs.2324$ID)
obs.2425$ID <- as.factor(obs.2425$ID)

random.2223$ID <- as.factor(random.2223$ID)
random.2324$ID <- as.factor(random.2324$ID)
random.2425$ID <- as.factor(random.2425$ID)

obs.2223$Winter <- as.factor(obs.2223$Winter)
obs.2324$Winter <- as.factor(obs.2324$Winter)
obs.2425$Winter <- as.factor(obs.2425$Winter)

random.2223$Winter <- as.factor(random.2223$Winter)
random.2324$Winter <- as.factor(random.2324$Winter)
random.2425$Winter <- as.factor(random.2425$Winter)

obs.2223$Site <- as.factor(obs.2223$Site)
obs.2324$Site <- as.factor(obs.2324$Site)
obs.2425$Site <- as.factor(obs.2425$Site)

random.2223$Site <- as.factor(random.2223$Site)
random.2324$Site <- as.factor(random.2324$Site)
random.2425$Site <- as.factor(random.2425$Site)


rsf_data <- bind_rows(obs.2223, obs.2324, obs.2425, random.2223, random.2324, random.2425)

# =============================================================================================== #
# 1b. TESTING FOR CORRELATION BETWEEN BUFFER-SCALE VARIABLES                                     #
# =============================================================================================== #

# Select the key variables at each scale for correlation testing
# Do this for the full dataset and per site

buffer_vars <- rsf_data %>%
  dplyr::select(
    # Vegetation height mean
    Mean_Veg_Height_25_m, Mean_Veg_Height_50_m, Mean_Veg_Height_100_m, Mean_Veg_Height_200_m,
    # Vegetation height SD
    Sd_Veg_Height_25_m, Sd_Veg_Height_50_m, Sd_Veg_Height_100_m, Sd_Veg_Height_200_m,
    # Population density
    mean_Pop_Dens_25_m, mean_Pop_Dens_50_m, mean_Pop_Dens_100_m, mean_Pop_Dens_200_m,
    # Road length
    Road_Length_25_m, Road_Length_50_m, Road_Length_100_m, Road_Length_200_m
  )

# 1. Pearson correlation matrix
cor_matrix <- cor(buffer_vars, use = "complete.obs", method = "pearson")
print(round(cor_matrix, 2))

# 2. Visual correlation plot — look for dark red/blue blocks across scales
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         tl.cex = 0.7,
         addCoef.col = "black",
         number.cex = 0.5,
         title = "Correlation between buffer-scale variables",
         mar = c(0, 0, 2, 0))

# 3. Focus on cross-scale correlations for each variable type
# This is the most targeted way to answer the reviewer

cross_scale_cor <- function(data, var_prefix, scales = c("25_m", "50_m", "100_m", "200_m"), label) {
  vars <- paste0(var_prefix, scales)
  mat <- cor(data[, vars], use = "complete.obs")
  cat("\n--- Cross-scale correlations for", label, "---\n")
  print(round(mat, 3))
}

cross_scale_cor(rsf_data, "Mean_Veg_Height_", label = "Mean Vegetation Height")
cross_scale_cor(rsf_data, "Sd_Veg_Height_",   label = "SD Vegetation Height")
cross_scale_cor(rsf_data, "mean_Pop_Dens_",    label = "Population Density")
cross_scale_cor(rsf_data, "Road_Length_",      label = "Road Length")

# 4. Repeat per site for completeness (since your models are site-specific)
for (site_name in c("MBO", "BDU", "CON")) {
  site_data <- switch(site_name,
                      MBO = rsf_data_MBO,
                      BDU = rsf_data_BDU,
                      CON = rsf_data_CON
  )
  cat("\n========== Site:", site_name, "==========\n")
  cross_scale_cor(site_data, "Mean_Veg_Height_", label = "Mean Veg Height")
  cross_scale_cor(site_data, "Sd_Veg_Height_",   label = "SD Veg Height")
  cross_scale_cor(site_data, "mean_Pop_Dens_",    label = "Population Density")
  cross_scale_cor(site_data, "Road_Length_",      label = "Road Length")
}


# # Let's first try with all the data for all the years.
# rsf_model_25 <- glm(Used ~ Mean_Veg_Height_25_m + 
#                       Sd_Veg_Height_25_m + 
#                       mean_Pop_Dens_25_m + 
#                       Road_Length_25_m + 
#                       Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data)
# rsf_model_50 <- glm(Used ~ Mean_Veg_Height_50_m + 
#                       Sd_Veg_Height_50_m + 
#                       mean_Pop_Dens_50_m + 
#                       Road_Length_50_m + 
#                       Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data)
# rsf_model_100 <- glm(Used ~ Mean_Veg_Height_100_m + 
#                        Sd_Veg_Height_100_m + 
#                        mean_Pop_Dens_100_m + 
#                        Road_Length_100_m + 
#                        Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data)
# rsf_model_200 <- glm(Used ~ Mean_Veg_Height_200_m + 
#                        Sd_Veg_Height_200_m + 
#                        mean_Pop_Dens_200_m + 
#                        Road_Length_200_m + 
#                        Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data)
# 
# summary(rsf_model_25) # AIC = 2472.4 <---
# summary(rsf_model_50) # AIC = 2523.3
# summary(rsf_model_100) # AIC = 2557
# summary(rsf_model_200) # AIC = 2552
# # Best model is the 25 meters model.
# vif(rsf_model_25) 
# vif(rsf_model_50)
# vif(rsf_model_100)
# vif(rsf_model_200) # All VIF < 10.
# 
# 
# 
# 
# ## prepare data by scaling
# noca_data_scaled.df <- rsf_data
# noca_data_scaled.df[,c(8:76)]  <-  scale(noca_data_scaled.df[,c(8:76)], center = TRUE, scale = TRUE)
# noca_data_scaled.df$ID <- as.factor(noca_data_scaled.df$ID)
# noca_data_scaled.df$Winter <- as.factor(noca_data_scaled.df$Winter)
# noca_data_scaled.df$Site <- as.factor(noca_data_scaled.df$Site)
# 
# # Candidate models for all of the data, using 25m buffer.
# noca_data_full_25 <- glmer(Used ~ Mean_Veg_Height_25_m + 
#                              Sd_Veg_Height_25_m + 
#                              mean_Pop_Dens_25_m + 
#                              Road_Length_25_m + 
#                              Dist_to_Nearest_Feeder_m +
#                              (1|Winter) + 
#                              (1|ID) + 
#                              (1|Site), family = binomial(link = "logit"), data = noca_data_scaled.df)
# 
# noca_tree_height <- glmer(Used ~ Mean_Veg_Height_25_m + 
#                             (1|Winter) + 
#                             (1|ID) + 
#                             (1|Site), family = binomial(link = "logit"), data = noca_data_scaled.df)
# 
# noca_tree_div <- glmer(Used ~ Sd_Veg_Height_25_m + 
#                          (1|Winter) + 
#                          (1|ID) + 
#                          (1|Site), family = binomial(link = "logit"), data = noca_data_scaled.df)
# 
# noca_human_model <- glmer(Used ~ mean_Pop_Dens_25_m + 
#                             Road_Length_25_m + 
#                             Dist_to_Nearest_Feeder_m +
#                             (1|Winter) + 
#                             (1|ID) + 
#                             (1|Site), family = binomial(link = "logit"), data = noca_data_scaled.df)
# 
# noca_data_null_25 <- glmer(Used ~ (1|Winter) + 
#                              (1|ID) + 
#                              (1|Site), family = binomial(link = "logit"), data = noca_data_scaled.df)
# 
# # Comparing candidate set - the full model is the best one.
# noca_data_AIC <- AIC(noca_data_full_25, noca_tree_height, noca_tree_div, noca_human_model, noca_data_null_25)
# noca_data_AIC
# # Full model is the best with AIC = 2451.912
# summary(noca_data_full_25)
# 
# 
# # summary of each model
# summary(noca_data_full_25)
# summary(noca_tree_height)
# summary(noca_tree_div)
# summary(noca_data_null_25)
# 
# # compare marginal vs conditional
# r.squaredGLMM(noca_data_full_25, nullfx =noca_data_null_25)
# 
# #H-S GOF test 
# hoslem.test(rsf_data$Used, fitted(noca_data_full_25), g = 10)
# hoslem.test(rsf_data$Used, fitted(noca_tree_height), g = 10)
# hoslem.test(rsf_data$Used, fitted(noca_tree_div), g = 10)
# hoslem.test(rsf_data$Used, fitted(noca_human_model), g = 10)
# # Every p-value < 0.05 meaning no model fits the data.


# =============================================================================================== #
# 2. INVESTIGATING HABITAT SELECTION SCALE FOR 3 SITES                                            #
# =============================================================================================== #
# Prepare 4 models for each of the scale and each site including the vegetation variables.

rsf_data_MBO <- rsf_data %>% 
  filter(Site %in% c(1,2))

rsf_data_BDU <- rsf_data %>% 
  filter(Site %in% c(3,4))

rsf_data_CON <- rsf_data %>% 
  filter(Site %in% c(5,6))

# MBO.
# Creating models including all the variables, at the 4 different scales.

rsf_model_MBO_25 <- glm(Used ~ Mean_Veg_Height_25_m + 
                          Sd_Veg_Height_25_m + 
                          mean_Pop_Dens_25_m + 
                          Road_Length_25_m +
                          Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_MBO)
rsf_model_MBO_50 <- glm(Used ~ Mean_Veg_Height_50_m + 
                          Sd_Veg_Height_50_m + 
                          mean_Pop_Dens_50_m + 
                          Road_Length_50_m +
                          Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_MBO)
rsf_model_MBO_100 <- glm(Used ~ Mean_Veg_Height_100_m + 
                           Sd_Veg_Height_100_m + 
                           mean_Pop_Dens_100_m + 
                           Road_Length_100_m +
                           Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_MBO)
rsf_model_MBO_200 <- glm(Used ~ Mean_Veg_Height_200_m + 
                           Sd_Veg_Height_200_m + 
                           mean_Pop_Dens_200_m + 
                           Road_Length_200_m +
                           Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_MBO)

# Comparing AICs for all of the models.
noca_MBO_AIC <- AIC(rsf_model_MBO_25, rsf_model_MBO_50, rsf_model_MBO_100, rsf_model_MBO_200)
noca_MBO_AIC
# 25 meters = 1456.957
# 50 meters = 1439.941 <--- Over 10 AIC point difference.
# 100 meters = 1493.194
# 200 meters = 1486.899

library(AICcmodavg) # more info AICc

models <- list(rsf_model_MBO_25, rsf_model_MBO_50, rsf_model_MBO_100, rsf_model_MBO_200) #define list of models
mod.names <- c('25 m', '50 m', '100 m', '200 m') #specify model names

aictab(cand.set = models, modnames = mod.names)

model_list <- list(rsf_model_MBO_25, rsf_model_MBO_50, rsf_model_MBO_100, rsf_model_MBO_200)
model_names <- c('25 m', '50 m', '100 m', '200 m')

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))

# Build AIC table
aic_table <- data.frame(
  Model = model_names,
  AIC = sapply(model_list, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Substantial = Delta_AIC <= 2
  )

# Preserve original plotting order
aic_table$Model <- factor(aic_table$Model, levels = model_names)

# Identify best model
min_aic <- min(aic_table$AIC)
best_model <- aic_table$Model[which.min(aic_table$AIC)]

# Plot
aic_scale_mbo <- ggplot(aic_table, aes(x = Model, y = AIC)) +
  # Grey shaded substantial-support area
  annotate("rect",
           xmin = 0.5, xmax = length(model_names) + 0.5,
           ymin = min_aic, ymax = min_aic + 2,
           alpha = 0.2, fill = "grey") +
  # Black line and points
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  # Labels
  labs(
    title = "AIC Comparison of Habitat Association Scale Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_scale_mbo

# summary(rsf_model_MBO_25)
summary(rsf_model_MBO_50) # <---
# summary(rsf_model_MBO_100)
# summary(rsf_model_MBO_200)
# The best model is the 50 meters one.

# vif(rsf_model_MBO_25)
vif(rsf_model_MBO_50) # For the best model, the 50 meters one, all VIFs < 10.
# vif(rsf_model_MBO_100)
# vif(rsf_model_MBO_200)


# BDU
rsf_model_BDU_25 <- glm(Used ~ Mean_Veg_Height_25_m + 
                          Sd_Veg_Height_25_m + 
                          mean_Pop_Dens_25_m + 
                          Road_Length_25_m +
                          Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_BDU)
rsf_model_BDU_50 <- glm(Used ~ Mean_Veg_Height_50_m + 
                          Sd_Veg_Height_50_m + 
                          mean_Pop_Dens_50_m + 
                          Road_Length_50_m +
                          Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_BDU)
rsf_model_BDU_100 <- glm(Used ~ Mean_Veg_Height_100_m + 
                           Sd_Veg_Height_100_m + 
                           mean_Pop_Dens_100_m + 
                           Road_Length_100_m +
                           Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_BDU)
rsf_model_BDU_200 <- glm(Used ~ Mean_Veg_Height_200_m + 
                           Sd_Veg_Height_200_m + 
                           mean_Pop_Dens_200_m + 
                           Road_Length_200_m +
                           Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_BDU)

# Comparing AICs for all of the models.
noca_BDU_AIC <- AIC(rsf_model_BDU_25, rsf_model_BDU_50, rsf_model_BDU_100, rsf_model_BDU_200)
noca_BDU_AIC
# 25 meters = 604.1074 <--- Over 10 AIC point difference.
# 50 meters = 688.6257 
# 100 meters = 700.4956
# 200 meters = 711.8367

models <- list(rsf_model_BDU_25, rsf_model_BDU_50, rsf_model_BDU_100, rsf_model_BDU_200) #define list of models
mod.names <- c('25 m', '50 m', '100 m', '200 m') #specify model names

aictab(cand.set = models, modnames = mod.names)

model_list <- list(rsf_model_BDU_25, rsf_model_BDU_50, rsf_model_BDU_100, rsf_model_BDU_200)
model_names <- c('25 m', '50 m', '100 m', '200 m')

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))

# Build AIC table
aic_table <- data.frame(
  Model = model_names,
  AIC = sapply(model_list, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Substantial = Delta_AIC <= 2
  )

# Preserve original plotting order
aic_table$Model <- factor(aic_table$Model, levels = model_names)

# Identify best model
min_aic <- min(aic_table$AIC)
best_model <- aic_table$Model[which.min(aic_table$AIC)]

# Plot
aic_scale_bdu <- ggplot(aic_table, aes(x = Model, y = AIC)) +
  # Grey shaded substantial-support area
  annotate("rect",
           xmin = 0.5, xmax = length(model_names) + 0.5,
           ymin = min_aic, ymax = min_aic + 2,
           alpha = 0.2, fill = "grey") +
  # Black line and points
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  # Labels
  labs(
    title = "AIC Comparison of Habitat Association Scale Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_scale_bdu

summary(rsf_model_BDU_25) # <---
# summary(rsf_model_BDU_50) 
# summary(rsf_model_BDU_100) 
# summary(rsf_model_BDU_200) 


vif(rsf_model_BDU_25) # For the best model, the 25 meters one, all VIFs < 10.
# vif(rsf_model_BDU_50)
# vif(rsf_model_BDU_100)
# vif(rsf_model_BDU_200)


# CON
rsf_model_CON_25 <- glm(Used ~ Mean_Veg_Height_25_m + 
                          Sd_Veg_Height_25_m + 
                          mean_Pop_Dens_25_m + 
                          Road_Length_25_m +
                          Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_CON)
rsf_model_CON_50 <- glm(Used ~ Mean_Veg_Height_50_m + 
                          Sd_Veg_Height_50_m + 
                          mean_Pop_Dens_50_m + 
                          Road_Length_50_m +
                          Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_CON)
rsf_model_CON_100 <- glm(Used ~ Mean_Veg_Height_100_m + 
                           Sd_Veg_Height_100_m + 
                           mean_Pop_Dens_100_m + 
                           Road_Length_100_m +
                           Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_CON)
rsf_model_CON_200 <- glm(Used ~ Mean_Veg_Height_200_m + 
                           Sd_Veg_Height_200_m + 
                           mean_Pop_Dens_200_m + 
                           Road_Length_200_m +
                           Dist_to_Nearest_Feeder_m, family = binomial(link = "logit"), data = rsf_data_CON)

# Comparing AICs for all of the models.
noca_CON_AIC <- AIC(rsf_model_CON_25, rsf_model_CON_50, rsf_model_CON_100, rsf_model_CON_200)
noca_CON_AIC
# 25 meters = 434.1834 <--- Over 10 AIC point difference.
# 50 meters = 463.2899 
# 100 meters = 460.0303
# 200 meters = 469.5799

models <- list(rsf_model_CON_25, rsf_model_CON_50, rsf_model_CON_100, rsf_model_CON_200) #define list of models
mod.names <- c('25 m', '50 m', '100 m', '200 m') #specify model names

aictab(cand.set = models, modnames = mod.names)

model_list <- list(rsf_model_CON_25, rsf_model_CON_50, rsf_model_CON_100, rsf_model_CON_200)
model_names <- c('25 m', '50 m', '100 m', '200 m')

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))

# Build AIC table
aic_table <- data.frame(
  Model = model_names,
  AIC = sapply(model_list, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Substantial = Delta_AIC <= 2
  )

# Preserve original plotting order
aic_table$Model <- factor(aic_table$Model, levels = model_names)

# Identify best model
min_aic <- min(aic_table$AIC)
best_model <- aic_table$Model[which.min(aic_table$AIC)]

# Plot
aic_scale_con <- ggplot(aic_table, aes(x = Model, y = AIC)) +
  # Grey shaded substantial-support area
  annotate("rect",
           xmin = 0.5, xmax = length(model_names) + 0.5,
           ymin = min_aic, ymax = min_aic + 2,
           alpha = 0.2, fill = "grey") +
  # Black line and points
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  # Labels
  labs(
    title = "AIC Comparison of Habitat Association Scale Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_scale_con

summary(rsf_model_CON_25) # <---
# summary(rsf_model_CON_50)
# summary(rsf_model_CON_100)
# summary(rsf_model_CON_200)

vif(rsf_model_CON_25) # For the best model, the 25 meters one, all VIFs < 10.
# vif(rsf_model_CON_50)
# vif(rsf_model_CON_100)
# vif(rsf_model_CON_200)

# |
# | In summary, for MBO the best scale is 50 meters, and for both BDU and CON the best scale is
# | 25 meters.
# |

# =============================================================================================== #
# 3. INVESTIGATING CANDIDATE MODELS FOR EACH LOCATIONS                                            #
# =============================================================================================== #

## prepare data by scaling
noca_MBO_scaled.df <- rsf_data_MBO
noca_MBO_scaled.df[,c(8:76)]  <-  scale(noca_MBO_scaled.df[,c(8:76)], center = TRUE, scale = TRUE)
noca_MBO_scaled.df$ID <- as.factor(noca_MBO_scaled.df$ID)

noca_BDU_scaled.df <- rsf_data_BDU
noca_BDU_scaled.df[,c(8:76)]  <-  scale(noca_BDU_scaled.df[,c(8:76)], center = TRUE, scale = TRUE)
noca_BDU_scaled.df$ID <- as.factor(noca_BDU_scaled.df$ID)

noca_CON_scaled.df <- rsf_data_CON
noca_CON_scaled.df[,c(8:76)]  <-  scale(noca_CON_scaled.df[,c(8:76)], center = TRUE, scale = TRUE)
noca_CON_scaled.df$ID <- as.factor(noca_CON_scaled.df$ID)

# Candidate models for MBO using 50m buffer:
noca_MBO_full_50 <- glmer(Used ~ Mean_Veg_Height_50_m + 
                            Sd_Veg_Height_50_m +
                            mean_Pop_Dens_50_m + 
                            Road_Length_50_m + 
                            Dist_to_Nearest_Feeder_m +
                            (1|Winter) + 
                            (1|ID) + 
                            (1|Site), family = binomial(link = "logit"), data = noca_MBO_scaled.df)

noca_tree_height_MBO <- glmer( Used ~ Mean_Veg_Height_50_m + 
                                 (1|Winter) + 
                                 (1|ID) + 
                                 (1|Site), family = binomial(link = "logit"), data = noca_MBO_scaled.df)

noca_tree_div_MBO <- glmer( Used ~ Sd_Veg_Height_50_m + 
                              (1|Winter) + 
                              (1|ID) + 
                              (1|Site), family = binomial(link = "logit"), data = noca_MBO_scaled.df)

noca_human_MBO <- glmer(Used ~ mean_Pop_Dens_50_m + 
                          Road_Length_50_m + 
                          Dist_to_Nearest_Feeder_m +
                          (1|Winter) + 
                          (1|ID) + 
                          (1|Site), family = binomial(link = "logit"), data = noca_MBO_scaled.df)

noca_MBO_null_50 <- glmer(Used ~ (1|Winter) + (1|ID) + (1|Site), data = noca_MBO_scaled.df, family = binomial)

# Comparing candidate set - the full model is the best one.
noca_MBO_AIC <- AIC(noca_MBO_full_50, noca_tree_height_MBO, noca_tree_div_MBO, noca_human_MBO, noca_MBO_null_50)
noca_MBO_AIC
# Full = 1445.941 <---
# Tree height = 1584.254
# Tree div = 1598.544
# Human = 1586.486
# Null = 1616.101

model_list <- list(noca_MBO_full_50, noca_tree_height_MBO, noca_tree_div_MBO, noca_human_MBO, noca_MBO_null_50)
model_names <- c('Full_50m', 'Tree_Height_50m', 'Tree_Diversity_50m', 'Human_50m', 'Null_50m')

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))

# Build AIC table
aic_table <- data.frame(
  Model = model_names,
  AIC = sapply(model_list, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Substantial = Delta_AIC <= 2
  )

# Preserve original plotting order
aic_table$Model <- factor(aic_table$Model, levels = model_names)

# Identify best model
min_aic <- min(aic_table$AIC)
best_model <- aic_table$Model[which.min(aic_table$AIC)]

# Plot
aic_candidate_mbo <- ggplot(aic_table, aes(x = Model, y = AIC)) +
  # Grey shaded substantial-support area
  annotate("rect",
           xmin = 0.5, xmax = length(model_names) + 0.5,
           ymin = min_aic, ymax = min_aic + 2,
           alpha = 0.2, fill = "grey") +
  # Black line and points
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  # Labels
  labs(
    title = "AIC Comparison of Candidate Models MBO",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_candidate_mbo

# summary of each model
summary(noca_MBO_full_50)
# summary(noca_tree_height_MBO)
# summary(noca_tree_div_MBO)
# summary(noca_MBO_null_50)

# Get a null model to calculate deviance explained.
noca_MBO_null <- glmer(
  Used ~ 1 + (1 | Winter) + (1 | ID) + (1 | Site),
  family = binomial,
  data = noca_MBO_scaled.df
)
1-deviance(noca_MBO_full_50)/deviance(noca_MBO_null)

# compare marginal vs conditional
r.squaredGLMM(noca_MBO_full_50, nullfx =noca_MBO_null_50)

#H-S GOF test 
hoslem.test(rsf_data_MBO$Used, fitted(noca_MBO_full_50), g = 10) # p-value < 0.05, doesn't fit data.
# hoslem.test(rsf_data_MBO$Used, fitted(noca_tree_height_MBO), g = 10)
# hoslem.test(rsf_data_MBO$Used, fitted(noca_tree_div_MBO), g = 10)
# hoslem.test(rsf_data_MBO$Used, fitted(noca_human_MBO), g = 10)

# Candidate models for BDU using 25 m buffer:
noca_BDU_full_25 <- glmer(Used ~ Mean_Veg_Height_25_m + 
                            Sd_Veg_Height_25_m +
                            mean_Pop_Dens_25_m + 
                            Road_Length_25_m + 
                            Dist_to_Nearest_Feeder_m +
                            (1|Winter) + 
                            (1|ID) + 
                            (1|Site),
                          family = binomial(link = "logit"), data = noca_BDU_scaled.df)

noca_tree_height_BDU <- glmer(Used ~ Mean_Veg_Height_25_m + 
                                (1|Winter) + 
                                (1|ID) + 
                                (1|Site), family = binomial(link = "logit"), data = noca_BDU_scaled.df)

noca_tree_div_BDU <- glmer( Used ~ Sd_Veg_Height_25_m + 
                              (1|Winter) + 
                              (1|ID) + 
                              (1|Site), family = binomial(link = "logit"), data = noca_BDU_scaled.df)

noca_human_BDU <- glmer(Used ~ mean_Pop_Dens_25_m + 
                          Road_Length_25_m + 
                          Dist_to_Nearest_Feeder_m +
                          (1|Winter) + 
                          (1|ID) + 
                          (1|Site), family = binomial(link = "logit"), data = noca_BDU_scaled.df)

noca_BDU_null_25 <- glmer(Used ~ (1|Winter) + (1|ID) + (1|Site), data = noca_BDU_scaled.df, family = binomial)

# Comparing candidate set - the tree height model is the best one.
noca_BDU_AIC <- AIC(noca_BDU_full_25, noca_tree_height_BDU, noca_tree_div_BDU, noca_human_BDU, noca_BDU_null_25)
noca_BDU_AIC
# Full = 596.6455 <-
# Tree height = 727.5336
# Tree div = 724.9683
# Human = 596.0550 <---
# Null = 726.1005

model_list <- list(noca_BDU_full_25, noca_tree_height_BDU, noca_tree_div_BDU, noca_human_BDU, noca_BDU_null_25)
model_names <- c('Full_25m', 'Tree_Height_25m', 'Tree_Diversity_25m', 'Human_25m', 'Null_25m')

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))

# Build AIC table
aic_table <- data.frame(
  Model = model_names,
  AIC = sapply(model_list, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Substantial = Delta_AIC <= 2
  )

# Preserve original plotting order
aic_table$Model <- factor(aic_table$Model, levels = model_names)

# Identify best model
min_aic <- min(aic_table$AIC)
best_model <- aic_table$Model[which.min(aic_table$AIC)]

# Plot
aic_candidate_bdu <- ggplot(aic_table, aes(x = Model, y = AIC)) +
  # Grey shaded substantial-support area
  annotate("rect",
           xmin = 0.5, xmax = length(model_names) + 0.5,
           ymin = min_aic, ymax = min_aic + 2,
           alpha = 0.2, fill = "grey") +
  # Black line and points
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  # Labels
  labs(
    title = "AIC Comparison of Candidate Models BDU",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_candidate_bdu

# Less than 1 AICc point difference, we may consider doing a model average.

# Average parameter estimates across top models
avg_model <- model.avg(noca_BDU_full_25, noca_human_BDU)

# summary of each model
summary(noca_BDU_full_25)
# summary(noca_tree_height_BDU)
# summary(noca_tree_div_BDU)
summary(noca_human_BDU)

# Get a null model to calculate deviance explained.
noca_BDU_null <- glmer(
  Used ~ 1 + (1 | Winter) + (1 | ID) + (1 | Site),
  family = binomial,
  data = noca_BDU_scaled.df
)
1-deviance(noca_human_BDU)/deviance(noca_BDU_null)

# summary(noca_BDU_null_25)
summary(avg_model)

# compare marginal vs conditional
r.squaredGLMM(noca_human_BDU, nullfx =noca_BDU_null_25)

#H-S GOF test 
hoslem.test(rsf_data_BDU$Used, fitted(noca_BDU_full_25), g = 10)
hoslem.test(rsf_data_BDU$Used, fitted(noca_tree_height_BDU), g = 10)
hoslem.test(rsf_data_BDU$Used, fitted(noca_tree_div_BDU), g = 10)
hoslem.test(rsf_data_BDU$Used, fitted(noca_human_BDU), g = 10) 
# Human model has p-value = 0.5791 > 0.05 so it does fit the data.

# Candidate models for CON using 25m buffer:
noca_CON_full_25 <- glmer(Used ~ Mean_Veg_Height_25_m + 
                            Sd_Veg_Height_25_m +
                            mean_Pop_Dens_25_m + 
                            Road_Length_25_m + 
                            Dist_to_Nearest_Feeder_m +
                            (1|Winter) + 
                            (1|ID) + 
                            (1|Site), family = binomial(link = "logit"), data = noca_CON_scaled.df)

noca_tree_height_CON <- glmer( Used ~ Mean_Veg_Height_25_m + 
                                 (1|Winter) + 
                                 (1|ID) + 
                                 (1|Site), family = binomial(link = "logit"), data = noca_CON_scaled.df)

noca_tree_div_CON <- glmer( Used ~ Sd_Veg_Height_25_m + 
                              (1|Winter) + 
                              (1|ID) + 
                              (1|Site), family = binomial(link = "logit"), data = noca_CON_scaled.df)

noca_human_CON <- glmer(Used ~ mean_Pop_Dens_25_m + 
                          Road_Length_25_m + 
                          Dist_to_Nearest_Feeder_m + 
                          (1|Winter) + 
                          (1|ID) + 
                          (1|Site), family = binomial(link = "logit"), data = noca_CON_scaled.df)

noca_CON_null_25 <- glmer(Used ~ (1|Winter) + (1|ID) + (1|Site), data = noca_CON_scaled.df, family = binomial)

# Comparing candidate set - the full model is the best one.
noca_CON_AIC <- AIC(noca_CON_full_25, noca_tree_height_CON, noca_tree_div_CON, noca_human_CON, noca_CON_null_25)
noca_CON_AIC
# Full = 440.1834 <--- Over 2 AICc point difference.
# Tree height = 486.8834
# Tree div = 485.9206
# Human = 442.6474 <-
# Null = 484.8853

model_list <- list(noca_CON_full_25, noca_tree_height_CON, noca_tree_div_CON, noca_human_CON, noca_CON_null_25)
model_names <- c('Full_25m', 'Tree_Height_25m', 'Tree_Diversity_25m', 'Human_25m', 'Null_25m')

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))

# Build AIC table
aic_table <- data.frame(
  Model = model_names,
  AIC = sapply(model_list, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Substantial = Delta_AIC <= 2
  )

# Preserve original plotting order
aic_table$Model <- factor(aic_table$Model, levels = model_names)

# Identify best model
min_aic <- min(aic_table$AIC)
best_model <- aic_table$Model[which.min(aic_table$AIC)]

# Plot
aic_candidate_con <- ggplot(aic_table, aes(x = Model, y = AIC)) +
  # Grey shaded substantial-support area
  annotate("rect",
           xmin = 0.5, xmax = length(model_names) + 0.5,
           ymin = min_aic, ymax = min_aic + 2,
           alpha = 0.2, fill = "grey") +
  # Black line and points
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  # Labels
  labs(
    title = "AIC Comparison of Candidate Models CON",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_candidate_con

# Plot all AIC graphs
(aic_scale_mbo|aic_scale_bdu|aic_scale_con)

(aic_candidate_mbo|aic_candidate_bdu|aic_candidate_con)

# summary of each model
summary(noca_CON_full_25)
# Get a null model to calculate deviance explained.
noca_CON_null <- glmer(
  Used ~ 1 + (1 | Winter) + (1 | ID) + (1 | Site),
  family = binomial,
  data = noca_CON_scaled.df
)
1-deviance(noca_CON_full_25)/deviance(noca_CON_null)

# summary(noca_tree_height_CON)
# summary(noca_tree_div_CON)
# summary(noca_human_CON)
# summary(noca_CON_null_25)

# compare marginal vs conditional
r.squaredGLMM(noca_CON_full_25, nullfx =noca_CON_null_25)

#H-S GOF test 
hoslem.test(rsf_data_CON$Used, fitted(noca_CON_full_25), g = 10) # p-value < 0.05, doesn't fit data.
# hoslem.test(rsf_data_CON$Used, fitted(noca_tree_height_CON), g = 10)
# hoslem.test(rsf_data_CON$Used, fitted(noca_tree_div_CON), g = 10)
# hoslem.test(rsf_data_CON$Used, fitted(noca_human_CON), g = 10)

# =============================================================================================== #
# 4. CREATING PREDICTION PLOTS FOR EACH MODEL AT EACH LOCATION                                    #
# =============================================================================================== #

## --- 4.1. Creating prediction plots for MBO --- ##

# Creating databases to predict from for each variable.
# Starting with the model for MBO: noca_MBO_full_50 with the data from
# noca_MBO_scaled.df

pred.df.mbo <- data.frame(     
  Mean_Veg_Height_50_m = rep(0, 500), 
  Sd_Veg_Height_50_m = 0,
  mean_Pop_Dens_50_m = 0,
  Road_Length_50_m = 0,
  Dist_to_Nearest_Feeder_m = 0,
  Site = sample(noca_MBO_scaled.df$Site,500),
  ID = sample(noca_MBO_scaled.df$ID,500))

# We use 2.33 because it represents the 99% quantile on scaled variables. This is to deal with
# outliers, anything above or below 2.33 will be ignored.

# cut_min <- function(x){ifelse(min(x) < -2.33, -2.33, min(x))}
# cut_max <- function(x){ifelse(max(x) > 2.33, 2.33, max(x))}
# Barbara suggested not doing using this for now.

# The uPredictions is what is used on the x axis for plotting (non-logged, non-scaled values).

Mean_Veg_pred.df.MBO <- pred.df.mbo %>% mutate(
  Mean_Veg_Height_50_m = seq(min(noca_MBO_scaled.df$Mean_Veg_Height_50_m), max(noca_MBO_scaled.df$Mean_Veg_Height_50_m), length.out=500),
  uPredictions = (Mean_Veg_Height_50_m * sd(rsf_data_MBO$Mean_Veg_Height_50_m) + mean(rsf_data_MBO$Mean_Veg_Height_50_m)))

Sd_Veg_Height_pred.df.MBO <- pred.df.mbo %>% mutate(
  Sd_Veg_Height_50_m = seq(min(noca_MBO_scaled.df$Sd_Veg_Height_50_m), max(noca_MBO_scaled.df$Sd_Veg_Height_50_m), length.out=500),
  uPredictions = (Sd_Veg_Height_50_m * sd(rsf_data_MBO$Sd_Veg_Height_50_m) + mean(rsf_data_MBO$Sd_Veg_Height_50_m)))

mean_Pop_Dens_pred.df.MBO <- pred.df.mbo %>% mutate(
 mean_Pop_Dens_50_m = seq(min(noca_MBO_scaled.df$mean_Pop_Dens_50_m), max(noca_MBO_scaled.df$mean_Pop_Dens_50_m), length.out=500),
 uPredictions = (mean_Pop_Dens_50_m * sd(rsf_data_MBO$mean_Pop_Dens_50_m) + mean(rsf_data_MBO$mean_Pop_Dens_50_m)))

Road_Length_pred.df.MBO <- pred.df.mbo %>% mutate(
  Road_Length_50_m = seq(min(noca_MBO_scaled.df$Road_Length_50_m), max(noca_MBO_scaled.df$Road_Length_50_m), length.out=500),
  uPredictions = (Road_Length_50_m * sd(rsf_data_MBO$Road_Length_50_m) + mean(rsf_data_MBO$Road_Length_50_m)))

Feeders_pred.df.MBO <- pred.df.mbo %>% mutate(
  Dist_to_Nearest_Feeder_m = seq(min(noca_MBO_scaled.df$Dist_to_Nearest_Feeder_m), max(noca_MBO_scaled.df$Dist_to_Nearest_Feeder_m), length.out=500),
  uPredictions = (Dist_to_Nearest_Feeder_m * sd(rsf_data_MBO$Dist_to_Nearest_Feeder_m) + mean(rsf_data_MBO$Dist_to_Nearest_Feeder_m)))


# Creating the prediction for Type (1 = Used, 0 = Random).
# We're adding se.fit to get the standard errors and create confidence intervals on the graphs
# later on.

Mean_Veg_pred2.df.MBO <-  predict(noca_MBO_full_50, newdata = Mean_Veg_pred.df.MBO, type="response", re.form = NA, se.fit = TRUE)

Sd_Veg_Height_pred2.df.MBO <-  predict(noca_MBO_full_50, newdata = Sd_Veg_Height_pred.df.MBO, type="response", re.form = NA, se.fit = TRUE)

mean_Pop_Dens_pred2.df.MBO <-  predict(noca_MBO_full_50, newdata = mean_Pop_Dens_pred.df.MBO, type="response", re.form = NA, se.fit = TRUE)

Road_Length_pred2.df.MBO <-  predict(noca_MBO_full_50, newdata = Road_Length_pred.df.MBO, type="response", re.form = NA, se.fit = TRUE)

Feeders_pred2.df.MBO <-  predict(noca_MBO_full_50, newdata = Feeders_pred.df.MBO, type="response", re.form = NA, se.fit = TRUE)

# Adding in both the prediction (fit) and the SE.

Mean_Veg_pred.df.MBO$Presence <- Mean_Veg_pred2.df.MBO$fit
Sd_Veg_Height_pred.df.MBO$Presence <- Sd_Veg_Height_pred2.df.MBO$fit
mean_Pop_Dens_pred.df.MBO$Presence <- mean_Pop_Dens_pred2.df.MBO$fit
Road_Length_pred.df.MBO$Presence <- Road_Length_pred2.df.MBO$fit
Feeders_pred.df.MBO$Presence <- Feeders_pred2.df.MBO$fit

Mean_Veg_pred.df.MBO$SE <- Mean_Veg_pred2.df.MBO$se.fit
Sd_Veg_Height_pred.df.MBO$SE <- Sd_Veg_Height_pred2.df.MBO$se.fit
mean_Pop_Dens_pred.df.MBO$SE <- mean_Pop_Dens_pred2.df.MBO$se.fit
Road_Length_pred.df.MBO$SE <- Road_Length_pred2.df.MBO$se.fit
Feeders_pred.df.MBO$SE <- Feeders_pred2.df.MBO$se.fit

# Using the SE to calculate the 95% CI from the geom_ribbon of the ggplot.
Mean_Veg_pred.df.MBO <- Mean_Veg_pred.df.MBO%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Sd_Veg_Height_pred.df.MBO <- Sd_Veg_Height_pred.df.MBO%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

mean_Pop_Dens_pred.df.MBO <- mean_Pop_Dens_pred.df.MBO%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Road_Length_pred.df.MBO <- Road_Length_pred.df.MBO%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Feeders_pred.df.MBO <- Feeders_pred.df.MBO%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

# Using the SE to calculate the 67% CI from the geom_ribbon of the ggplpot.
# Tried this after Kyle suggested it since the 95% ICs are very large.
Mean_Veg_pred.df.MBO <- Mean_Veg_pred.df.MBO %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Sd_Veg_Height_pred.df.MBO <- Sd_Veg_Height_pred.df.MBO %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

mean_Pop_Dens_pred.df.MBO <- mean_Pop_Dens_pred.df.MBO %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Road_Length_pred.df.MBO <- Road_Length_pred.df.MBO %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Feeders_pred.df.MBO <- Feeders_pred.df.MBO %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

# Plotting parameter effects of the full model at 50 meters for MBO.
summod<- summary(noca_MBO_full_50)

signifVariables <- rownames(summod$coefficients)[summod$coefficients[,4]<0.05]

dataList <- list(Mean_Veg_pred.df.MBO,  Sd_Veg_Height_pred.df.MBO, Feeders_pred.df.MBO) # only signif variables

Variables <- c("Mean Vegetation Height", "Sd Vegetation Height", "Distance to Nearest Feeder")

axisLabels <- c("Mean Vegetation Height (m)", "Standard Deviation of Vegetation Height", "Distance to Nearest Feeder (m)")

plotList <- list()

for(i in 1:length(dataList)){
  significant <- any(grepl(Variables[i], signifVariables))
  dat <- dataList[[i]]
  plotList[[i]] <- ggplot(data=dat, aes(x=uPredictions,y=Presence)) +
    geom_ribbon(data=dat, aes(x=uPredictions, ymin = LCL67, ymax = UCL67), 
                color = NA, fill = ifelse(significant, "#44AA99", "#44AA99"), alpha = 0.3) + 
    geom_line(colour = "white", size = 1.5)  + 
    scale_x_continuous(axisLabels[i]) +
    # scale_y_continuous(limits=c(0,1)) +
    theme(axis.title.x = element_text(size = 8), 
          axis.title.y = element_text(size = 8),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}

plot_grid(plotlist=plotList, ncol=3, labels=c("A)", "B)", "C)"))


### with observed points
library(ggpubr)

# Create 1 big dataframe using all the observation tables.
obs <- bind_rows(obs.2223, obs.2324, obs.2425)
random <- bind_rows(random.2223, random.2324, random.2425)

mean_veg_plot_mbo <- ggplot(data= Mean_Veg_pred.df.MBO, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#44AA99", alpha = 0.3) +  # confidence ribbon
  geom_line(colour = "#44AA99", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[1], n.breaks=4, limits = c(0,25)) +
  scale_y_continuous(n.breaks = 4, limits=c(-1.5,1.5), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = Mean_Veg_Height_50_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = Mean_Veg_Height_50_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

sd_veg_plot_mbo <- ggplot(data= Sd_Veg_Height_pred.df.MBO, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#44AA99", alpha = 0.3) +
  geom_line(colour = "#44AA99", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[2], n.breaks=4, limits = c(0,13)) +
  scale_y_continuous(n.breaks = 4, limits=c(-1.5,1.5), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = Sd_Veg_Height_50_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = Sd_Veg_Height_50_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

feeder_plot_mbo <- ggplot(data= Feeders_pred.df.MBO, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#44AA99", alpha = 0.3) +
  geom_line(colour = "#44AA99", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[3], n.breaks=4, limits = c(0,700)) +
  scale_y_continuous(n.breaks = 4, limits=c(-1.5,1.5), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = Dist_to_Nearest_Feeder_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = Dist_to_Nearest_Feeder_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

library(patchwork)

# Assigning plots tag manually before plotting together.

mean_veg_plot_mbo <- mean_veg_plot_mbo + ggtitle("A)")
sd_veg_plot_mbo <- sd_veg_plot_mbo + ggtitle("B)")
feeder_plot_mbo <- feeder_plot_mbo + ggtitle("C)")

# Then combining the 3 plots into one with a title.

prediction_plots_mbo <- wrap_elements(mean_veg_plot_mbo) + 
  wrap_elements(sd_veg_plot_mbo) + 
  wrap_elements(feeder_plot_mbo) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title = "McGill Bird Observatory Prediction Plots",
    theme = theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

## --- 4.2. Creating prediction plots for BDU --- ##

# Creating databases to predict from for each variable.
# Starting with the model for BDU: noca_human_BDU at 25 meters with the data from
# noca_BDU_scaled.df

pred.df.BDU <- data.frame(     
  Mean_Veg_Height_25_m = rep(0, 500), 
  Sd_Veg_Height_25_m = 0,
  mean_Pop_Dens_25_m = 0,
  Road_Length_25_m = 0,
  Dist_to_Nearest_Feeder_m = 0,
  Site = sample(noca_BDU_scaled.df$Site,500, replace = TRUE),
  ID = sample(noca_BDU_scaled.df$ID,500, replace = TRUE)) # Adding replace = TRUE shouldn't
                                                          # change results because we're not looking
                                                          # at their effects.

# We use 2.33 because it represents the 99% quantile on scaled variables. This is to deal with
# outliers, anything above or below 2.33 will be ignored.

# cut_min <- function(x){ifelse(min(x) < -2.33, -2.33, min(x))}
# cut_max <- function(x){ifelse(max(x) > 2.33, 2.33, max(x))}

# The uPredictions is what is used on the x axis for plotting (non-logged, non-scaled values).

Mean_Veg_pred.df.BDU <- pred.df.BDU %>% mutate(
  Mean_Veg_Height_25_m = seq(min(noca_BDU_scaled.df$Mean_Veg_Height_25_m), max(noca_BDU_scaled.df$Mean_Veg_Height_25_m), length.out=500),
  uPredictions = (Mean_Veg_Height_25_m * sd(rsf_data_BDU$Mean_Veg_Height_25_m) + mean(rsf_data_BDU$Mean_Veg_Height_25_m)))

Sd_Veg_Height_pred.df.BDU <- pred.df.BDU %>% mutate(
  Sd_Veg_Height_25_m = seq(min(noca_BDU_scaled.df$Sd_Veg_Height_25_m), max(noca_BDU_scaled.df$Sd_Veg_Height_25_m), length.out=500),
  uPredictions = (Sd_Veg_Height_25_m * sd(rsf_data_BDU$Sd_Veg_Height_25_m) + mean(rsf_data_BDU$Sd_Veg_Height_25_m)))

mean_Pop_Dens_pred.df.BDU <- pred.df.BDU %>% mutate(
  mean_Pop_Dens_25_m = seq(min(noca_BDU_scaled.df$mean_Pop_Dens_25_m), max(noca_BDU_scaled.df$mean_Pop_Dens_25_m), length.out=500),
  uPredictions = (mean_Pop_Dens_25_m * sd(rsf_data_BDU$mean_Pop_Dens_25_m) + mean(rsf_data_BDU$mean_Pop_Dens_25_m)))

Road_Length_pred.df.BDU <- pred.df.BDU %>% mutate(
  Road_Length_25_m = seq(min(noca_BDU_scaled.df$Road_Length_25_m), max(noca_BDU_scaled.df$Road_Length_25_m), length.out=500),
  uPredictions = (Road_Length_25_m * sd(rsf_data_BDU$Road_Length_25_m) + mean(rsf_data_BDU$Road_Length_25_m)))

Feeders_pred.df.BDU <- pred.df.BDU %>% mutate(
  Dist_to_Nearest_Feeder_m = seq(min(noca_BDU_scaled.df$Dist_to_Nearest_Feeder_m), max(noca_BDU_scaled.df$Dist_to_Nearest_Feeder_m), length.out=500),
  uPredictions = (Dist_to_Nearest_Feeder_m * sd(rsf_data_BDU$Dist_to_Nearest_Feeder_m) + mean(rsf_data_BDU$Dist_to_Nearest_Feeder_m)))


# Creating the prediction for Type (1 = Used, 0 = Random).

Mean_Veg_pred2.df.BDU <-  predict(noca_human_BDU, newdata = Mean_Veg_pred.df.BDU, type="response", re.form = NA, se.fit = TRUE)

Sd_Veg_Height_pred2.df.BDU <-  predict(noca_human_BDU, newdata = Sd_Veg_Height_pred.df.BDU, type="response", re.form = NA, se.fit = TRUE)

mean_Pop_Dens_pred2.df.BDU <-  predict(noca_human_BDU, newdata = mean_Pop_Dens_pred.df.BDU, type="response", re.form = NA, se.fit = TRUE)

Road_Length_pred2.df.BDU <-  predict(noca_human_BDU, newdata = Road_Length_pred.df.BDU, type="response", re.form = NA, se.fit = TRUE)

Feeders_pred2.df.BDU <-  predict(noca_human_BDU, newdata = Feeders_pred.df.BDU, type="response", re.form = NA, se.fit = TRUE)

# Adding in both the prediction (fit) and the SE.

Mean_Veg_pred.df.BDU$Presence <- Mean_Veg_pred2.df.BDU$fit
Sd_Veg_Height_pred.df.BDU$Presence <- Sd_Veg_Height_pred2.df.BDU$fit
mean_Pop_Dens_pred.df.BDU$Presence <- mean_Pop_Dens_pred2.df.BDU$fit
Road_Length_pred.df.BDU$Presence <- Road_Length_pred2.df.BDU$fit
Feeders_pred.df.BDU$Presence <- Feeders_pred2.df.BDU$fit

Mean_Veg_pred.df.BDU$SE <- Mean_Veg_pred2.df.BDU$se.fit
Sd_Veg_Height_pred.df.BDU$SE <- Sd_Veg_Height_pred2.df.BDU$se.fit
mean_Pop_Dens_pred.df.BDU$SE <- mean_Pop_Dens_pred2.df.BDU$se.fit
Road_Length_pred.df.BDU$SE <- Road_Length_pred2.df.BDU$se.fit
Feeders_pred.df.BDU$SE <- Feeders_pred2.df.BDU$se.fit

# Using the SE to calculate the 95% CI from the geom_ribbon of the ggplot.
Mean_Veg_pred.df.BDU <- Mean_Veg_pred.df.BDU %>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Sd_Veg_Height_pred.df.BDU <- Sd_Veg_Height_pred.df.BDU %>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

mean_Pop_Dens_pred.df.BDU <- mean_Pop_Dens_pred.df.BDU %>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Road_Length_pred.df.BDU <- Road_Length_pred.df.BDU %>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Feeders_pred.df.BDU <- Feeders_pred.df.BDU %>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

# Using the SE to calculate the 67% CI from the geom_ribbon of the ggplot.
Mean_Veg_pred.df.BDU <- Mean_Veg_pred.df.BDU %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Sd_Veg_Height_pred.df.BDU <- Sd_Veg_Height_pred.df.BDU %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

mean_Pop_Dens_pred.df.BDU <- mean_Pop_Dens_pred.df.BDU %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Road_Length_pred.df.BDU <- Road_Length_pred.df.BDU %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Feeders_pred.df.BDU <- Feeders_pred.df.BDU %>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

# Plotting parameter effects of the human model at 25 meters for BDU.

summod<- summary(noca_human_BDU)

signifVariables <- rownames(summod$coefficients)[summod$coefficients[,4]<0.05]

dataList <- list(mean_Pop_Dens_pred.df.BDU,  Road_Length_pred.df.BDU, Feeders_pred.df.BDU) # only signif variables

Variables <- c("Mean Population Density", "Mean Road Length", "Distance to Nearest Feeder")

axisLabels <- c("Mean Population Density (people per km2)", "Mean Road Length (m)", "Distance to Nearest Feeder (m)")

plotList <- list()

for(i in 1:length(dataList)){
  significant <- any(grepl(Variables[i], signifVariables))
  dat <- dataList[[i]]
  plotList[[i]] <- ggplot(data=dat, aes(x=uPredictions,y=Presence)) +
    geom_ribbon(data=dat, aes(x=uPredictions, ymin = LCL67, ymax = UCL67), 
                color = NA, fill = ifelse(significant, "#DDCC77", "#DDCC77"), alpha = 0.3) + 
    geom_line(colour = "white", size = 1.5)  + 
    scale_x_continuous(axisLabels[i]) +
    # scale_y_continuous(limits=c(0,1)) +
    theme(axis.title.x = element_text(size = 8), 
          axis.title.y = element_text(size = 8),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}
plot_grid(plotlist=plotList, ncol=3, labels=c("D)", "E)", "F)"))


### with observed points
library(ggpubr)

# Create 1 big dataframe using all the observation tables.
obs <- bind_rows(obs.2223, obs.2324, obs.2425)
random <- bind_rows(random.2223, random.2324, random.2425)

mean_pop_dens_plot_bdu <- ggplot(data= mean_Pop_Dens_pred.df.BDU, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#DDCC77", alpha = 0.3) +
  geom_line(colour = "#DDCC77", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[1], n.breaks=4, limits = c(0,1250)) +
  scale_y_continuous(n.breaks = 4, limits=c(-1.2,1.2), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = mean_Pop_Dens_25_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = mean_Pop_Dens_25_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

mean_road_plot_bdu <- ggplot(data= Road_Length_pred.df.BDU, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#DDCC77", alpha = 0.3) +
  geom_line(colour = "#DDCC77", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[2], n.breaks=4, limits = c(0,800)) +
  scale_y_continuous(n.breaks = 4, limits=c(-1,2.5), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = Road_Length_25_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = Road_Length_25_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

feeder_plot_bdu <- ggplot(data= Feeders_pred.df.BDU, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#DDCC77", alpha = 0.3) +
  geom_line(colour = "#DDCC77", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[3], n.breaks=4, limits = c(0,700)) +
  scale_y_continuous(n.breaks = 4, limits=c(-2,2), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = Dist_to_Nearest_Feeder_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = Dist_to_Nearest_Feeder_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

library(patchwork)

# Assigning plots tag manually before plotting together.

mean_pop_dens_plot_bdu <- mean_pop_dens_plot_bdu + ggtitle("D)")
mean_road_plot_bdu <- mean_road_plot_bdu + ggtitle("E)")
feeder_plot_bdu <- feeder_plot_bdu + ggtitle("F)")

# Then combining the 3 plots into one with a title.

prediction_plots_bdu <- wrap_elements(mean_pop_dens_plot_bdu) + 
  wrap_elements(mean_road_plot_bdu) + 
  wrap_elements(feeder_plot_bdu) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Town of Baie d'Urfé Prediction Plots",
    theme = theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

## --- 4.3. Creating prediction plots for CON --- ##

# Creating databases to predict from for each variable.
# Starting with the model for CON: noca_CON_full_25 with the data from noca_CON_scaled.df

pred.df.CON <- data.frame(     
  Mean_Veg_Height_25_m = rep(0, 500), 
  Sd_Veg_Height_25_m = 0,
  mean_Pop_Dens_25_m = 0,
  Road_Length_25_m = 0,
  Dist_to_Nearest_Feeder_m = 0,
  Site = sample(noca_BDU_scaled.df$Site,500, replace = TRUE),
  ID = sample(noca_BDU_scaled.df$ID,500, replace = TRUE)) # Adding replace = TRUE shouldn't
# change results because we're not looking
# at their effects.

# We use 2.33 because it represents the 99% quantile on scaled variables. This is to deal with
# outliers, anything above or below 2.33 will be ignored.

# cut_min <- function(x){ifelse(min(x) < -2.33, -2.33, min(x))}
# cut_max <- function(x){ifelse(max(x) > 2.33, 2.33, max(x))}

# The uPredictions is what is used on the x axis for plotting (non-logged, non-scaled values).

Mean_Veg_pred.df.CON <- pred.df.CON %>% mutate(
  Mean_Veg_Height_25_m = seq(min(noca_CON_scaled.df$Mean_Veg_Height_25_m), max(noca_CON_scaled.df$Mean_Veg_Height_25_m), length.out=500),
  uPredictions = (Mean_Veg_Height_25_m * sd(rsf_data_CON$Mean_Veg_Height_25_m) + mean(rsf_data_CON$Mean_Veg_Height_25_m)))

Sd_Veg_Height_pred.df.CON <- pred.df.CON %>% mutate(
  Sd_Veg_Height_25_m = seq(min(noca_CON_scaled.df$Sd_Veg_Height_25_m), max(noca_CON_scaled.df$Sd_Veg_Height_25_m), length.out=500),
  uPredictions = (Sd_Veg_Height_25_m * sd(rsf_data_CON$Sd_Veg_Height_25_m) + mean(rsf_data_CON$Sd_Veg_Height_25_m)))

mean_Pop_Dens_pred.df.CON <- pred.df.CON %>% mutate(
  mean_Pop_Dens_25_m = seq(min(noca_CON_scaled.df$mean_Pop_Dens_25_m), max(noca_CON_scaled.df$mean_Pop_Dens_25_m), length.out=500),
  uPredictions = (mean_Pop_Dens_25_m * sd(rsf_data_CON$mean_Pop_Dens_25_m) + mean(rsf_data_CON$mean_Pop_Dens_25_m)))

Road_Length_pred.df.CON <- pred.df.CON %>% mutate(
  Road_Length_25_m = seq(min(noca_CON_scaled.df$Road_Length_25_m), max(noca_CON_scaled.df$Road_Length_25_m), length.out=500),
  uPredictions = (Road_Length_25_m * sd(rsf_data_CON$Road_Length_25_m) + mean(rsf_data_CON$Road_Length_25_m)))

Feeders_pred.df.CON <- pred.df.CON %>% mutate(
  Dist_to_Nearest_Feeder_m = seq(min(noca_CON_scaled.df$Dist_to_Nearest_Feeder_m), max(noca_CON_scaled.df$Dist_to_Nearest_Feeder_m), length.out=500),
  uPredictions = (Dist_to_Nearest_Feeder_m * sd(rsf_data_CON$Dist_to_Nearest_Feeder_m) + mean(rsf_data_CON$Dist_to_Nearest_Feeder_m)))


# Creating the prediction for Type (1 = Used, 0 = Random).

Mean_Veg_pred2.df.CON <-  predict(noca_CON_full_25, newdata = Mean_Veg_pred.df.CON, type="response", re.form = NA, se.fit = TRUE)

Sd_Veg_Height_pred2.df.CON <-  predict(noca_CON_full_25, newdata = Sd_Veg_Height_pred.df.CON, type="response", re.form = NA, se.fit = TRUE)

mean_Pop_Dens_pred2.df.CON <-  predict(noca_CON_full_25, newdata = mean_Pop_Dens_pred.df.CON, type="response", re.form = NA, se.fit = TRUE)

Road_Length_pred2.df.CON <-  predict(noca_CON_full_25, newdata = Road_Length_pred.df.CON, type="response", re.form = NA, se.fit = TRUE)

Feeders_pred2.df.CON <-  predict(noca_CON_full_25, newdata = Feeders_pred.df.CON, type="response", re.form = NA, se.fit = TRUE)

# Adding in both the prediction (fit) and the SE.

Mean_Veg_pred.df.CON$Presence <- Mean_Veg_pred2.df.CON$fit
Sd_Veg_Height_pred.df.CON$Presence <- Sd_Veg_Height_pred2.df.CON$fit
mean_Pop_Dens_pred.df.CON$Presence <- mean_Pop_Dens_pred2.df.CON$fit
Road_Length_pred.df.CON$Presence <- Road_Length_pred2.df.CON$fit
Feeders_pred.df.CON$Presence <- Feeders_pred2.df.CON$fit

Mean_Veg_pred.df.CON$SE <- Mean_Veg_pred2.df.CON$se.fit
Sd_Veg_Height_pred.df.CON$SE <- Sd_Veg_Height_pred2.df.CON$se.fit
mean_Pop_Dens_pred.df.CON$SE <- mean_Pop_Dens_pred2.df.CON$se.fit
Road_Length_pred.df.CON$SE <- Road_Length_pred2.df.CON$se.fit
Feeders_pred.df.CON$SE <- Feeders_pred2.df.CON$se.fit

# Using the SE to calculate the 95% CI from the geom_ribbon of the ggplot.
Mean_Veg_pred.df.CON <- Mean_Veg_pred.df.CON%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Sd_Veg_Height_pred.df.CON <- Sd_Veg_Height_pred.df.CON%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

mean_Pop_Dens_pred.df.CON <- mean_Pop_Dens_pred.df.CON%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Road_Length_pred.df.CON <- Road_Length_pred.df.CON%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

Feeders_pred.df.CON <- Feeders_pred.df.CON%>% mutate(
  LCL95= Presence - 1.96 * SE, 
  UCL95= Presence + 1.96 * SE)

# Using the SE to calculate the 67% CI from the geom_ribbon of the ggplot.
Mean_Veg_pred.df.CON <- Mean_Veg_pred.df.CON%>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Sd_Veg_Height_pred.df.CON <- Sd_Veg_Height_pred.df.CON%>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

mean_Pop_Dens_pred.df.CON <- mean_Pop_Dens_pred.df.CON%>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Road_Length_pred.df.CON <- Road_Length_pred.df.CON%>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

Feeders_pred.df.CON <- Feeders_pred.df.CON%>% mutate(
  LCL67= Presence - 1.00 * SE, 
  UCL67= Presence + 1.00 * SE)

# Plotting parameter effects of the full model at 25 meters for CON.
summod<- summary(noca_CON_full_25)

signifVariables <- rownames(summod$coefficients)[summod$coefficients[,4]<0.05]

dataList <- list(Sd_Veg_Height_pred.df.CON,  Road_Length_pred.df.CON, Feeders_pred.df.CON) # only signif variables

Variables <- c("Sd Vegetation Height", "Mean Road Length", "Distance to Nearest Feeder")

axisLabels <- c("Standard Deviation of Vegetation Height", "Mean Road Length (m)", "Distance to Nearest Feeder (m)")

plotList <- list()

for(i in 1:length(dataList)){
  significant <- any(grepl(Variables[i], signifVariables))
  dat <- dataList[[i]]
  plotList[[i]] <- ggplot(data=dat, aes(x=uPredictions,y=Presence)) +
    geom_ribbon(data=dat, aes(x=uPredictions, ymin = LCL67, ymax = UCL67), 
                color = NA, fill = ifelse(significant, "#ed6859", "#ed9c93"), alpha = 0.3) + 
    geom_line(colour = "white", size = 1.5)  + 
    scale_x_continuous(axisLabels[i]) +
    # scale_y_continuous(limits=c(0,1)) +
    theme(axis.title.x = element_text(size = 8), 
          axis.title.y = element_text(size = 8),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}

plot_grid(plotlist=plotList, ncol=3, labels=c("G)", "H)", "I)"))


### with observed points
library(ggpubr)

# Create 1 big dataframe using all the observation tables.
obs <- bind_rows(obs.2223, obs.2324, obs.2425)
random <- bind_rows(random.2223, random.2324, random.2425)

sd_veg_plot_con <- ggplot(data= Sd_Veg_Height_pred.df.CON, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#ed9c93", alpha = 0.3) +
  geom_line(colour = "#ed6859", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[2], n.breaks= 4, limits = c(0,12)) +
  scale_y_continuous(n.breaks = 4, limits=c(-2,2.5), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = Sd_Veg_Height_50_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = Sd_Veg_Height_50_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

mean_road_plot_con <- ggplot(data= Road_Length_pred.df.CON, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#ed9c93", alpha = 0.3) +
  geom_line(colour = "#ed6859", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[2], n.breaks=4, limits = c(0,900)) +
  scale_y_continuous(n.breaks = 4, limits=c(-0.5,2), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = Road_Length_25_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = Road_Length_25_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

feeder_plot_con <- ggplot(data= Feeders_pred.df.CON, aes(x=uPredictions,y=Presence)) +
  geom_ribbon(aes(ymin = LCL67, ymax = UCL67), fill = "#ed9c93", alpha = 0.3) +
  geom_line(colour = "#ed6859", size = 1.5) + # or colour = ifelse(significant, "dodgerblue", "grey")
  scale_x_continuous(axisLabels[3], n.breaks=4, limits = c(0,400)) +
  scale_y_continuous(n.breaks = 4, limits=c(-1,2), labels = scales::number_format(accuracy = 1))+
  geom_point(data = obs, mapping = aes(x = Dist_to_Nearest_Feeder_m, y = 1), color = "gray40", alpha=0.5) + # observed points
  geom_point(data = random, mapping = aes(x = Dist_to_Nearest_Feeder_m, y = 0), color = "gray40", alpha=0.5) + # random points
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

library(patchwork)

# Assigning plots tag manually before plotting together.

sd_veg_plot_con <- sd_veg_plot_con + ggtitle("G)")
mean_road_plot_con <- mean_road_plot_con + ggtitle("H)")
feeder_plot_con <- feeder_plot_con + ggtitle("I)")

# Then combining the 3 plots into one with a title.

prediction_plots_con <- wrap_elements(sd_veg_plot_con) + 
  wrap_elements(mean_road_plot_con) + 
  wrap_elements(feeder_plot_con) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Concordia University Loyola Campus Prediction Plots",
    theme = theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

# Merging all the prediction plots into 1 figure.

all_prediction_plots <- wrap_elements(prediction_plots_mbo) + wrap_elements(prediction_plots_bdu) + 
  wrap_elements(prediction_plots_con) + 
  plot_layout(ncol = 1, guides = "collect")
all_prediction_plots

#ggsave("~/Library/CloudStorage/OneDrive-McGillUniversity/Urban-NOCA-MTL/05_figures/combined_all_prediction_plots.png",
#       plot = all_prediction_plots,
#       bg = "transparent",
#       width = 12, height = 12, dpi = 300)

