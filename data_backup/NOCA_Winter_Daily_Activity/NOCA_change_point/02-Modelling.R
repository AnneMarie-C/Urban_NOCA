## June 20, 2025
## Anne-Marie Cousineau

# =============================================================================================== #
# 1. CREATING MODELS TO INVESTIGATE ACTIVITY AND TEMPERATURE LAG                                  #
# =============================================================================================== #

# Loading appropriate packages.
library(tidyverse)
library(dplyr)
library(lubridate)
library(patchwork)
library(tidyr)
library(mgcv)
library(gratia)
library(ggplot2)
library(AICcmodavg) # more info AICc


## --- 1.1. Preparing data for temperature lag (onset) --- ##
cp_onset <- read.csv("02_cleaned_data/cp_onset_cleaned.csv")

# Creating the lag columns for 5 days prior to the initial Min_temp recorded.
cp_onset <- cp_onset %>%
  arrange(Date) %>%  # Ensure the data is in chronological order
  mutate(
    Min_temp_lag1 = lag(Min_temp, 1),
    Min_temp_lag2 = lag(Min_temp, 2),
    Min_temp_lag3 = lag(Min_temp, 3),
    Min_temp_lag4 = lag(Min_temp, 4),
    Min_temp_lag5 = lag(Min_temp, 5)
  )

# Filtering out rows that don't have 5 days of lag (otherwise NA present).
cp_onset <- cp_onset %>%
  filter(!is.na(Min_temp_lag5))  # removes rows that don't have 5-day lag

# Adding columns that are the mean temperatures for the amount of days the lag is calculated.
# Will try to phrase this better...
cp_onset <- cp_onset %>%
  mutate(
    Mean_lag1 = rowMeans(across(c(Min_temp, Min_temp_lag1)), na.rm = TRUE),
    Mean_lag2 = rowMeans(across(c(Min_temp, Min_temp_lag1, Min_temp_lag2)), na.rm = TRUE),
    Mean_lag3 = rowMeans(across(c(Min_temp, Min_temp_lag1, Min_temp_lag2, Min_temp_lag3)), na.rm = TRUE),
    Mean_lag4 = rowMeans(across(c(Min_temp, Min_temp_lag1, Min_temp_lag2, Min_temp_lag3, Min_temp_lag4)), na.rm = TRUE),
    Mean_lag5 = rowMeans(across(c(Min_temp, Min_temp_lag1, Min_temp_lag2, Min_temp_lag3, Min_temp_lag4, Min_temp_lag5)), na.rm = TRUE)
  )

# We'll use the log of rain because the effect of going from 0 to 10 mm is probably stronger
# than going from 30 to 40 mm.
cp_onset <- cp_onset %>%
  mutate(log_precip = log(Total_precip + 1)) # Adding 1 so zeros can be included.

# Making sure ID, Winter and Site are all factors.
cp_onset$ID <- as.factor(cp_onset$ID)
cp_onset$Winter <- as.factor(cp_onset$Winter)
cp_onset$Site <- as.factor(cp_onset$Site)

## --- 1.2. Creating models for the onset of activity --- ##
# model.nolag <- gam(cpt.sr ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) +
#                      s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") +
#                      s(Site, bs = "re"), data = cp_onset)
# 
# model.lag1 <- gam(cpt.sr ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag1) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") +
#                     s(Site, bs = "re"), data = cp_onset)
# 
# model.lag2 <- gam(cpt.sr ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag2) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") +
#                     s(Site, bs = "re"), data = cp_onset)
# 
# model.lag3 <- gam(cpt.sr ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag3) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") +
#                     s(Site, bs = "re"), data = cp_onset)
# 
# model.lag4 <- gam(cpt.sr ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag4) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") + 
#                     s(Site, bs = "re"), data = cp_onset)
# 
# model.lag5 <- gam(cpt.sr ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag5) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") + 
#                     s(Site, bs = "re"), data = cp_onset)
# 
# model.null <- gam(cpt.sr ~ s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
#                   data = cp_onset)
# 
# # Concurvity test (done by Barbara).
# C0 <- concurvity(model.nolag, full = FALSE)
# C0 <- as.data.frame(C0$observed)
# C0$para <- NULL
# C0 = C0[-1,]
# 
# C1 <- concurvity(model.lag1, full = FALSE)
# C1 <- as.data.frame(C1$observed)
# C1$para <- NULL
# C1 = C1[-1,]
# 
# C2 <- concurvity(model.lag2, full = FALSE)
# C2 <- as.data.frame(C2$observed)
# C2$para <- NULL
# C2 = C2[-1,]
# 
# C3 <- concurvity(model.lag3, full = FALSE)
# C3 <- as.data.frame(C3$observed)
# C3$para <- NULL
# C3 = C3[-1,]
# 
# C4 <- concurvity(model.lag4, full = FALSE)
# C4 <- as.data.frame(C4$observed)
# C4$para <- NULL
# C4 = C4[-1,]
# 
# C5 <- concurvity(model.lag5, full = FALSE)
# C5 <- as.data.frame(C5$observed)
# C5$para <- NULL
# C5 = C5[-1,]
# 
# C6 <- concurvity(model.null, full = FALSE)
# C6 <- as.data.frame(C6$observed)
# C6$para <- NULL
# C6 = C6[-1,]

# After looking at concurvity, the temperature variables are highly correlated so we're only
# going to use Min_temp.

# New models:
model.nolag <- gam(cpt.sr ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_onset)

model.lag1 <- gam(cpt.sr ~ s(Min_temp)  + s(Mean_lag1) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_onset)

model.lag2 <- gam(cpt.sr ~ s(Min_temp) + s(Mean_lag2) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_onset)

model.lag3 <- gam(cpt.sr ~ s(Min_temp) + s(Mean_lag3) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_onset)

model.lag4 <- gam(cpt.sr ~ s(Min_temp) + s(Mean_lag4) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_onset)

model.lag5 <- gam(cpt.sr ~ s(Min_temp) + s(Mean_lag5) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_onset)

model.null <- gam(cpt.sr ~ s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"), data = cp_onset)

lag.AIC <- AIC(model.nolag, model.lag1, model.lag2, model.lag3, model.lag4, model.lag5, model.null)
lag.AIC
# model.nolag = 26753.43 AIC
# model.lag1 = 26749.10 AIC
# model.lag2 = 26741.93 AIC 
# model.lag3 = 26736.89 AIC
# model.lag4 = 26734.24 AIC
# model.lag5 = 26733.69 AIC <---
# model.null = 26867.95 AIC

model_list <- list(model.nolag, model.lag1, model.lag2, model.lag3, model.lag4, model.lag5, model.null)
model_names <- c("No Lag", "1-Day", "2-Day", "3-Day", "4-Day", "5-Day", "Null")

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))
aic_table

summary(model.lag5)
draw(model.lag5)

# Plotting AIC graph.
# Calculate AIC table for new models
aic_table <- data.frame(
  Model = model_names,
  AIC = sapply(model_list, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Substantial = Delta_AIC <= 2
  )

# Keep the original order for plotting
aic_table$Model <- factor(aic_table$Model, levels = model_names)

# Identify best model
min_aic <- min(aic_table$AIC)
best_model <- aic_table$Model[which.min(aic_table$AIC)]

# Plot
aic_lag_sr <- ggplot(aic_table, aes(x = Model, y = AIC)) +
  # Shaded substantial-support area
  annotate("rect",
           xmin = 0.5, xmax = length(model_names) + 0.5,
           ymin = min_aic, ymax = min_aic + 2,
           alpha = 0.2, fill = "grey") +
  # Black line and points
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  # Labels
  labs(
    title = "AIC Comparison of Lag Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_lag_sr


## --- 1.3. Preparing data for temperature lag (end) --- ##
cp_end <- read.csv("02_cleaned_data/cp_end_cleaned.csv")

# Creating the lag columns for 5 days prior to the initial Min_temp recorded.
cp_end <- cp_end %>%
  arrange(Date) %>%  # Ensure the data is in chronological order
  mutate(
    Min_temp_lag1 = lag(Min_temp, 1),
    Min_temp_lag2 = lag(Min_temp, 2),
    Min_temp_lag3 = lag(Min_temp, 3),
    Min_temp_lag4 = lag(Min_temp, 4),
    Min_temp_lag5 = lag(Min_temp, 5)
  )

# Filtering out rows that don't have 5 days of lag (otherwise NA present).
cp_end <- cp_end %>%
  filter(!is.na(Min_temp_lag5))  # removes rows that don't have 5-day lag

# Adding columns that are the mean temperatures for the amount of days the lag is calculated.
# Will try to phrase this better...
cp_end <- cp_end %>%
  mutate(
    Mean_lag1 = rowMeans(across(c(Min_temp, Min_temp_lag1)), na.rm = TRUE),
    Mean_lag2 = rowMeans(across(c(Min_temp, Min_temp_lag1, Min_temp_lag2)), na.rm = TRUE),
    Mean_lag3 = rowMeans(across(c(Min_temp, Min_temp_lag1, Min_temp_lag2, Min_temp_lag3)), na.rm = TRUE),
    Mean_lag4 = rowMeans(across(c(Min_temp, Min_temp_lag1, Min_temp_lag2, Min_temp_lag3, Min_temp_lag4)), na.rm = TRUE),
    Mean_lag5 = rowMeans(across(c(Min_temp, Min_temp_lag1, Min_temp_lag2, Min_temp_lag3, Min_temp_lag4, Min_temp_lag5)), na.rm = TRUE)
  )

# We'll use the log of rain because the effect of going from 0 to 10 mm is probably stronger
# than going from 30 to 40 mm.
cp_end <- cp_end %>%
  mutate(log_precip = log(Total_precip + 1)) # Adding 1 so zeros can be included.

# Making sure ID, Winter and Site are all factors.
cp_end$ID <- as.factor(cp_end$ID)
cp_end$Winter <- as.factor(cp_end$Winter)
cp_end$Site <- as.factor(cp_end$Site)

## --- 1.2. Creating models for the end of activity --- ##
# model.nolag <- gam(cpt.ss ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Speed_max_gust_wind) +
#                      s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") +
#                      s(Site, bs = "re"), data = cp_end)
# 
# model.lag1 <- gam(cpt.ss ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag1) + s(Speed_max_gust_wind) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") +
#                     s(Site, bs = "re"), data = cp_end)
# 
# model.lag2 <- gam(cpt.ss ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag2) + s(Speed_max_gust_wind) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") +
#                     s(Site, bs = "re"), data = cp_end)
# 
# model.lag3 <- gam(cpt.ss ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag3) + s(Speed_max_gust_wind) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") +
#                     s(Site, bs = "re"), data = cp_end)
# 
# model.lag4 <- gam(cpt.ss ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag4) + s(Speed_max_gust_wind) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") + 
#                     s(Site, bs = "re"), data = cp_end)
# 
# model.lag5 <- gam(cpt.ss ~ s(Min_temp) + s(Max_temp) + s(Mean_temp) + s(Mean_lag5) + s(Speed_max_gust_wind) +
#                     s(log_precip) + s(Snow_on_grnd) + s(Winter, bs = "re") + s(ID, bs = "re") + 
#                     s(Site, bs = "re"), data = cp_end)
# 
# model.null <- gam(cpt.ss ~ s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
#                   data = cp_end)
# 
# # Concurvity test (done by Barbara).
# C0 <- concurvity(model.nolag, full = FALSE)
# C0 <- as.data.frame(C0$observed)
# C0$para <- NULL
# C0 = C0[-1,]
# 
# C1 <- concurvity(model.lag1, full = FALSE)
# C1 <- as.data.frame(C1$observed)
# C1$para <- NULL
# C1 = C1[-1,]
# 
# C2 <- concurvity(model.lag2, full = FALSE)
# C2 <- as.data.frame(C2$observed)
# C2$para <- NULL
# C2 = C2[-1,]
# 
# C3 <- concurvity(model.lag3, full = FALSE)
# C3 <- as.data.frame(C3$observed)
# C3$para <- NULL
# C3 = C3[-1,]
# 
# C4 <- concurvity(model.lag4, full = FALSE)
# C4 <- as.data.frame(C4$observed)
# C4$para <- NULL
# C4 = C4[-1,]
# 
# C5 <- concurvity(model.lag5, full = FALSE)
# C5 <- as.data.frame(C5$observed)
# C5$para <- NULL
# C5 = C5[-1,]
# 
# C6 <- concurvity(model.null, full = FALSE)
# C6 <- as.data.frame(C6$observed)
# C6$para <- NULL
# C6 = C6[-1,]

# After looking at concurvity, the temperature variables are highly correlated so we're only
# going to use Min_temp.

# New models:
model.nolag <- gam(cpt.ss ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                   + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                   data = cp_end)

model.lag1 <- gam(cpt.ss ~ s(Min_temp)  + s(Mean_lag1) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_end)

model.lag2 <- gam(cpt.ss ~ s(Min_temp) + s(Mean_lag2) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_end)

model.lag3 <- gam(cpt.ss ~ s(Min_temp) + s(Mean_lag3) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_end)

model.lag4 <- gam(cpt.ss ~ s(Min_temp) + s(Mean_lag4) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_end)

model.lag5 <- gam(cpt.ss ~ s(Min_temp) + s(Mean_lag5) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind)
                  + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                  data = cp_end)

model.null <- gam(cpt.ss ~ s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"), data = cp_end)

lag.AIC <- AIC(model.nolag, model.lag1, model.lag2, model.lag3, model.lag4, model.lag5, model.null)
lag.AIC
# model.nolag = 24024.54 AIC
# model.lag1 = 24005.69 AIC
# model.lag2 = 23965.94 AIC <---
# model.lag3 = 24007.50 AIC 
# model.lag4 = 24018.26 AIC 
# model.lag5 = 24019.82 AIC
# model.null = 24090.85 AIC

model_list <- list(model.nolag, model.lag1, model.lag2, model.lag3, model.lag4, model.lag5, model.null)
model_names <- c("No Lag", "1-Day", "2-Day", "3-Day", "4-Day", "5-Day", "Null")

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)
aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))

# Calculate AIC table
aic_table <- data.frame(
  Model = model_names,
  AIC = sapply(model_list, AIC)
) %>%
  arrange(AIC) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Substantial = Delta_AIC <= 2
  )

# Preserve original order for plotting
aic_table$Model <- factor(aic_table$Model, levels = model_names)

# Identify best model
min_aic <- min(aic_table$AIC)
best_model <- aic_table$Model[which.min(aic_table$AIC)]

# Plot
aic_lag_ss <- ggplot(aic_table, aes(x = Model, y = AIC)) +
  # Shaded substantial-support area
  annotate("rect",
           xmin = 0.5, xmax = length(model_names) + 0.5,
           ymin = min_aic, ymax = min_aic + 2,
           alpha = 0.2, fill = "grey") +
  # Black line and points
  geom_line(group = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  # Labels
  labs(
    title = "AIC Comparison of Lag Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_lag_ss

summary(model.lag2)
draw(model.lag2)

# =============================================================================================== #
# 2. CREATING MODELS TO INVESTIGATE EFFECT OF COLD SNAPS                                          #
# =============================================================================================== #
# Identify and define the sudden temperature drops as "cold snap".

cp_onset <- cp_onset %>%
  mutate(
    temp_drop = Min_temp - Min_temp_lag1,
    cold_snap_5 = if_else(temp_drop <= -5, 1, 0),
    cold_snap_7 = if_else(temp_drop <= -7, 1, 0),
    cold_snap_10 = if_else(temp_drop <= -10, 1, 0)
  )
cp_onset <- cp_onset %>%
  mutate(
    temp_drop = Min_temp - Min_temp_lag1,
    heat_snap_5 = if_else(temp_drop >= 5, 1, 0),
    heat_snap_7 = if_else(temp_drop >= 7, 1, 0),
    heat_snap_10 = if_else(temp_drop >= 10, 1, 0)
  )

cp_end <- cp_end %>%
  mutate(
    temp_drop = Min_temp - Min_temp_lag1,
    cold_snap_5 = if_else(temp_drop <= -5, 1, 0),
    cold_snap_7 = if_else(temp_drop <= -7, 1, 0),
    cold_snap_10 = if_else(temp_drop <= -10, 1, 0)
  )
cp_end <- cp_end %>%
  mutate(
    temp_drop = Min_temp - Min_temp_lag1,
    heat_snap_5 = if_else(temp_drop >= 5, 1, 0),
    heat_snap_7 = if_else(temp_drop >= 7, 1, 0),
    heat_snap_10 = if_else(temp_drop >= 10, 1, 0)
  )

# Quickly plot the patterns.

# Step 1: Pivot to long format for cold snap thresholds.

cp_onset_long <- cp_onset %>%
  pivot_longer(
    cols = starts_with("cold_snap_"),
    names_to = "snap_type",
    values_to = "cold_snap"
  ) %>%
  mutate(
    snap_type = recode(snap_type,
                       "cold_snap_5" = "Cold Snap ≥5°C",
                       "cold_snap_7" = "Cold Snap ≥7°C",
                       "cold_snap_10" = "Cold Snap ≥10°C")
  )

cp_end_long <- cp_end %>%
  pivot_longer(
    cols = starts_with("cold_snap_"),
    names_to = "snap_type",
    values_to = "cold_snap"
  ) %>%
  mutate(
    snap_type = recode(snap_type,
                       "cold_snap_5" = "Cold Snap ≥5°C",
                       "cold_snap_7" = "Cold Snap ≥7°C",
                       "cold_snap_10" = "Cold Snap ≥10°C")
  )

# Build GAM model to test for the effects of temperature drop.

# Onset.

# model.drop.onset <- gam(cpt.sr ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
#                          s(Site, bs = "re") + s(ID, bs = "re"), data = cp_onset, method = "REML")
# model.drop.onset.null <- gam(cpt.sr ~ s(Site, bs = "re") + s(ID, bs = "re"), data = cp_onset, method = "REML")

# C0 <- concurvity(model.drop.onset, full = FALSE)
# C0 <- as.data.frame(C0$observed)
# C0$para <- NULL
# C0 = C0[-1,]

# C1 <- concurvity(model.drop.onset.null, full = FALSE)
# C1 <- as.data.frame(C1$observed)
# C1$para <- NULL
# C1 = C1[-1,]

# No issues with concurvity, even between temp_drop and Min_temp.

# AIC.table.drop <- AIC(model.drop.onset, model.drop.onset.null)
# AIC.table.drop
# model.drop.onset = 26789.96 AIC <---
# model.drop.onset.null = 26873.30 AIC
# summary(model.drop.onset)
# draw(model.drop.onset)

# End.

# model.drop.end <- gam(cpt.ss ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
#                          s(Site, bs = "re") + s(ID, bs = "re"), data = cp_end, method = "REML")
# model.drop.end.null <- gam(cpt.ss ~ s(Site, bs = "re") + s(ID, bs = "re"), data = cp_end, method = "REML")

# C0 <- concurvity(model.drop.end, full = FALSE)
# C0 <- as.data.frame(C0$observed)
# C0$para <- NULL
# C0 = C0[-1,]

# C1 <- concurvity(model.drop.end.null, full = FALSE)
# C1 <- as.data.frame(C1$observed)
# C1$para <- NULL
# C1 = C1[-1,]

# Concurvity still good.

# AIC.table.drop <- AIC(model.drop.end, model.drop.end.null)
# AIC.table.drop
# model.drop.end = 24095.15 AIC <---
# model.drop.end.null = 24106.95 AIC
# summary(model.drop.end)
# draw(model.drop.end)

# Build GAM model to investigate the effect of the cold snaps.
# First making sure the cold snaps are factors.
cp_onset <- cp_onset %>%
  mutate(
    cold_snap_5 = as.factor(cold_snap_5),
    cold_snap_7 = as.factor(cold_snap_7),
    cold_snap_10 = as.factor(cold_snap_10)
  )
cp_onset <- cp_onset %>%
  mutate(
    heat_snap_5 = as.factor(heat_snap_5),
    heat_snap_7 = as.factor(heat_snap_7),
    heat_snap_10 = as.factor(heat_snap_10)
  )

cp_end <- cp_end %>%
  mutate(
    cold_snap_5 = as.factor(cold_snap_5),
    cold_snap_7 = as.factor(cold_snap_7),
    cold_snap_10 = as.factor(cold_snap_10)
  )
cp_end <- cp_end %>%
  mutate(
    heat_snap_5 = as.factor(heat_snap_5),
    heat_snap_7 = as.factor(heat_snap_7),
    heat_snap_10 = as.factor(heat_snap_10)
  )

# In the models, we're not smoothing the cold snaps because they're binary variables.
model.nosnap.onset <- gam(cpt.sr ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
                           s(Site, bs = "re") + s(ID, bs = "re"),
                         data = cp_onset, method = "REML")

model.cold5.onset <- gam(cpt.sr ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
                         cold_snap_5 +
                         s(Site, bs = "re") + s(ID, bs = "re"),
                       data = cp_onset, method = "REML")

model.cold7.onset <- gam(cpt.sr ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
                         cold_snap_7 +
                         s(Site, bs = "re") + s(ID, bs = "re"),
                       data = cp_onset, method = "REML")

model.cold10.onset <- gam(cpt.sr ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
                          cold_snap_10 +
                          s(Site, bs = "re") + s(ID, bs = "re"),
                        data = cp_onset, method = "REML")

model.null.onset <- gam(cpt.sr ~ s(Site, bs = "re") + s(ID, bs = "re"),
                          data = cp_onset, method = "REML")

# C1 <- concurvity(model.cold5.onset, full = FALSE)
# C1 <- as.data.frame(C1$observed)
# C1$para <- NULL
# C1 = C1[-1,]
# 
# C2 <- concurvity(model.cold7.onset, full = FALSE)
# C2 <- as.data.frame(C2$observed)
# C2$para <- NULL
# C2 = C2[-1,]
# 
# C3 <- concurvity(model.cold10.onset, full = FALSE)
# C3 <- as.data.frame(C3$observed)
# C3$para <- NULL
# C3 = C3[-1,]

AIC(model.nosnap.onset, model.cold5.onset, model.cold7.onset, model.cold10.onset, model.null.onset)
# model.drop.onset = 26768.80 AIC
# model.cold5.onset = 26761.39 AIC <---
# model.cold7.onset = 26766.92 AIC 
# model.cold10.onset = 26770.65 AIC
# model.drop.onset.null = 26873.30 AIC

model_list <- list(model.nosnap.onset, model.cold5.onset, model.cold7.onset, model.cold10.onset, model.null.onset)
model_names <- c("No Snap", "-5˚C", "-7˚C", "-10˚C", "Null")

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
aic_cs_sr <- ggplot(aic_table, aes(x = Model, y = AIC)) +
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
    title = "AIC Comparison of Cold Snap Onset Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_cs_sr

summary(model.cold5.onset)
draw(model.cold5.onset)

model.nosnap.end <- gam(cpt.ss ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
                            s(Site, bs = "re") + s(ID, bs = "re"),
                          data = cp_end, method = "REML")

model.cold5.end <- gam(cpt.ss ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
                     cold_snap_5 +
                     s(Site, bs = "re") + s(ID, bs = "re"),
                   data = cp_end, method = "REML")

model.cold7.end <- gam(cpt.ss ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
                     cold_snap_7 +
                     s(Site, bs = "re") + s(ID, bs = "re"),
                   data = cp_end, method = "REML")

model.cold10.end <- gam(cpt.ss ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Speed_max_gust_wind) +
                      cold_snap_10 +
                      s(Site, bs = "re") + s(ID, bs = "re"),
                    data = cp_end, method = "REML")

model.null.end <- gam(cpt.ss ~ s(Site, bs = "re") + s(ID, bs = "re"),
                       data = cp_end, method = "REML")

# C1 <- concurvity(model.cold5.end, full = FALSE)
# C1 <- as.data.frame(C1$observed)
# C1$para <- NULL
# C1 = C1[-1,]
# 
# C2 <- concurvity(model.cold7.end, full = FALSE)
# C2 <- as.data.frame(C2$observed)
# C2$para <- NULL
# C2 = C2[-1,]
# 
# C3 <- concurvity(model.cold10.end, full = FALSE)
# C3 <- as.data.frame(C3$observed)
# C3$para <- NULL
# C3 = C3[-1,]

AIC(model.nosnap.end, model.cold5.end, model.cold7.end, model.cold10.end, model.null.end)
# model.drop.end = 24095.99 AIC
# model.cold5.end = 24088.47 AIC
# model.cold7.end = 24097.29 AIC 
# model.cold10.end = 24085.50 AIC <---
# model.drop.end.null = 24106.95 AIC

model_list <- list(model.nosnap.end, model.cold5.end, model.cold7.end, model.cold10.end, model.null.end)
model_names <- c("No Snap", "-5˚C", "-7˚C", "-10˚C", "Null")

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
aic_cs_ss <- ggplot(aic_table, aes(x = Model, y = AIC)) +
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
    title = "AIC Comparison of Cold Snap End Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_cs_ss

summary(model.cold10.end)
draw(model.cold10.end)

# Final plot to select best model according to AIC.
(aic_lag_sr|aic_cs_sr) /
  (aic_lag_ss|aic_cs_ss)

# Building global model for both Onset and End using cold snaps and temp lag.

# Onset
full.model.onset <- gam(cpt.sr ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Mean_lag5) + cold_snap_5
                   + s(Speed_max_gust_wind) + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                   data = cp_onset, method = "REML")
precip.model.onset <- gam(cpt.sr ~ s(log_precip) + s(Snow_on_grnd) + s(Mean_lag5) + cold_snap_5
                          + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                          data = cp_onset, method = "REML")
temp.model.onset <- gam(cpt.sr ~ s(Min_temp) + s(Mean_lag5) + cold_snap_5
                        + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                        data = cp_onset, method = "REML")
variation.model.onset <- gam(cpt.sr ~ s(Mean_lag5) + cold_snap_5
                             + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                             data = cp_onset, method = "REML")
null.model.onset <- gam(cpt.sr ~ s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"), 
                        data = cp_onset)
AIC(full.model.onset, precip.model.onset, temp.model.onset, variation.model.onset, null.model.onset)
# full.model.onset = 26749.97 AIC <---
# precip.model.onset = 26758.09 AIC
# temp.model.onset = 26762.13 AIC 
# variation.model.onset = 26784.44 AIC
# null.model.onset = 26867.95 AIC

model_list <- list(full.model.onset, precip.model.onset, temp.model.onset, variation.model.onset, null.model.onset)
model_names <- c("Full", "Precipitation", "Temperature", "Extreme", "Null")

aic_values <- sapply(model_list, AIC)

aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

aic_table <- aic_table %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = AIC - min(AIC))
aic_table

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
aic_total_model_sr <- ggplot(aic_table, aes(x = Model, y = AIC)) +
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
    title = "AIC Comparison of Total Onset Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_total_model_sr

summary(full.model.onset)
draw(full.model.onset)

# We can remove winter and site since they're not significant, along with log_precip.
full.model.onset <- gam(cpt.sr ~ s(Min_temp) + s(Snow_on_grnd) + s(Mean_lag5) + cold_snap_5 + s(Speed_max_gust_wind)
                        + s(ID, bs = "re"),
                        data = cp_onset, method = "REML")
summary(full.model.onset)
draw(full.model.onset)

# End
full.model.end <- gam(cpt.ss ~ s(Min_temp) + s(log_precip) + s(Snow_on_grnd) + s(Mean_lag2) + cold_snap_10
                        + s(Speed_max_gust_wind) + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                        data = cp_end, method = "REML")
precip.model.end <- gam(cpt.ss ~ s(log_precip) + s(Snow_on_grnd) + s(Mean_lag2) + cold_snap_10
                          + s(Speed_max_gust_wind) + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                          data = cp_end, method = "REML")
temp.model.end <- gam(cpt.ss ~ s(Min_temp) + s(Mean_lag2) + cold_snap_10
                        + s(Speed_max_gust_wind) + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                        data = cp_end, method = "REML")
variation.model.end <- gam(cpt.ss ~ s(Mean_lag2) + cold_snap_10
                             + s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"),
                             data = cp_end, method = "REML")
null.model.end <- gam(cpt.ss ~ s(Winter, bs = "re") + s(ID, bs = "re") + s(Site, bs = "re"), 
                        data = cp_end)
AIC(full.model.end, precip.model.end, temp.model.end, variation.model.end, null.model.end)
# full.model.end = 24035.18 AIC 
# precip.model.end = 24074.17 AIC
# temp.model.end = 24032.54 AIC <---
# variation.model.end = 24077.56 AIC
# null.model.end = 24090.85 AIC

model_list <- list(full.model.end, precip.model.end, temp.model.end, variation.model.end, null.model.end)
model_names <- c("Full", "Precipitation", "Temperature", "Extreme", "Null")

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
aic_total_model_ss <- ggplot(aic_table, aes(x = Model, y = AIC)) +
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
    title = "AIC Comparison of Total End Models",
    y = "AIC",
    x = "Model"
  ) +
  theme_minimal(base_size = 14)
aic_total_model_ss

summary(temp.model.end)
draw(temp.model.end)

# Plotting the two graphs together.
(aic_total_model_sr|aic_total_model_ss)

# We can remove winter and site since they're not significant, along with wind.
temp.model.end <- gam(cpt.ss ~ s(Min_temp) + s(Mean_lag2) + cold_snap_10
                      + s(ID, bs = "re"),
                      data = cp_end, method = "REML")
summary(temp.model.end)
draw(temp.model.end)

# # Just looking at what we get when we consider ID as a variable.
# full.model.onset.id <- gam(cpt.sr ~ s(Min_temp) + s(Snow_on_grnd) + s(Mean_lag5) + cold_snap_5 + s(Speed_max_gust_wind)
#                         + (ID),
#                         data = cp_onset, method = "REML")
# summary(full.model.onset.id)
# draw(full.model.onset.id)
# 
# temp.model.end.id <- gam(cpt.ss ~ s(Min_temp) + s(Mean_lag2) + cold_snap_10
#                       + (ID),
#                       data = cp_end, method = "REML")
# summary(temp.model.end.id)
# draw(temp.model.end.id)

# FIGURES

#-------------------#
# Helper Function to Create Smooth Plots
#-------------------#
make_smooth_plot <- function(model, term, x_label, y_label) {
  draw(model, select = term) +
    theme_minimal(base_size = 13) +
    labs(
      title = NULL,
      x = x_label,
      y = y_label
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
}

#-------------------#
# Helper Function to Plot Cold Snap Based on Model Predictions
#-------------------#
make_cold_snap_boxplot_from_model <- function(model, data, x_var, x_label, y_label) {
  # Predict fitted values from the model
  data$fitted <- predict(model, newdata = data, type = "response")
  
  # Handle case where x_var is not in data
  if (!x_var %in% names(data)) {
    stop(paste("Column", x_var, "not found in data"))
  }
  
  # Get p-value from parametric table
  cold_snap_term <- grep("cold_snap", rownames(summary(model)$p.table), value = TRUE)
  p_value <- summary(model)$p.table[cold_snap_term, "Pr(>|t|)"]
  p_label <- paste0("p = ", format.pval(p_value, digits = 3, eps = .001))
  
  # Compute y max for annotation
  y_max <- max(data$fitted, na.rm = TRUE)
  y_label_pos <- y_max + (0.05 * y_max)  # 5% above max
  
  # Create plot
  ggplot(data, aes(x = factor(.data[[x_var]]), y = fitted)) +
    geom_boxplot(fill = "steelblue", alpha = 0.6) +
    annotate("text", x = 1.5, y = y_label_pos, label = p_label,
             size = 5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
}
#===========================#
# FIGURE 1: ONSET TIME
#===========================#

# Smooth terms & labels
onset_terms <- c("s(Min_temp)", "s(Snow_on_grnd)", "s(Mean_lag5)", "s(Speed_max_gust_wind)")
onset_labels <- c(
  "s(Min_temp)" = "Minimum temperature (°C)",
  "s(Snow_on_grnd)" = "Snow accumulation on the ground (cm)",
  "s(Mean_lag5)" = "5-Day mean temperature lag (°C)",
  "s(Speed_max_gust_wind)" = "Maximum wind gust speed (km/h)"
)

# Make smooth plots
onset_smooths <- lapply(onset_terms, function(term) {
  make_smooth_plot(full.model.onset, term, onset_labels[[term]], "Partial effect on onset time")
})

# Cold snap boxplot based on GAM
onset_box <- make_cold_snap_boxplot_from_model(
  model = full.model.onset,
  data = cp_onset,
  x_var = "cold_snap_5",
  x_label = "Cold Snap -5 °C (0 = No, 1 = Yes)",
  y_label = "Predicted onset time"
)

# Combine onset plots
onset_fig <- wrap_plots(c(onset_smooths, list(onset_box)), ncol = 3) +
  plot_annotation(
    title = "Effects of environmental variables on the northern cardinal activity onset time relative to sunrise",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

#===========================#
# FIGURE 2: END TIME
#===========================#

# Smooth terms & labels
end_terms <- c("s(Min_temp)", "s(Mean_lag2)")
end_labels <- c(
  "s(Min_temp)" = "Minimum Temperature (°C)",
  "s(Mean_lag2)" = "2-Day mean temperature lag (°C)"
)

# Make smooth plots
end_smooths <- lapply(end_terms, function(term) {
  make_smooth_plot(temp.model.end, term, end_labels[[term]], "Partial effect on end time")
})

# Cold snap boxplot (with significance label now)
end_box <- make_cold_snap_boxplot_from_model(
  model = temp.model.end,
  data = cp_onset,
  x_var = "cold_snap_10",
  x_label = "Cold Snap -10 °C (0 = No, 1 = Yes)",
  y_label = "Predicted end time"
)

# Combine end plots
end_fig <- wrap_plots(c(end_smooths, list(end_box)), ncol = 3) +
  plot_annotation(
    title = "Effects of environmental variables on the northern cardinal activity end time relative to sunset",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

#===========================#
# Display or Save
#===========================#
onset_fig
end_fig

ggsave("onset_plot.png", plot = onset_fig, width = 12, height = 9, dpi = 300, bg = "transparent")
ggsave("end_plot.png", plot = end_fig, width = 12, height = 6, dpi = 300, bg = "transparent")