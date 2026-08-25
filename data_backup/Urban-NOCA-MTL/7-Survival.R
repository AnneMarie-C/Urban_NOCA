# Survival Analysis for Northern Cardinals (NOCA)
# Apparent and Fated Survival Analysis

# Load required packages
library(tidyverse)
library(lubridate)
library(RMark)  # For mark-recapture analysis
library(survival)  # For survival curves

# Read the data
df <- read.csv("~/Library/CloudStorage/OneDrive-McGillUniversity/Urban-NOCA-MTL/02_cleaned_data/noca_2425_cleaned.csv",
               stringsAsFactors = FALSE)

# Data preparation
df <- df %>%
  mutate(
    Date = as.Date(Date),
    Month = month(Date),
    Year = year(Date),
    Season = case_when(
      Month %in% c(10, 11, 12, 1, 2, 3) ~ "Winter",
      Month %in% c(4, 5) ~ "Spring",
      TRUE ~ "Other"
    ),
    Week = floor_date(Date, "week")
  )

# ============================================================
# 1. APPARENT SURVIVAL ANALYSIS (Cormack-Jolly-Seber Model)
# ============================================================

# Create encounter history for CJS models
# First, let's create a weekly encounter matrix

create_encounter_history <- function(data) {
  # Get unique weeks and individuals
  weeks <- sort(unique(data$Week))
  birds <- unique(data$ID)
  
  # Create encounter matrix
  encounter_mat <- matrix(0, nrow = length(birds), ncol = length(weeks))
  rownames(encounter_mat) <- birds
  colnames(encounter_mat) <- as.character(weeks)
  
  # Fill in encounters (1 = detected, 0 = not detected)
  for (i in 1:nrow(data)) {
    bird <- data$ID[i]
    week <- as.character(data$Week[i])
    encounter_mat[bird, week] <- 1
  }
  
  # Convert to character string format for RMark
  encounter_history <- apply(encounter_mat, 1, paste, collapse = "")
  
  return(data.frame(
    ch = encounter_history,
    ID = birds
  ))
}

enc_history <- create_encounter_history(df)

# Add covariates (sex, site, etc.)
bird_info <- df %>%
  group_by(ID) %>%
  summarise(
    Sex = first(Sex),
    Site = first(Site),
    First_Date = min(Date),
    Last_Date = max(Date)
  )

enc_data <- left_join(enc_history, bird_info, by = "ID")

# Run CJS model using RMark
# Model 1: Constant survival and detection
cjs_model_constant <- function() {
  # Process data
  processed <- process.data(enc_data, model = "CJS")
  
  # Create design data
  ddl <- make.design.data(processed)
  
  # Define models
  # Phi = survival probability
  # p = detection probability
  
  # Constant survival and detection
  Phi_dot <- list(formula = ~1)
  p_dot <- list(formula = ~1)
  
  # Time-varying survival
  Phi_time <- list(formula = ~time)
  
  # Sex effect on survival
  Phi_sex <- list(formula = ~Sex)
  
  # Run models
  model_list <- create.model.list("CJS")
  results <- mark.wrapper(model_list, data = processed, ddl = ddl)
  
  return(results)
}

# Model 2: Time-varying survival and constant detection
cjs_model_time <- function() {
  processed <- process.data(enc_data, model = "CJS")
  ddl <- make.design.data(processed)
  
  Phi_time <- list(formula = ~time)
  p_dot <- list(formula = ~1)
  
  model <- mark(processed, ddl, 
                model.parameters = list(Phi = Phi_time, p = p_dot))
  
  return(model)
}

# ============================================================
# 2. FATED SURVIVAL ANALYSIS (Known-fate Model)
# ============================================================

# For fated survival, we need to know which birds are confirmed dead
# vs. those that are just not detected

# Create known-fate encounter history
create_known_fate_history <- function(data) {
  bird_summary <- data %>%
    group_by(ID) %>%
    summarise(
      First_Date = min(Date),
      Last_Date = max(Date),
      N_Detections = n(),
      Sex = first(Sex),
      Site = first(Site),
      Study_End = max(data$Date)
    ) %>%
    mutate(
      # Days from first detection to last detection
      Days_Tracked = as.numeric(Last_Date - First_Date),
      # Fate: 1 = survived study, 0 = presumed dead (not seen for >30 days before study end)
      Days_Since_Last = as.numeric(Study_End - Last_Date),
      Fate = ifelse(Days_Since_Last < 30, 1, 0)
    )
  
  return(bird_summary)
}

fate_data <- create_known_fate_history(df)

# Kaplan-Meier survival curves
library(survminer)

# Create survival object
surv_obj <- Surv(time = fate_data$Days_Tracked, 
                 event = 1 - fate_data$Fate)  # event = 1 means death

# Overall survival curve
km_fit <- survfit(surv_obj ~ 1, data = fate_data)

# Plot Kaplan-Meier curve
plot(km_fit, 
     xlab = "Days Since First Detection",
     ylab = "Survival Probability",
     main = "Kaplan-Meier Survival Curve for Northern Cardinals")

# Survival by sex
km_sex <- survfit(surv_obj ~ Sex, data = fate_data)
plot(km_sex, 
     col = c("blue", "red"),
     xlab = "Days Since First Detection",
     ylab = "Survival Probability",
     main = "Survival by Sex")
legend("topright", c("Female", "Male"), col = c("blue", "red"), lty = 1)

# Log-rank test for difference between sexes
survdiff(surv_obj ~ Sex, data = fate_data)

# Cox proportional hazards model
cox_model <- coxph(surv_obj ~ Sex + Site, data = fate_data)
summary(cox_model)

# ============================================================
# 3. DETECTION PROBABILITY ANALYSIS
# ============================================================

# Calculate detection probability by period
detection_analysis <- df %>%
  group_by(Week, ID) %>%
  summarise(Detected = 1, .groups = "drop") %>%
  complete(Week, ID, fill = list(Detected = 0)) %>%
  group_by(Week) %>%
  summarise(
    N_Individuals = n_distinct(ID),
    N_Detected = sum(Detected),
    Detection_Rate = N_Detected / N_Individuals
  )

# Plot detection rates over time
ggplot(detection_analysis, aes(x = Week, y = Detection_Rate)) +
  geom_line() +
  geom_point() +
  labs(title = "Weekly Detection Rates",
       x = "Week",
       y = "Detection Probability") +
  theme_minimal()

# ============================================================
# 4. SUMMARY STATISTICS
# ============================================================

# Individual bird summary
bird_summary <- df %>%
  group_by(ID, Sex, Site) %>%
  summarise(
    First_Seen = min(Date),
    Last_Seen = max(Date),
    Days_Tracked = as.numeric(max(Date) - min(Date)),
    N_Detections = n(),
    .groups = "drop"
  )

print("Summary Statistics by Sex:")
bird_summary %>%
  group_by(Sex) %>%
  summarise(
    N_Birds = n(),
    Mean_Days_Tracked = mean(Days_Tracked),
    SD_Days_Tracked = sd(Days_Tracked),
    Mean_Detections = mean(N_Detections)
  )

# ============================================================
# 5. EXPORT RESULTS
# ============================================================

# Save encounter history for external analysis
write.csv(enc_data, "encounter_history.csv", row.names = FALSE)

# Save fate data
write.csv(fate_data, "fate_data.csv", row.names = FALSE)

# Save bird summary
write.csv(bird_summary, "bird_summary.csv", row.names = FALSE)

print("Analysis complete!")