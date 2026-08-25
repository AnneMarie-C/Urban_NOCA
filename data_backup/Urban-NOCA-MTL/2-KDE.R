# =============================================================================================== #
# CALCULATING KERNEL DENSITY ESTIMATOR FOR NORTHERN CARDINALS                                     #
# =============================================================================================== #

# Load the appropriate packages.
library(sp)
# library(tidyverse)
library(ggplot2)
# library(ggpubr)
library(dplyr)
library(sf)
library(adehabitatHR)
library(units)
# library(car) -> doesn't work anymore, replaced by rstatix.
library(rstatix)
library(rcompanion)
library(patchwork)
library(ggsignif)
library(purrr)
library(multcomp)
library(multcompView)


# =============================================================================================== #
# 1. KDE FOR DATA COLLECTED IN 2022 & 2023                                                        #
# =============================================================================================== #

## --- 1.1. Load & prepare the data --- ##

noca_2223_kde <- read.csv("02_cleaned_data/noca_2223_cleaned.csv")

# Add a site column so it can run future functions.
noca_2223_kde <- noca_2223_kde %>% 
  mutate(Site = 1)

# Convert the coordinates to an sf object R can use.
noca_2223_sf <- st_as_sf(noca_2223_kde, coords = c("Lon", "Lat"), crs = 4326)
# Convert sf object back to spatial points for KDE
noca_2223_sp <- as(noca_2223_sf, "Spatial")

# Split into a list by NOCA_ID
noca_2223_list <- split(noca_2223_sf, noca_2223_sf$ID)


## --- 1.2. Applying kernelUD to the observations with enough data --- ##
kde_results_noca_2223 <- lapply(noca_2223_list, function(bird_data) {
  # Convert the filtered data (per individual) to a SpatialPointsDataFrame
  noca_2223_sp <- as(bird_data, "Spatial")  # Convert only the filtered individual
  
  # Run kernelUD with a larger extent (increase extent if necessary)
  kde <- kernelUD(noca_2223_sp, h = "href", extent =  2)  # Set extent for grid size
  return(kde)
})


## --- 1.3. Visualizing the KDE before saving --- ##

# Extract 50% and 95% UDs for each bird
noca_2223_kde_50 <- lapply(kde_results_noca_2223, getverticeshr, percent = 50)
noca_2223_kde_95 <- lapply(kde_results_noca_2223, getverticeshr, percent = 95)

# Convert them to sf objects for plotting & add bird IDs to KDE 50% and 95% UDs
noca_2223_kde_50_sf <- do.call(rbind, lapply(names(noca_2223_kde_50), function(ID) {
  sf_obj <- st_as_sf(noca_2223_kde_50[[ID]])
  sf_obj$ID <- ID  # Add bird ID as a column
  return(sf_obj)
}))

noca_2223_kde_95_sf <- do.call(rbind, lapply(names(noca_2223_kde_95), function(ID) {
  sf_obj <- st_as_sf(noca_2223_kde_95[[ID]])
  sf_obj$ID <- ID  # Add bird ID as a column
  return(sf_obj)
}))

noca_2223_sf$ID <- as.factor(noca_2223_sf$ID)
noca_2223_kde_50_sf$ID <- as.factor(noca_2223_kde_50_sf$ID)
noca_2223_kde_95_sf$ID <- as.factor(noca_2223_kde_95_sf$ID)

## --- 1.4. Plotting the results --- ##
ggplot() +
  geom_sf(data = noca_2223_kde_95_sf, aes(fill = ID), alpha = 0.1, color = NA) +  # 95% UD with transparent fill
  geom_sf(data = noca_2223_kde_50_sf, aes(fill = ID), alpha = 0.4, color = NA) +  # 50% UD with more opaque fill
  ggtitle("Kernel Density Estimators (KDEs) for every Northern Cardinal in 2022-2023") +
  theme_minimal() + # GPS points with ID-based colors
  geom_sf(data = noca_2223_sf, aes(color = ID), size = 0.8) +
  theme(legend.title = element_blank())  # Optional: Remove legend title

## --- 1.5. Export the resulting KDE for the 50% and 95% as a shapefile --- ##

folder_path <- "03_processed_data"

noca_2223_kde_50_sf_clean <- noca_2223_kde_50_sf[, c("area", "ID", "geometry")]
st_write(noca_2223_kde_50_sf_clean, paste0(folder_path, "/noca_2223_kde_50.shp"), delete_layer = TRUE)

noca_2223_kde_95_sf_clean <- noca_2223_kde_95_sf[, c("area", "ID", "geometry")]
st_write(noca_2223_kde_95_sf_clean, paste0(folder_path, "/noca_2223_kde_95.shp"), delete_layer = TRUE)


# =============================================================================================== #
# 2. KDE FOR DATA COLLECTED IN 2023 & 2024                                                        #
# =============================================================================================== #

## --- 2.1. Load & prepare the data --- ##

noca_2324_kde <- read.csv("02_cleaned_data/noca_2324_cleaned.csv")

# Remove data points in BDU for birds that switched to MBO so that the KDE only calcunew_lates
# the home range once they are at MBO.
# NOCA 93b.
bird_id <- "93b"
dates_to_remove <- as.Date(c("2023-11-02"))
noca_2324_kde <- noca_2324_kde %>%
  filter(!(ID == bird_id & Date %in% dates_to_remove))

# NOCA 130.
bird_id <- "130"
dates_to_remove <- as.Date(c("2023-10-23", "2023-10-25", "2023-10-27", "2023-11-06", "2023-11-15"))
noca_2324_kde <- noca_2324_kde %>%
  filter(!(ID == bird_id & Date %in% dates_to_remove))

# NOCA 121.
bird_id <- "121"
dates_to_remove <- as.Date(c("2023-10-18", "2023-10-20", "2023-10-25"))
noca_2324_kde <- noca_2324_kde %>%
  filter(!(ID == bird_id & Date %in% dates_to_remove))

# Convert the coordinates to an sf object R can use.
noca_2324_sf <- st_as_sf(noca_2324_kde, coords = c("Lon", "Lat"), crs = 4326)
# Convert sf object back to spatial points for KDE
noca_2324_sp <- as(noca_2324_sf, "Spatial")

# Split into a list by NOCA_ID
noca_2324_list <- split(noca_2324_sf, noca_2324_sf$ID)


## --- 2.2. Applying kernelUD to the observations with enough data --- ##
kde_results_noca_2324 <- lapply(noca_2324_list, function(bird_data) {
  # Convert the filtered data (per individual) to a SpatialPointsDataFrame
  noca_2324_sp <- as(bird_data, "Spatial")  # Convert only the filtered individual
  
  # Run kernelUD with a larger extent (increase extent if necessary)
  kde <- kernelUD(noca_2324_sp, h = "href", extent =  2)  # Set extent for grid size
  return(kde)
})


## --- 2.3. Visualizing the KDE before saving --- ##

# Extract 50% and 95% UDs for each bird
noca_2324_kde_50 <- lapply(kde_results_noca_2324, getverticeshr, percent = 50)
noca_2324_kde_95 <- lapply(kde_results_noca_2324, getverticeshr, percent = 95)

# Convert them to sf objects for plotting & add bird IDs to KDE 50% and 95% UDs
noca_2324_kde_50_sf <- do.call(rbind, lapply(names(noca_2324_kde_50), function(ID) {
  sf_obj <- st_as_sf(noca_2324_kde_50[[ID]])
  sf_obj$ID <- ID  # Add bird ID as a column
  return(sf_obj)
}))

noca_2324_kde_95_sf <- do.call(rbind, lapply(names(noca_2324_kde_95), function(ID) {
  sf_obj <- st_as_sf(noca_2324_kde_95[[ID]])
  sf_obj$ID <- ID  # Add bird ID as a column
  return(sf_obj)
}))

noca_2324_sf$ID <- as.factor(noca_2324_sf$ID)
noca_2324_kde_50_sf$ID <- as.factor(noca_2324_kde_50_sf$ID)
noca_2324_kde_95_sf$ID <- as.factor(noca_2324_kde_95_sf$ID)

## --- 2.4. Plotting the results --- ##
ggplot() +
  geom_sf(data = noca_2324_kde_95_sf, aes(fill = ID), alpha = 0.1, color = NA) +  # 95% UD with transparent fill
  geom_sf(data = noca_2324_kde_50_sf, aes(fill = ID), alpha = 0.4, color = NA) +  # 50% UD with more opaque fill
  ggtitle("Kernel Density Estimators (KDEs) for every Northern Cardinal in 2023-2024") +
  theme_minimal() + # GPS points with ID-based colors
  geom_sf(data = noca_2324_sf, aes(color = ID), size = 0.8) +
  theme(legend.title = element_blank())  # Optional: Remove legend title

## --- 2.5. Export the resulting KDE for the 50% and 95% as a shapefile --- ##

folder_path <- "03_processed_data"

noca_2324_kde_50_sf_clean <- noca_2324_kde_50_sf[, c("area", "ID", "geometry")]
st_write(noca_2324_kde_50_sf_clean, paste0(folder_path, "/noca_2324_kde_50.shp"), delete_layer = TRUE)

noca_2324_kde_95_sf_clean <- noca_2324_kde_95_sf[, c("area", "ID", "geometry")]
st_write(noca_2324_kde_95_sf_clean, paste0(folder_path, "/noca_2324_kde_95.shp"), delete_layer = TRUE)


# =============================================================================================== #
# 3. KDE FOR DATA COLLECTED IN 2024 & 2025                                                        #
# =============================================================================================== #

## --- 3.1. Load & prepare the data --- ##

noca_2425_kde <- read.csv("02_cleaned_data/noca_2425_cleaned.csv")

# Remove data points in BDU for birds that switched to MBO so that the KDE only calcunew_lates
# the home range once they are at MBO.
# NOCA 627.
bird_id <- "627"
dates_to_remove <- as.Date(c("2024-11-08", "2024-11-10"))
noca_2425_kde <- noca_2425_kde %>%
  filter(!(ID == bird_id & Date %in% dates_to_remove))

# Convert the coordinates to an sf object R can use.
noca_2425_sf <- st_as_sf(noca_2425_kde, coords = c("Lon", "Lat"), crs = 4326)
# Convert sf object back to spatial points for KDE
noca_2425_sp <- as(noca_2425_sf, "Spatial")

# Split into a list by NOCA_ID
noca_2425_list <- split(noca_2425_sf, noca_2425_sf$ID)


## --- 3.2. Applying kernelUD to the observations with enough data --- ##
kde_results_noca_2425 <- lapply(noca_2425_list, function(bird_data) {
  # Convert the filtered data (per individual) to a SpatialPointsDataFrame
  noca_2425_sp <- as(bird_data, "Spatial")  # Convert only the filtered individual
  
  # Run kernelUD with a larger extent (increase extent if necessary)
  kde <- kernelUD(noca_2425_sp, h = "href", extent =  2.5)  # Set extent for grid size
  return(kde)
})


## --- 3.3. Visualizing the KDE before saving --- ##

# Extract 50% and 95% UDs for each bird
noca_2425_kde_50 <- lapply(kde_results_noca_2425, getverticeshr, percent = 50)
noca_2425_kde_95 <- lapply(kde_results_noca_2425, getverticeshr, percent = 95)

# Convert them to sf objects for plotting & add bird IDs to KDE 50% and 95% UDs
noca_2425_kde_50_sf <- do.call(rbind, lapply(names(noca_2425_kde_50), function(ID) {
  sf_obj <- st_as_sf(noca_2425_kde_50[[ID]])
  sf_obj$ID <- ID  # Add bird ID as a column
  return(sf_obj)
}))

noca_2425_kde_95_sf <- do.call(rbind, lapply(names(noca_2425_kde_95), function(ID) {
  sf_obj <- st_as_sf(noca_2425_kde_95[[ID]])
  sf_obj$ID <- ID  # Add bird ID as a column
  return(sf_obj)
}))

noca_2425_sf$ID <- as.factor(noca_2425_sf$ID)
noca_2425_kde_50_sf$ID <- as.factor(noca_2425_kde_50_sf$ID)
noca_2425_kde_95_sf$ID <- as.factor(noca_2425_kde_95_sf$ID)

## --- 3.4. Plotting the results --- ##
ggplot() +
  geom_sf(data = noca_2425_kde_95_sf, aes(fill = ID), alpha = 0.1, color = NA) +  # 95% UD with transparent fill
  geom_sf(data = noca_2425_kde_50_sf, aes(fill = ID), alpha = 0.4, color = NA) +  # 50% UD with more opaque fill
  ggtitle("Kernel Density Estimators (KDEs) for every Northern Cardinal in 2024-2025") +
  theme_minimal() + # GPS points with ID-based colors
  geom_sf(data = noca_2425_sf, aes(color = ID), size = 0.8) +
  theme(legend.title = element_blank())  # Optional: Remove legend title

## --- 3.5. Export the resulting KDE for the 50% and 95% as a shapefile --- ##

folder_path <- "03_processed_data"

noca_2425_kde_50_sf_clean <- noca_2425_kde_50_sf[, c("area", "ID", "geometry")]
st_write(noca_2425_kde_50_sf_clean, paste0(folder_path, "/noca_2425_kde_50.shp"), delete_layer = TRUE)

noca_2425_kde_95_sf_clean <- noca_2425_kde_95_sf[, c("area", "ID", "geometry")]
st_write(noca_2425_kde_95_sf_clean, paste0(folder_path, "/noca_2425_kde_95.shp"), delete_layer = TRUE)


# =============================================================================================== #
# 4. PLOTTING KDE AREAS FOR 2024-2025                                                                       #
# =============================================================================================== #

## --- 4.1. Function to process the kde data and extract the area. --- ##
process_kde_data <- function(kde_50, kde_95, noca_data, year) {
  # Validate geometries
  kde_50 <- if (any(!st_is_valid(kde_50))) st_make_valid(kde_50) else kde_50
  kde_95 <- if (any(!st_is_valid(kde_95))) st_make_valid(kde_95) else kde_95
  
  # Compute areas
  kde_50$Area_50 <- st_area(kde_50)
  kde_95$Area_95 <- st_area(kde_95)
  
  # Add site information
  kde_50 <- kde_50 %>%
    left_join(noca_data %>% dplyr::select(ID, Site) %>% distinct(), by = "ID")
  kde_95 <- kde_95 %>%
    left_join(noca_data %>% dplyr::select(ID, Site) %>% distinct(), by = "ID")
  
  # Combine KDE areas and add year
  kde_data <- bind_rows(
    kde_50 %>%
      dplyr::mutate(KDE = "50", Area = as.numeric(Area_50)) %>%
      dplyr::select(ID, Site, KDE, Area),
    kde_95 %>%
      dplyr::mutate(KDE = "95", Area = as.numeric(Area_95)) %>%
      dplyr::select(ID, Site, KDE, Area)
  ) %>%
    mutate(Year = year)  # Add year column
  
  return(kde_data)
}

## --- 4.2. Prepare data for plotting --- ##

# Make sure the ID column is a factor.
noca_2223_kde_50_sf$ID <- as.factor(noca_2223_kde_50_sf$ID)
noca_2223_kde_95_sf$ID <- as.factor(noca_2223_kde_95_sf$ID)
noca_2223_kde$ID <- as.factor(noca_2223_kde$ID)
noca_2324_kde_50_sf$ID <- as.factor(noca_2324_kde_50_sf$ID)
noca_2324_kde_95_sf$ID <- as.factor(noca_2324_kde_95_sf$ID)
noca_2324_kde$ID <- as.factor(noca_2324_kde$ID)
noca_2425_kde_50_sf$ID <- as.factor(noca_2425_kde_50_sf$ID)
noca_2425_kde_95_sf$ID <- as.factor(noca_2425_kde_95_sf$ID)
noca_2425_kde$ID <- as.factor(noca_2425_kde$ID)


# Run the function created in 4.1.
kde_2223_data <- process_kde_data(noca_2223_kde_50_sf, noca_2223_kde_95_sf, noca_2223_kde, "2022-2023")

kde_2324_data <- process_kde_data(noca_2324_kde_50_sf, noca_2324_kde_95_sf, noca_2324_kde, "2024-2025")

kde_2425_data <- process_kde_data(noca_2425_kde_50_sf, noca_2425_kde_95_sf, noca_2425_kde, "2024-2025")


# Combine all years.
kde_all_data <- bind_rows(kde_2223_data, kde_2324_data, kde_2425_data)

# Convert Area to square kilometers (1 km² = 1,000,000 m²).
kde_all_data <- kde_all_data %>%
  mutate(Area_km2 = Area / 1e6)

# Add a new grouping column based on Site values.
kde_all_data <- kde_all_data %>%
  mutate(Group = case_when(
    Site %in% c(1, 2) ~ "MBO",
    Site %in% c(3, 4) ~ "BDU",
    Site %in% c(5, 6) ~ "CON",
    TRUE ~ "Other"
  ))

# Order the Group factor to ensure correct facet ordering.
kde_all_data$Group <- factor(kde_all_data$Group, levels = c("MBO", "BDU", "CON"))

# Create a color scale for the sites with KDE intensity.
kde_color_scale <- c("MBO_50" = "#4B8378", "MBO_95" = "#44AA99", 
                     "BDU_50" = "#887F51", "BDU_95" = "#DDCC77", 
                     "CON_50" = "#874A54", "CON_95" = "#CC6677")

# Generate a grouping variable for color based on site and KDE.
kde_all_data <- kde_all_data %>%
  mutate(Site_KDE = paste0(Group, "_", KDE))

# Summarize mean area for each Year, Group, and KDE.
kde_summary <- kde_all_data %>%
  group_by(Year, Group, KDE) %>%
  summarise(Mean_Area = mean(Area_km2, na.rm = TRUE))

kde_summary_all_years <- kde_all_data %>%
  group_by(Group, KDE) %>%
  summarise(
    Mean_Area = mean(Area_km2, na.rm = TRUE),
    SD_Area = sd(Area_km2, na.rm = TRUE),
    SE_Area = SD_Area / sqrt(n())  # Standard error
  )

## --- 4.3. Plot the mean area for each site and each year --- ##

# Filter data for KDE 50 and KDE 95 separately
kde_summary_50 <- kde_summary_all_years %>% filter(KDE == 50)
kde_summary_95 <- kde_summary_all_years %>% filter(KDE == 95)

# Define custom legend labels
custom_legend_labels <- c(
  "MBO_50" = "KDE 50% for MBO", "MBO_95" = "KDE 95% for MBO",
  "BDU_50" = "KDE 50% for BDU", "BDU_95" = "KDE 95% for BDU",
  "CON_50" = "KDE 50% for CON", "CON_95" = "KDE 95% for CON"
)

# Define the order of groups for the legend
legend_order <- c("MBO_50", "MBO_95", "BDU_50", "BDU_95", "CON_50", "CON_95")

# KDE 50 Plot
p1 <- ggplot(kde_summary_50, aes(x = Group, y = Mean_Area, fill = paste0(Group, "_", KDE))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = Mean_Area - SE_Area, ymax = Mean_Area + SE_Area), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(
    title = "Mean Home Range Size (KDE 50) per Site (all years combined)",
    x = "Site",
    y = "Mean Home Range Area (km²)",
    fill = "KDE Type"
  ) +
  scale_fill_manual(
    values = kde_color_scale, 
    labels = custom_legend_labels, 
    breaks = legend_order  # Reorder legend items
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),  # Adjust legend text size
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "black", linewidth = 0.2)
  )

# KDE 95 Plot
p2 <- ggplot(kde_summary_95, aes(x = Group, y = Mean_Area, fill = paste0(Group, "_", KDE))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = Mean_Area - SE_Area, ymax = Mean_Area + SE_Area), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(
    title = "Mean Home Range Size (KDE 95) per Site (all years combined)",
    x = "Site",
    y = "Mean Home Range Area (km²)",
    fill = "KDE Type"
  ) +
  scale_fill_manual(
    values = kde_color_scale, 
    labels = custom_legend_labels, 
    breaks = legend_order  # Reorder legend items
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "black", linewidth = 0.2)
  )

# Print both plots
print(p1)
print(p2)

# Recreate the Site_KDE variable to control fill color
kde_summary_all_years <- kde_summary_all_years %>%
  mutate(
    Site_KDE = paste0(Group, "_", KDE),
    KDE_label = factor(KDE, levels = c(50, 95), labels = c("KDE 50%", "KDE 95%"))
  )

# Plot: Side-by-side grouped bars for KDE 50% and 95%
ggplot(kde_summary_all_years, aes(x = Group, y = Mean_Area, fill = Site_KDE)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), color = "black") +
  geom_errorbar(aes(ymin = Mean_Area - SE_Area, ymax = Mean_Area + SE_Area),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = kde_color_scale) +
  labs(
    title = "Average Home Range Area (km²) by Site and KDE Level",
    x = "Site",
    y = "Mean Home Range Area (km²)",
    fill = "Site & KDE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey80"),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )

# Save the plot.
ggsave("kde_home_range_barplot.png",
       width = 8, height = 6, dpi = 300, bg = "transparent")

# Testing difference across site.
kde_50 <- kde_all_data %>%
  filter(KDE == 50) %>%        # numeric, not "50"
  st_drop_geometry()           # IMPORTANT for rstatix

kde_95 <- kde_all_data %>%
  filter(KDE == 95) %>%
  st_drop_geometry()

# Check normality for each KDE.
shapiro.test(kde_50$Area_km2)  # Core
shapiro.test(kde_95$Area_km2)  # Total
# Both not normal, let's try to log the values before moving on to Kruskal-Wallis.
kde_50 <- kde_50 %>% 
  mutate(log_area = log(Area)) %>% 
  mutate(log_area_km2 = log(Area_km2))

kde_95 <- kde_95 %>% 
  mutate(log_area = log(Area)) %>% 
  mutate(log_area_km2 = log(Area_km2))

shapiro.test(kde_50$log_area_km2)
shapiro.test(kde_95$log_area_km2)
# Now p-values are > 0.05 meaning the data is normal.

# Doing a Levene test for variance.
levene_test(kde_50, log_area ~ Group)
levene_test(kde_95, log_area ~ Group)

# Doign a one way anova.
kde_50$Group <- factor(kde_50$Group, levels = c("MBO","BDU","CON"))
kde_95$Group <- factor(kde_95$Group, levels = c("MBO","BDU","CON"))
anova_50 <- aov(log_area_km2 ~ Group, data = kde_50)
summary(anova_50)
# Not significant, but approaching a trend.
tukey_50 <- TukeyHSD(anova_50)

anova_95 <- aov(log_area_km2 ~ Group, data = kde_95)
summary(anova_95)
# Significant, let's run post-hoc tests.
tukey_95 <- TukeyHSD(anova_95)

# Prepare for plotting.

tukey_letters_50 <- multcompLetters(tukey_50$Group[,"p adj"])
tukey_letters_95 <- multcompLetters(tukey_95$Group[,"p adj"])

# Add group letters to the datasets
group_letters_50 <- data.frame(
  Group = names(tukey_letters_50$Letters),
  Letters = tukey_letters_50$Letters
)

group_letters_95 <- data.frame(
  Group = names(tukey_letters_95$Letters),
  Letters = tukey_letters_95$Letters
)

# Add means and letters to data.
means_50 <- kde_50 %>%
  group_by(Group) %>%
  summarize(mean_log_area = mean(log_area_km2)) %>%
  left_join(group_letters_50, by = "Group")

means_95 <- kde_95 %>%
  group_by(Group) %>%
  summarize(mean_log_area = mean(log_area_km2)) %>%
  left_join(group_letters_95, by = "Group")

# 1. Extract Tukey results to a data frame
tukey_95_df <- as.data.frame(tukey_95$Group)
tukey_95_df$comparison <- rownames(tukey_95_df)

# 2. Convert to tibble before using dplyr verbs
tukey_95_df <- as_tibble(tukey_95_df)

# 3. Separate and filter significant comparisons
sig_95 <- tukey_95_df %>%
  filter(`p adj` < 0.05) %>%
  separate(comparison, into = c("group1", "group2"), sep = "-")

# 4. Now safely extract pairs for ggsignif
comparisons_95 <- sig_95 %>%
  mutate(pair = purrr::map2(group1, group2, ~ sort(c(.x, .y)))) %>%
  pull(pair)

# Define a named vector of colors
site_colors <- c(
  "MBO" = "#C4B14D",  # Red
  "BDU" = "#DF646D",  # Blue
  "CON" = "#973769"   # Green
)

base_theme <- theme_minimal(base_size = 14) +  # Set base text size
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),  # Bold + centered title
    axis.title = element_text(size = 16),       # Axis title font size
    axis.text = element_text(size = 14),        # Axis tick label font size
    legend.position = "none"
  )
# Home Ranges plot

y_limits <- range(
  c(kde_50$log_area_km2, kde_95$log_area_km2),
  na.rm = TRUE
)

y_limits[2] <- y_limits[2] + 1

y_pos <- max(kde_95$log_area_km2, na.rm = TRUE) + 0.5

p1 <- ggplot(kde_50, aes(x = Group, y = log_area_km2, fill = Group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, color = "black") +
  scale_fill_manual(values = site_colors) +
  scale_x_discrete(labels = c("MBO\n(peri-urban)", "BDU\n(suburban)", "CON\n(dense urban)")) +
  coord_cartesian(ylim = y_limits) +
  labs(
    title = "Core home range size (KDE 50%)",
    y = "Log(area in km²)",
    x = "Site"
  ) +
  base_theme


# Total Home Range Plot (95%) – show significant brackets only

p2 <- ggplot(kde_95, aes(x = Group, y = log_area_km2, fill = Group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, color = "black") +
  scale_fill_manual(values = site_colors) +
  scale_x_discrete(labels = c("MBO\n(peri-urban)", "BDU\n(suburban)", "CON\n(dense urban)")) +
  coord_cartesian(ylim = y_limits) +
  labs(
    title = "Total home range size (KDE 95%)",
    y = "Log(area in km²)",
    x = "Site"
  ) +
  base_theme +
  geom_signif(
    comparisons = comparisons_95,
    annotations = rep("*", length(comparisons_95)),
    y_position = seq(4.5, by = 0.6, length.out = length(comparisons_95)),
    tip_length = 0.02
  )


# Combine plots

full_plot <- p1 + p2 + plot_layout(ncol = 2)

full_plot

ggsave("home_range_plot_50.png", plot = p1, 
       width = 8, height = 6, dpi = 300, bg = "transparent")
ggsave("home_range_plot_95.png", plot = p2, 
       width = 8, height = 6, dpi = 300, bg = "transparent")
ggsave("home_ranges_plot.png", plot = full_plot,
       width = 14, height = 10, dpi = 300, bg = "transparent")





