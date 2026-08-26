
## March 18, 2025
## Anne-Marie Cousineau

# =============================================================================================== #
# DATA CLEANING FOR URBAN NORTHERN CARDINAL PROJECT                                               #
# =============================================================================================== #

# Loading appropriate packages.
library(tidyverse)
library(dplyr)
library(geosphere)

# =============================================================================================== #
# DATA COLLECTED DURING WINTER 2022 & 2023 ONLY AT MBO                                            #
# =============================================================================================== #
# Load the raw data file.
noca_2223 <- read.csv("00_raw_data/noca_2022-2023_raw.csv")
print(noca_2223)

# Fix the date format for R to recognize it as dd/mm/yyyy.
noca_2223$Date <- as.Date(noca_2223$Date, format = "%d/%m/%Y")

# Rename the Tag.ID column for ID.
noca_2223 <- noca_2223 %>% 
  dplyr::rename(ID = Tag.ID)

# Fix the ID column for it to be factors and not numbers.
noca_2223$ID <- as.factor(noca_2223$ID)

# The column Lon doesn't have negatives so we'll add them to all the values.
noca_2223 <- noca_2223 %>%
  mutate(Lon = ifelse(Lon > 0, -Lon, Lon))

# Removing data where either lat or lon is NA.
noca_2223 <- noca_2223 %>%
  filter(!is.na(Lat) & !is.na(Lon))

# Removing data where the date is NA.
noca_2223 <- noca_2223 %>%
  filter(!is.na(Date))

# Removing data where the lat or lon is 0.
noca_2223 <- noca_2223 %>%
  filter(Lat != 0 & Lon != 0)

# Removing all data points that have a gain superior to 40 and a signal of less than 140 for
# older receivers, and a signal of less than 190 for newer receivers.
# Here we keep the NAs as they're the release points.
noca_2223 <- noca_2223 %>%
  filter(Gain <= 40)   # This is to remove data points with gains over 40.
noca_2223 <- noca_2223 %>%
  filter(!(Receiver == "New" & Signal <=189))   # This is to remove data points with signals lower
# than 190 for the new receivers.
noca_2223 <- noca_2223 %>%
  filter(!(Receiver == "Old" & Signal <= 139))   # This is to remove data points with signals lower
# than 140 for the old receivers.

# Removing data with less than 5 observations.
noca_2223 <- noca_2223 %>% 
  group_by(ID) %>%    # Group by the individual identifier "ID".
  filter(n() >= 5) %>%   # Keep only groups with 5 or more rows.
  ungroup()

# After inspection of the results, there's a point with a Lat of around 45.48, it should be removed.
noca_2223 <- noca_2223 %>% 
  filter(Lat != 45.4873)

# =============================================================================================== #
# DATA COLLECTED DURING WINTER 2023 & 2024 AT MBO, BDU AND CON                                    #
# =============================================================================================== #
# Load the raw data file.
noca_2324 <- read.csv("00_raw_data/noca_2023-2024_raw.csv")

# Fix the date format for R to recognize it as dd/mm/yyyy.
noca_2324$Date <- as.Date(noca_2324$Date, format = "%d/%m/%Y")

# Fix the ID column for it to be factors and not numbers.
noca_2324$ID <- as.factor(noca_2324$ID)

# Removing data where either lat or lon is NA.
noca_2324 <- noca_2324 %>%
  filter(!is.na(Lat) & !is.na(Lon))

# Removing data where the date is NA.
noca_2324 <- noca_2324 %>%
  filter(!is.na(Date))

# Removing data where the lat or lon is 0.
noca_2324 <- noca_2324 %>%
  filter(Lat != 0 & Lon != 0)

# Removing all data points that have a gain superior to 40 and a signal of less than 140 for
# older receivers, and a signal of less than 190 for newer receivers.
# Here we keep the NAs as they're the release points.
noca_2324 <- noca_2324 %>%
  filter(Gain <= 40)   # This is to remove data points with gains over 40.
noca_2324 <- noca_2324 %>%
  filter(!(Receiver == "NEW" & Signal <=189))   # This is to remove data points with signals lower
#  than 190 for the new receivers.
noca_2324 <- noca_2324 %>%
  filter(!(Receiver == "OLD" & Signal <= 139))   # This is to remove data points with signals lower
#  than 140 for the old receivers.

# Removing data with less than 5 observations.
noca_2324 <- noca_2324 %>% 
  group_by(ID) %>%    # Group by the individual identifier "ID".
  filter(n() >= 5) %>%   # Keep only groups with 5 or more rows.
  ungroup()

# After inspection of the results, there are 3 points that should be removed.
# 1 point at Lon -73.98851; another at Lon -73.84356 and the last one at Lat 45.40657.
noca_2324 <- noca_2324 %>% 
  filter(Lat != 45.40657) %>% 
  filter(Lon != -73.98851) %>% 
  filter(Lon != -73.84356)

# After discussion with Kyle & Barbara, we decided that the birds with points in both MBO and BDU
# will have the points assigned to both locations independently of the release site.
# All BDU points will be assigned to MBO site 1.

# NOCA 93b.
bird_id <- "93b"  # Select the bird ID
dates_to_keep <- as.Date(c("2023-11-02"))  # Dates that should NOT change.
# Reassign the site for the specific ID: keep original site for selected dates, change for all others

# Remove the selected dates from the dataset
noca_2324 <- noca_2324 %>%
  mutate(Site = case_when(
    ID == bird_id & !Date %in% dates_to_keep ~ 1,# Change to new site for all other dates
                          TRUE ~ Site))  # Keep the original site for the specified dates & all other birds

# NOCA 130.
bird_id <- "130"  # Select the bird ID
dates_to_keep <- as.Date(c("2023-10-23", "2023-10-25", "2023-10-27", "2023-11-06", "2023-11-15"))  # Dates that should NOT change.

# Remove the selected dates from the dataset
noca_2324 <- noca_2324 %>%
  mutate(Site = case_when(
    ID == bird_id & !Date %in% dates_to_keep ~ 1,# Change to new site for all other dates
    TRUE ~ Site))

# NOCA 121.
bird_id <- "121"  # Select the bird ID
dates_to_keep <- as.Date(c("2023-10-18", "2023-10-20", "2023-10-25"))  # Dates that should NOT change.

# Remove the selected dates from the dataset
noca_2324 <- noca_2324 %>%
  mutate(Site = case_when(
    ID == bird_id & !Date %in% dates_to_keep ~ 1,# Change to new site for all other dates
    TRUE ~ Site))

# =============================================================================================== #
# DATA COLLECTED DURING WINTER 2024 & 2025 AT MBO, BDU AND CON                                    #
# =============================================================================================== #
# Load the raw data file.
noca_2425 <- read.csv("00_raw_data/noca_2024-2025_raw.csv")

# Fix the date format for R to recognize it as dd/mm/yyyy.
noca_2425$Date <- as.Date(noca_2425$Date, format = "%m-%d-%Y")

# Rename the NOCA_ID column for ID.
noca_2425 <- noca_2425 %>% 
  dplyr::rename(ID = NOCA_ID)

# Fix the ID column for it to be factors and not numbers.
noca_2425$ID <- as.factor(noca_2425$ID)

# For some reason it seems that Lat & Lon are characters and not numbers. Run this to fix the
# issue.
noca_2425$Lon <- gsub("[^0-9.-]", "", noca_2425$Lon)  # Keep only numbers, dots, and dashes
noca_2425$Lat <- gsub("[^0-9.-]", "", noca_2425$Lat)
noca_2425$Lon <- as.numeric(noca_2425$Lon)  # Convert again
noca_2425$Lat <- as.numeric(noca_2425$Lat)

# Removing data where either lat or lon is NA.
noca_2425 <- noca_2425 %>%
  filter(!is.na(Lat) & !is.na(Lon))

# Removing data where the date is NA.
noca_2425 <- noca_2425 %>%
  filter(!is.na(Date))

# Removing data where the lat or lon is 0.
noca_2425 <- noca_2425 %>%
  filter(Lat != 0 & Lon != 0)

# Removing all data points that have a gain superior to 40 and a signal of less than 140 for
# older receivers, and a signal of less than 190 for newer receivers.
# Here we keep the NAs as they're the release points.
noca_2425 <- noca_2425 %>%
  filter(Gain <= 40)   # This is to remove data points with gains over 40.
noca_2425 <- noca_2425 %>%
  filter(!(Receiver == "NEW" & Signal <=189))   # This is to remove data points with signals lower
#  than 190 for the new receivers.
noca_2425 <- noca_2425 %>%
  filter(!(Receiver == "OLD" & Signal <= 139))   # This is to remove data points with signals lower
#  than 140 for the old receivers.

# Removing data with less than 5 observations.
noca_2425 <- noca_2425 %>% 
  group_by(ID) %>%    # Group by the individual identifier "ID".
  filter(n() >= 5) %>%   # Keep only groups with 5 or more rows.
  ungroup()

# After inspection of the results, there is a point that should be removed with Lon -73.83896.
# In addition, we make sure Lat and Lon are numerical values.
noca_2425 <- noca_2425 %>% 
  filter(Lon != -73.83896)%>% 
  filter(Lon != -73.93052)%>% 
  filter(Lon != -73.02889)

# After discussion with Kyle & Barbara, we decided that the birds with points in both MBO and BDU
# will have the points assigned to both locations independently of the release site.
# All BDU points will be assigned to MBO site 1.

# MODIFIED THE CODE TO REMOVE THE DATES IN BDU.

# NOCA 627.
bird_id <- "627"  # Select the bird ID
dates_to_keep <- as.Date(c("2024-11-10", "2024-11-08"))  # Dates that should NOT change.

# Remove the selected dates from the dataset
noca_2324 <- noca_2324 %>%
  mutate(Site = case_when(
    ID == bird_id & !Date %in% dates_to_keep ~ 1,# Change to new site for all other dates
    TRUE ~ Site))

# =============================================================================================== #
# ADJUSTING EVERY GPS LOCATION TO TAKE IN ACCOUNT ANGLE AND DISTANCE                                         #
# =============================================================================================== #
# Assigning new names for the most recently cleaned versions.
noca_2223_cleaned <- noca_2223
noca_2324_cleaned <- noca_2324
noca_2425_cleaned <- noca_2425

# After discussion with Kyle and Barbara I'll try to adjust the GPS locations to take in account
# the angle and gain.
# This is because we initially saw a very strong correlation of habitat selection with road length
# but feared that that was only the case because all the GPS locations were taken from the street.

# Assign distances for the gains. After investigating in QGIS, here are the assigned distances.
# Gain of 10 is about 10 meters away.
# Gain of 20 is about 15 meters away.
# Gain of 30 is about 25 meters away.
# Gain of 40 is about 35 meters away.

noca_2223_cleaned <- noca_2223_cleaned %>%
  mutate(Distance_new = case_when(
    Gain <= 10 ~ 10,
    Gain <= 20 ~ 15,
    Gain <= 30 ~ 25,
    Gain <= 40 ~ 35,
    TRUE ~ 0
  ))

# Apply destPoint() row-wise to generate the new lat lon according to Angle and Distance_new.
noca_2223_cleaned <- noca_2223_cleaned %>%
  rowwise() %>%
  mutate(
    new_coords = list(
      if (!is.na(Angle)) {
        destPoint(c(Lon, Lat), Angle, Distance_new)
      } else {
        c(Lon, Lat)
      }
    ),
    new_lon = new_coords[[1]],
    new_lat = new_coords[[2]]
  ) %>%
  ungroup() %>%
  select(-new_coords)  # Optional: remove list column if no longer needed

# Plot to verify that it ran properly.
ggplot(noca_2223_cleaned) +
  geom_point(aes(x = Lon, y = Lat), color = "blue", size = 3) +
  geom_point(aes(x = new_lon, y = new_lat), color = "red", size = 3) +
  geom_segment(aes(x = Lon, y = Lat, xend = new_lon, yend = new_lat),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray50") +
  coord_fixed() +
  labs(title = "Original (blue) and Shifted (red) Points",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# Repeat for 2324.
noca_2324_cleaned <- noca_2324_cleaned %>%
  mutate(Distance_new = case_when(
    Gain <= 10 ~ 10,
    Gain <= 20 ~ 15,
    Gain <= 30 ~ 25,
    Gain <= 40 ~ 35,
    TRUE ~ 0
  ))

# Apply destPoint() row-wise to generate the new lat lon according to Angle and Distance_new.
noca_2324_cleaned <- noca_2324_cleaned %>%
  rowwise() %>%
  mutate(
    new_coords = list(
      if (!is.na(Angle)) {
        destPoint(c(Lon, Lat), Angle, Distance_new)
      } else {
        c(Lon, Lat)
      }
    ),
    new_lon = new_coords[[1]],
    new_lat = new_coords[[2]]
  ) %>%
  ungroup() %>%
  select(-new_coords)  # Optional: remove list column if no longer needed

# Plot to verify that it ran properly.
noca_2324_cleaned_plot <- noca_2324_cleaned %>% 
  filter( Site %in% c(5, 6))

ggplot(noca_2324_cleaned_plot) +
  geom_point(aes(x = Lon, y = Lat), color = "blue", size = 3) +
  geom_point(aes(x = new_lon, y = new_lat), color = "red", size = 3) +
  geom_segment(aes(x = Lon, y = Lat, xend = new_lon, yend = new_lat),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray50") +
  coord_fixed() +
  labs(title = "Original (blue) and Shifted (red) Points",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# Repeat for 2425.
noca_2425_cleaned <- noca_2425_cleaned %>%
  mutate(Distance_new = case_when(
    Gain <= 10 ~ 10,
    Gain <= 20 ~ 15,
    Gain <= 30 ~ 25,
    Gain <= 40 ~ 35,
    TRUE ~ 0
  ))

# Apply destPoint() row-wise to generate the new lat lon according to Angle and Distance_new.
noca_2425_cleaned <- noca_2425_cleaned %>%
  rowwise() %>%
  mutate(
    new_coords = list(
      if (!is.na(Angle)) {
        destPoint(c(Lon, Lat), Angle, Distance_new)
      } else {
        c(Lon, Lat)
      }
    ),
    new_lon = new_coords[[1]],
    new_lat = new_coords[[2]]
  ) %>%
  ungroup() %>%
  select(-new_coords)  # Optional: remove list column if no longer needed

# Plot to verify that it ran properly.
noca_2425_cleaned_plot <- noca_2425_cleaned %>% 
  filter( Site %in% c(5, 6))

ggplot(noca_2425_cleaned_plot) +
  geom_point(aes(x = Lon, y = Lat), color = "blue", size = 3) +
  geom_point(aes(x = new_lon, y = new_lat), color = "red", size = 3) +
  geom_segment(aes(x = Lon, y = Lat, xend = new_lon, yend = new_lat),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray50") +
  coord_fixed() +
  labs(title = "Original (blue) and Shifted (red) Points",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# =============================================================================================== #
# SAVING CLEANED DATA AS CSV FILES FOR USE IN SCRIPT NO.2                                         #
# =============================================================================================== #

# Writing the csv files in the folder "data-processed".
write.csv(noca_2223_cleaned,"02_cleaned_data/noca_2223_cleaned.csv", row.names = FALSE)
write.csv(noca_2324_cleaned,"02_cleaned_data/noca_2324_cleaned.csv", row.names = FALSE)
write.csv(noca_2425_cleaned,"02_cleaned_data/noca_2425_cleaned.csv", row.names = FALSE)