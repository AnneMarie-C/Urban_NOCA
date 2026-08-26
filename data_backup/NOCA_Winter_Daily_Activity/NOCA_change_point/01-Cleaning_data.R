## June 19, 2025
## Anne-Marie Cousineau

# =============================================================================================== #
# DATA CLEANING FOR NORTHERN CARDINAL CHANGEPOINT                                                 #
# =============================================================================================== #

# Loading appropriate packages.
library(tidyverse)
library(dplyr)
library(lubridate)
library(patchwork)
library(tidyr)


# Load the raw data files needed.
cp_onset <- read.csv("00_raw_data/change_point_onset.csv")
cp_end <- read.csv("00_raw_data/change_point_end.csv")
temp_data <- read.csv("00_raw_data/weather_data.csv")

# Remove extra rows not needed in the temp_data.
temp_data <- temp_data %>% 
  select(-X, -X.1, -X.2, -X.3, -X.4, -X.5, -X.6, -X.7, -X.8, -X.9)

# Adding the a column for the onset and end datasets to later combine them.
cp_onset <- cp_onset %>% 
  mutate(status = "onset") %>% 
  select(status, everything())

cp_end <- cp_end %>% 
  mutate(status = "end") %>% 
  select(status, everything())

# Adding a column with the corresponding birds' ID.
lookup <- data.frame(
  motusTagID = c(68906,68907,68908,68909,68910,68911,68912,68913,68914,68915,68916,
                 68917,68918,68919,68920,68921,68922,68923,68924,68925,77867,77868,77869,
                 77870,77871,77872,77873,77874,77875,77876,77877,77878,77879,77880,77881,
                 77882,77883,77884,77885,77886,77887,77888,77889,77890,77891,77892,77893,77894,
                 77895,77896,77897,77898,77899,77900,77901,77902,77903,77904,77905,
                 77906,77907,77908,90326,90327,90328,90329,90330,90331,90332,90333,90334,
                 90335,90336,90337,90338,90339,90340,90341,90342,90343,90344,90345,
                 90346,90347,90348,90349,90350,90351,90352,90353,90354,90355,90356,90357,
                 90358,90359,90360,90361,90362,90363),
  ID = c("206","207","208","209","210","211","212","213","214","215","216","217","218","219","220",
         "221","222","223","224","225","90","91","92","93","94","95","96","97","98","99","100","101",
         "102","103","104","105","106","107","108","109","110","111","112","113","114","115","116",
         "117","118","119","120","121","123","124","125","126","127","129","130","131","133","134",
         "608","609","610","611","612","613","614","615","616","617","618","619","620","621","622",
         "623","624","625","626","627","628","629","630","631","632","633","634","635","636","637",
         "638","639","640","641","642","643","644","645"))

cp_onset <- cp_onset %>% 
  left_join(lookup, by = "motusTagID") %>% 
  select(ID, everything())
cp_end <- cp_end %>% 
  left_join(lookup, by = "motusTagID") %>% 
  select(ID, everything())

# Adding site to the data frame.
site.df <- read.csv("00_raw_data/Site.csv")
site.df$ID <- as.factor(site.df$ID)
site.df$Site <- as.factor(site.df$Site)

cp_onset <- cp_onset %>%
  left_join(site.df, by = "ID")

cp_end <- cp_end %>%
  left_join(site.df, by = "ID")

# Filter out NAs.
cp_onset <- cp_onset %>%
  filter(!is.na(Site))

cp_end <- cp_end %>%
  filter(!is.na(Site))

# Rename "Dates" column to "Date".
cp_onset <- cp_onset %>% 
  rename(Date = Dates)

cp_end <- cp_end %>% 
  rename(Date = Dates)

# Adding weather data to the corresponding dates.
cp_onset$Date <- as.Date(cp_onset$Date)
cp_end$Date <- as.Date(cp_end$Date)
temp_data$Date <- as.Date(temp_data$Date)

# Add a winter number for an easier time when plotting.
cp_onset <- cp_onset %>%
  mutate(Winter = case_when(
    Date >= as.Date("2022-09-01") & Date < as.Date("2023-07-01") ~ 1,
    Date >= as.Date("2023-10-01") & Date < as.Date("2024-08-03") ~ 2,
    Date >= as.Date("2024-10-01") & Date < as.Date("2025-07-01") ~ 3
  ))

cp_end <- cp_end %>%
  mutate(Winter = case_when(
    Date >= as.Date("2022-09-01") & Date < as.Date("2023-07-01") ~ 1,
    Date >= as.Date("2023-10-01") & Date < as.Date("2024-08-03") ~ 2,
    Date >= as.Date("2024-10-01") & Date < as.Date("2025-07-01") ~ 3
  ))

cp_onset$Winter <- as.factor(cp_onset$Winter)
cp_end$Winter <- as.factor(cp_end$Winter)

# Keeping only temp data from PETI.
temp_data <- temp_data %>% 
  filter(Station_ID == "PETI")

cp_onset <- cp_onset %>% 
  left_join(temp_data, by = "Date")

cp_end <- cp_end %>% 
  left_join(temp_data, by = "Date")

# Filter out "bad" quality data.
cp_onset <- cp_onset %>%
  filter(postfilter == "good")

cp_end <- cp_end %>%
  filter(postfilter == "good")

cp_onset <- cp_onset %>%
  mutate(
    # Combine date + time (in decimal hours) into POSIXct UTC datetime
    datetime_utc = as.POSIXct(Date, tz = "UTC") + 
      hours(floor(time)) + 
      minutes(round((time %% 1) * 60)),
    
    # Convert to local Montreal time (handles DST automatically)
    datetime_mtl = with_tz(datetime_utc, tzone = "America/Montreal")
  )

cp_end <- cp_end %>%
  mutate(
    # Combine date + time (in decimal hours) into POSIXct UTC datetime
    datetime_utc = as.POSIXct(Date, tz = "UTC") + 
      hours(floor(time)) + 
      minutes(round((time %% 1) * 60)),
    
    # Convert to local Montreal time (handles DST automatically)
    datetime_mtl = with_tz(datetime_utc, tzone = "America/Montreal")
  )

# Filtering out NAs one last time before saving.
cp_onset <- na.omit(cp_onset)
cp_end <- na.omit(cp_end)

write.csv(cp_onset, "~/Library/CloudStorage/OneDrive-McGillUniversity/NOCA_Winter_Activity/NOCA_change_point/02_cleaned_data/cp_onset_cleaned.csv", row.names = FALSE)
write.csv(cp_end, "~/Library/CloudStorage/OneDrive-McGillUniversity/NOCA_Winter_Activity/NOCA_change_point/02_cleaned_data/cp_end_cleaned.csv", row.names = FALSE)

# Plotting just to have an idea if the transformations make sense...
# Prepare the data
cp_onset <- cp_onset %>%
  mutate(
    Date = as.Date(datetime_mtl),
    Hour = hour(datetime_mtl) + minute(datetime_mtl) / 60
  )

ggplot(cp_onset, aes(x = Date)) +
  # Temperature line with gradient color
  geom_line(aes(y = Min_temp, color = Min_temp), size = 0.8) +
  
  # Change point time as orange points
  geom_point(aes(y = Hour), color = "#D55E00", size = 1.5, alpha = 0.6) +
  
  # Color scale for temperature
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Min Temp (°C)"
  ) +
  
  # Facet by Winter
  facet_wrap(~ Winter, scales = "free_x") +
  
  scale_y_continuous(
    name = "Hour of Change Point / Temperature (°C)",
    breaks = seq(-30, 30, by = 5),
    limits = c(-30, 30)
  ) +
  scale_x_date(name = "Date") +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Morning Change Points and Temperature by Winter",
    subtitle = "Points = Change Point Hour; Line = Temperature"
  )

# Plotting just to have an idea if the transformations make sense...
# Prepare the data
cp_end <- cp_end %>%
  mutate(
    Date = as.Date(datetime_mtl),
    Hour = hour(datetime_mtl) + minute(datetime_mtl) / 60
  )

ggplot(cp_end, aes(x = Date)) +
  # Temperature line with gradient color
  geom_line(aes(y = Min_temp, color = Min_temp), size = 0.8) +
  
  # Change point time as orange points
  geom_point(aes(y = Hour), color = "#0072B2", size = 1.5, alpha = 0.6) +
  
  # Color scale for temperature
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Min Temp (°C)"
  ) +
  
  # Facet by Winter
  facet_wrap(~ Winter, scales = "free_x") +
  
  scale_y_continuous(
    name = "Hour of Change Point / Temperature (°C)",
    breaks = seq(-30, 30, by = 5),
    limits = c(-30, 30)
  ) +
  scale_x_date(name = "Date") +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Evening Change Points and Temperature by Winter",
    subtitle = "Points = Change Point Hour; Line = Temperature"
  )