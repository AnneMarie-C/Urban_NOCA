# =============================================================================================== #
# GENERATING RANDOM POINTS FOR THE RSF ANALYSES                                                   #
# =============================================================================================== #

# Installing and loading relevant packages
library(dplyr)
library(sf)
library(sp)
library(geosphere)
library(ggplot2)
library(mapview)

# Code adapted from Random points_from Kristen scprit. It was created by Kristen Lalla and Alison.

# Setting a seed makes it so R always "spits" out random points at the same locations.
# For this to work, you need to add set.seed(232409) before any randomization process, whether a 
# function or a command.

set.seed <- 232409

# =============================================================================================== #
# 1. CUSTOM FUNCTION TO GENERATE A RANDOM DISTANCE AND ANGLE                                      #
# =============================================================================================== #

## --- 1.1. Function to generate random Distance --- #

set.seed(232409)
randomDist <- function(max_distance, n) {
  # Generate random distances between 0 and max_distance
  runif(n, min = 0, max = max_distance)
}

## --- 1.2. Function to generate random Angle --- #

set.seed(232409)
randomDistAngle <- function(n) {
  randang <- runif(n, 0, 360) # Generates random angles between 0 and 360.
  # Output vector of random angles.
  randang
}

## --- 1.3. Function to calculate XY coordinates from distance, angle, easting, and northing --- ##

set.seed(232409)
calculateXY <- function(dist, angle, col_x, col_y) {
  # Calculate distance in meters between origin (0,0) (i.e. the habitat centroid) and new point.
  # The calculation is simpler if you use angles rather than radians.
  xdist <- dist * cos(angle)  
  ydist <- dist * sin(angle)
  X <- col_x + xdist  # Add distance to origin.
  Y <- col_y + ydist
  data.frame(  # Output data frame with x and y.
    randX = X, randY =Y
  )
}

# =============================================================================================== #
# 2. LOAD & PREPARE THE DATA FOR THE CARDINALS TRACKED IN 2023-2024                               #
# =============================================================================================== #
noca_2324 <- read.csv("02_cleaned_data/noca_2324_cleaned.csv")

# Make sure ID is a factor.
noca_2324$ID <- as.factor(noca_2324$ID)
noca <- noca_2324

# Add column of distance maximums since that what I want to use to create random points
# In my case, for the cardinals, I want to create the points within 300m.
noca <- noca %>% # Add column of maxs.
  group_by(ID) %>% 
  mutate(rp.distance.m = 300)

# Remove the old Lat and Lon, replace the new ones so that the function doesn't have to change.
#noca <- noca %>% 
  #dplyr::select(-Lat, -Lon)
#noca <- noca %>% 
  #rename(Lat = new_lat,
         #Lon = new_lon)

# =============================================================================================== #
# 3. GENERATING RANDOM POINTS WITH A 1:1 RATIO                                                    #
# =============================================================================================== #
# 1:1 RATIO for the Random Points (means that for each "real" observation, 1 random point will be
# generated).

set.seed(232409)

## --- 3.1. Function to generate random latitude and longitude --- ##

calculateLatLon <- function(lat, lon, distance, bearing) {
  # Use the destPoint function from geosphere to calculate new coordinates.
  # Lat, lon are the starting point.
  # Distance is in meters.
  # Bearing is in degrees.
  new_coords <- destPoint(c(lon, lat), bearing, distance)
  # Return new latitude and longitude.
  return(new_coords)
}

set.seed(232409)

# Now apply it in `mutate()` function
randomPts <- replicate(n = 1, 
                       noca %>% 
                         dplyr::select(ID, Date, Site, Lat, Lon) %>%
                         group_by(ID) %>% 
                         dplyr::mutate(
                           Type = "random",
                           randist = randomDist(300, n()),  # Pass '300' and the length of the group.
                           randang = randomDistAngle(n())  # Generate random angles for each entry.
                         ), simplify = FALSE)

RP1 <- do.call(rbind, randomPts)
RP1 <- RP1 %>% 
  arrange(ID)
randomPts <- as.data.frame(RP1)

# Generate new latitudes and longitudes (the code above generates X and Y points).
set.seed(232409)
randomLatLon <- mapply(calculateLatLon, randomPts$Lat, randomPts$Lon, randomPts$randist, randomPts$randang)
randomPts$randLat <- randomLatLon[1, ]
randomPts$randLon <- randomLatLon[2, ]

# View the results
head(randomPts)


# =============================================================================================== #
# 4. CYCLE THROUGH EACH INDIVIDUAL AND CREATE A MAP                                               #
# =============================================================================================== #

set.seed(232409)

## --- 4.1. Function to generate a random distance constrained to 300m --- ##

randomDist <- function(max_distance, n) {
  # Generate `n` random values between 0 and `max_distance`
  runif(n, min = 0, max = max_distance)
}

set.seed(232409)

## --- 4.2. Function to calculate new lat-lon based on distance and angle --- ##

calculateLatLon <- function(lat, lon, dist, angle) {
  # Radius of Earth in meters
  R <- 6371000
  set.seed(232409)
  # Convert latitude and longitude to radians
  lat_rad <- lat * pi / 180
  lon_rad <- lon * pi / 180
  set.seed(232409)
  # Convert distance and angle to radians
  delta_lat <- dist / R * cos(angle * pi / 180)
  delta_lon <- dist / (R * cos(lat_rad)) * sin(angle * pi / 180)
  set.seed(232409)
  # Calculate new latitude and longitude
  new_lat <- lat_rad + delta_lat
  new_lon <- lon_rad + delta_lon
  set.seed(232409)
  # Convert radians back to degrees
  c(new_lat * 180 / pi, new_lon * 180 / pi)
}

## --- 4.3.Loop through each NOCA ID and plot the resulting points --- ##

set.seed(232409)

for (id in unique(noca$ID)) {
  # Subset data for the current ID
  temp <- subset(noca, noca$ID == id)
  set.seed(232409)
  # Debug: Check the structure of `temp`
  print(paste("Processing ID:", id, "with", nrow(temp), "rows"))
  set.seed(232409)
  # Convert observation points to an sf object
  obs <- st_as_sf(temp, coords = c("Lon", "Lat"), crs = "EPSG:4326")
  set.seed(232409)
  # Generate random distances and angles
  rand_distances <- randomDist(300, nrow(temp))  # Pass max_distance and number of rows
  rand_angles <- randomDistAngle(nrow(temp))    # Generate random angles
  set.seed(232409)
  # Debug: Check generated distances and angles
  print(head(rand_distances))
  print(head(rand_angles))
  set.seed(232409)
  # Add random distances and angles to `temp`
  rand_points <- temp %>%
    mutate(
      randist = rand_distances,
      randang = rand_angles
    )
  set.seed(232409)
  
  # Calculate new random latitudes and longitudes
  
  randomLatLon <- mapply(
    calculateLatLon,
    rand_points$Lat,
    rand_points$Lon,
    rand_points$randist,
    rand_points$randang
  )
  
  set.seed(232409)
  # Debug: Check randomLatLon
  print(head(randomLatLon))
  
  set.seed(232409)
  # Assign calculated coordinates
  rand_points$randLat <- randomLatLon[1, ]
  rand_points$randLon <- randomLatLon[2, ]
  
  set.seed(232409)
  # Convert random points to sf object
  rand_sf <- st_as_sf(rand_points, coords = c("randLon", "randLat"), crs = "EPSG:4326")
  set.seed(232409)
  # Plot observation and random points
  print(
    mapview(obs, col.regions = 'blue') + 
      mapview(rand_sf, col.regions = 'red')
  )
}

## --- 4.5. Save the random points in a dataframe --- ##
write.csv(randomPts, "03_processed_data/noca_2324_ran_points.csv")

# =============================================================================================== #
# 5. LOAD & PREPARE THE DATA FOR THE CARDINALS TRACKED IN 2022-2023                               #
# =============================================================================================== #
noca_2223 <- read.csv("02_cleaned_data/noca_2223_cleaned.csv")

# Make sure ID is a factor.
noca_2223$ID <- as.factor(noca_2223$ID)
noca <- noca_2223

# Add column of distance maximums since that what I want to use to create random points
# In my case, for the cardinals, I want to create the points within 300m.
noca <- noca %>% # Add column of maxs.
  group_by(ID) %>% 
  mutate(rp.distance.m = 300)

# Add a column with the site being 1 for MBO.
noca <- noca %>% 
  mutate(Site = 1)

# Remove the old Lat and Lon, replace the new ones so that the function doesn't have to change.
#noca <- noca %>% 
  #dplyr::select(-Lat, -Lon)
#noca <- noca %>% 
  #rename(Lat = new_lat,
         #Lon = new_lon)

# =============================================================================================== #
# 6. GENERATING RANDOM POINTS WITH A 1:1 RATIO                                                    #
# =============================================================================================== #
# 1:1 RATIO for the Random Points (means that for each "real" observation, 1 random point will be
# generated).

set.seed(232409)

## --- 6.1. Function to generate random latitude and longitude --- ##

calculateLatLon <- function(lat, lon, distance, bearing) {
  # Use the destPoint function from geosphere to calculate new coordinates.
  # Lat, lon are the starting point.
  # Distance is in meters.
  # Bearing is in degrees.
  new_coords <- destPoint(c(lon, lat), bearing, distance)
  # Return new latitude and longitude.
  return(new_coords)
}

set.seed(232409)

# Now apply it in `mutate()` function
randomPts <- replicate(n = 1, 
                       noca %>% 
                         dplyr::select(ID, Date, Site, Lat, Lon) %>%
                         group_by(ID) %>% 
                         dplyr::mutate(
                           Type = "random",
                           randist = randomDist(300, n()),  # Pass '300' and the length of the group.
                           randang = randomDistAngle(n())  # Generate random angles for each entry.
                         ), simplify = FALSE)

RP1 <- do.call(rbind, randomPts)
RP1 <- RP1 %>% 
  arrange(ID)
randomPts <- as.data.frame(RP1)

# Generate new latitudes and longitudes (the code above generates X and Y points).
set.seed(232409)
randomLatLon <- mapply(calculateLatLon, randomPts$Lat, randomPts$Lon, randomPts$randist, randomPts$randang)
randomPts$randLat <- randomLatLon[1, ]
randomPts$randLon <- randomLatLon[2, ]

# View the results
head(randomPts)


# =============================================================================================== #
# 7. CYCLE THROUGH EACH INDIVIDUAL AND CREATE A MAP                                               #
# =============================================================================================== #

set.seed(232409)

## --- 7.1. Function to generate a random distance constrained to 300m --- ##

randomDist <- function(max_distance, n) {
  # Generate `n` random values between 0 and `max_distance`
  runif(n, min = 0, max = max_distance)
}

set.seed(232409)

## --- 7.2. Function to calculate new lat-lon based on distance and angle --- ##

calculateLatLon <- function(lat, lon, dist, angle) {
  # Radius of Earth in meters
  R <- 6371000
  set.seed(232409)
  # Convert latitude and longitude to radians
  lat_rad <- lat * pi / 180
  lon_rad <- lon * pi / 180
  set.seed(232409)
  # Convert distance and angle to radians
  delta_lat <- dist / R * cos(angle * pi / 180)
  delta_lon <- dist / (R * cos(lat_rad)) * sin(angle * pi / 180)
  set.seed(232409)
  # Calculate new latitude and longitude
  new_lat <- lat_rad + delta_lat
  new_lon <- lon_rad + delta_lon
  set.seed(232409)
  # Convert radians back to degrees
  c(new_lat * 180 / pi, new_lon * 180 / pi)
}

## --- 7.3.Loop through each NOCA ID and plot the resulting points --- ##

set.seed(232409)

for (id in unique(noca$ID)) {
  # Subset data for the current ID
  temp <- subset(noca, noca$ID == id)
  set.seed(232409)
  # Debug: Check the structure of `temp`
  print(paste("Processing ID:", id, "with", nrow(temp), "rows"))
  set.seed(232409)
  # Convert observation points to an sf object
  obs <- st_as_sf(temp, coords = c("Lon", "Lat"), crs = "EPSG:4326")
  set.seed(232409)
  # Generate random distances and angles
  rand_distances <- randomDist(300, nrow(temp))  # Pass max_distance and number of rows
  rand_angles <- randomDistAngle(nrow(temp))    # Generate random angles
  set.seed(232409)
  # Debug: Check generated distances and angles
  print(head(rand_distances))
  print(head(rand_angles))
  set.seed(232409)
  # Add random distances and angles to `temp`
  rand_points <- temp %>%
    mutate(
      randist = rand_distances,
      randang = rand_angles
    )
  set.seed(232409)
  
  # Calculate new random latitudes and longitudes
  
  randomLatLon <- mapply(
    calculateLatLon,
    rand_points$Lat,
    rand_points$Lon,
    rand_points$randist,
    rand_points$randang
  )
  
  set.seed(232409)
  # Debug: Check randomLatLon
  print(head(randomLatLon))
  
  set.seed(232409)
  # Assign calculated coordinates
  rand_points$randLat <- randomLatLon[1, ]
  rand_points$randLon <- randomLatLon[2, ]
  
  set.seed(232409)
  # Convert random points to sf object
  rand_sf <- st_as_sf(rand_points, coords = c("randLon", "randLat"), crs = "EPSG:4326")
  set.seed(232409)
  # Plot observation and random points
  print(
    mapview(obs, col.regions = 'blue') + 
      mapview(rand_sf, col.regions = 'red')
  )
}

## --- 7.5. Save the random points in a dataframe --- ##
write.csv(randomPts, "03_processed_data/noca_2223_ran_points.csv")

# =============================================================================================== #
# 8. LOAD & PREPARE THE DATA FOR THE CARDINALS TRACKED IN 2024-2025                               #
# =============================================================================================== #
noca_2425 <- read.csv("02_cleaned_data/noca_2425_cleaned.csv")

# Make sure ID is a factor.
noca_2425$ID <- as.factor(noca_2425$ID)
noca <- noca_2425

# Add column of distance maximums since that what I want to use to create random points
# In my case, for the cardinals, I want to create the points within 300m.
noca <- noca %>% # Add column of maxs.
  group_by(ID) %>% 
  mutate(rp.distance.m = 300)

# Remove the old Lat and Lon, replace the new ones so that the function doesn't have to change.
#noca <- noca %>% 
  #dplyr::select(-Lat, -Lon)
#noca <- noca %>% 
  #rename(Lat = new_lat,
         #Lon = new_lon)

# =============================================================================================== #
# 9. GENERATING RANDOM POINTS WITH A 1:1 RATIO                                                    #
# =============================================================================================== #
# 1:1 RATIO for the Random Points (means that for each "real" observation, 1 random point will be
# generated).

set.seed(232409)

## --- 9.1. Function to generate random latitude and longitude --- ##

calculateLatLon <- function(lat, lon, distance, bearing) {
  # Use the destPoint function from geosphere to calculate new coordinates.
  # Lat, lon are the starting point.
  # Distance is in meters.
  # Bearing is in degrees.
  new_coords <- destPoint(c(lon, lat), bearing, distance)
  # Return new latitude and longitude.
  return(new_coords)
}

set.seed(232409)

# Now apply it in `mutate()` function
randomPts <- replicate(n = 1, 
                       noca %>% 
                         dplyr::select(ID, Date, Site, Lat, Lon) %>%
                         group_by(ID) %>% 
                         dplyr::mutate(
                           Type = "random",
                           randist = randomDist(300, n()),  # Pass '300' and the length of the group.
                           randang = randomDistAngle(n())  # Generate random angles for each entry.
                         ), simplify = FALSE)

RP1 <- do.call(rbind, randomPts)
RP1 <- RP1 %>% 
  arrange(ID)
randomPts <- as.data.frame(RP1)

# Generate new latitudes and longitudes (the code above generates X and Y points).
set.seed(232409)
randomLatLon <- mapply(calculateLatLon, randomPts$Lat, randomPts$Lon, randomPts$randist, randomPts$randang)
randomPts$randLat <- randomLatLon[1, ]
randomPts$randLon <- randomLatLon[2, ]

# View the results
head(randomPts)


# =============================================================================================== #
# 10. CYCLE THROUGH EACH INDIVIDUAL AND CREATE A MAP                                               #
# =============================================================================================== #

set.seed(232409)

## --- 10.1. Function to generate a random distance constrained to 300m --- ##

randomDist <- function(max_distance, n) {
  # Generate `n` random values between 0 and `max_distance`
  runif(n, min = 0, max = max_distance)
}

set.seed(232409)

## --- 10.2. Function to calculate new lat-lon based on distance and angle --- ##

calculateLatLon <- function(lat, lon, dist, angle) {
  # Radius of Earth in meters
  R <- 6371000
  set.seed(232409)
  # Convert latitude and longitude to radians
  lat_rad <- lat * pi / 180
  lon_rad <- lon * pi / 180
  set.seed(232409)
  # Convert distance and angle to radians
  delta_lat <- dist / R * cos(angle * pi / 180)
  delta_lon <- dist / (R * cos(lat_rad)) * sin(angle * pi / 180)
  set.seed(232409)
  # Calculate new latitude and longitude
  new_lat <- lat_rad + delta_lat
  new_lon <- lon_rad + delta_lon
  set.seed(232409)
  # Convert radians back to degrees
  c(new_lat * 180 / pi, new_lon * 180 / pi)
}

## --- 7.3.Loop through each NOCA ID and plot the resulting points --- ##

set.seed(232409)

for (id in unique(noca$ID)) {
  # Subset data for the current ID
  temp <- subset(noca, noca$ID == id)
  set.seed(232409)
  # Debug: Check the structure of `temp`
  print(paste("Processing ID:", id, "with", nrow(temp), "rows"))
  set.seed(232409)
  # Convert observation points to an sf object
  obs <- st_as_sf(temp, coords = c("Lon", "Lat"), crs = "EPSG:4326")
  set.seed(232409)
  # Generate random distances and angles
  rand_distances <- randomDist(300, nrow(temp))  # Pass max_distance and number of rows
  rand_angles <- randomDistAngle(nrow(temp))    # Generate random angles
  set.seed(232409)
  # Debug: Check generated distances and angles
  print(head(rand_distances))
  print(head(rand_angles))
  set.seed(232409)
  # Add random distances and angles to `temp`
  rand_points <- temp %>%
    mutate(
      randist = rand_distances,
      randang = rand_angles
    )
  set.seed(232409)
  
  # Calculate new random latitudes and longitudes
  
  randomLatLon <- mapply(
    calculateLatLon,
    rand_points$Lat,
    rand_points$Lon,
    rand_points$randist,
    rand_points$randang
  )
  
  set.seed(232409)
  # Debug: Check randomLatLon
  print(head(randomLatLon))
  
  set.seed(232409)
  # Assign calculated coordinates
  rand_points$randLat <- randomLatLon[1, ]
  rand_points$randLon <- randomLatLon[2, ]
  
  set.seed(232409)
  # Convert random points to sf object
  rand_sf <- st_as_sf(rand_points, coords = c("randLon", "randLat"), crs = "EPSG:4326")
  set.seed(232409)
  # Plot observation and random points
  print(
    mapview(obs, col.regions = 'blue') + 
      mapview(rand_sf, col.regions = 'red')
  )
}

## --- 10.5. Save the random points in a dataframe --- ##
write.csv(randomPts, "03_processed_data/noca_2425_ran_points.csv")
