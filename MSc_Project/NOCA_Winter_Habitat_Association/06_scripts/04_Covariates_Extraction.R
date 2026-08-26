# =============================================================================================== #
# USED VS AVAILABLE HABITAT ANALYSIS FOR NORTHERN CARDINALS                                       #
# =============================================================================================== #
# Code adapted from Vanessa Poirier, MSc candidate, McGill University &
# Barbara Frei, ECCC (barbara.frei@ec.gc.ca)

# Last updated on March 20th 2025.

## --- TABLE OF CONTENT --- ##

# 1. LOADING FILES FOR TRUE AND RANDOM OBSERVATIONS, ALONG WITH RASTERS AND DESCRIPTION DATAFRAME

# 2. EXTRACTING TREE HEIGHT COVARIATES FOR THE YEAR 2023-2024

# 3. EXTRACTING COVARIATES FOR POPULATION DATA FOR 2023-2024

# 4. DATA PREPARATION FOR THE YEAR 2022-2023

# 5. EXTRACTING COVARIATES FOR TREE HEIGHT DATA FOR 2022-2023

# 6. EXTRACTING COVARIATES FOR POPULATION DATA FOR 2022-2023

# 7. DATA PREPARATION FOR THE YEAR 2024-2025

# 8. EXTRACTING COVARIATES FOR TREE HEIGHT DATA FOR 2024-2025

# 9. EXTRACTING COVARIATES FOR POPULATION DATA FOR 2024-2025

# 10. EXTRACTING ROADS LENGTHS COVARIATES FOR THE YEAR 2022-2023

# 11. EXTRACTING ROADS LENGTHS COVARIATES FOR THE YEAR 2023-2024

# 12. EXTRACTING ROADS LENGTHS COVARIATES FOR THE YEAR 2024-2025

# 13. EXTRACTING DISTANCE TO CLOSEST FEEDER FOR THE YEAR 2022-2023

# 14. EXTRACTING DISTANCE TO CLOSEST FEEDER FOR THE YEAR 2023-2024

# 15.  EXTRACTING DISTANCE TO CLOSEST FEEDER FOR THE YEAR 2024-2025

# 16. INVESTIGATING SITE DIFFERENCES WITH POPULATION DENSITY


## --- Packages Used --- ##

# library(rgdal) # rgdal isn't available anymore, the suggested replacement are the sf or terra
# packages.
# library(rgeos) -> not available anymore, integrated in the sf and terra packages.
# library(landscapetools) -> not available anymore, integrated in the sf and terra packages.

library(cancensus)
library(sf)
library(terra)
library(raster)
library(plyr)
library(ggplot2)
library(amt)
library(sp)
library(landscapemetrics)
library(purrr)
library(dplyr)
library(pdftools)
library(tidyr)
library(readr)
library(lubridate)
library(progressr)

rasterOptions(progress = "text")

# =============================================================================================== #
# 1. LOADING FILES FOR TRUE AND RANDOM OBSERVATIONS, ALONG WITH RASTERS AND DESCRIPTION DATAFRAME #
# =============================================================================================== #

## --- 1.1. RASTERS --- ##

# First, loading the raster file for the canopy height extracted from OpenForest.
canopy.ras.1 <- rast("01_downloaded_data/raster_canopee_2950_zone1.tif")
canopy.ras.2 <- rast("01_downloaded_data/raster_canopee_2950_zone2.tif")

canopy.ras <- merge(canopy.ras.1, canopy.ras.2)
canopy.ras <- raster(canopy.ras)

# Second, loading the population data from cancensus package.
# To use the data from CanCensus you need a personal API key. Enter it here after creating an
# account.
options(cancensus.api_key='CensusMapper_6206ab5a091356ca39a9dc940abfc253')

# PR = 24 for quebec province (diff numbers for diff provinces)
# DA is smallest census level <- I'll use this one.

# v_CA_21_1/2 - 2021, 2016 census data, if you just want the most recent census data, 
# erase "v_CA_21_2"; you can also just get total private dwellings (v_CA21_4), 
# pop density (v_CA21_6) and land area (v_CA21_7) as vectors.

# level CD = census division
# level = CSD = census subdivision
# level = CT = census tract
# level = DA = dissemination area

# Query Census Data for Montreal (Dissemination Area Level)
montreal_data <- get_census(
  dataset = "CA21",
  regions = list(CMA = "24462"),  # Montreal GeoUID
  vectors = c("v_CA21_6", "v_CA21_435", "v_CA21_436", "v_CA21_437", # All those extras ar for housing.
              "v_CA21_438", "v_CA21_439", "v_CA21_440", "v_CA21_441",
              "v_CA21_442"),  # Population density
  labels = "detailed",
  level = "DA",  # Use DA instead of DB
  geo_format = "sf"  # Get spatial data
)

# Define function to get census data and rasterize.
get_raster <- function(vector_code, field_name) {
  data <- get_census(
    dataset = "CA21",
    regions = list(CMA = "24462"),
    vectors = c(vector_code),
    level = "DA",
    geo_format = "sf"
  )
  
  template_raster <- raster(ncol = 1000, nrow = 1000, ext = extent(c(-74.2, -73.5, 45.4, 45.8)))
  points <- as(data, "Spatial")
  rasterized <- rasterize(points, template_raster, field = field_name)
  
  return(rasterized)
}

# List of variables to extract.
variables <- list(
  pop_density = list(code = "v_CA21_6", field = "v_CA21_6..Population.density.per.square.kilometre"),
  single_detached = list(code = "v_CA21_435", field = "v_CA21_435..Single.detached.house"),
  semi_detached = list(code = "v_CA21_436", field = "v_CA21_436..Semi.detached.house"),
  row_house = list(code = "v_CA21_437", field = "v_CA21_437..Row.house"),
  app_duplex = list(code = "v_CA21_438", field = "v_CA21_438..Apartment.or.flat.in.a.duplex"),
  app_less_5 = list(code = "v_CA21_439", field = "v_CA21_439..Apartment.in.a.building.that.has.fewer.than.five.storeys"),
  app_over_5 = list(code = "v_CA21_440", field = "v_CA21_440..Apartment.in.a.building.that.has.five.or.more.storeys")
)

# Apply function to all variables.
raster_results <- lapply(variables, function(v) get_raster(v$code, v$field))

# Assign individual raster objects.
pop_density_raster <- raster_results$pop_density
single_detached_raster <- raster_results$single_detached
semi_detached_raster <- raster_results$semi_detached
row_house_raster <- raster_results$row_house
app_duplex_raster <- raster_results$app_duplex
app_less_5_raster <- raster_results$app_less_5
app_over_5_raster <- raster_results$app_over_5

# Third, get streets and roads data.
# Load road polygons
roads <- st_read("01_downloaded_data/streets_22/VOI_VOIRIE_S_V22.shp")

## --- 1.2. NOCA POINTS & RANDOM POINTS --- ##

# Loading the csv file for the true observations, filtered in the project noca_2324_random_points_creation.R.
noca.obs2324.df <- read.csv("02_cleaned_data/noca_2324_cleaned.csv")

# Loading the csv file for the randomly generated points with a 1:1 ratio from noca_2324_random_points_creation.R.
noca.ran2324.df <- read.csv("03_processed_data/noca_2324_ran_points.csv")

## --- 1.3. DATA PREPARATION --- ##

# Filter the true points dataframe to only contain the randlat, randlon and ID.
noca.obs2324.df <- noca.obs2324.df %>% 
  dplyr::select(Date, Site, ID, Lat, Lon)
head(noca.obs2324.df)

# Assign a number for the winter of sampling to eventually compare years.
noca.obs2324.df <- noca.obs2324.df %>%
  mutate(Winter = 2)

# Filter the random points dataframe to only contain the randlat, randlon and ID.
noca.ran2324.df <- noca.ran2324.df %>% 
  dplyr::select(ID, Site, randLat, randLon)
head(noca.ran2324.df)

# In this dataframe (noca.ran2324.df), the column names for lat and lon are reversed.
# We'll switch it back to normal.
noca.ran2324.df <- noca.ran2324.df %>%
  mutate(temp = randLat,  # Store randLat in a temporary column
         randLat = randLon, 
         randLon = temp) %>%
  dplyr::select(-temp)  # Remove the temporary column
head(noca.ran2324.df)

# Assign a number for the winter of sampling to eventually compare years.
noca.ran2324.df <- noca.ran2324.df %>%
  mutate(Winter = 2)

## =============================================================================================== #
# 2.  EXTRACTING TREE HEIGHT COVARIATES FOR THE YEAR 2023-2024                                     #
# ================================================================================================ #

noca.obs2324.sf <- st_as_sf(noca.obs2324.df, coords = c("Lon", "Lat"), crs = 4326)
# Ensure the CRS matches the CRS of the raster before buffering.
noca.obs2324.sf <- st_transform(noca.obs2324.sf, crs = st_crs(canopy.ras))

# Ensure the CRS matches the CRS of the raster before buffering.
noca.ran2324.sf <- st_as_sf(noca.ran2324.df, coords = c("randLon", "randLat"), crs = 4326)
noca.ran2324.sf <- st_transform(noca.ran2324.sf, crs = st_crs(canopy.ras))

# Define buffer distances
buffers <- c(25, 50, 100, 200)

# Extract mean values for true observations.
obs_mean_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Mean_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.obs2324.sf, buffer = .x, fun = "mean")), col_name)
})

# Extract standard deviation values for random observations.
obs_sd_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Sd_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.obs2324.sf, buffer = .x, fun = "sd")), col_name)
})

# Add to original dataframe
noca.obs2324.df <- bind_cols(noca.obs2324.df, obs_mean_values, obs_sd_values)

# Extract mean values for random observations.
ran_mean_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Mean_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.ran2324.sf, buffer = .x, fun = "mean")), col_name)
})

# Extract standard deviation values for random observations.
ran_sd_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Sd_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.ran2324.sf, buffer = .x, fun = "sd")), col_name)
})

# Add to original dataframe
noca.ran2324.df <- bind_cols(noca.ran2324.df, ran_mean_values, ran_sd_values)

# =============================================================================================== #
# 3. EXTRACTING COVARIATES FOR POPULATION DATA FOR 2023-2024                                      #
# =============================================================================================== #

# Population data covariates extraction.
# Define buffers and rasters.
buffers <- c(25, 50, 100, 200)
rasters <- list(
  Pop_Dens = pop_density_raster,
  Single_House = single_detached_raster,
  Semi_Detached_House = semi_detached_raster,
  Row_House = row_house_raster,
  App_Duplex = app_duplex_raster,
  App_Less_5_Storeys = app_less_5_raster,
  App_Over_5_Storeys = app_over_5_raster
)

# Function to extract mean and SD for each raster and buffer.
extract_covariates <- function(df, sf_object, rasters, buffers, stat) {
  for (raster_name in names(rasters)) {
    raster_obj <- rasters[[raster_name]]
    for (buffer in buffers) {
      column_name <- paste0(stat, "_", raster_name, "_", buffer, "_m")
      df[[column_name]] <- raster::extract(raster_obj, sf_object, buffer = buffer, fun = stat)
    }
  }
  return(df)
}

# Extract mean values for observed and random observations.
# We do get some warnings, but those are just because the crs of the df is adjusted to match
# the crs of the rasters.
noca.obs2324.df <- extract_covariates(noca.obs2324.df, noca.obs2324.sf, rasters, buffers, 'mean')
noca.ran2324.df <- extract_covariates(noca.ran2324.df, noca.ran2324.sf, rasters, buffers, 'mean')

# Extract standard deviation values for observed and random observations.
# We do get some warnings, but those are just because the crs of the df is adjusted to match
# the crs of the rasters.
noca.obs2324.df <- extract_covariates(noca.obs2324.df, noca.obs2324.sf, rasters, buffers, 'sd')
noca.ran2324.df <- extract_covariates(noca.ran2324.df, noca.ran2324.sf, rasters, buffers, 'sd')

# Replace NAs with 0.
noca.obs2324.df <- noca.obs2324.df %>% mutate(across(everything(), ~ replace_na(.x, 0)))
noca.ran2324.df <- noca.ran2324.df %>% mutate(across(everything(), ~ replace_na(.x, 0)))

write.csv(noca.obs2324.df, "03_processed_data/noca_obs_df_mea_sd_2324.csv")
write.csv(noca.ran2324.df, "03_processed_data/noca_ran_df_mea_sd_2324.csv")

# ================================================================================================ #
# 4. DATA PREPARATION FOR THE YEAR 2022-2023                                                       #
# ================================================================================================ #

## --- 4.2. NOCA POINTS & RANDOM POINTS --- ##

# Loading the csv file for the true observations, filtered in the project noca_2324_random_points_creation.R.
noca.obs2223.df <- read.csv("02_cleaned_data/noca_2223_cleaned.csv")

# Loading the csv file for the randomly generated points with a 1:1 ratio from noca_2324_random_points_creation.R.
noca.ran2223.df <- read.csv("03_processed_data/noca_2223_ran_points.csv")

## --- 4.3. DATA PREPARATION --- ##

# Filter the true points dataframe to only contain the lat, lon and ID.
noca.obs2223.df <- noca.obs2223.df %>% 
  mutate(Site = 1)
noca.obs2223.df <- noca.obs2223.df %>% 
  dplyr::select(Date, Site, ID, Lat, Lon)
head(noca.obs2223.df)
noca.obs2223.df$Date <- as.Date(noca.obs2223.df$Date, format = "%Y-%m-%d") 

# Assign a number for the winter of sampling to eventually compare years.
noca.obs2223.df <- noca.obs2223.df %>%
  mutate(Winter = 1)

# Filter the random points dataframe to only contain the randlat, randlon and ID and site.
noca.ran2223.df <- noca.ran2223.df %>% 
  dplyr::select(ID, Site, randLat, randLon)
head(noca.ran2223.df)

# In this dataframe (noca.ran2223.df), the column names for lat and lon are reversed.
# We'll switch it back to normal.
noca.ran2223.df <- noca.ran2223.df %>%
  mutate(temp = randLat,  # Store randLat in a temporary column
         randLat = randLon, 
         randLon = temp) %>%
  dplyr::select(-temp)  # Remove the temporary column
head(noca.ran2223.df)

# Assign a number for the winter of sampling to eventually compare years.
noca.ran2223.df <- noca.ran2223.df %>%
  mutate(Winter = 1)

## =============================================================================================== #
# 5.  EXTRACTING COVARIATES FOR TREE HEIGHT DATA FOR 2022-2023                                     #
# ================================================================================================ #

noca.obs2223.sf <- st_as_sf(noca.obs2223.df, coords = c("Lon", "Lat"), crs = 4326)
# Ensure the CRS matches the CRS of the raster before buffering.
noca.obs2223.sf <- st_transform(noca.obs2223.sf, crs = st_crs(canopy.ras))

# Ensure the CRS matches the CRS of the raster before buffering.
noca.ran2223.sf <- st_as_sf(noca.ran2223.df, coords = c("randLon", "randLat"), crs = 4326)
noca.ran2223.sf <- st_transform(noca.ran2223.sf, crs = st_crs(canopy.ras))

# Define buffer distances
buffers <- c(25, 50, 100, 200)

# Extract mean values for true observations.
obs_mean_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Mean_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.obs2223.sf, buffer = .x, fun = "mean")), col_name)
})

# Extract standard deviation values for random observations.
obs_sd_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Sd_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.obs2223.sf, buffer = .x, fun = "sd")), col_name)
})

# Add to original dataframe
noca.obs2223.df <- bind_cols(noca.obs2223.df, obs_mean_values, obs_sd_values)

# Extract mean values for random observations.
ran_mean_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Mean_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.ran2223.sf, buffer = .x, fun = "mean")), col_name)
})

# Extract standard deviation values for random observations.
ran_sd_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Sd_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.ran2223.sf, buffer = .x, fun = "sd")), col_name)
})

# Add to original dataframe
noca.ran2223.df <- bind_cols(noca.ran2223.df, ran_mean_values, ran_sd_values)

# =============================================================================================== #
# 6. EXTRACTING COVARIATES FOR POPULATION DATA FOR 2022-2023                                     #
# =============================================================================================== #

# Covariates extraction.
# Define buffers and rasters.
buffers <- c(25, 50, 100, 200)
rasters <- list(
  Pop_Dens = pop_density_raster,
  Single_House = single_detached_raster,
  Semi_Detached_House = semi_detached_raster,
  Row_House = row_house_raster,
  App_Duplex = app_duplex_raster,
  App_Less_5_Storeys = app_less_5_raster,
  App_Over_5_Storeys = app_over_5_raster
)

# Function to extract mean and SD for each raster and buffer.
extract_covariates <- function(df, sf_object, rasters, buffers, stat) {
  for (raster_name in names(rasters)) {
    raster_obj <- rasters[[raster_name]]
    for (buffer in buffers) {
      column_name <- paste0(stat, "_", raster_name, "_", buffer, "_m")
      df[[column_name]] <- raster::extract(raster_obj, sf_object, buffer = buffer, fun = stat)
    }
  }
  return(df)
}

# Extract mean values for observed and random observations.
# We do get some warnings, but those are just because the crs of the df is adjusted to match
# the crs of the rasters.
noca.obs2223.df <- extract_covariates(noca.obs2223.df, noca.obs2223.sf, rasters, buffers, 'mean')
noca.ran2223.df <- extract_covariates(noca.ran2223.df, noca.ran2223.sf, rasters, buffers, 'mean')

# Extract standard deviation values for observed and random observations.
# We do get some warnings, but those are just because the crs of the df is adjusted to match
# the crs of the rasters.
noca.obs2223.df <- extract_covariates(noca.obs2223.df, noca.obs2223.sf, rasters, buffers, 'sd')
noca.ran2223.df <- extract_covariates(noca.ran2223.df, noca.ran2223.sf, rasters, buffers, 'sd')

# Replace NAs with 0.
noca.obs2223.df <- noca.obs2223.df %>% mutate(across(everything(), ~ replace_na(.x, 0)))
noca.ran2223.df <- noca.ran2223.df %>% mutate(across(everything(), ~ replace_na(.x, 0)))

write.csv(noca.obs2223.df, "03_processed_data/noca_obs_df_mea_sd_2223.csv")
write.csv(noca.ran2223.df, "03_processed_data/noca_ran_df_mea_sd_2223.csv")

# ================================================================================================ #
# 7. DATA PREPARATION FOR THE YEAR 2024-2025                                                       #
# ================================================================================================ #

## --- 7.1. NOCA POINTS & RANDOM POINTS --- ##

# Loading the csv file for the true observations, filtered in the project noca_2324_random_points_creation.R.
noca.obs2425.df <- read.csv("02_cleaned_data/noca_2425_cleaned.csv")

# Loading the csv file for the randomly generated points with a 1:1 ratio from noca_2324_random_points_creation.R.
noca.ran2425.df <- read.csv("03_processed_data/noca_2425_ran_points.csv")

## --- 7.2. DATA PREPARATION --- ##

# Filter the true points dataframe to only contain the lat, lon and ID.
noca.obs2425.df <- noca.obs2425.df %>% 
  dplyr::select(Date, Site, ID, Lat, Lon)
head(noca.obs2425.df)
noca.obs2425.df$Date <- as.Date(noca.obs2425.df$Date, format = "%Y-%m-%d") 

# Assign a number for the winter of sampling to eventually compare years.
noca.obs2425.df <- noca.obs2425.df %>%
  mutate(Winter = 3)

# Filter the random points dataframe to only contain the randlat, randlon and ID and site.
noca.ran2425.df <- noca.ran2425.df %>% 
  dplyr::select(Date, ID, Site, randLat, randLon)
head(noca.ran2425.df)

# In this dataframe (noca.ran2425.df), the column names for lat and lon are reversed.
# We'll switch it back to normal.
noca.ran2425.df <- noca.ran2425.df %>%
  mutate(temp = randLat,  # Store randLat in a temporary column
         randLat = randLon, 
         randLon = temp) %>%
  dplyr::select(-temp)  # Remove the temporary column
head(noca.ran2425.df)

# Assign a number for the winter of sampling to eventually compare years.
noca.ran2425.df <- noca.ran2425.df %>%
  mutate(Winter = 3)

## =============================================================================================== #
# 8.  EXTRACTING COVARIATES FOR TREE HEIGHT DATA FOR 2024-2025                                     #
# ================================================================================================ #

noca.obs2425.sf <- st_as_sf(noca.obs2425.df, coords = c("Lon", "Lat"), crs = 4326)
# Ensure the CRS matches the CRS of the raster before buffering.
noca.obs2425.sf <- st_transform(noca.obs2425.sf, crs = st_crs(canopy.ras))

# Ensure the CRS matches the CRS of the raster before buffering.
noca.ran2425.sf <- st_as_sf(noca.ran2425.df, coords = c("randLon", "randLat"), crs = 4326)
noca.ran2425.sf <- st_transform(noca.ran2425.sf, crs = st_crs(canopy.ras))

# Define buffer distances
buffers <- c(25, 50, 100, 200)

# Extract mean values for true observations.
obs_mean_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Mean_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.obs2425.sf, buffer = .x, fun = "mean")), col_name)
})

# Extract standard deviation values for random observations.
obs_sd_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Sd_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.obs2425.sf, buffer = .x, fun = "sd")), col_name)
})

# Add to original dataframe
noca.obs2425.df <- bind_cols(noca.obs2425.df, obs_mean_values, obs_sd_values)

# Extract mean values for random observations.
ran_mean_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Mean_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.ran2425.sf, buffer = .x, fun = "mean")), col_name)
})

# Extract standard deviation values for random observations.
ran_sd_values <- map_dfc(buffers, ~ {
  col_name <- paste0("Sd_Veg_Height_", .x, "_m")
  setNames(data.frame(raster::extract(canopy.ras, noca.ran2425.sf, buffer = .x, fun = "sd")), col_name)
})

# Add to original dataframe
noca.ran2425.df <- bind_cols(noca.ran2425.df, ran_mean_values, ran_sd_values)

# =============================================================================================== #
# 9. EXTRACTING COVARIATES FOR POPULATION DATA FOR 2024-2025                                     #
# =============================================================================================== #

# Covariates extraction.
# Define buffers and rasters.
buffers <- c(25, 50, 100, 200)
rasters <- list(
  Pop_Dens = pop_density_raster,
  Single_House = single_detached_raster,
  Semi_Detached_House = semi_detached_raster,
  Row_House = row_house_raster,
  App_Duplex = app_duplex_raster,
  App_Less_5_Storeys = app_less_5_raster,
  App_Over_5_Storeys = app_over_5_raster
)

# Function to extract mean and SD for each raster and buffer.
extract_covariates <- function(df, sf_object, rasters, buffers, stat) {
  for (raster_name in names(rasters)) {
    raster_obj <- rasters[[raster_name]]
    for (buffer in buffers) {
      column_name <- paste0(stat, "_", raster_name, "_", buffer, "_m")
      df[[column_name]] <- raster::extract(raster_obj, sf_object, buffer = buffer, fun = stat)
    }
  }
  return(df)
}

# Extract mean values for observed and random observations.
# We do get some warnings, but those are just because the crs of the df is adjusted to match
# the crs of the rasters.
noca.obs2425.df <- extract_covariates(noca.obs2425.df, noca.obs2425.sf, rasters, buffers, 'mean')
noca.ran2425.df <- extract_covariates(noca.ran2425.df, noca.ran2425.sf, rasters, buffers, 'mean')

# Extract standard deviation values for observed and random observations.
# We do get some warnings, but those are just because the crs of the df is adjusted to match
# the crs of the rasters.
noca.obs2425.df <- extract_covariates(noca.obs2425.df, noca.obs2425.sf, rasters, buffers, 'sd')
noca.ran2425.df <- extract_covariates(noca.ran2425.df, noca.ran2425.sf, rasters, buffers, 'sd')

# Replace NAs with 0.
noca.obs2425.df <- noca.obs2425.df %>% mutate(across(everything(), ~ replace_na(.x, 0)))
noca.ran2425.df <- noca.ran2425.df %>% mutate(across(everything(), ~ replace_na(.x, 0)))

write.csv(noca.obs2425.df, "03_processed_data/noca_obs_df_mea_sd_2425.csv")
write.csv(noca.ran2425.df, "03_processed_data/noca_ran_df_mea_sd_2425.csv")

## =============================================================================================== #
# 10.  EXTRACTING ROADS LENGTHS COVARIATES FOR THE YEAR 2022-2023                                  #
# ================================================================================================ #

## --- 10.1. True NOCA observations. --- ##

noca.obs2223.sf <- st_as_sf(noca.obs2223.df, coords = c("Lon", "Lat"), crs = 4326)
# Ensure the CRS matches the CRS of the roads polygon before buffering.
noca.obs2223.sf <- st_transform(noca.obs2223.sf, crs = st_crs(roads))

# Convert roads from polygon to MULTILINESTRING so that we are working with lengths.
roads_lines <- st_cast(roads, "MULTILINESTRING")

buffers <- c(25, 50, 100, 200)

# Create unique ID to track row positions
noca.obs2223.sf$id <- seq_len(nrow(noca.obs2223.sf))

# Total steps = number of points × number of buffer sizes
total_steps <- length(buffers) * nrow(noca.obs2223.sf)

# Start progress-aware operation
obs_road_lengths <- with_progress({
  p <- progressor(steps = total_steps)
  
  map_dfc(buffers, function(buffer_size) {
    col_name <- paste0("Road_Length_", buffer_size, "_m")  # New column name reflects length
    
    buffer_lengths <- map_dfr(1:nrow(noca.obs2223.sf), function(i) {
      p(message = paste("Point", i, "Buffer", buffer_size, "m"))  # Progress update
      
      point_buffer <- st_buffer(noca.obs2223.sf[i, ], dist = buffer_size)
      intersected <- tryCatch(
        st_intersection(point_buffer, roads_lines),  # <- updated from `roads` to `roads_lines`
        error = function(e) NULL
      )
      
      if (!is.null(intersected) && nrow(intersected) > 0) {
        road_length <- sum(as.numeric(st_length(intersected)))  # <- updated from `st_area()` to `st_length()`
      } else {
        road_length <- NA
      }
      
      tibble(!!col_name := road_length)
    })
    
    return(buffer_lengths)
  })
})

# Add to original dataframe
noca.obs2223.df <- bind_cols(noca.obs2223.df, obs_road_lengths)

# There are a bunch of NAs from areas without streets, replace them with 0.
noca.obs2223.df[names(obs_road_lengths)] <- obs_road_lengths %>%
  mutate(across(everything(), ~replace_na(., 0)))

## --- 10.2. Random NOCA observations. --- ##

noca.ran2223.sf <- st_as_sf(noca.ran2223.df, coords = c("randLon", "randLat"), crs = 4326)
# Ensure the CRS matches the CRS of the roads polygon before buffering.
noca.ran2223.sf <- st_transform(noca.ran2223.sf, crs = st_crs(roads))

# Convert roads from polygon to MULTILINESTRING so that we are working with lengths.
roads_lines <- st_cast(roads, "MULTILINESTRING")

buffers <- c(25, 50, 100, 200)

# Create unique ID to track row positions
noca.ran2223.sf$id <- seq_len(nrow(noca.ran2223.sf))

# Total steps = number of points × number of buffer sizes
total_steps <- length(buffers) * nrow(noca.ran2223.sf)

# Start progress-aware operation
ran_road_lengths <- with_progress({
  p <- progressor(steps = total_steps)
  
  map_dfc(buffers, function(buffer_size) {
    col_name <- paste0("Road_Length_", buffer_size, "_m")  # New column name reflects length
    
    buffer_lengths <- map_dfr(1:nrow(noca.ran2223.sf), function(i) {
      p(message = paste("Point", i, "Buffer", buffer_size, "m"))  # Progress update
      
      point_buffer <- st_buffer(noca.ran2223.sf[i, ], dist = buffer_size)
      intersected <- tryCatch(
        st_intersection(point_buffer, roads_lines),  # <- updated from `roads` to `roads_lines`
        error = function(e) NULL
      )
      
      if (!is.null(intersected) && nrow(intersected) > 0) {
        road_length <- sum(as.numeric(st_length(intersected)))  # <- updated from `st_area()` to `st_length()`
      } else {
        road_length <- NA
      }
      
      tibble(!!col_name := road_length)
    })
    
    return(buffer_lengths)
  })
})

# Add to original dataframe
noca.ran2223.df <- bind_cols(noca.ran2223.df, ran_road_lengths)

# There are a bunch of NAs from areas without streets, replace them with 0.
noca.ran2223.df[names(ran_road_lengths)] <- ran_road_lengths %>%
  mutate(across(everything(), ~replace_na(., 0)))

write.csv(noca.obs2223.df, "03_processed_data/noca_obs_df_mea_sd_2223.csv")
write.csv(noca.ran2223.df, "03_processed_data/noca_ran_df_mea_sd_2223.csv")

## =============================================================================================== #
# 11.  EXTRACTING ROADS LENGTHS COVARIATES FOR THE YEAR 2023-2024                                  #
# ================================================================================================ #

## --- 11.1. True NOCA observations. --- ##

noca.obs2324.sf <- st_as_sf(noca.obs2324.df, coords = c("Lon", "Lat"), crs = 4326)
# Ensure the CRS matches the CRS of the roads polygon before buffering.
noca.obs2324.sf <- st_transform(noca.obs2324.sf, crs = st_crs(roads))

# Convert roads from polygon to MULTILINESTRING so that we are working with lengths.
roads_lines <- st_cast(roads, "MULTILINESTRING")

buffers <- c(25, 50, 100, 200)

# Create unique ID to track row positions
noca.obs2324.sf$id <- seq_len(nrow(noca.obs2324.sf))

# Total steps = number of points × number of buffer sizes
total_steps <- length(buffers) * nrow(noca.obs2324.sf)

# Start progress-aware operation
obs_road_lengths <- with_progress({
  p <- progressor(steps = total_steps)
  
  map_dfc(buffers, function(buffer_size) {
    col_name <- paste0("Road_Length_", buffer_size, "_m")  # New column name reflects length
    
    buffer_lengths <- map_dfr(1:nrow(noca.obs2324.sf), function(i) {
      p(message = paste("Point", i, "Buffer", buffer_size, "m"))  # Progress update
      
      point_buffer <- st_buffer(noca.obs2324.sf[i, ], dist = buffer_size)
      intersected <- tryCatch(
        st_intersection(point_buffer, roads_lines),  # <- updated from `roads` to `roads_lines`
        error = function(e) NULL
      )
      
      if (!is.null(intersected) && nrow(intersected) > 0) {
        road_length <- sum(as.numeric(st_length(intersected)))  # <- updated from `st_area()` to `st_length()`
      } else {
        road_length <- NA
      }
      
      tibble(!!col_name := road_length)
    })
    
    return(buffer_lengths)
  })
})

# Add to original dataframe
noca.obs2324.df <- bind_cols(noca.obs2324.df, obs_road_lengths)

# There are a bunch of NAs from areas without streets, replace them with 0.
noca.obs2324.df[names(obs_road_lengths)] <- obs_road_lengths %>%
  mutate(across(everything(), ~replace_na(., 0)))

## --- 11.2. Random NOCA observations. --- ##

noca.ran2324.sf <- st_as_sf(noca.ran2324.df, coords = c("randLon", "randLat"), crs = 4326)
# Ensure the CRS matches the CRS of the roads polygon before buffering.
noca.ran2324.sf <- st_transform(noca.ran2324.sf, crs = st_crs(roads))

# Convert roads from polygon to MULTILINESTRING so that we are working with lengths.
roads_lines <- st_cast(roads, "MULTILINESTRING")

buffers <- c(25, 50, 100, 200)

# Create unique ID to track row positions
noca.ran2324.sf$id <- seq_len(nrow(noca.ran2324.sf))

# Total steps = number of points × number of buffer sizes
total_steps <- length(buffers) * nrow(noca.ran2324.sf)

# Start progress-aware operation
ran_road_lengths <- with_progress({
  p <- progressor(steps = total_steps)
  
  map_dfc(buffers, function(buffer_size) {
    col_name <- paste0("Road_Length_", buffer_size, "_m")  # New column name reflects length
    
    buffer_lengths <- map_dfr(1:nrow(noca.ran2324.sf), function(i) {
      p(message = paste("Point", i, "Buffer", buffer_size, "m"))  # Progress update
      
      point_buffer <- st_buffer(noca.ran2324.sf[i, ], dist = buffer_size)
      intersected <- tryCatch(
        st_intersection(point_buffer, roads_lines),  # <- updated from `roads` to `roads_lines`
        error = function(e) NULL
      )
      
      if (!is.null(intersected) && nrow(intersected) > 0) {
        road_length <- sum(as.numeric(st_length(intersected)))  # <- updated from `st_area()` to `st_length()`
      } else {
        road_length <- NA
      }
      
      tibble(!!col_name := road_length)
    })
    
    return(buffer_lengths)
  })
})

# Add to original dataframe
noca.ran2324.df <- bind_cols(noca.ran2324.df, ran_road_lengths)

# There are a bunch of NAs from areas without streets, replace them with 0.
noca.ran2324.df[names(ran_road_lengths)] <- ran_road_lengths %>%
  mutate(across(everything(), ~replace_na(., 0)))

write.csv(noca.obs2324.df, "03_processed_data/noca_obs_df_mea_sd_2324.csv")
write.csv(noca.ran2324.df, "03_processed_data/noca_ran_df_mea_sd_2324.csv")

## =============================================================================================== #
# 12.  EXTRACTING ROADS LENGTHS COVARIATES FOR THE YEAR 2024-2025                                  #
# ================================================================================================ #

## --- 12.1. True NOCA observations. --- ##

noca.obs2425.sf <- st_as_sf(noca.obs2425.df, coords = c("Lon", "Lat"), crs = 4326)
# Ensure the CRS matches the CRS of the roads polygon before buffering.
noca.obs2425.sf <- st_transform(noca.obs2425.sf, crs = st_crs(roads))

# Convert roads from polygon to MULTILINESTRING so that we are working with lengths.
roads_lines <- st_cast(roads, "MULTILINESTRING")

buffers <- c(25, 50, 100, 200)

# Create unique ID to track row positions
noca.obs2425.sf$id <- seq_len(nrow(noca.obs2425.sf))

# Total steps = number of points × number of buffer sizes
total_steps <- length(buffers) * nrow(noca.obs2425.sf)

# Start progress-aware operation
obs_road_lengths <- with_progress({
  p <- progressor(steps = total_steps)
  
  map_dfc(buffers, function(buffer_size) {
    col_name <- paste0("Road_Length_", buffer_size, "_m")  # New column name reflects length
    
    buffer_lengths <- map_dfr(1:nrow(noca.obs2425.sf), function(i) {
      p(message = paste("Point", i, "Buffer", buffer_size, "m"))  # Progress update
      
      point_buffer <- st_buffer(noca.obs2425.sf[i, ], dist = buffer_size)
      intersected <- tryCatch(
        st_intersection(point_buffer, roads_lines),  # <- updated from `roads` to `roads_lines`
        error = function(e) NULL
      )
      
      if (!is.null(intersected) && nrow(intersected) > 0) {
        road_length <- sum(as.numeric(st_length(intersected)))  # <- updated from `st_area()` to `st_length()`
      } else {
        road_length <- NA
      }
      
      tibble(!!col_name := road_length)
    })
    
    return(buffer_lengths)
  })
})

# Add to original dataframe
noca.obs2425.df <- bind_cols(noca.obs2425.df, obs_road_lengths)

# There are a bunch of NAs from areas without streets, replace them with 0.
noca.obs2425.df[names(obs_road_lengths)] <- obs_road_lengths %>%
  mutate(across(everything(), ~replace_na(., 0)))

## --- 12.2. Random NOCA observations. --- ##

noca.ran2425.sf <- st_as_sf(noca.ran2425.df, coords = c("randLon", "randLat"), crs = 4326)
# Ensure the CRS matches the CRS of the roads polygon before buffering.
noca.ran2425.sf <- st_transform(noca.ran2425.sf, crs = st_crs(roads))

# Convert roads from polygon to MULTILINESTRING so that we are working with lengths.
roads_lines <- st_cast(roads, "MULTILINESTRING")

buffers <- c(25, 50, 100, 200)

# Create unique ID to track row positions
noca.ran2425.sf$id <- seq_len(nrow(noca.ran2425.sf))

# Total steps = number of points × number of buffer sizes
total_steps <- length(buffers) * nrow(noca.ran2425.sf)

# Start progress-aware operation
ran_road_lengths <- with_progress({
  p <- progressor(steps = total_steps)
  
  map_dfc(buffers, function(buffer_size) {
    col_name <- paste0("Road_Length_", buffer_size, "_m")  # New column name reflects length
    
    buffer_lengths <- map_dfr(1:nrow(noca.ran2425.sf), function(i) {
      p(message = paste("Point", i, "Buffer", buffer_size, "m"))  # Progress update
      
      point_buffer <- st_buffer(noca.ran2425.sf[i, ], dist = buffer_size)
      intersected <- tryCatch(
        st_intersection(point_buffer, roads_lines),  # <- updated from `roads` to `roads_lines`
        error = function(e) NULL
      )
      
      if (!is.null(intersected) && nrow(intersected) > 0) {
        road_length <- sum(as.numeric(st_length(intersected)))  # <- updated from `st_area()` to `st_length()`
      } else {
        road_length <- NA
      }
      
      tibble(!!col_name := road_length)
    })
    
    return(buffer_lengths)
  })
})

# Add to original dataframe
noca.ran2425.df <- bind_cols(noca.ran2425.df, ran_road_lengths)

# There are a bunch of NAs from areas without streets, replace them with 0.
noca.ran2425.df[names(ran_road_lengths)] <- ran_road_lengths %>%
  mutate(across(everything(), ~replace_na(., 0)))

write.csv(noca.obs2425.df, "03_processed_data/noca_obs_df_mea_sd_2425.csv")
write.csv(noca.ran2425.df, "03_processed_data/noca_ran_df_mea_sd_2425.csv")

## =============================================================================================== #
# 13.  EXTRACTING DISTANCE TO CLOSEST FEEDER FOR THE YEAR 2022-2023                                #
# ================================================================================================ #

# Load feeders data.
feeders_df <- read.csv("00_raw_data/noca_feeders.csv")

# Transform observations, random points and feeders to sf objects.
noca.obs2223.sf <- st_as_sf(noca.obs2223.df, coords = c("Lon", "Lat"), crs = 4326)
noca.ran2223.sf <- st_as_sf(noca.ran2223.df, coords = c("randLon", "randLat"), crs = 4326)
feeders_sf <- st_as_sf(feeders_df, coords = c("Lon", "Lat"), crs = 4326)

# Reproject to a projected CRS for accurate distance (UTM zone zone 18N, EPSG 32618)
noca.obs2223.sf <- st_transform(noca.obs2223.sf, crs = 32618)
noca.ran2223.sf <- st_transform(noca.ran2223.sf, crs = 32618)
feeders_sf <- st_transform(feeders_sf, crs = 32618)

# Compute distance to the nearest feeder for true and random points
nearest_distances_obs <- st_distance(noca.obs2223.sf, feeders_sf) %>%
  apply(1, min)  # find the minimum distance for each point

nearest_distances_ran <- st_distance(noca.ran2223.sf, feeders_sf) %>%
  apply(1, min)

# Add to non-spatial data frame
noca.obs2223.df$Dist_to_Nearest_Feeder_m <- as.numeric(nearest_distances_obs)
noca.ran2223.df$Dist_to_Nearest_Feeder_m <- as.numeric(nearest_distances_ran)

## =============================================================================================== #
# 14.  EXTRACTING DISTANCE TO CLOSEST FEEDER FOR THE YEAR 2023-2024                                #
# ================================================================================================ #

# Load feeders data.
feeders_df <- read.csv("00_raw_data/noca_feeders.csv")

# Transform observations, random points and feeders to sf objects.
noca.obs2324.sf <- st_as_sf(noca.obs2324.df, coords = c("Lon", "Lat"), crs = 4326)
noca.ran2324.sf <- st_as_sf(noca.ran2324.df, coords = c("randLon", "randLat"), crs = 4326)
feeders_sf <- st_as_sf(feeders_df, coords = c("Lon", "Lat"), crs = 4326)

# Reproject to a projected CRS for accurate distance (UTM zone zone 18N, EPSG 32618)
noca.obs2324.sf <- st_transform(noca.obs2324.sf, crs = 32618)
noca.ran2324.sf <- st_transform(noca.ran2324.sf, crs = 32618)
feeders_sf <- st_transform(feeders_sf, crs = 32618)

# Compute distance to the nearest feeder for true and random points
nearest_distances_obs <- st_distance(noca.obs2324.sf, feeders_sf) %>%
  apply(1, min)  # find the minimum distance for each point

nearest_distances_ran <- st_distance(noca.ran2324.sf, feeders_sf) %>%
  apply(1, min)

# Add to non-spatial data frame
noca.obs2324.df$Dist_to_Nearest_Feeder_m <- as.numeric(nearest_distances_obs)
noca.ran2324.df$Dist_to_Nearest_Feeder_m <- as.numeric(nearest_distances_ran)

## =============================================================================================== #
# 15.  EXTRACTING DISTANCE TO CLOSEST FEEDER FOR THE YEAR 2024-2025                                #
# ================================================================================================ #

# Load feeders data.
feeders_df <- read.csv("00_raw_data/noca_feeders.csv")

# Transform observations, random points and feeders to sf objects.
noca.obs2425.sf <- st_as_sf(noca.obs2425.df, coords = c("Lon", "Lat"), crs = 4326)
noca.ran2425.sf <- st_as_sf(noca.ran2425.df, coords = c("randLon", "randLat"), crs = 4326)
feeders_sf <- st_as_sf(feeders_df, coords = c("Lon", "Lat"), crs = 4326)

# Reproject to a projected CRS for accurate distance (UTM zone zone 18N, EPSG 32618)
noca.obs2425.sf <- st_transform(noca.obs2425.sf, crs = 32618)
noca.ran2425.sf <- st_transform(noca.ran2425.sf, crs = 32618)
feeders_sf <- st_transform(feeders_sf, crs = 32618)

# Compute distance to the nearest feeder for true and random points
nearest_distances_obs <- st_distance(noca.obs2425.sf, feeders_sf) %>%
  apply(1, min)  # find the minimum distance for each point

nearest_distances_ran <- st_distance(noca.ran2425.sf, feeders_sf) %>%
  apply(1, min)

# Add to non-spatial data frame
noca.obs2425.df$Dist_to_Nearest_Feeder_m <- as.numeric(nearest_distances_obs)
noca.ran2425.df$Dist_to_Nearest_Feeder_m <- as.numeric(nearest_distances_ran)

write.csv(noca.obs2223.df, "03_processed_data/noca_obs_df_mea_sd_2223.csv")
write.csv(noca.ran2223.df, "03_processed_data/noca_ran_df_mea_sd_2223.csv")

write.csv(noca.obs2324.df, "03_processed_data/noca_obs_df_mea_sd_2324.csv")
write.csv(noca.ran2324.df, "03_processed_data/noca_ran_df_mea_sd_2324.csv")

write.csv(noca.obs2425.df, "03_processed_data/noca_obs_df_mea_sd_2425.csv")
write.csv(noca.ran2425.df, "03_processed_data/noca_ran_df_mea_sd_2425.csv")


## =============================================================================================== #
# 16. INVESTIGATING SITE DIFFERENCES WITH POPULATION DENSITY                                       #
# ================================================================================================ #

# Crop a bounding box for each site and extract cell values for all of them.
bbox_MBO <- extent(c(-73.947554, -73.922805, 45.421547, 45.435010))
pop_MBO <- crop(pop_density_raster, bbox_MBO)
vals_MBO <- getValues(pop_MBO)

bbox_BDU <- extent(c(-73.932519, -73.899103, 45.405677, 45.416159))
pop_BDU <- crop(pop_density_raster, bbox_BDU)
vals_BDU <- getValues(pop_BDU)

bbox_CON <- extent(c(-73.659980, -73.618136, 45.441673, 45.475709))
pop_CON <- crop(pop_density_raster, bbox_CON)
vals_CON <- getValues(pop_CON)

# Clean up and combine in one dataframe.
df <- data.frame(
  value = c(vals_MBO, vals_BDU, vals_CON),
  site = factor(c(
    rep("MBO", length(vals_MBO)),
    rep("BDU", length(vals_BDU)),
    rep("CON", length(vals_CON))
  ))
)

# Remove NAs.
df <- na.omit(df)

# Look at normality and variance homogeneity.
shapiro.test(df$value[df$site == "MBO"])  # check normality
leveneTest(value ~ site, data = df)       # check variance homogeneity
shapiro.test(df$value[df$site == "BDU"])  # check normality
leveneTest(value ~ site, data = df)       # check variance homogeneity
set.seed(123)  # So that the random sampling is always the same.
sample_CON <- sample(df$value[df$site == "CON"], 5000)
shapiro.test(sample_CON)
# All p-values are < 0.05 meaning that nothing is normal or homogenous. Moving on with a
# Kruskal test instead of an ANOVA.

kruskal.test(value ~ site, data = df)
pairwise.wilcox.test(df$value, df$site, p.adjust.method = "bonferroni")
# All p-values are < 0.05, meaning that all the sites are indeed different from one another on the
# basis of population density.

# Function to extract mean and standard error
get_pop_density_stats <- function(raster_layer, bbox) {
  cropped <- crop(raster_layer, bbox)
  values <- getValues(cropped)
  values <- values[!is.na(values)]  # remove NA values
  mean_val <- mean(values)
  se_val <- sd(values) / sqrt(length(values))  # standard error
  return(c(mean = mean_val, se = se_val))
}

# Define bounding boxes for each site
bbox_MBO <- extent(c(-73.947554, -73.922805, 45.421547, 45.435010))
bbox_BDU <- extent(c(-73.932519, -73.899103, 45.405677, 45.416159))
bbox_CON <- extent(c(-73.659980, -73.618136, 45.441673, 45.475709))

# Apply the function to each site
MBO_stats <- get_pop_density_stats(pop_density_raster, bbox_MBO)
BDU_stats <- get_pop_density_stats(pop_density_raster, bbox_BDU)
CON_stats <- get_pop_density_stats(pop_density_raster, bbox_CON)

# Combine results into a data frame
pop_density_summary <- data.frame(
  Site = c("MBO", "BDU", "CON"),
  Mean = c(MBO_stats["mean"], BDU_stats["mean"], CON_stats["mean"]),
  SE = c(MBO_stats["se"], BDU_stats["se"], CON_stats["se"])
)

print(pop_density_summary)