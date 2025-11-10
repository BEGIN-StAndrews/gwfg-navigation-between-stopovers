# ==============================================================================
# Title: Initial GPS preprocessing & aggregation for GWfG datasets
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Remove low-satellite fixes and duplicate timestamps per individual
#   2) Filter dense GPS bursts (keep first + gaps > 10 s)
#   3) Apply geographic bounds and exclude problematic IDs
#   4) Select common columns and aggregate cleaned files
#   5) Optionally filter by species column and location error
#   6) Remove speed outliers using {move} (threshold 45 m/s)
# Inputs:
#   - Raw CSVs (see vector in script): *.csv
# Outputs:
#   - cleaned_<original>.csv
#   - gwfg_cleaned_aggregated.csv
#   - gwfg_outlier_removed.csv
# Usage:
#   Rscript R/01_preprocess_gps_cleaning.R
# Notes:
#   - This repository accompanies a manuscript submitted to a peer-reviewed journal.
#   - Raw GPS data may be restricted; see README and Availability statement.
# ==============================================================================


# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)
library(sp)

# ------------------------------------------------------------------------------
# Define Custom Functions
# ------------------------------------------------------------------------------
remove_low_fixes <- function(data) {
  data %>% 
    filter(gps.satellite.count > 3 | is.na(gps.satellite.count)) %>% 
    group_by(individual.local.identifier) %>%  
    distinct(timestamp, .keep_all = TRUE) %>%  
    ungroup()
}


filter_gps_bursts <- function(data) {
  data %>%
    mutate(timestamp_dt = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
    arrange(individual.local.identifier, timestamp_dt) %>%  
    group_by(individual.local.identifier) %>%  
    mutate(
      time_diff = as.numeric(timestamp_dt - lag(timestamp_dt))
    ) %>%  
    filter(time_diff > 10 | is.na(time_diff)) %>%  
    ungroup() %>%                     
    dplyr::select(-timestamp_dt, -time_diff)
}


remove_inconsistent_locations <- function(data) {
  data %>%
    filter(between(location.long, -2, 120) & between(location.lat, 40, 82))
}


remove_all_problematic_ids <- function(data) {
  problematic_ids <- c(
    "474_KOL_F", "A25_KOL_F", "62_Wolodja_J", "20_Aleska_F", 
    "GWFG_2015_424", "85_KOL_F", "Lily_17", "16_Clementina_F", 
    "15_Kamilla_F", "22_Alissa_F", "21_Kelemen_M",
    "44_Adriana_F", "45_Adele_J", "47_Adam_J", "50_Tina_J", 
    "53_Eva_F", "56_Evita_J", "59_Wolka_F", "67_Jasminka_J", 
    "70_Janika_F", "83_Charlotta_J", "84_Charly_J", "48_Tineke_F", 
    "46_Adriaan_J", "82_Chris_M", 
    
    "Hector_3555_NK7_J","Heidi_3553_NK5_J","Heleen_3557_NK8_J",
    "Henk_4011_NK6_J",  "Hennie_3554_NK4_F","Jonathan_3920_BG7_J",
    "Josanne_3991_BG4_J", "Jouke_3556_BG5_J", "Jonas_4191_BG3_J",
    "Josephine_4192_BG6_J","GWFG_2015_442","64_Jakob_J", 
    "65_Jasper_J","GWFG_2018_462","428_KOL_F"
      ) %>% unique()
  
  data %>% filter(!individual.local.identifier %in% problematic_ids)
}





# ------------------------------------------------------------------------------
# Processing Pipeline
# ------------------------------------------------------------------------------
process_dataset <- function(dataset) {
  read.csv(dataset) %>%
    remove_all_problematic_ids() %>%
    remove_low_fixes() %>%
    filter_gps_bursts() %>%
    remove_inconsistent_locations()
}


# ------------------------------------------------------------------------------
# Main Execution - Process & Save Each Dataset
# ------------------------------------------------------------------------------
datasets <- c(
  "Disturbance of GWFG by IFV and IWWR 2017-2018.csv",
  "Disturbance of GWFG by IFV and IWWR .csv",
  "Geese_MPIAB_Kolguev2016.csv",
  "LifeTrack Geese IG-RAS MPIAB ICARUS 2.csv",
  "Geese_MPIAB_IFV_IWWR_Kolguev_HUN2018.csv"
)


lapply(datasets, function(dataset) {
  df <- process_dataset(dataset)
    df2 <- df %>% 
    dplyr::select(timestamp, location.long, location.lat, ground.speed, heading, height.above.msl, individual.local.identifier, individual.taxon.canonical.name, gps.fix.type, location.error.numerical)
  
  write.csv(df2, paste0("cleaned_", dataset), row.names = FALSE)
})


################################################################################

# ------------------------------------------------------------------------------
# Aggregation with Common Columns AFTER Processing
# ------------------------------------------------------------------------------
select_common_columns <- function(data_list) {
  common_cols <- Reduce(intersect, lapply(data_list, colnames))
  lapply(data_list, function(df) df[, common_cols, drop = FALSE])
}


processed_files <- list.files(pattern = "^cleaned_.*\\.csv$")
processed_data_list <- lapply(processed_files, read.csv)
processed_data_common <- select_common_columns(processed_data_list)

final_aggregated_data <- bind_rows(processed_data_common)

if ("individual.taxon.canonical.name" %in% colnames(final_aggregated_data)) {
  final_aggregated_data <- final_aggregated_data %>%
    filter(individual.taxon.canonical.name != "Anser fabalis" | is.na(individual.taxon.canonical.name))             
}


if ("location.error.numerical" %in% colnames(final_aggregated_data)) {   
  final_aggregated_data2 <- final_aggregated_data %>%
    filter(location.error.numerical <= 30 | is.na(location.error.numerical))
}

write.csv(final_aggregated_data2, "gwfg_cleaned_aggregated.csv", row.names = FALSE) 


################################################################################
# ------------------------------------------------------------------------------
# Step 1: Remove Outlier AFTER Aggregation
# ------------------------------------------------------------------------------
rm(list = ls())

library(move)
library(dplyr)      
library(lubridate)

# Reload the aggregated dataset.
filtered_data <- read.csv("gwfg_cleaned_aggregated.csv")
filtered_data$individual.local.identifier <- paste0("id_", filtered_data$individual.local.identifier)

track_data <- move(x = filtered_data$location.long, y = filtered_data$location.lat, 
                   time = as.POSIXct(filtered_data$timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"), 
                   proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"), 
                   data = filtered_data, animal = filtered_data$individual.local.identifier)

track_data$Calspeed <- unlist(lapply(speed(track_data), c, NA))
track_df <- as.data.frame(track_data)

track_df <- track_df %>%
  filter(Calspeed < 45)

track_df$individual.local.identifier <- track_df$trackId

track_df2 <- track_df %>%
  dplyr::select(timestamp, location.long, location.lat, ground.speed, heading, height.above.msl, individual.local.identifier)



write.csv(track_df2, "gwfg_outlier_removed.csv", row.names = FALSE)



