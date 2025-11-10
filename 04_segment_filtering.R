# ==============================================================================
# Title: Segment filtering
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Compute per-segment metrics (distance, duration, speeds, gaps)
#   2) Keep segments with great-circle ≥150 km and acceptable time gaps (≤2 h logic)
#   3) Export the filtered segment table and corresponding point data
# Inputs:  points_all.csv
# Outputs: segments_migratory.csv, points_migratory.csv
# Usage:   Rscript R/04_segment_filtering.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(lubridate)
library(geosphere)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
TIME_GAP_THRESHOLD_MIN <- 2*60   # 2 hours in minutes
DISTANCE_THRESHOLD_KM  <- 30     # For identifying large displacement
MIN_DISTANCE           <- 150    # Minimum great-circle distance for Step 2

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------
compute_segment_metrics <- function(df) {
  n_full <- nrow(df)
  
  if (n_full < 4) {
    return(data.frame(
      NumberOFpoints = 0,
      cumulative_distance = 0,
      great_circle_distance = 0,
      median_speed = NA_real_,
      median_time_interval = NA_real_,
      max_time_interval = NA_real_,
      time_gap_for_long_displacement = NA_real_,
      duration = 0,
      first_timestamp = NA, first_long = NA, first_lat = NA,
      last_timestamp  = NA, last_long  = NA, last_lat  = NA
    ))
  }
  
  df_sub <- df[2:(n_full - 1), ] %>%
    arrange(timestamp) %>%           
    mutate(row_idx = row_number())   
  
  n_sub  <- nrow(df_sub)
  
    metrics <- list(
    NumberOFpoints = n_sub,
    cumulative_distance = 0,                 
    great_circle_distance = 0,               
    median_speed = NA_real_,                 
    median_time_interval = NA_real_,         
    max_time_interval = NA_real_,            
    time_gap_for_long_displacement = NA_real_,
    duration = as.numeric(difftime(df_sub$timestamp[n_sub], df_sub$timestamp[1], units = "hours")),  # overall hours from raw data
    first_timestamp = df_sub$timestamp[1],
    first_long      = df_sub$location.long[1],
    first_lat       = df_sub$location.lat[1],
    last_timestamp  = df_sub$timestamp[n_sub],
    last_long       = df_sub$location.long[n_sub],
    last_lat        = df_sub$location.lat[n_sub]
  )
  
  if (n_sub >= 2) {
    next_pts <- lead(df_sub[, c("location.long", "location.lat")])
    dist_m   <- distHaversine(
      as.matrix(df_sub[, c("location.long", "location.lat")]),
      as.matrix(next_pts)
    )
    
    dist_m <- dist_m[!is.na(dist_m)]
    dist_km <- dist_m / 1000
    
    time_diffs_sec <- as.numeric(difftime(lead(df_sub$timestamp),df_sub$timestamp,units = "secs"))
    time_diffs_sec <- time_diffs_sec[!is.na(time_diffs_sec)]
    
    metrics$cumulative_distance <- sum(dist_km, na.rm = TRUE)
    metrics$great_circle_distance <- distHaversine(
      c(df_sub$location.long[1], df_sub$location.lat[1]),
      c(df_sub$location.long[n_sub], df_sub$location.lat[n_sub])
  ) / 1000
    
    time_diffs_min <- time_diffs_sec / 60
    
    if (length(time_diffs_min) > 0) {
      
      metrics$median_time_interval <- median(time_diffs_min, na.rm = TRUE)
      metrics$max_time_interval    <- max(time_diffs_min,    na.rm = TRUE)
      speeds_m_s <- dist_m / time_diffs_sec
      metrics$median_speed <- median(speeds_m_s, na.rm = TRUE)
      
      long_Time_gap_idx <- which(time_diffs_min > TIME_GAP_THRESHOLD_MIN )
      
      if (length(long_Time_gap_idx)> 0){
        
        long_Time_dist_gap_idx <- which(time_diffs_min > TIME_GAP_THRESHOLD_MIN & dist_km > DISTANCE_THRESHOLD_KM)
        
        if (length(long_Time_dist_gap_idx)< 1){
          
          metrics$time_gap_for_long_displacement<- -1 
          
        }else {
          metrics$time_gap_for_long_displacement<- max(time_diffs_min[long_Time_dist_gap_idx])
        }
      }
    }
  }
  
  out <- as.data.frame(metrics)
  out$cumulative_distance   <- round(out$cumulative_distance, 2)
  out$great_circle_distance <- round(out$great_circle_distance, 2)
  out$median_speed          <- round(out$median_speed, 6)
  out$median_time_interval  <- round(out$median_time_interval, 2)
  out$max_time_interval     <- round(out$max_time_interval, 2)
  out$time_gap_for_long_displacement   <- round(out$time_gap_for_long_displacement, 2)
  out$duration              <- round(out$duration, 2)
  
  return(out)
}

filter_valid_migratory_segments <- function(df,distance_cutoff = MIN_DISTANCE,time_gap_cutoff = TIME_GAP_THRESHOLD_MIN) {
  
  
  cat("=== Starting Data ===\n")
  cat("Number of segments at start:", nrow(df), "\n\n")

  n_before_1 <- nrow(df)
  df_removed_step1 <- df %>% filter(NumberOFpoints == 0)
  
  df_step1_passed <- df %>% filter(NumberOFpoints > 0)
  
  n_after_1 <- nrow(df_step1_passed)
  cat("Step 1: Remove NumberOFpoints == 0\n")
  cat("   Before:", n_before_1, "rows\n")
  cat("   After: ", n_after_1, "rows\n")
  cat("   Lost:  ", (n_before_1 - n_after_1),
      sprintf(" (%.1f%%)\n\n", 
              100 * (n_before_1 - n_after_1) / n_before_1))
  
  n_before_2 <- n_after_1
  
  df_removed_step2 <- df_step1_passed %>% 
    filter(great_circle_distance < distance_cutoff)
  
  df_step2_passed <- df_step1_passed %>%
    filter(great_circle_distance >= distance_cutoff)
  
  n_after_2 <- nrow(df_step2_passed)
  cat(sprintf("Step 2: Remove great_circle_distance < %s\n", distance_cutoff))
  cat("   Before:", n_before_2, "rows\n")
  cat("   After: ", n_after_2, "rows\n")
  cat("   Lost:  ", (n_before_2 - n_after_2),
      sprintf(" (%.1f%%)\n\n", 
              100 * (n_before_2 - n_after_2) / n_before_2))
  
  n_before_3 <- n_after_2
  
  df_step3_logic <- df_step2_passed %>%
    mutate(keep_it = case_when(
      max_time_interval <= time_gap_cutoff ~ TRUE,
      max_time_interval > time_gap_cutoff & time_gap_for_long_displacement == -1 ~ TRUE,
      max_time_interval > time_gap_cutoff & is.na(time_gap_for_long_displacement) ~ FALSE,
      max_time_interval > time_gap_cutoff &
        !is.na(time_gap_for_long_displacement) &
        time_gap_for_long_displacement != -1 ~ 
        (time_gap_for_long_displacement < time_gap_cutoff),
          TRUE ~ FALSE
    ))
  
  df_removed_step3 <- df_step3_logic %>%
    filter(keep_it == FALSE)
  
  df_step3_passed <- df_step3_logic %>%
    filter(keep_it == TRUE) %>%
    dplyr::select(-keep_it)
  
  n_after_3 <- nrow(df_step3_passed)
  cat("Step 3: Time-gap criteria\n")
  cat("   Before:", n_before_3, "rows\n")
  cat("   After: ", n_after_3, "rows\n")
  cat("   Lost:  ", (n_before_3 - n_after_3),
      sprintf(" (%.1f%%)\n\n", 
              100 * (n_before_3 - n_after_3) / n_before_3))
  
  cat("=== Done! ===\n")
  
  df_removed_step1  <- df_removed_step1  %>% mutate(removed_step = 1)
  df_removed_step2  <- df_removed_step2  %>% mutate(removed_step = 2)
  df_removed_step3  <- df_removed_step3  %>% mutate(removed_step = 3)
  
  df_non_migratory <- bind_rows(df_removed_step1,
                                df_removed_step2,
                                df_removed_step3)
  
  df_migratory     <- df_step3_passed 
  
  return(list(
    migratory     = df_migratory,
    non_migratory = df_non_migratory
  ))
}

# ------------------------------------------------------------------------------
# Pipeline: compute metrics per segment, select migratory, export
# ------------------------------------------------------------------------------
data <- read.csv("points_all.csv", stringsAsFactors = FALSE) %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  arrange(individual.local.identifier, id_year, unique_seg, timestamp)

season_metrics <- data %>%
  group_by(individual.local.identifier, id_year,season,segment_id, unique_seg) %>%
  group_modify(~ compute_segment_metrics(.x)) %>%
  ungroup()

season_metrics <- season_metrics %>% dplyr::select(individual.local.identifier, everything())
result  <- filter_valid_migratory_segments(season_metrics, distance_cutoff = 150, time_gap_cutoff = 120)
write.csv(result$migratory,     "segments_migratory.csv",     row.names = FALSE)

# ------------------------------------------------------------------------------
# Post-processing: export points for migratory segments only
# ------------------------------------------------------------------------------
migratory_ids <- read.csv("segments_migratory.csv", stringsAsFactors = FALSE)
point_data <- read.csv("points_all.csv", stringsAsFactors = FALSE)
point_data_migratory <- point_data %>%
  semi_join(migratory_ids, by = c("unique_seg"))
migratory_cleaned <- point_data_migratory %>% filter(is.na(stopover_status))
migratory_cleaned$stopover_status <- NULL

write.csv(migratory_cleaned, "points_migratory.csv", row.names = FALSE)


