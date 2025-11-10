# ==============================================================================
# Title: Detecting long/short stopovers & pruning stopover points
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Detect long stopovers (>= 48h, radius <= 30 km)
#   2) Detect short stopovers (>= 4h, radius <= 10 km)
#   3) Remove short stopovers overlapping long stopovers
#   4) Remove points within long stopovers (keep first & last point)
# Inputs:
#   - gwfg_outlier_removed.csv
# Outputs:
#   - stopovers_long.csv
#   - stopovers_short_raw.csv
#   - stopovers_short_filtered.csv
#   - gwfg_no_stops.csv
# Usage:
#   Rscript R/02_detect_stopovers.R

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
# Long stopovers
stopover_params_long <- list(
  min_time  = 48 * 3600,  # 48 hours (seconds)
  max_radius = 30000,     # 30 km (meters)
  coverage  = 0.95,
  step_size = 1
)

# Short stopovers 
stopover_params_short <- list(
  min_time  = 4 * 3600,   # 4 hours (seconds)
  max_radius = 10000,     # 10 km (meters)
  coverage  = 0.95,
  step_size = 1
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
robustCircle95 <- function(lon, lat, coverage) {
  n <- length(lon)
  if (n == 0) return(NULL)
  if (n == 1) {
    return(list(
      centerLon = lon[1],
      centerLat = lat[1],
      radius = 0,
      inlierIndex = 1
    ))
  }
  
  center0_lon <- median(lon)
  center0_lat <- median(lat)
  
  dists <- distHaversine(
    cbind(center0_lon, center0_lat),
    cbind(lon, lat)
  )
  
  idx_sorted <- order(dists)
  keep_n <- max(1, floor(coverage * n))
  inliers_local <- idx_sorted[seq_len(keep_n)]
  
  centerLon <- median(lon[inliers_local])
  centerLat <- median(lat[inliers_local])
  
  final_dists <- distHaversine(
    cbind(centerLon, centerLat),
    cbind(lon[inliers_local], lat[inliers_local])
  )
  
  list(
    centerLon = centerLon,
    centerLat = centerLat,
    radius = max(final_dists),
    inlierIndex = inliers_local
  )
}


isValidStop <- function(i, j, lon, lat, timestamp_num, min_time, max_radius, coverage) {
  if (j <= i || j > length(lon)) return(FALSE)
  if ((timestamp_num[j] - timestamp_num[i]) < min_time) return(FALSE)
  
  sub_lon <- lon[i:j]
  sub_lat <- lat[i:j]
  
  circle <- robustCircle95(sub_lon, sub_lat, coverage)
  if (is.null(circle) || circle$radius > max_radius) return(FALSE)
  
  last_dist <- distHaversine(
    cbind(circle$centerLon, circle$centerLat),
    cbind(lon[j], lat[j])
  )
  
  last_dist <= max_radius
}


detect_stopovers <- function(lon, lat, timestamp, min_time, max_radius, coverage, step_size) {
  n <- length(lon)
  if (n == 0) return(data.frame())
  timestamp_num <- as.numeric(timestamp)
  
  stopovers <- list()
  stop_count <- 0
  i <- 1
  
  while (i < n) {
    jLow <- i + 1
    while (jLow <= n && (timestamp_num[jLow] - timestamp_num[i]) < min_time) {
      jLow <- jLow + step_size
    }
    if (jLow > n) break
    
    if (distHaversine(cbind(lon[i], lat[i]), cbind(lon[jLow], lat[jLow])) > 2 * max_radius) {
      i <- i + 1
      next
    }
    
    validJ <- -1
    lo <- jLow
    hi <- n
    while (lo <= hi) {
      mid <- (lo + hi) %/% 2
      if (isValidStop(i, mid, lon, lat, timestamp_num, min_time, max_radius, coverage)) {
        validJ <- mid
        lo <- mid + 1
      } else {
        hi <- mid - 1
      }
    }
    
    if (validJ > i) {
      sub_lon <- lon[i:validJ]
      sub_lat <- lat[i:validJ]
      circle <- robustCircle95(sub_lon, sub_lat, coverage)
      
      global_idx <- i - 1 + circle$inlierIndex
      start_idx <- min(global_idx)
      
      stop_count <- stop_count + 1
      stopovers[[stop_count]] <- data.frame(
        arrival = timestamp[start_idx],
        departure = timestamp[validJ],
        centerLon = circle$centerLon,
        centerLat = circle$centerLat,
        duration = as.numeric(difftime(timestamp[validJ], timestamp[start_idx], units = "days")),
        
        radius = circle$radius
      )
      i <- validJ + 1  
    } else {
      i <- i + 1
    }
  }
  
  if (stop_count == 0) return(data.frame())
  do.call(rbind, stopovers)
}


# ------------------------------------------------------------------------------
# Pipeline A: Detect long stopovers
# ------------------------------------------------------------------------------
tracking_data <- read.csv("gwfg_outlier_removed.csv") %>% 
  mutate(
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  ) %>% 
  arrange(individual.local.identifier, timestamp)

results <- tracking_data %>% 
  group_by(individual.local.identifier) %>% 
  group_modify(~ {
    cat("Processing", .y$individual.local.identifier, "\n")
    detect_stopovers(
      .x$location.long, 
      .x$location.lat,
      .x$timestamp,
      stopover_params_long$min_time,
      stopover_params_long$max_radius,
      stopover_params_long$coverage,
      stopover_params_long$step_size
    )
  })

write.csv(results, "stopovers_long.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Pipeline B: Detect short stopovers
# ------------------------------------------------------------------------------
tracking_data <- read.csv("gwfg_outlier_removed.csv") %>% 
  mutate(
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  ) %>% 
  arrange(individual.local.identifier, timestamp)

results2 <- tracking_data %>% 
  group_by(individual.local.identifier) %>% 
  group_modify(~ {
    cat("Processing", .y$individual.local.identifier, "\n")
    detect_stopovers(
      .x$location.long, 
      .x$location.lat,
      .x$timestamp,
      stopover_params_short$min_time,
      stopover_params_short$max_radius,
      stopover_params_short$coverage,
      stopover_params_short$step_size
    )
  })

write.csv(results2, "stopovers_short_raw.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Pipeline C: Remove short stopovers overlapping long stopovers
# ------------------------------------------------------------------------------
long_stops <- read.csv("stopovers_long.csv") %>%
  mutate(
    arrival = as.POSIXct(arrival, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    departure = as.POSIXct(departure, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

short_stops <- read.csv("stopovers_short_raw.csv") %>%
  mutate(
    arrival = as.POSIXct(arrival, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    departure = as.POSIXct(departure, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

is_overlapping <- function(short_arrival, short_departure, long_arrival, long_departure) {
  return((short_arrival < long_departure) & (short_departure > long_arrival))
}

filtered_short_stops <- short_stops %>%
  group_by(individual.local.identifier) %>%
  group_modify(~ {
    current_id <- .y$individual.local.identifier
    
    individual_long_stops <- long_stops %>% 
      filter(individual.local.identifier == current_id)
    
    if (nrow(individual_long_stops) == 0) return(.x)
    
    .x %>%
      rowwise() %>%
      mutate(overlap = any(is_overlapping(arrival, departure,
                                          individual_long_stops$arrival,
                                          individual_long_stops$departure))) %>%
      ungroup() %>%
      filter(!overlap) %>%
      dplyr::select(-overlap)
  }) %>%
  ungroup()

write.csv(filtered_short_stops, "stopovers_short_filtered.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Pipeline D: Remove points within long stopovers (keep first & last)
# ------------------------------------------------------------------------------
tracking_data <- read.csv("gwfg_outlier_removed.csv", stringsAsFactors = FALSE) %>% 
  mutate(
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  ) %>%   arrange(individual.local.identifier, timestamp)


stopover_data <- read.csv("stopovers_long.csv", stringsAsFactors = FALSE)
stopover_data$arrival <- as.POSIXct(stopover_data$arrival, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
stopover_data$departure <- as.POSIXct(stopover_data$departure, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

tracking_data <- tracking_data %>%
  mutate(
    to_remove = 1,  # 1 = keep by default, 0 = remove
    stopover_status = NA_character_
  )


for (i in seq_len(nrow(stopover_data))) {
    individual_id <- stopover_data$individual.local.identifier[i]
  arrival_time <- stopover_data$arrival[i]
  departure_time <- stopover_data$departure[i]
  
  idx_stop <- which(tracking_data$individual.local.identifier == individual_id & 
                      tracking_data$timestamp >= arrival_time & 
                      tracking_data$timestamp <= departure_time)
  
  if (length(idx_stop) == 0) {
    next  
  }
  
  stop_points <- tracking_data[idx_stop, ] %>%
    arrange(timestamp)
  
  n_stops <- nrow(stop_points)
  
  if (n_stops > 2) {
    tracking_data$to_remove[idx_stop[1]] <- 1  # Keep first point
    tracking_data$stopover_status[idx_stop[1]] <- "entering stopover"
    
    tracking_data$to_remove[idx_stop[n_stops]] <- 1  # Keep last point
    tracking_data$stopover_status[idx_stop[n_stops]] <- "exiting stopover"
    
    if (n_stops > 2) {
      unique_idx <- unique(idx_stop[2:(n_stops - 1)])
      tracking_data$to_remove[unique_idx] <- 0
    }
    
  } else {
    tracking_data$to_remove[idx_stop] <- 1
    
    if (n_stops == 2) {
      tracking_data$stopover_status[idx_stop[1]] <- "entering stopover"
      tracking_data$stopover_status[idx_stop[2]] <- "exiting stopover"
    } else {
      tracking_data$stopover_status[idx_stop] <- "short stopover"
    }
  }
}

cleaned_data <- tracking_data[tracking_data$to_remove == 1, ]
cleaned_data$to_remove <- NULL  

write.csv(cleaned_data, "gwfg_no_stops.csv", row.names = FALSE)

