# ==============================================================================
# Title: Flight-segment delineation between consecutive stopovers
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Define migration windows: spring (1 Mar–31 May), autumn (15 Aug–15 Nov)
#   2) For each id_year&season: identify origin and destination stopovers and label
#      intermediates; assign per-point segment_id and segment_type
#   3) Attach previous/next stopover timestamps, centres, and durations to points
#   4) Flag segments that intersect any short stopover
# Inputs:
#   - gwfg_no_stops.csv
#   - stopovers_long.csv
#   - stopovers_short_filtered.csv
# Outputs:
#   - points_all.csv
# Usage:
#   Rscript R/03_flight_segment_delineation.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(lubridate)


# ------------------------------------------------------------------------------
# 1) Load data
# ------------------------------------------------------------------------------
tracking_data <- read.csv("gwfg_no_stops.csv", stringsAsFactors = FALSE)
stopover_data <- read.csv("stopovers_long.csv", stringsAsFactors = FALSE)

tracking_data <- tracking_data %>%
  mutate(
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    year      = year(timestamp),
    id_year   = paste(individual.local.identifier, year, sep = "_")
  ) %>%  arrange(individual.local.identifier, timestamp)

stopover_data$arrival <- as.POSIXct(   stopover_data$arrival, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
stopover_data$departure <- as.POSIXct( stopover_data$departure, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# ------------------------------------------------------------------------------
# 2) Define migration periods (spring/autumn)
# ------------------------------------------------------------------------------
migration_periods <- list(
  spring = list(start = "-03-01", end = "-05-31"),
  autumn = list(start = "-08-15", end = "-11-15")
)

# ------------------------------------------------------------------------------
# 3) process_migration(): lookup of previous/next stopovers
# ------------------------------------------------------------------------------
process_migration <- function(tracking_data, stopover_data, current_id_year, 
                              period_name, period_dates) {
  
 tracking_data_period <- tracking_data %>%
    filter(id_year == current_id_year)
  
  if (nrow(tracking_data_period) == 0) {
    warning(paste("Warning:", current_id_year, 
                  "has no tracking data rows; returning NULL."))
    return(NULL)
  }
  
  year       <- unique(tracking_data_period$year)
  individual <- unique(tracking_data_period$individual.local.identifier)
  
  stopover_data_period <- stopover_data %>%
    filter(
      individual.local.identifier == individual,
      year(arrival)   <= year,
      year(departure) >= year
    )
  
  period_start <- as.POSIXct(paste0(year, period_dates$start), tz = "UTC")
  period_end   <- as.POSIXct(paste0(year, period_dates$end),   tz = "UTC")
  
  stopover_data_period <- stopover_data_period %>%
    filter(
      (departure >= period_start & departure <= period_end) |
        (arrival   >= period_start & arrival   <= period_end)
    )
  
  if (nrow(stopover_data_period) == 0) {
    warning(paste("Warning: No stopovers for", period_name, 
                  "in id_year", current_id_year))
    return(NULL)
  }
  
  origin_stopover <- stopover_data_period %>%
    filter(departure >= period_start) %>%
    slice_min(order_by = departure, n = 1)
  
  destination_stopover <- stopover_data_period %>%
    filter(arrival <= period_end) %>%
    slice_max(order_by = arrival, n = 1)
  
  if (nrow(origin_stopover) == 0 || nrow(destination_stopover) == 0) {
    warning(paste("Warning: Could not identify origin or destination stopover for",
                  period_name, "in id_year", current_id_year))
    return(NULL)
  }
  
  stopover_data_period <- stopover_data_period %>%
    mutate(
      stopover_type = case_when(
        departure == origin_stopover$departure    ~ "origin",
        arrival   == destination_stopover$arrival ~ "destination",
        TRUE                                      ~ "intermediate"
      )
    )
  
  tracking_data_period <- tracking_data_period %>%
    filter(timestamp >= period_start & timestamp <= period_end) %>%
    mutate(
      segment_type                  = NA_character_,
      previous_stopover_depar_time        = NA_character_,
      previous_stopover_center_lon  = NA_real_,
      previous_stopover_center_lat  = NA_real_,
      previous_stopover_duration    = NA_real_,
      next_stopover_arriv_time            = NA_character_,
      next_stopover_center_lon      = NA_real_,
      next_stopover_center_lat      = NA_real_,
      next_stopover_duration        = NA_real_,
      segment_id                    = NA_integer_
      )
  
  if (nrow(tracking_data_period) == 0) {
    return(NULL)
  }
  
 stopover_data_period <- stopover_data_period %>%
    mutate(original_row = row_number())  # track original row inside this subset
  
  stopover_dep_sorted <- stopover_data_period %>%
    arrange(departure, original_row)
  dep_vec <- stopover_dep_sorted$departure
  
  stopover_arr_sorted <- stopover_data_period %>%
    arrange(arrival, original_row)
  arr_vec <- stopover_arr_sorted$arrival
  

  get_closest_previous <- function(current_time) {
    idx <- findInterval(current_time, dep_vec)
    if (idx == 0) {
      return(stopover_dep_sorted[0, ])
    }

    chosen_dep_time <- dep_vec[idx]
    matches <- which(stopover_dep_sorted$departure == chosen_dep_time)
    row_index <- matches[length(matches)]
    return(stopover_dep_sorted[row_index, ])
  }
  
  get_closest_next_by_departure <- function(current_time) {
    idx <- findInterval(current_time, dep_vec)
    candidate <- idx + 1
    if (candidate > length(dep_vec)) {
      return(stopover_dep_sorted[0, ])
    }
    chosen_dep_time <- dep_vec[candidate]
    matches <- which(stopover_dep_sorted$departure == chosen_dep_time)
    row_index <- matches[1]  # the first row among those with that dep
    return(stopover_dep_sorted[row_index, ])
  }
  
    get_closest_next_by_arrival <- function(current_time) {
    idx <- findInterval(current_time, arr_vec)
    candidate <- idx + 1
    if (candidate > length(arr_vec)) {
      return(stopover_arr_sorted[0, ])
    }
    chosen_arr_time <- arr_vec[candidate]
    matches <- which(stopover_arr_sorted$arrival == chosen_arr_time)
    row_index <- matches[1]
    return(stopover_arr_sorted[row_index, ])
  }
  
 tracking_data_period <- tracking_data_period %>%
    arrange(timestamp)
  
  segment_counter <- 1L
  
  for (i in seq_len(nrow(tracking_data_period))) {
    current_time <- tracking_data_period$timestamp[i]
    
    previous_stop <- get_closest_previous(current_time)
    
    if (!is.na(tracking_data_period$stopover_status[i]) &&
        tracking_data_period$stopover_status[i] == "entering stopover") {
      next_stop <- get_closest_next_by_departure(current_time)
    } else {
      next_stop <- get_closest_next_by_arrival(current_time)
    }
    
    if (nrow(previous_stop) > 0 && nrow(next_stop) > 0) {
      
      if (previous_stop$stopover_type == "origin") {
        tracking_data_period$segment_type[i] <- "origin to stopover"
      } else if (next_stop$stopover_type == "destination") {
        tracking_data_period$segment_type[i] <- "stopover to destination"
      } else {
        tracking_data_period$segment_type[i] <- "stopover to stopover"
      }
      
      if (is.na(tracking_data_period$segment_id[i])) {
        tracking_data_period$segment_id[i] <- segment_counter
      }
      
      if (!is.na(tracking_data_period$stopover_status[i]) &&
          tracking_data_period$stopover_status[i] == "entering stopover") {
        if (i < nrow(tracking_data_period)) {
          segment_counter <- segment_counter + 1L
        }
      }
      
      tracking_data_period$previous_stopover_depar_time[i] <- as.character(previous_stop$departure)
      tracking_data_period$next_stopover_arriv_time[i]     <- as.character(next_stop$arrival)
      
      tracking_data_period$previous_stopover_center_lon[i] <- previous_stop$centerLon
      tracking_data_period$previous_stopover_center_lat[i] <- previous_stop$centerLat
      tracking_data_period$next_stopover_center_lon[i]     <- next_stop$centerLon
      tracking_data_period$next_stopover_center_lat[i]     <- next_stop$centerLat
      
      tracking_data_period$previous_stopover_duration[i] <- previous_stop$duration
      tracking_data_period$next_stopover_duration[i]     <- next_stop$duration
      
    }
  }

  return(tracking_data_period)
}

# ------------------------------------------------------------------------------
# 4) Process all id_year across seasons; build master table
# ------------------------------------------------------------------------------
all_segments_list <- list()
idx <- 1L

for(season in names(migration_periods)) {
  
  period_dates <- migration_periods[[season]]
  
  for(current_id_year in unique(tracking_data$id_year)) {
    
    seg <- process_migration(
      tracking_data, stopover_data,
      current_id_year, season, period_dates
    )
    if (!is.null(seg)) {
      seg$season <- season
      all_segments_list[[idx]] <- seg
      idx <- idx + 1L
    }
  }
}

all_segments <- do.call(rbind, all_segments_list)

all_segments <- all_segments %>%
  filter(!is.na(segment_id)) %>%
  mutate(unique_seg = paste(id_year, segment_id, season, sep = "_"))


# ------------------------------------------------------------------------------
# 5) Short stop flag function
# ------------------------------------------------------------------------------
add_short_stop_flag <- function(season_data, stopover_data) {
  
  stopover_data <- stopover_data %>%
    mutate(year    = year(arrival),
           id_year = paste(individual.local.identifier, year, sep = "_"))
  
  segments <- season_data %>%
    group_by(individual.local.identifier, id_year, segment_id,unique_seg) %>%
    summarise(
      segment_start = min(timestamp, na.rm = TRUE),
      segment_end   = max(timestamp, na.rm = TRUE),
      .groups       = "drop"
    )
  
  segments_with_flag <- segments %>%
    left_join(stopover_data, by = c("individual.local.identifier", "id_year")) %>%
    mutate(
      overlap = (arrival < segment_end) & (departure > segment_start)
    ) %>%
    group_by(individual.local.identifier, id_year, segment_id,unique_seg) %>%
    summarise(
      contain_short_stop = any(overlap, na.rm = TRUE),
      .groups            = "drop"
    )
  
  season_data_flagged <- season_data %>%
    left_join(segments_with_flag, 
              by = c("individual.local.identifier", "id_year", "segment_id", "unique_seg")) %>%
    mutate(contain_short_stop = if_else(is.na(contain_short_stop), FALSE, contain_short_stop))
  
  return(season_data_flagged)
}


# ------------------------------------------------------------------------------
# 6) Apply short stop flag & export
# ------------------------------------------------------------------------------
Short_stopover_data <- read.csv("stopovers_short_filtered.csv", stringsAsFactors = FALSE)
Short_stopover_data$arrival <- as.POSIXct(Short_stopover_data$arrival, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
Short_stopover_data$departure <- as.POSIXct(Short_stopover_data$departure, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
 
flagged_all <- add_short_stop_flag(all_segments, Short_stopover_data)

flagged_all_sorted <- flagged_all %>%
  arrange(individual.local.identifier, unique_seg, timestamp)

write.csv(flagged_all_sorted, "points_all",row.names = FALSE)


