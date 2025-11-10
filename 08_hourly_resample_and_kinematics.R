# ==============================================================================
# Title: Hourly resampling & kinematics for migratory segments
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read Kp-annotated, heading-ready points
#   2) Resample each segment onto an hourly grid via linear interpolation
#   3) Reattach immutable segment metadata
#   4) Compute CalSpeed / CalTurnAngle / CalHeading using {move}
#   5) Export resampled points and refresh cumulative_distance for segments
# Inputs:
#   - points_with_initial_heading.csv
#   - segments_with_initial_heading.csv
# Outputs:
#   - points_resampled.csv
#   - segments_resample.csv
# Usage:
#   Rscript R/08_hourly_resample_and_kinematics.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(lubridate)
library(readr)
library(move)
library(geosphere)

# ------------------------------------------------------------------------------
# Meta-column names
# ------------------------------------------------------------------------------
meta_cols <- c(
  "individual.local.identifier", "id_year", "segment_id","season",
  "previous_stopover_depar_time", "previous_stopover_center_lon", "previous_stopover_center_lat","previous_stopover_duration",
  "next_stopover_arriv_time", "next_stopover_center_lon","next_stopover_center_lat","next_stopover_duration",
  "segment_type", "contain_short_stop", "MaxKp")

# ------------------------------------------------------------------------------
# Linear resampler (hourly)
# ------------------------------------------------------------------------------
resample_segment_hourly_lin <- function(segment_df, alt_col = "height.above.msl") {
  segment_df <- segment_df[order(segment_df$timestamp), ]
  if (nrow(segment_df) < 2) return(segment_df)
  
  t_grid <- seq(from = min(segment_df$timestamp),
                to   = max(segment_df$timestamp),
                by   = "hour")
  if (last(t_grid) < max(segment_df$timestamp))
    t_grid <- c(t_grid, max(segment_df$timestamp))
  
  ts_num   <- as.numeric(segment_df$timestamp)
  grid_num <- as.numeric(t_grid)
  
  lon_out <- approx(x    = ts_num,
                    y    = segment_df$location.long,
                    xout = grid_num)$y
  lat_out <- approx(x    = ts_num,
                    y    = segment_df$location.lat,
                    xout = grid_num)$y
  
  valid_alts <- segment_df[[alt_col]]
  n_valid    <- sum(!is.na(valid_alts))
  
  if (n_valid == 0) {
    alt_out <- rep(NA_real_, length(t_grid))
  } else if (n_valid == 1) {
    alt_out <- rep(valid_alts[!is.na(valid_alts)], length(t_grid))
  } else {
    alt_out <- approx(x    = ts_num,
                      y    = valid_alts,
                      xout = grid_num,
                      method = "linear",
                      ties   = mean 
    )$y
  }
  
  data.frame(
    unique_seg       = segment_df$unique_seg[1],
    timestamp        = t_grid,
    location.long    = lon_out,
    location.lat     = lat_out,
    height.above.msl = alt_out,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Read input points (with initial headings)
# ------------------------------------------------------------------------------
raw_df <- read.csv("points_with_initial_heading.csv", stringsAsFactors = FALSE) %>%
  mutate(timestamp = as.POSIXct(timestamp,format = "%Y-%m-%d %H:%M:%S",tz = "UTC"))%>%
  arrange(unique_seg, timestamp)

# ------------------------------------------------------------------------------
# Resample per segment
# ------------------------------------------------------------------------------
seg_ids     <- unique(raw_df$unique_seg)
resampled_lst <- setNames(vector("list", length(seg_ids)), seg_ids)

for (seg in seg_ids) {
  
  seg_df <- raw_df %>% filter(unique_seg == seg)
  resampled_lst[[seg]] <- resample_segment_hourly_lin(seg_df)
}

# ------------------------------------------------------------------------------
# Reattach segment metadata
# ------------------------------------------------------------------------------
segment_meta <- raw_df %>%
  dplyr::select(unique_seg, all_of(meta_cols)) %>%
  distinct()

combined_resampled <- bind_rows(resampled_lst) %>%
  left_join(segment_meta, by = "unique_seg") %>%
  arrange(unique_seg, timestamp)

# ------------------------------------------------------------------------------
# Kinematics via {move}
# ------------------------------------------------------------------------------
track_mv <- move(
  x    = combined_resampled$location.long,
  y    = combined_resampled$location.lat,
  time = combined_resampled$timestamp,
  proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
  data = combined_resampled,
  animal = combined_resampled$unique_seg
)

track_mv$CalSpeed     <- unlist(lapply(speed(track_mv),      c, NA))
track_mv$CalTurnAngle <- unlist(lapply(turnAngleGc(track_mv), function(x) c(NA, x, NA)))
track_mv$CalHeading   <- unlist(lapply(angle(track_mv),       c, NA))

track_df <- as.data.frame(track_mv)

convert_heading <- function(heading) {
  ifelse(heading < 0, heading + 360, heading)
}

track_df$CalHeading <- convert_heading(track_df$CalHeading)

desired_cols <- c(
  "timestamp", "location.long", "location.lat","unique_seg",
  meta_cols,
  "CalSpeed", "CalTurnAngle", "CalHeading"
)

track_df <- track_df[, intersect(names(track_df), desired_cols)]
write.csv(track_df, "points_resampled.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Refresh cumulative_distance in segments table
# ------------------------------------------------------------------------------
compute_cumulative_distance <- function(df) {
  n <- nrow(df)
  if (n < 2) return(0)
  coords1  <- cbind(df$location.long[1:(n-1)], df$location.lat[1:(n-1)])
  coords2  <- cbind(df$location.long[2:n],   df$location.lat[2:n])
  dists_m  <- distHaversine(coords1, coords2)
  sum(dists_m, na.rm = TRUE) / 1000
}

segments_df    <- read.csv("segments_with_initial_heading.csv",stringsAsFactors = FALSE) %>%
  dplyr::select(-max_time_interval, -time_gap_for_long_displacement)%>%
  arrange(unique_seg)

resamp_points  <- read.csv("points_resampled.csv" ,stringsAsFactors = FALSE) %>%
  mutate(timestamp = as.POSIXct(timestamp,format = "%Y-%m-%d %H:%M:%S",tz = "UTC")) %>%
  arrange(unique_seg, timestamp)

new_cumdist <- resamp_points %>%
  group_by(unique_seg) %>%
  summarise(
    new_cumulative_distance = compute_cumulative_distance(cur_data()),
    .groups = "drop")%>%
  arrange(unique_seg)

segments_updated <- segments_df %>%
  left_join(new_cumdist, by = "unique_seg") %>%
  mutate(cumulative_distance = new_cumulative_distance) %>%
  dplyr::select(-new_cumulative_distance)

write.csv(segments_updated, "segments_resample.csv", row.names = FALSE)
