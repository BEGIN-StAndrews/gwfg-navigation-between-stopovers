# ==============================================================================
# Title: Initial heading extraction
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Remove points falling in short stopovers
#   2) Filter slow points (<10 m/s)
#   3) Compute distance to previous stopover centre and step headings
#   4) Derive initial mean/median heading within 100 km (fallback 120 km)
#   5) Join headings to segment table and export points/segments
# Inputs:  segments_kpinfo.csv, points_kpinfo.csv, stopovers_short_filtered.csv
# Outputs: segments_with_initial_heading.csv, points_with_initial_heading.csv
# Usage:   Rscript R/06_initial_heading.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(lubridate)
library(geosphere)    
library(circular)     
library(move)         

# ------------------------------------------------------------------------------
# 1) Load data
# ------------------------------------------------------------------------------
segments     <- read.csv("segments_kpinfo.csv", stringsAsFactors = FALSE)

tracking_data <- read.csv("points_kpinfo.csv", stringsAsFactors = FALSE) %>% 
  mutate(
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  ) %>%   arrange(individual.local.identifier, timestamp)

# ------------------------------------------------------------------------------
# 2) Remove points inside short stopovers (drop all interior points)
# ------------------------------------------------------------------------------
stopover_data <- read.csv("stopovers_short_filtered.csv", stringsAsFactors = FALSE)
stopover_data$arrival <- as.POSIXct(stopover_data$arrival, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
stopover_data$departure <- as.POSIXct(stopover_data$departure, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

tracking_data <- tracking_data %>%
  mutate(
    to_remove = 1  # 1 = keep by default, 0 = remove
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
  
  tracking_data$to_remove[idx_stop] <- 0
  
}
points_clean <- tracking_data[tracking_data$to_remove == 1, ]
points_clean$to_remove <- NULL  # Remove helper column
points_clean <- points_clean %>% arrange(unique_seg, timestamp)
write.csv(points_clean,    "points_shortstops_removed.csv",   row.names = FALSE)

# ------------------------------------------------------------------------------
# 3) Filter slow points (<10 m/s)
# ------------------------------------------------------------------------------
mv <- move(
  x     = points_clean$location.long,
  y     = points_clean$location.lat,
  time  = as.POSIXct(points_clean$timestamp,format = "%Y-%m-%d %H:%M:%S", tz = "UTC"), 
  proj  = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
  data  = points_clean,
  animal= points_clean$unique_seg
)

mv$Calspeed <- unlist(lapply(speed(mv), c, NA))
points_speed <- as.data.frame(mv) %>%
  filter(Calspeed >= 10) %>%
  dplyr::select(names(points_clean), Calspeed)

# ------------------------------------------------------------------------------
# 4) Distance to previous stopover & step headings
# ------------------------------------------------------------------------------
pts2 <- points_speed %>%
  rowwise() %>%
  mutate(
    dist_to_prev_km = distHaversine(
      c(previous_stopover_center_lon, previous_stopover_center_lat),
      c(location.long,               location.lat)
    ) / 1000
  ) %>%
  ungroup() %>%
  arrange(individual.local.identifier, id_year, unique_seg, timestamp) %>%
  group_by(individual.local.identifier, id_year, unique_seg) %>%
  mutate(
    lon_next    = lead(location.long),
    lat_next    = lead(location.lat),
    raw_bearing = bearing(
      p1 = cbind(location.long, location.lat),
      p2 = cbind(lon_next,       lat_next)
    )
  ) %>%
  ungroup() %>%
  mutate(
    heading_deg = (raw_bearing + 360) %% 360
  ) %>%  dplyr::select(-lon_next, -lat_next, -raw_bearing)

# ------------------------------------------------------------------------------
# 5) Initial heading within 100 km (fallback 120 km)
# ------------------------------------------------------------------------------
heading_summary <- pts2 %>%
  filter(!is.na(heading_deg), dist_to_prev_km < 100) %>%
  group_by(unique_seg) %>%
  summarise(
    n_pts_for_heading = n(),
    init_mean_heading  = as.numeric(mean.circular(
      circular(heading_deg, units="degrees", template="geographics"),
      na.rm = TRUE
    )),
    init_median_heading = as.numeric(median.circular(
      circular(heading_deg, units="degrees", template="geographics"),
      na.rm = TRUE
    )),
    .groups = "drop"
  )%>%
  mutate(
    init_mean_heading   = (init_mean_heading   + 360) %% 360,
    init_median_heading = (init_median_heading + 360) %% 360
  )

missing_segs <- setdiff(
  segments$unique_seg,
  heading_summary$unique_seg)

heading_summary_120 <- pts2 %>%
  filter(unique_seg %in% missing_segs, !is.na(heading_deg), dist_to_prev_km < 120) %>%
  group_by(unique_seg) %>%
  summarise(
    n_pts_for_heading  = n(),
    init_mean_heading   = as.numeric(mean.circular(
      circular(heading_deg, units="degrees", template="geographics"),
      na.rm = TRUE
    )),
    init_median_heading = as.numeric(median.circular(
      circular(heading_deg, units="degrees", template="geographics"),
      na.rm = TRUE
    )),
    .groups = "drop"
  ) %>%
  mutate(
    init_mean_heading   = (init_mean_heading   + 360) %% 360,
    init_median_heading = (init_median_heading + 360) %% 360
  )  

heading_summary <- bind_rows(heading_summary, heading_summary_120) %>%
  distinct(unique_seg, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 6) Join to segments, filter missing headings, export
# ------------------------------------------------------------------------------
segments_joint <- segments %>%left_join(heading_summary, by = "unique_seg")

segments_joint %>% 
  summarise(
    n_segments        = n(),
    na_mean_heading   = sum(is.na(init_mean_heading)),
    na_median_heading = sum(is.na(init_median_heading))
  ) %>% print()

valid_segs <- segments_joint %>%  filter(!is.na(init_mean_heading)) %>% pull(unique_seg)

segments_clean <- segments_joint %>%  filter(unique_seg %in% valid_segs)
write.csv(segments_clean, "segments_with_initial_heading.csv",row.names = FALSE)

tracking_data_clean <- tracking_data %>%  filter(unique_seg %in% valid_segs)
write.csv(tracking_data_clean, "points_with_initial_heading.csv", row.names = FALSE)
