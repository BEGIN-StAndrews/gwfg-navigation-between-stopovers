# ==============================================================================
# Title: Kp annotation
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Join hourly Kp to all migratory points (by season, id_year, timestamp)
#   2) Summarise per segment: first time Kp > threshold; MaxKp
#   3) Compute hours from segment start until first high-Kp
#   4) Export point- and segment-level tables with Kp info
# Inputs:  kp_raster.csv (derived via MagGeo), points_migratory.csv, segments_migratory.csv
# Outputs: segments_kpinfo.csv, points_kpinfo.csv
# Usage:   Rscript R/05_kp_annotation.R

# Citation:
#   - Benitez-Paez, F., Brum-Bastos, V., Beggan, C.D., Long, J.A., Demšar, U. (2021).
#     Fusion of wildlife tracking and satellite geomagnetic data for the study of animal migration.
#     Movement Ecology 9:31.
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(lubridate)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
high_kp_thresh <- 5

# ------------------------------------------------------------------------------
# 1) Load data
# ------------------------------------------------------------------------------
kp_raster        <- read.csv("kp_raster.csv", stringsAsFactors = FALSE)
points_combined  <- read.csv("points_migratory.csv",  stringsAsFactors = FALSE)
segments_combined <- read.csv("segments_migratory.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 2) Ensure POSIXct timestamps and order
# ------------------------------------------------------------------------------
kp_raster <- kp_raster %>%
  mutate(timestamp = as.POSIXct(timestamp,format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))%>%
  arrange(individual.local.identifier, timestamp)

points_combined <- points_combined %>%
  mutate(timestamp = as.POSIXct(timestamp,format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))%>%
  arrange(individual.local.identifier, timestamp)

segments_combined <- segments_combined %>%
  mutate(first_timestamp = as.POSIXct(first_timestamp,format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))%>%
  arrange(individual.local.identifier, first_timestamp)

# ------------------------------------------------------------------------------
# 3) Point-level join: attach Kp to every point
# ------------------------------------------------------------------------------
points_with_kp <- points_combined %>%
  left_join(kp_raster %>% 
              dplyr::select(season, id_year, timestamp, Kp),by = c("season", "id_year", "timestamp"))

# ------------------------------------------------------------------------------
# 4) Segment-level summaries (first high-Kp time; MaxKp)
# ------------------------------------------------------------------------------
high_kp_info<- points_with_kp %>%
  filter(!is.na(Kp) & Kp > high_kp_thresh) %>%
  group_by(season, id_year, segment_id, unique_seg) %>%
  summarise(first_high_kp_time = min(timestamp, na.rm = TRUE), .groups = "drop")

max_kp_info <- points_with_kp %>%
  group_by(season, id_year, segment_id, unique_seg) %>%
  summarise(MaxKp = ifelse(all(is.na(Kp)), NA, max(Kp, na.rm = TRUE)),
            .groups = "drop")

# ------------------------------------------------------------------------------
# 5) Merge summaries onto the segment table
# ------------------------------------------------------------------------------
segments_with_kp <- segments_combined %>%
  left_join(high_kp_info,
            by = c("season", "id_year", "segment_id", "unique_seg")) %>%
  left_join(max_kp_info,
            by = c("season", "id_year", "segment_id", "unique_seg")) %>%
  mutate(flight_duration_until_high_kp = as.numeric(difftime(first_high_kp_time, first_timestamp, units = "hours")))

# ------------------------------------------------------------------------------
# 6) put MaxKp on every point
# ------------------------------------------------------------------------------
points_with_kp <- points_with_kp %>%
  left_join(max_kp_info,
            by = c("season", "id_year", "segment_id", "unique_seg")) %>%
  arrange(individual.local.identifier, timestamp)

# ------------------------------------------------------------------------------
# 7) Quick summary report
# ------------------------------------------------------------------------------
total_segments <- nrow(max_kp_info)
high_kp_count  <- sum(!is.na(max_kp_info$MaxKp) & max_kp_info$MaxKp > high_kp_thresh)
low_kp_count   <- sum(!is.na(max_kp_info$MaxKp) & max_kp_info$MaxKp <=  high_kp_thresh)
na_kp_count    <- sum(is.na(max_kp_info$MaxKp))

# ------------------------------------------------------------------------------
# 8) Export
# ------------------------------------------------------------------------------
write.csv(segments_with_kp, "segments_kpinfo.csv", row.names = FALSE)
write.csv(points_with_kp,    "points_kpinfo.csv",   row.names = FALSE)
