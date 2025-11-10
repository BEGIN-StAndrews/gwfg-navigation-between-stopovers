# ==============================================================================
# Title: Assemble segment metadata + attach best-route labels (per metric/level)
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   PART A — Build a unified segment table:
#     1) Read base segments
#     2) Add selected per-point meta (stopover fields) collapsed to one row/segment
#     3) Join wind summaries and start-light classes
#     4) Drop unwanted columns and write segments_all_meta.csv
#   PART B — Add best_route labels:
#     5) Read dominance_by_segment.csv (unique_seg × metric × level)
#     6) For each metric (med, dtw, dir), join best_route separately for:
#        (i) efficiency level → segments_all_meta_efficiency_<metric>.csv
#        (ii) compass level   → segments_all_meta_compass_<metric>.csv
# Inputs:
#   - segments_final.csv
#   - points_wind_metrics.csv
#   - segments_wind_summary.csv
#   - segments_light_start.csv
#   - dominance_by_segment.csv
# Outputs:
#   - segments_all_meta.csv
#   - segments_all_meta_efficiency_med.csv
#   - segments_all_meta_efficiency_dtw.csv
#   - segments_all_meta_efficiency_dir.csv
#   - segments_all_meta_compass_med.csv
#   - segments_all_meta_compass_dtw.csv
#   - segments_all_meta_compass_dir.csv
# Usage:
#   Rscript R/14_3_segments_all_meta_and_labels.R
# ==============================================================================

# ──────────────────────────────────────────────────────────────────────────────
# PART A — Build unified segment table (segments_all_meta.csv)
# ──────────────────────────────────────────────────────────────────────────────

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------
library(dplyr)

# ------------------------------------------------------------------------------
# Read base segments
# ------------------------------------------------------------------------------
segments <- read.csv("segments_final.csv", stringsAsFactors = FALSE) %>%
  mutate(
    first_timestamp     = as.POSIXct(first_timestamp,format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    first_lat           = as.numeric(first_lat),
    first_long          = as.numeric(first_long),
    last_lat            = as.numeric(last_lat),
    last_long           = as.numeric(last_long),
    duration            = as.numeric(duration),
    cumulative_distance = as.numeric(cumulative_distance),
    great_circle_distance = as.numeric(great_circle_distance),
    MaxKp               = as.numeric(MaxKp),
    unique_seg          = as.character(unique_seg),
    straightness = great_circle_distance / cumulative_distance,  
    geomag_storm = ifelse(!is.na(MaxKp) & MaxKp > 5, TRUE, FALSE))

# ------------------------------------------------------------------------------
# Read per-point file and collapse to one row per segment (selected columns)
# ------------------------------------------------------------------------------
points <- read.csv("points_wind_metrics.csv", stringsAsFactors = FALSE) %>%
  mutate(unique_seg = as.character(unique_seg))

cols_to_keep <- c("unique_seg","previous_stopover_duration", "next_stopover_duration", "segment_type", "contain_short_stop")
points_meta <- points %>% dplyr::select(all_of(cols_to_keep)) %>% distinct()

# ------------------------------------------------------------------------------
# Read wind summary and start-light labels
# ------------------------------------------------------------------------------
wind_summary <- read.csv("segments_wind_summary.csv", stringsAsFactors = FALSE) %>%
  mutate(unique_seg = as.character(unique_seg))

light_start <- read.csv("segments_light_start.csv",stringsAsFactors = FALSE) %>%
  mutate(unique_seg = as.character(unique_seg))

# ------------------------------------------------------------------------------
# Join all to segments and preview
# ------------------------------------------------------------------------------
segments_all <- segments %>%
  left_join(points_meta,   by = "unique_seg") %>%
  left_join(wind_summary,  by = "unique_seg") %>%
  left_join(light_start,   by = "unique_seg")

# ------------------------------------------------------------------------------
# Drop unwanted columns
# ------------------------------------------------------------------------------
segments_all_clean <- segments_all %>%
  select(
    -id_year,
    -segment_id,
    -NumberOFpoints,
    -median_time_interval,
    -first_timestamp,
    -first_long,
    -first_lat,
    -last_timestamp,
    -last_long,
    -last_lat,
    -first_high_kp_time,
    -flight_duration_until_high_kp,
    -n_pts_for_heading,
    -init_mean_heading,
    -init_median_heading,
    -MaxKp,
    -median_speed,
    -straightness)

# ------------------------------------------------------------------------------
# Write unified meta table
# ------------------------------------------------------------------------------
write.csv(segments_all_clean, file      = "segments_all_meta.csv",row.names = FALSE)



# ──────────────────────────────────────────────────────────────────────────────
# PART B — Add best_route labels from dominance_by_segment.csv
# ──────────────────────────────────────────────────────────────────────────────

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# Read unified meta and dominance tables
# ------------------------------------------------------------------------------
segments_all <- read_csv("segments_all_meta.csv", show_col_types = FALSE) %>%
  mutate(unique_seg = as.character(unique_seg))

dominance <- read_csv("dominance_by_segment.csv", show_col_types = FALSE) %>%
  mutate(unique_seg = as.character(unique_seg))

# ------------------------------------------------------------------------------
# Efficiency level: write per-metric labeled segments
# ------------------------------------------------------------------------------
for (m in c("med", "dtw", "dir")) {
  best_m <- dominance %>%
    filter(metric == m, level=="efficiency") %>%
    select(unique_seg, best_route)
    segments_m <- segments_all %>%
    left_join(best_m, by = "unique_seg")
    out_fname <- paste0("segments_all_meta_efficiency_", m, ".csv")
  write_csv(segments_m, out_fname)
  
  message(sprintf("✔ Wrote %s (added best_route for metric '%s')", out_fname, m))
}

# ------------------------------------------------------------------------------
# Compass level: write per-metric labeled segments
# ------------------------------------------------------------------------------
for (m in c("med", "dtw", "dir")) {
    best_m <- dominance %>%
    filter(metric == m, level=="compass") %>%
    select(unique_seg, best_route)
    segments_m <- segments_all %>%
    left_join(best_m, by = "unique_seg")
  out_fname <- paste0("segments_all_meta_compass_", m, ".csv")
  write_csv(segments_m, out_fname)
  
  message(sprintf("✔ Wrote %s (added best_route for metric '%s')", out_fname, m))
}
