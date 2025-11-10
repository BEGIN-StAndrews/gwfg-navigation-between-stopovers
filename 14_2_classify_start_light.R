# ==============================================================================
# Title: Start-light classification for segment departures
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Define classify_start_light(lat, lon, t0) returning one of:
#      "polar_night" | "night" | "day" | "polar_day"
#   2) Apply to the departure fix (first_lat, first_long, first_timestamp)
#      for every segment
#   3) Summarise counts and percentages by class
#   4) Save per-segment labels to CSV
# Inputs:
#   - segments_final.csv  (expects columns: unique_seg, first_timestamp,
#                          first_lat, first_long, plus other metadata)
# Outputs:
#   - segments_light_start.csv  (unique_seg × light_start)
# Usage:
#   Rscript R/14_2_classify_start_light.R
# Notes:
#   - Sun elevations from {suncalc}; polar day/night handled via sunrise/sunset
# ==============================================================================

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------
library(suncalc)
library(lubridate)
library(dplyr)
library(scales)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------
segments <- read.csv("segments_final.csv", stringsAsFactors=FALSE) %>%
  mutate(
    first_timestamp     = as.POSIXct(first_timestamp,format = "%Y-%m-%d %H:%M:%OS", tz="UTC"),
    first_high_kp_time  = as.POSIXct(first_high_kp_time,format = "%Y-%m-%d %H:%M:%OS", tz="UTC"),
    first_lat           = as.numeric(first_lat),
    first_long          = as.numeric(first_long),
    last_lat            = as.numeric(last_lat),
    last_long           = as.numeric(last_long),
    duration            = as.numeric(duration),
    cumulative_distance = as.numeric(cumulative_distance),
    MaxKp               = as.numeric(MaxKp),
    unique_seg          = as.character(unique_seg) )

# ------------------------------------------------------------------------------
# Define classifier
# ------------------------------------------------------------------------------
classify_start_light <- function(lat, lon, t0) {
  st <- getSunlightTimes(
    date   = as_date(t0),
    lat    = lat,
    lon    = lon,
    keep   = c("sunrise", "sunset", "solarNoon", "nadir")
  )
  alt_now  <- getSunlightPosition(t0, lat, lon)$altitude
  
  if (!is.na(st$sunrise) && !is.na(st$sunset)) {
    if (alt_now < -6 * pi/180) {
      "night"
    } else {
      "day"
    }
  } else {
    # polar day/night
    alt_noon <- getSunlightPosition(st$solarNoon, lat, lon)$altitude
    if (alt_noon > 0) {
      "polar_day"
    } else {
      "polar_night"
    }
  }
}

# ------------------------------------------------------------------------------
# Apply per segment
# ------------------------------------------------------------------------------
segments <- segments %>%
  rowwise() %>%
  mutate(
    light_start = classify_start_light(
      first_lat, first_long, first_timestamp
    )
  ) %>%  ungroup()

# ------------------------------------------------------------------------------
# Summarise results (counts + %)
# ------------------------------------------------------------------------------
light_summary <- segments %>%
  count(light_start) %>%
  mutate(
    pct = scales::percent(n / sum(n))
  )

print(light_summary)

# ------------------------------------------------------------------------------
# Save per-segment labels
# ------------------------------------------------------------------------------
unique_light <- segments %>%distinct(unique_seg, light_start)
write.csv(unique_light,file = "segments_light_start.csv", row.names = FALSE)
