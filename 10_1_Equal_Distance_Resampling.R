# ==============================================================================
# Title: Equal-interval resampling (Actual, GC, Wind) + Common-segment filter
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Load Actual hourly points, Great-circle (GC) nodes, Wind-optimal DP nodes
#   2) Load other route tables (Geomagnetic, Geographic, Magnetoclinic, Sun, Local-wind)
#   3) Compute the intersection of segment IDs present in *all* routes (common set)
#   4) Filter the master segments file to this common set and save as segments_final.csv
#   5) For each common segment, resample Actual/GC/Wind to the same number of
#      equally spaced points as the Actual track (equal distance along polyline)
#   6) Export harmonised CSVs (actual_equal_interval.csv, GreatCircle_equal_interval.csv,
#      wind_optimal_equal_interval.csv) with identical per-segment indexing (idx)
# Inputs:
#   - points_resampled.csv               (Actual hourly points)
#   - GreatCircle_routes.csv             (GC route nodes)
#   - wind_optimal_routes_raw.csv        (DP nodes per slice)
#   - Geomagnetic_routes.csv             (Geomagnetic route nodes)
#   - Geographic_routes.csv              (Geographic loxodrome nodes)
#   - Magnetoclinic_routes.csv           (Magnetoclinic route nodes)
#   - SunCompass_classic_routes.csv      (Sun-compass route nodes)
#   - local_wind_routes.csv              (Local wind-aligned route nodes)
#   - segments_resample.csv              (Master segment table)
# Outputs:
#   - segments_final.csv                 (Segments present in *all* routes)
#   - actual_equal_interval.csv
#   - GreatCircle_equal_interval.csv
#   - wind_optimal_equal_interval.csv
# Usage:
#   Rscript R/10_1_equal_interval_resampling.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(readr)
library(geosphere)
library(tibble)

# ------------------------------------------------------------------------------
# Parameters (file paths)
# ------------------------------------------------------------------------------
actual_file <- "points_resampled.csv"
gc_file     <- "GreatCircle_routes.csv"
wind_file   <- "wind_optimal_routes_raw.csv"

geomagnetic_file   <- "Geomagnetic_routes.csv"
geographic_file    <- "Geographic_routes.csv"
magnetoclinic_file <- "Magnetoclinic_routes.csv"
sun_file           <- "SunCompass_classic_routes.csv"
local_wind_file    <- "local_wind_routes.csv"

segments_master_in  <- "segments_resample.csv"
segments_final_out  <- "segments_final.csv"

out_actual <- "actual_equal_interval.csv"
out_gc     <- "GreatCircle_equal_interval.csv"
out_wind   <- "wind_optimal_equal_interval.csv"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
# Robust route ordering by any known index column (falls back to as-is)
arrange_route <- function(df) {
  order_cols <- c("idx", "step", "time", "timestamp", "layer")
  present <- intersect(order_cols, names(df))
  if (length(present) > 0) {
    df %>% arrange(unique_seg, !!!rlang::syms(present))
  } else {
    df %>% arrange(unique_seg)
  }
}

# Equal-distance resampler for a route data.frame and target N
resample_equal <- function(route_df, N) {
  stopifnot(all(c("lon","lat") %in% names(route_df)))
  route_df <- route_df %>% dplyr::select(lon, lat)
  if (N < 2L || nrow(route_df) < 2L) {
    # if too short, repeat first point to length N to keep shape consistent
    return(tibble(lon = rep(route_df$lon[1], max(N, 1)),
                  lat = rep(route_df$lat[1], max(N, 1))) %>%
             slice(1:N))
  }
  coords  <- as.matrix(route_df[, c("lon","lat")])
  dstep   <- distHaversine(coords[-nrow(coords), ], coords[-1, ]) / 1000
  cumd    <- c(0, cumsum(dstep))
  targets <- seq(0, tail(cumd, 1), length.out = N)
  lon_eq  <- approx(cumd, route_df$lon, xout = targets, rule = 2)$y
  lat_eq  <- approx(cumd, route_df$lat, xout = targets, rule = 2)$y
  tibble(lon = lon_eq, lat = lat_eq)
}

# ------------------------------------------------------------------------------
# Load inputs (Actual/GC/Wind) with robust ordering
# ------------------------------------------------------------------------------
actual_raw <- read.csv(actual_file, stringsAsFactors = FALSE) %>% arrange_route()
gc_raw     <- read.csv(gc_file,     stringsAsFactors = FALSE) %>% arrange_route()
wind_raw   <- read.csv(wind_file,   stringsAsFactors = FALSE) %>% arrange_route()

# ------------------------------------------------------------------------------
# Load other route tables (IDs only needed, ordering not required here)
# ------------------------------------------------------------------------------
geomagnetic_raw   <- read.csv(geomagnetic_file,   stringsAsFactors = FALSE)
geographic_raw    <- read.csv(geographic_file,    stringsAsFactors = FALSE)
magnetoclinic_raw <- read.csv(magnetoclinic_file, stringsAsFactors = FALSE)
sun_raw           <- read.csv(sun_file,           stringsAsFactors = FALSE)
local_wind_raw    <- read.csv(local_wind_file,    stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Common set of segments present in *all* routes
# ------------------------------------------------------------------------------
segs_all <- Reduce(intersect, list(
  unique(actual_raw$unique_seg),
  unique(gc_raw$unique_seg),
  unique(wind_raw$unique_seg),
  unique(geomagnetic_raw$unique_seg),
  unique(geographic_raw$unique_seg),
  unique(magnetoclinic_raw$unique_seg),
  unique(sun_raw$unique_seg),
  unique(local_wind_raw$unique_seg)
))

# ------------------------------------------------------------------------------
# Filter the master segments table and save as segments_final.csv
# ------------------------------------------------------------------------------
segments_master <- read.csv(segments_master_in, stringsAsFactors = FALSE)
segments_final  <- segments_master %>% dplyr::filter(unique_seg %in% segs_all)
write.csv(segments_final, segments_final_out, row.names = FALSE)

# ------------------------------------------------------------------------------
# Resample per segment (using only the common set)
# ------------------------------------------------------------------------------
segs <- segs_all

actual_list <- gc_list <- wind_list <- vector("list", length(segs))
names(actual_list) <- names(gc_list) <- names(wind_list) <- segs

for (seg in segs) {
  act <- dplyr::filter(actual_raw, unique_seg == seg)
  gc  <- dplyr::filter(gc_raw,     unique_seg == seg)
  wd  <- dplyr::filter(wind_raw,   unique_seg == seg)
  
  N <- nrow(act)  # target count taken from Actual
  
  actual_list[[seg]] <- resample_equal(act, N) %>%
    mutate(unique_seg = seg, idx = dplyr::row_number())
  
  gc_list[[seg]]     <- resample_equal(gc,  N) %>%
    mutate(unique_seg = seg, idx = dplyr::row_number())
  
  wind_list[[seg]]   <- resample_equal(wd,  N) %>%
    mutate(unique_seg = seg, idx = dplyr::row_number())
}

# ------------------------------------------------------------------------------
# Bind & export
# ------------------------------------------------------------------------------
bind_rows(actual_list) %>%
  dplyr::select(unique_seg, idx, lon, lat) %>%
  write.csv(out_actual, row.names = FALSE)

bind_rows(gc_list) %>%
  dplyr::select(unique_seg, idx, lon, lat) %>%
  write.csv(out_gc, row.names = FALSE)

bind_rows(wind_list) %>%
  dplyr::select(unique_seg, idx, lon, lat) %>%
  write.csv(out_wind, row.names = FALSE)
