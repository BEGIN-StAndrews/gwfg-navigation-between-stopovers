# ==============================================================================
# Title: Static geomagnetic corridor raster (low-Kp segments)
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Filter segments to low Kp (≤5)
#   2) Build one local LCC per segment (init/first/last bounding box)
#   3) Generate a widening corridor along initial heading with start/end buffers
#   4) Rasterize corridor to 10 km grid cells and take centroids
#   5) Export full raster for MagGeo annotation
# Inputs:  segments_resample.csv
# Outputs: geomag_corridor_static.csv
# Usage:   Rscript R/07_1_geomag_corridor_static.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(sf)
library(geosphere)

# ------------------------------------------------------------------------------
# 1) Read inputs
# ------------------------------------------------------------------------------
segments_file <- "segments_resample.csv"
segments <- read_csv(segments_file)%>%  filter(is.na(MaxKp) | MaxKp <= 5)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------
build_combined_lcc_crs <- function(first_lat, first_lon, last_lat, last_lon) {
  lat_min <- min(first_lat, last_lat)
  lat_max <- max(first_lat, last_lat)
  lon_min <- min(first_lon, last_lon)
  lon_max <- max(first_lon, last_lon)
  
  lat0 <- mean(c(lat_min, lat_max))
  lon0 <- mean(c(lon_min, lon_max))
  
  paste0(
    "+proj=lcc +lat_1=", lat_min,
    " +lat_2=", lat_max,
    " +lat_0=", lat0,
    " +lon_0=", lon0,
    " +datum=WGS84 +units=m +no_defs"
  )
}

corridor_generation <- function(start_lat, start_lon,end_lat,end_lon,heading_deg,lcc_crs,start_buffer_m,end_buffer_m) {
  
  if (any(is.na(c(start_lat, start_lon, end_lat, end_lon, heading_deg)))) {
    return(NULL)
  }
  pt_start <- st_sfc(st_point(c(start_lon, start_lat)), crs = 4326)
  pt_end   <- st_sfc(st_point(c(end_lon,   end_lat)),   crs = 4326)
  start_lcc <- st_transform(pt_start, lcc_crs)
  end_lcc   <- st_transform(pt_end,   lcc_crs)
  
  start_xy <- st_coordinates(start_lcc)
  end_xy   <- st_coordinates(end_lcc)
  
  base_dist_m <- as.numeric(st_distance(start_lcc, end_lcc))
  if (is.na(base_dist_m) || base_dist_m < 1) {
    return(NULL)
  }
  segment_distance <- start_buffer_m + base_dist_m + end_buffer_m
  
  W_start <- 20000
  W_end   <- 0.2 * segment_distance
  
  if (W_end < W_start) {
    W_end <- W_start
  }
  step_m <- 30000
  dist_seq <- seq(0, segment_distance, by = step_m)
  if (tail(dist_seq, 1) < segment_distance) {
    dist_seq <- c(dist_seq, segment_distance)
  }
  heading_rad <- heading_deg * pi / 180  # 0=North, 90=East
  half_width <- function(d) {
    (W_start + (d / segment_distance)*(W_end - W_start)) / 2
  }
  
  start_xy_buffered <- start_xy +
        c(  -start_buffer_m * sin(heading_rad),
            -start_buffer_m * cos(heading_rad))
  
    center_pts <- lapply(dist_seq, function(d) {
        dx <- d * sin(heading_rad)
        dy <- d * cos(heading_rad)
        c(start_xy_buffered[1] + dx, start_xy_buffered[2] + dy) })
    
  polygons <- list()
  for (i in seq_len(length(dist_seq) - 1)) {
    d1 <- dist_seq[i]
    d2 <- dist_seq[i + 1]
    p1 <- center_pts[[i]]
    p2 <- center_pts[[i+1]]
    w1 <- half_width(d1)
    w2 <- half_width(d2)
    
    left1  <- c(p1[1] + w1*sin(heading_rad - pi/2),
                p1[2] + w1*cos(heading_rad - pi/2))
    right1 <- c(p1[1] + w1*sin(heading_rad + pi/2),
                p1[2] + w1*cos(heading_rad + pi/2))
    left2  <- c(p2[1] + w2*sin(heading_rad - pi/2),
                p2[2] + w2*cos(heading_rad - pi/2))
    right2 <- c(p2[1] + w2*sin(heading_rad + pi/2),
                p2[2] + w2*cos(heading_rad + pi/2))
    
    coords <- rbind(left1, left2, right2, right1, left1)
    polygons[[i]] <- st_polygon(list(coords))
  }
  
  multi_sf       <- st_sfc(polygons, crs = lcc_crs)
  corridor_union <- suppressMessages(st_union(multi_sf))
  corridor_union <- st_make_valid(corridor_union)
  if (is.null(corridor_union) || length(corridor_union) == 0) {
    return(NULL)
  }
  st_sf(geometry = corridor_union)
}


corridor_initialisation <- function(segment_row) {
  first_lat <- segment_row$first_lat
  first_lon <- segment_row$first_long
  last_lat  <- segment_row$last_lat
  last_lon  <- segment_row$last_long
  init_hdg  <- segment_row$init_mean_heading
  seg_length <- segment_row$great_circle_distance
  
  if (any(is.na(c(init_hdg, first_lat, first_lon, last_lat, last_lon)))) { return(NULL) }
  
  lcc_crs <- build_combined_lcc_crs(first_lat, first_lon, last_lat, last_lon)
  
    corridor <- corridor_generation(
    start_lat      = first_lat,
    start_lon      = first_lon,
    end_lat        = last_lat,
    end_lon        = last_lon,
    heading_deg    = init_hdg,
    lcc_crs        = lcc_crs,
    start_buffer_m = 10000,                       
    end_buffer_m   = max(15000, seg_length*10)
  ) 
  
  corridor_bbox <- st_bbox(corridor)
  grid_size_m   <- 10000
  
  grid_lcc <- st_make_grid(
    st_as_sfc(corridor_bbox), 
    cellsize = c(grid_size_m, grid_size_m), 
    square   = TRUE )
  
  grid_lcc <- st_sf(geometry = grid_lcc, crs = lcc_crs)
  
  idx <- st_intersects(grid_lcc, corridor, sparse = FALSE)[,1]
  if (!any(idx)) return(NULL)
  
  selected_cells <- grid_lcc[idx,]
  
  selected_wgs84 <- st_transform(selected_cells, 4326)
  centroids      <- st_centroid(selected_wgs84)
  coords         <- st_coordinates(centroids)
  
  out <- segment_row %>%
    dplyr::select(
      unique_seg,
      first_timestamp,  
      last_timestamp,
      init_mean_heading)
  
  out_rep <- out[rep(1, nrow(coords)), ]
  out_rep$grid_lon <- coords[, "X"]
  out_rep$grid_lat <- coords[, "Y"]
  
  return(out_rep)
}

# ------------------------------------------------------------------------------
# Pipeline
# ------------------------------------------------------------------------------
segments_list <- split(segments, seq_len(nrow(segments)))
results_list  <- lapply(segments_list, corridor_initialisation)
final_result  <- bind_rows(results_list)

output_file <- "GeomagRaster_static.csv"
write.csv(final_result, output_file, row.names = FALSE)
