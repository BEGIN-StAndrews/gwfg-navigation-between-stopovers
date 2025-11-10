# ==============================================================================
# Title: Dynamic geomagnetic corridor raster (Kp > 5, hourly slices)
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Keep segments with Kp > 5 (storm exposure)
#   2) Use resampled track points to compute cumulative distance per segment
#   3) Build static pre-storm corridor + hourly post-storm dynamic slices (LCC)
#   4) Rasterise each slice to grid-cell centroids (WGS84) with valid time
#   5) Export a single CSV ready for MagGeo annotation
# Inputs:
#   - segments_with_initial_heading.csv
#     (must include: unique_seg, segment_id, individual.local.identifier,
#      first_timestamp, first_high_kp_time, last_timestamp,
#      first_lat/first_long, last_lat/last_long,
#      init_mean_heading, great_circle_distance, MaxKp, duration)
#   - points_with_initial_heading.csv
#     (must include: unique_seg, timestamp, location.long, location.lat)
# Output:
#   - geomag_corridor_dynamic_hourly.csv
# Usage:
#   Rscript R/07_2_geomag_corridor_dynamic_hourly.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(sf)
library(lubridate)    
library(geosphere)    


# ------------------------------------------------------------------------------
# 1) Read inputs
# ------------------------------------------------------------------------------
seg_file <- "segments_with_initial_heading.csv"
pts_file <- "points_with_initial_heading.csv"

segments <- readr::read_csv(seg_file, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(MaxKp) & MaxKp > 5) %>%
  dplyr::mutate(
    first_timestamp    = as.POSIXct(first_timestamp,    format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    first_high_kp_time = as.POSIXct(first_high_kp_time, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    last_timestamp     = as.POSIXct(last_timestamp,     format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    duration           = as.numeric(duration)
  )

points_df <- readr::read_csv(pts_file, show_col_types = FALSE) %>%
  dplyr::mutate(
    timestamp  = as.POSIXct(timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    unique_seg = as.character(unique_seg)
  ) %>%
  dplyr::arrange(unique_seg, timestamp)


# ----------------------------------------------------------------------
# 2) Distance along each segment
# ----------------------------------------------------------------------
points_df <- points_df %>%
  group_by(unique_seg) %>%
  arrange(timestamp) %>%
  mutate(
    prev_lon = lag(location.long, default = first(location.long)),
    prev_lat = lag(location.lat,  default = first(location.lat)),
    step_km  = distHaversine(
      cbind(prev_lon, prev_lat),
      cbind(location.long, location.lat)
    ) / 1000,
    step_km  = ifelse(is.na(step_km), 0, step_km),  
    cumulative_distance = cumsum(step_km)
  ) %>%
  ungroup() %>%
  dplyr::select(-prev_lon, -prev_lat, -step_km) %>%
  arrange(unique_seg, timestamp)


# ----------------------------------------------------------------------
# 3) Lambert CRS helper
# ----------------------------------------------------------------------
build_combined_lcc_crs <- function(lat1, lon1, lat2, lon2) {
  paste0(
    "+proj=lcc +lat_1=", min(lat1, lat2),
    " +lat_2=", max(lat1, lat2),
    " +lat_0=", mean(c(lat1, lat2)),
    " +lon_0=", mean(c(lon1, lon2)),
    " +datum=WGS84 +units=m +no_defs"
  )
}

# ----------------------------------------------------------------------
# 4) Corridor geometry helpers
# ----------------------------------------------------------------------
corridor_static <- function(
    start_lat, start_lon,
    end_lat,   end_lon,
    heading_deg,
    lcc_crs,
    start_buffer_m = 20000,
    end_buffer_m   = 50000
) {
  pt_start <- st_transform(st_sfc(st_point(c(start_lon, start_lat)), crs = 4326), lcc_crs)
  pt_end   <- st_transform(st_sfc(st_point(c(end_lon,   end_lat  )), crs = 4326), lcc_crs)
  
  base_dist_m <- as.numeric(st_distance(pt_start, pt_end))
  if (is.na(base_dist_m) || base_dist_m < 1) return(NULL)
  
  seg_dist <- start_buffer_m + base_dist_m + end_buffer_m
  W_start  <- 100000
  W_end    <- max(W_start, 0.3 * seg_dist)
  step_m   <- 30000
  dist_seq <- seq(0, seg_dist, by = step_m)
  if (tail(dist_seq, 1) < seg_dist) dist_seq <- c(dist_seq, seg_dist)
  
  heading_rad <- heading_deg * pi / 180
  half_width  <- function(d) (W_start + (d / seg_dist) * (W_end - W_start)) / 2
  
  start_xy     <- st_coordinates(pt_start)
  start_xy_buf <- start_xy + c(-start_buffer_m * sin(heading_rad),
                               -start_buffer_m * cos(heading_rad))
  
  polys <- vector("list", length(dist_seq) - 1)
  for (i in seq_len(length(dist_seq) - 1)) {
    d1 <- dist_seq[i]; d2 <- dist_seq[i + 1]
    p1 <- start_xy_buf + c(d1 * sin(heading_rad), d1 * cos(heading_rad))
    p2 <- start_xy_buf + c(d2 * sin(heading_rad), d2 * cos(heading_rad))
    w1 <- half_width(d1); w2 <- half_width(d2)
    l1 <- p1 + w1 * c(sin(heading_rad - pi / 2), cos(heading_rad - pi / 2))
    r1 <- p1 + w1 * c(sin(heading_rad + pi / 2), cos(heading_rad + pi / 2))
    l2 <- p2 + w2 * c(sin(heading_rad - pi / 2), cos(heading_rad - pi / 2))
    r2 <- p2 + w2 * c(sin(heading_rad + pi / 2), cos(heading_rad + pi / 2))
    polys[[i]] <- st_polygon(list(rbind(l1, l2, r2, r1, l1)))
  }
  suppressMessages(
    st_union(st_sfc(polys, crs = lcc_crs))
  ) %>% st_make_valid() %>% st_sf()
}

corridor_dynamic <- function(
    start_lat, start_lon,
    heading_deg,
    dist_min_m, dist_max_m,
    width_m,
    lcc_crs
) {
  if (dist_max_m <= dist_min_m || width_m < 1) return(NULL)
  
  start_xy <- st_coordinates(
    st_transform(
      st_sfc(st_point(c(start_lon, start_lat)), crs = 4326),
      lcc_crs)
  )[1, ]
  heading_rad <- heading_deg * pi / 180
  halfW <- width_m / 2
  
  pA <- start_xy + c(dist_min_m * sin(heading_rad),
                     dist_min_m * cos(heading_rad))
  pB <- start_xy + c(dist_max_m * sin(heading_rad),
                     dist_max_m * cos(heading_rad))
  lA <- pA + halfW * c(sin(heading_rad - pi / 2), cos(heading_rad - pi / 2))
  rA <- pA + halfW * c(sin(heading_rad + pi / 2), cos(heading_rad + pi / 2))
  lB <- pB + halfW * c(sin(heading_rad - pi / 2), cos(heading_rad - pi / 2))
  rB <- pB + halfW * c(sin(heading_rad + pi / 2), cos(heading_rad + pi / 2))
  
  st_union(
    st_sfc(st_polygon(list(rbind(lA, lB, rB, rA, lA))), crs = lcc_crs)
  ) %>% st_make_valid() %>% st_sf()
}

# ----------------------------------------------------------------------
# 5) Per-segment processing
# ----------------------------------------------------------------------
build_corridors_for_segment <- function(seg) {
    if (anyNA(
    c(seg$first_long,
      seg$first_lat,
      seg$last_long,
      seg$last_lat)
  )) {
    warning("Skipping segment ", seg$unique_seg, " – missing coordinates")
    return(NULL)
  }
  
  pts <- filter(points_df, unique_seg == seg$unique_seg)
  if (nrow(pts) < 2) return(NULL)
  
    cum_km_at <- function(t_abs) {
    approx(x = as.numeric(pts$timestamp),
           y = pts$cumulative_distance,
           xout = as.numeric(t_abs),
           rule = 2)$y
  }
  
  storm_time_abs <- seg$first_high_kp_time
  storm_km       <- cum_km_at(storm_time_abs)           
  storm_proj     <- geosphere::destPointRhumb(
    c(seg$first_long, seg$first_lat),
    seg$init_mean_heading,
    storm_km * 1000)                    
  storm_lat <- storm_proj[2]; storm_lon <- storm_proj[1]
  
  lcc_crs <- build_combined_lcc_crs(
    seg$first_lat, seg$first_long,
    seg$last_lat,  seg$last_long)
  
  static_poly <- corridor_static(
    start_lat = seg$first_lat, start_lon = seg$first_long,
    end_lat   = storm_lat,     end_lon   = storm_lon,
    heading_deg = seg$init_mean_heading,
    lcc_crs = lcc_crs)
  
  if (!is.null(static_poly)) {
    static_poly$part_type    <- "static"
    static_poly$polygon_time <- seg$first_timestamp
  }
  
  dur_hr   <- as.numeric(difftime(seg$last_timestamp,  seg$first_timestamp, units = "hours"))
  storm_hr <- as.numeric(difftime(storm_time_abs,      seg$first_timestamp, units = "hours"))
  
  dyn_list <- list()
  if (storm_hr < dur_hr) {
    
    t_seq <- seq(ceiling(storm_hr), dur_hr, by = 1)
    if (tail(t_seq, 1) < dur_hr) t_seq <- c(t_seq, dur_hr)
    
    final_static_w <- max(1.8 * storm_km * 1000 + 80000, 1.9 * storm_km * 1000)
    maxWidth    <- min(1.9 * seg$great_circle_distance * 1000, 1e9)
    margin_base <- max(seg$cumulative_distance / seg$duration * 2000, 200000)
    rangeTime   <- dur_hr - storm_hr
    
    dist_m_at <- function(t_hr) cum_km_at(seg$first_timestamp + hours(t_hr)) * 1000
    
    for (i in seq_len(length(t_seq) - 1)) {
      
      t1 <- t_seq[i]; t2 <- t_seq[i + 1]
      
      d1   <- dist_m_at(t1)
      frac <- if (rangeTime > 0) ((0.7 * (t1 + t2) - storm_hr) / rangeTime) else 1
      
      dyn_margin <- margin_base + 4 * frac * margin_base
      dist_min_m <- max(0, d1 - dyn_margin- 70000)
      dist_max_m <-  d1 + dyn_margin + 70000
      width_m    <- final_static_w + frac * (maxWidth - final_static_w)
      
      rect <- corridor_dynamic(
        start_lat  = seg$first_lat, start_lon = seg$first_long,
        heading_deg = seg$init_mean_heading,
        dist_min_m = dist_min_m, dist_max_m = dist_max_m,
        width_m    = width_m,
        lcc_crs    = lcc_crs)
      
      if (!is.null(rect)) {
        rect$part_type    <- "dynamic"
        rect$polygon_time <- seg$first_timestamp + hours(t1)
        dyn_list[[length(dyn_list) + 1]] <- rect
      }
    }
  }
  
  pieces <- Filter(Negate(is.null), c(list(static_poly), dyn_list))
  if (length(pieces) == 0) return(NULL)
  
  do.call(rbind, pieces) %>%
    mutate(
      segment_id                  = seg$segment_id,
      unique_seg                  = seg$unique_seg,
      individual.local.identifier = seg$individual.local.identifier
    ) %>%
    st_make_valid() %>%
    st_transform(4326)         
}

# ----------------------------------------------------------------------
# 6) Build all polygons
# ----------------------------------------------------------------------
polys_list <- lapply(
  split(segments, segments$unique_seg),
  build_corridors_for_segment
)
polys_list <- Filter(Negate(is.null), polys_list)
polys      <- do.call(rbind, polys_list)
cat("Polygon slices:", nrow(polys), "\n")

# ----------------------------------------------------------------------
# 7) Rasterise to 10 km grid
# ----------------------------------------------------------------------
gridify <- function(poly_sf, grid_size_m = 10000) {
  
  uni   <- st_union(poly_sf$geometry) %>% st_make_valid()
  bb    <- st_bbox(uni)
  lcc   <- build_combined_lcc_crs(bb["ymin"], bb["xmin"], bb["ymax"], bb["xmax"])
  uni_l <- st_transform(uni, lcc)
  
  grid  <- st_make_grid(st_as_sfc(st_bbox(uni_l)), cellsize = grid_size_m)
  sel   <- grid[st_intersects(st_sf(geometry = grid), uni_l, sparse = FALSE)[, 1]]
  ctr   <- st_centroid(sel) %>% st_transform(4326)
  xy    <- st_coordinates(ctr)
  
  data.frame(
    unique_seg                  = poly_sf$unique_seg[1],
    segment_id                  = poly_sf$segment_id[1],
    individual.local.identifier = poly_sf$individual.local.identifier[1],
    part_type                   = paste(unique(poly_sf$part_type), collapse = ","),
    cell_time                   = as.character(poly_sf$polygon_time[1]),
    grid_lon                    = xy[, 1],
    grid_lat                    = xy[, 2]
  )
}

points_out <- bind_rows(
  lapply(split(polys, seq_len(nrow(polys))), gridify)
)
cat("Grid cells:", nrow(points_out), "\n")

# ----------------------------------------------------------------------
# 8) Write output
# ----------------------------------------------------------------------
write.csv(points_out, "geomag_corridor_dynamic_hourly.csv", row.names = FALSE)
