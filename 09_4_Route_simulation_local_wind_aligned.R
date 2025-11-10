# ==============================================================================
# Title: Local wind-aligned route simulation (matched length, free goal)
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read segments (with initial headings) and hourly-resampled points
#   2) Open ERA5 wind per-year via a small NetCDF handle cache
#   3) Simulate Local wind-aligned route: each hour follow downwind bearing
#   4) Export all simulated tracks to local_wind_routes.csv
# Inputs:
#   - segments_resample.csv
#   - points_resampled.csv
#   - ERA5_wind_<YEAR>.nc   (vars: longitude, latitude, valid_time, u100, v100)
# Outputs:
#   - local_wind_routes.csv
# Usage:
#   Rscript R/09_4_Route_simulation_local_wind_aligned.R
#Note:
# - Wind data: ERA5 hourly single levels (u100, v100)
# - annual NetCDF files can be downloaded from the Copernicus Climate Data Store (CDS), DOI: 10.24381/cds.adbb2d47.

# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(purrr)
library(geosphere)
library(ncdf4)
library(lubridate)
library(tibble)
library(tidyr)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
Va_ms <- 15 # Bird air-speed (m/s)
cell_d <- 0.25 # ERA5 grid spacing (degrees)
D_min_km <- 5 # Minimum km to count as a flight hour

# ----------------------------------------------------------------------------
# Load DATA
# ----------------------------------------------------------------------------
segments <- read.csv("segments_resample.csv", stringsAsFactors = FALSE) %>%
  mutate(
    first_timestamp = as.POSIXct(first_timestamp,format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"),
    unique_seg      = as.character(unique_seg),
    across(c(first_lat, first_long, last_lat, last_long), as.numeric)) 
segments <- segments %>%
  mutate(wind_year = lubridate::year(first_timestamp))

points_df <- read.csv("points_resampled.csv", stringsAsFactors = FALSE) %>%
  mutate(
    timestamp  = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS",tz = "UTC"),
    unique_seg = as.character(unique_seg),
    lat        = as.numeric(location.lat),
    lon        = as.numeric(location.long)
  ) %>% arrange(unique_seg, timestamp)

# ------------------------------------------------------------------------------
# Define Custom Functions
# ------------------------------------------------------------------------------
find_corners <- function(lon_pt, lat_pt) {
  ix_raw <- floor((lon_pt - all_lons[1]) / cell_d) + 1
  ix     <- pmin(pmax(ix_raw, 1), length(all_lons) - 1)
  iy_raw <- floor((all_lats[1] - lat_pt) / cell_d) + 1
  iy     <- pmin(pmax(iy_raw, 1), length(all_lats) - 1)
  list(ix = ix, iy = iy)
}

nc_cache <- new.env(parent = emptyenv())

get_era5_for_year <- function(yr) {
  key <- as.character(yr)
  
  if (!exists(key, envir = nc_cache, inherits = FALSE)) {
    fpath <- sprintf("ERA5_wind_%d.nc", yr)
    if (!file.exists(fpath))
      stop("Cannot find ERA-5 file for year ", yr, "  (expected: ", fpath, ")")
    nc_obj <- nc_open(fpath)
    meta <- list(
      nc      = nc_obj,
      lons    = ncvar_get(nc_obj, "longitude"),
      lats    = ncvar_get(nc_obj, "latitude"),
      vt_posix= as.POSIXct(
        ncvar_get(nc_obj, "valid_time"),
        origin = "1970-01-01", tz = "UTC"
      )
    )
    assign(key, meta, envir = nc_cache)
  }
  get(key, envir = nc_cache, inherits = FALSE)
}

close_all_era5 <- function() {
  invisible(lapply(as.list(nc_cache),
                   function(x) try(nc_close(x$nc), silent = TRUE)))
}
on.exit(close_all_era5(), add = TRUE)


wind_uv <- function(lon, lat, t_posix) {
  cr <- find_corners(lon, lat)
  fx <- (lon - all_lons[cr$ix]) / cell_d
  fy <- (all_lats[cr$iy]   - lat)     / cell_d
  
  k  <- findInterval(t_posix, valid_time_posix)
  k  <- pmax(pmin(k, length(valid_time_posix) - 1), 1)
  k2 <- k + 1
  dt_secs <- as.numeric(difftime(valid_time_posix[k2], valid_time_posix[k], units="secs"))
  a <- if (dt_secs == 0) 0 else as.numeric(difftime(t_posix, valid_time_posix[k], "secs"))/dt_secs
  a <- min(max(a,0),1)
  read_slice <- function(var, idx) {
    ncvar_get(nc_data, var,
              start = c(cr$ix, cr$iy, idx),
              count = c(2,      2,      1))
  }
  
  U <- (1 - a) * read_slice("u100", k)  +  a * read_slice("u100", k2)
  V <- (1 - a) * read_slice("v100", k)  +  a * read_slice("v100", k2)
  
  # bilinear interpolation
  u <- (1-fx)*(1-fy)*U[1,1] + fx*(1-fy)*U[2,1] +
    (1-fx)*fy      *U[1,2] + fx*fy       *U[2,2]
  v <- (1-fx)*(1-fy)*V[1,1] + fx*(1-fy)*V[2,1] +
    (1-fx)*fy      *V[1,2] + fx*fy       *V[2,2]
  
  c(u = u, v = v)
}

# -- One-hour step: Local wind-aligned (downwind; matched length) ---------------
step_tailwind_fixed <- function(lon, lat, t_posix, dist_m, prev_brg_deg) {
  uv <- wind_uv(lon, lat, t_posix)
  if (sqrt(sum(uv^2)) < 1e-1) {
    brg <- prev_brg_deg
  } else {
    brg <- atan2(uv["u"], uv["v"]) * 180 / pi
  }
  dest <- destPoint(c(lon, lat), brg, dist_m)
  list(lon = dest[1], lat = dest[2], brg = brg %% 360)
}

# -- Segment simulator ---------------------------------------------------------
simulate_tailwind_matched <- function(seg) {
  pts <- points_df %>% filter(unique_seg == seg$unique_seg) %>% arrange(timestamp)
  step_km  <- distHaversine(
    cbind(pts$lon[-1], pts$lat[-1]),
    cbind(pts$lon[-nrow(pts)], pts$lat[-nrow(pts)])
  ) / 1000
  
  fly_mask <- step_km >= D_min_km
  n <- nrow(pts)
  track <- tibble(lon = numeric(n), lat = numeric(n), time = pts$timestamp, brg = NA_real_)
  track$lon[1] <- seg$first_long
  track$lat[1] <- seg$first_lat
  prev_brg <- seg$init_mean_heading
  
  for (i in 2:n) {
    if (fly_mask[i-1]) {
      st <- step_tailwind_fixed(track$lon[i-1], track$lat[i-1], track$time[i-1],
                                dist_m = step_km[i-1]*1000, prev_brg_deg = prev_brg)
      track$lon[i] <- st$lon
      track$lat[i] <- st$lat
      track$brg[i] <- st$brg
      prev_brg      <- st$brg
    } else {
      track$lon[i] <- track$lon[i-1]
      track$lat[i] <- track$lat[i-1]
    }
  }
  
  track %>% mutate(
    unique_seg   = seg$unique_seg,
    route_type   = "max_tailwind_matched",
    speed_mode   = "dyn",
    heading_type = "tailwind_matched",
    id_year      = seg$id_year,
    season       = seg$season,
    segment_id   = seg$segment_id
  )
}

# ------------------------------------------------------------------------------
# Run All Segments & Save
# ------------------------------------------------------------------------------
N <- nrow(segments)
routes_matched <- vector("list", N)

for (i in seq_len(N)) {
  seg <- segments[i, ]
  meta             <- get_era5_for_year(seg$wind_year)
  nc_data          <<- meta$nc
  all_lons         <<- meta$lons
  all_lats         <<- meta$lats
  valid_time_posix <<- meta$vt_posix
  routes_matched[[i]] <- simulate_tailwind_matched(seg)
  
}

bind_rows(routes_matched) %>% write.csv("local_wind_routes.csv", row.names = FALSE)

