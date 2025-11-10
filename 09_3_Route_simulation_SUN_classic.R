# ==============================================================================
# Title: Route simulation – Time-compensated Sun-Compass
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read resampled segments and hourly points with kinematics
#   2) Implement robust Sun-Compass (classic mode only):
#        • Handles polar day/night and night departures
#        • Uses one calibration event/instant per segment
#   3) Simulate routes with dynamic step lengths from hourly speeds
#   4) Export classic route table
# Inputs:
#   - segments_resample.csv
#   - points_resampled.csv
# Outputs:
#   - SunCompass_classic_routes.csv
# Usage:
#   Rscript R/09_3_Route_simulation_SUN_classic.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(lubridate)
library(geosphere)
library(suncalc)
library(data.table)
library(purrr)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------
segments <- read.csv("segments_resample.csv", stringsAsFactors=FALSE) %>%
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
    unique_seg          = as.character(unique_seg))

points_df <- read.csv("points_resampled.csv", stringsAsFactors=FALSE) %>%
  mutate(
    timestamp  = as.POSIXct(timestamp,format = "%Y-%m-%d %H:%M:%OS", tz="UTC"),
    speed_kmh  = as.numeric(CalSpeed) * 3.6,
    unique_seg = as.character(unique_seg)
  ) %>%  arrange(unique_seg, timestamp)

# ==============================================================================
# Helper functions
# ==============================================================================
sun_az_deg <- function(lat, lon, t){
  a <- getSunlightPosition(t, lat, lon)$azimuth
  if (is.na(a)) return(NA_real_)
  if (a <  pi) return(a * 180/pi + 180)
  if (a >  pi) return(a * 180/pi - 180)
  0
}

sun_speed_deg <- function(lat, lon, t){
  th1 <- sun_az_deg(lat, lon, t - 30)
  th2 <- sun_az_deg(lat, lon, t + 30)
  if (is.na(th1) | is.na(th2)) return(15 * sin(abs(lat) * pi/180))
  ((th2 - th1 + 540) %% 360 - 180) * 60       
}

circ_diff_deg <- function(a,b) ((a - b + 540) %% 360) - 180


safe_event_az_deg <- function(lat, lon, t, event = c("sunset","sunrise")){
  event <- match.arg(event)
  st <- getSunlightTimes(as_date(t), lat, lon,
                         keep = c(event,"solarNoon","nadir"))
  ev_time <- st[[event]]
  if (!is.na(ev_time))                          
    return(sun_az_deg(lat, lon, ev_time))
  
  alt_noon <- getSunlightPosition(st$solarNoon, lat, lon)$altitude
  if (alt_noon > 0)                             
    sun_az_deg(lat, lon, st$nadir)              
  else                                        
    sun_az_deg(lat, lon, st$solarNoon)          
}

last_sunrise_before <- function(lat, lon, t){
  for (d in c(0,-1)){                                 
    sr <- getSunlightTimes(as_date(t)+d, lat, lon,
                           keep="sunrise")$sunrise
    if (!is.na(sr) && sr <= t) return(sr)
  }
  NA                                                    
}

choose_cal_event <- function(lat, lon, t){
  st <- getSunlightTimes(as_date(t), lat, lon,
                         keep=c("sunrise","sunset","solarNoon","nadir"))
  
  if (!is.na(st$sunrise) && !is.na(st$sunset)) return("sunset")
  if (!is.na(st$sunrise))                         return("sunrise")
  if (getSunlightPosition(st$solarNoon, lat, lon)$altitude > 0)
    return("nadir")                                             # polar day
  "solarNoon"                                                      # polar night
}

cal_event_az <- function(lat, lon, date, cal_event){
  ev_time <- getSunlightTimes(date, lat, lon, keep = cal_event)[[cal_event]]
  if (is.na(ev_time))                               # event missing here
    ev_time <- attr(date, "ref_instant_utc")        # use ref-site instant
  sun_az_deg(lat, lon, ev_time)
}

## ----------------------------------------------------------------------
# Time-Compensated Sun-Compass stepper
## ----------------------------------------------------------------------
step_suncompass <- function(seg, beta0, times, step_km,mode = c("classic","proximate")){
  mode <- match.arg(mode)
  t0 <- times[1]
  if (getSunlightPosition(t0, seg$first_lat, seg$first_long)$altitude <
      -6 * pi/180){
    t_ref <- last_sunrise_before(seg$first_lat, seg$first_long, t0)
    if (is.na(t_ref)) return(NULL)                 
  } else {
    t_ref <- t0
  }
  
  cal_event <- choose_cal_event(seg$first_lat, seg$first_long, t_ref)
  ref_ev_time <- getSunlightTimes(as_date(t_ref),
                                  seg$first_lat, seg$first_long,
                                  keep = cal_event)[[cal_event]]
  lat_now  <- seg$first_lat
  lon_now  <- seg$first_long
  alpha_c  <- beta0                     
  lon_c    <- lon_now                   
  lat_sref <- lat_now                 
  omega_ref<- sun_speed_deg(lat_sref, lon_c, t_ref)
  n <- length(step_km) + 1
  rows <- vector("list", n)
  rows[[1]] <- list(lon = lon_now, lat = lat_now, time = t0)
  
  for (i in seq_along(step_km)){
    t_now <- times[i]
    date_i <- structure(as_date(t_now), ref_instant_utc = ref_ev_time)
    th_h <- cal_event_az(lat_now, lon_now, date_i, cal_event)
    th_r <- cal_event_az(lat_sref, lon_c , date_i, cal_event)
    if (is.na(th_h) | is.na(th_r)) return(NULL)     # should never happen
    
    dlon  <- (lon_now - lon_c + 540) %% 360 - 180
    omega <- if (mode == "proximate")
      sun_speed_deg(lat_now, lon_now, t_now) else omega_ref
    alpha_now <- (alpha_c +
                    circ_diff_deg(th_h, th_r) +
                    (dlon/15)*omega) %% 360
    
    mv <- destPointRhumb(c(lon_now, lat_now), alpha_now, step_km[i]*1000)
    lon_now <- mv[1];  lat_now <- mv[2]
    rows[[i+1]] <- list(lon = lon_now, lat = lat_now, time = times[i+1])
  }
  
  data.table::rbindlist(rows)
}

# ==============================================================================
# Simulation wrapper 
# ==============================================================================
simulate_path <- function(seg, heading_deg, heading_type,speed_mode, route_type){
  
  pts   <- filter(points_df, unique_seg == seg$unique_seg)
  times <- pts$timestamp
  dt_h  <- as.numeric(difftime(times[-1], times[-nrow(pts)], units = "hours"))
  
  step_km <- switch(speed_mode,
                    const = (seg$cumulative_distance / seg$duration) * dt_h,
                    dyn   = pmax(head(pts$speed_kmh, -1) * dt_h, 0))
  
   trk <- switch(route_type,
                magnetoclinic = step_magnetoclinic(seg, heading_deg, times, step_km, grid_static, grid_dyn),
                geomagnetic   = step_Geomag_loxodrome(seg, heading_deg, times, step_km,grid_static, grid_dyn),
                geographic    = step_geog_loxodrome(seg, heading_deg, times, step_km),
                sun_classic   = step_suncompass(seg, heading_deg, times, step_km,mode = "classic"),
                sun_proximate = step_suncompass(seg, heading_deg, times, step_km,mode = "proximate"),
                stop("Unknown route_type")
  )
  
  if (is.null(trk)) return(NULL)
  
  trk %>%
    mutate(unique_seg   = seg$unique_seg,
           heading_init = heading_deg,
           heading_type = heading_type,
           speed_mode   = speed_mode,
           route_type   = route_type,          
           id_year      = seg$id_year,
           season       = seg$season,
           segment_id   = seg$segment_id)
}

# Heading set (your established initial heading)
get_heading_set <- function(seg) {
  c(ini100km = seg$init_mean_heading)
}

# ==============================================================================
# Run simulations (classic × headings × speed modes)
# ==============================================================================
route_types <- c("sun_classic")
spead_modes<- c("dyn")

all_routes <- purrr::map_dfr(
  split(segments, segments$unique_seg),
  function(seg) {
    heads <- get_heading_set(seg)
    sim_routes <- purrr::map_dfr(names(heads), function(ht) {
      deg <- heads[[ht]]
      if (is.na(deg)) return(NULL)
      purrr::map_dfr(route_types, function(rt) {
        purrr::map_dfr(spead_modes, function(sm) {
          simulate_path(seg, deg, ht, sm, rt)
        })
      })
    })
    sim_routes
  }
)

write.csv(filter(all_routes, route_type=="sun_classic"),   "SunCompass_classic_routes.csv", row.names=FALSE)
