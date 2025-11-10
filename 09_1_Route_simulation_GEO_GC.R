# ==============================================================================
# Title: Route simulation – Geographic loxodrome & Great-circle
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read resampled segments and hourly points with kinematics
#   2) Simulate paths for:
#        - Geographic loxodrome (constant bearing)
#        - Great-circle (orthodrome)
#      using dynamic step lengths from hourly speeds
#   3) Export route tables by route type
# Inputs:
#   - segments_resample.csv
#   - points_resampled.csv
# Outputs:
#   - Geographic_routes.csv
#   - GreatCircle_routes.csv
# Usage:
#   Rscript R/09_1_Route_simulation_GEO_GC.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(geosphere)
library(RANN)
library(data.table)
library(lubridate)
library(suncalc)
library(purrr)
library(readr)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------
segments <- read.csv("segments_resample.csv", stringsAsFactors=FALSE) %>%
  mutate(
    first_timestamp     = as.POSIXct(first_timestamp, tz="UTC"),
    first_high_kp_time  = as.POSIXct(first_high_kp_time, tz="UTC"),
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
    timestamp  = as.POSIXct(timestamp, tz="UTC"),
    speed_kmh  = as.numeric(CalSpeed) * 3.6,
    unique_seg = as.character(unique_seg)
  ) %>% arrange(unique_seg, timestamp)

# ------------------------------------------------------------------------------
# Route steppers
# ------------------------------------------------------------------------------
# Geographic loxodrome (constant bearing, rhumb-line)
step_geog_loxodrome <- function(seg, beta0, times, step_km) {
  lon<-seg$first_long; lat<-seg$first_lat
  rows<-list(list(lon=lon,lat=lat,time=times[1]))
  for(i in seq_along(step_km)){
    mv  <- geosphere::destPointRhumb(c(lon,lat), beta0, step_km[i] * 1000)
    lon<-mv[1]; lat<-mv[2]
    rows[[i+1]]<-list(lon=lon,lat=lat,time=times[i+1])
  }
  trk<-rbindlist(rows);
  setnames(trk,c("lon","lat","time"))
  trk
}

# Great-circle (orthodrome) with N points (match observed count)
step_greatcircle <- function(seg, N) {
  gc_pts <- geosphere::gcIntermediate(
    p1            = c(seg$first_long, seg$first_lat),
    p2            = c(seg$last_long,  seg$last_lat),
    n             = N - 2,           # intervals = N-1 → points = N
    addStartEnd   = TRUE,
    breakAtDateLine = FALSE
  )
  data.table::data.table(
    lon = gc_pts[,1],
    lat = gc_pts[,2]
  )
}

# ------------------------------------------------------------------------------
# Path simulation wrapper
# ------------------------------------------------------------------------------
simulate_path <- function(seg, heading_deg, heading_type,speed_mode, route_type){
  
  pts   <- filter(points_df, unique_seg == seg$unique_seg)
  Npnt <- nrow( pts)
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
                great_circle = step_greatcircle(seg, Npnt),
                  
            
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

# ------------------------------------------------------------------------------
# Heading set
# ------------------------------------------------------------------------------
get_heading_set <- function(seg) {
  c(ini100km = seg$init_mean_heading)
}

# ------------------------------------------------------------------------------
# Run simulations (routes × headings × speed modes)
# ------------------------------------------------------------------------------
route_types <- c("geographic","great_circle")
speed_modes<- c("dyn")

all_routes <- purrr::map_dfr(
  split(segments, segments$unique_seg),
  function(seg) {
    heads <- get_heading_set(seg)
    sim_routes <- purrr::map_dfr(names(heads), function(ht) {
      deg <- heads[[ht]]
      if (is.na(deg)) return(NULL)
      purrr::map_dfr(route_types, function(rt) {
        purrr::map_dfr(speed_modes, function(sm) {
          simulate_path(seg, deg, ht, sm, rt)
        })
      })
    })

    sim_routes
  }
)

write.csv(filter(all_routes,route_type=="geographic"),   "Geographic_routes.csv",   row.names=FALSE)
write.csv(filter(all_routes, route_type=="great_circle"),"GreatCircle_routes.csv", row.names = FALSE)
