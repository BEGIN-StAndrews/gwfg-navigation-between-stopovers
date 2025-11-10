# ==============================================================================
# Title: Unified route simulation – Magnetoclinic & Geomagnetic
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read resampled segments and hourly points (+ geomagnetic grid)
#   2) Interpolate declination D and inclination I (static/dynamic) safely
#   3) Simulate:
#        - Magnetoclinic (apparent dip solver, heading updated each step)
#        - Geomagnetic loxodrome (constant geomagnetic bearing)
#      with dynamic step lengths from hourly speeds
#   4) Export routes and failure logs
# Inputs:
#   - segments_resample.csv
#   - points_resampled.csv
#   - GeomagRaster_complete_final.csv   # columns: unique_seg, cell_time, grid_lon, grid_lat, D, I
# Outputs:
#   - Magnetoclinic_routes.csv
#   - Geomagnetic_routes.csv
#   - declination_failures_static.csv
#   - geomag_problem_segments_with_metadata_all.csv
# Usage:
#   Rscript R/09_2_Route_simulation_MAG_GMAG.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(geosphere)
library(RANN)
library(data.table)
library(purrr)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------
segments <- read.csv("segments_resample.csv", stringsAsFactors=FALSE) %>%
  mutate(
    first_timestamp     = as.POSIXct(first_timestamp,format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    first_high_kp_time  = as.POSIXct(first_high_kp_time,format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    first_lat           = as.numeric(first_lat),
    first_long          = as.numeric(first_long),
    last_lat            = as.numeric(last_lat),
    last_long           = as.numeric(last_long),
    duration            = as.numeric(duration),
    cumulative_distance = as.numeric(cumulative_distance),
    MaxKp               = as.numeric(MaxKp),
    unique_seg          = as.character(unique_seg) )

points_df <- read.csv("points_resampled.csv", stringsAsFactors=FALSE) %>%
  mutate(
    timestamp  = as.POSIXct(timestamp,format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    speed_kmh  = as.numeric(CalSpeed) * 3.6,
    unique_seg = as.character(unique_seg)
  ) %>%  arrange(unique_seg, timestamp)

GeoMag_data <- read.csv("GeomagRaster_complete_final.csv") %>%
  mutate(
    cell_time  = as.POSIXct(cell_time,format="%Y-%m-%d %H:%M:%S", tz="UTC"),
    unique_seg = as.character(unique_seg)
  ) %>% dplyr::select(cell_time, grid_lon, grid_lat, D, I, unique_seg)


# ------------------------------------------------------------------------------
# Failure log
# ------------------------------------------------------------------------------
decl_failures <- data.frame(
  unique_seg  = character(),
  query_lon   = numeric(),
  query_lat   = numeric(),
  message     = character(),
  stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Spatial / spatio-temporal interpolation (safe)
# ------------------------------------------------------------------------------
interp_space_safe <- function(lon, lat, sub_t, seg_id, field, k = 4, dist_max_km = 500) {
  sub_t <- sub_t[complete.cases(sub_t[,c("grid_lon","grid_lat")]), ] 
  if (nrow(sub_t) < 3) {
    msg <- sprintf("only %d rows for %s", nrow(sub_t), field)
    decl_failures <<- rbind(decl_failures,
                            data.frame(unique_seg = seg_id, query_lon = lon, query_lat = lat, message = msg,
                                       stringsAsFactors = FALSE))
    return(NA_real_)
  }
  
  nn <- RANN::nn2(sub_t[, c("grid_lon", "grid_lat")],query = matrix(c(lon, lat), 1, 2), k = min(k, nrow(sub_t)))
  nbrs <- sub_t[nn$nn.idx[1,], ]
  nbrs <- nbrs[!is.na(nbrs[[field]]), ]
  
  if (nrow(nbrs) < 3) {
    msg <- sprintf("only %d valid neighbors for %s", nrow(nbrs), field)
    decl_failures <<- rbind(decl_failures,
                            data.frame(unique_seg = seg_id, query_lon = lon, query_lat = lat, message = msg,
                                       stringsAsFactors = FALSE))
    return(NA_real_)
  }
    dists <- geosphere::distHaversine(
    matrix(c(lon, lat), 1, 2),
    as.matrix(nbrs[, c("grid_lon", "grid_lat")])
  ) / 1000
  if (max(dists, na.rm = TRUE) > dist_max_km) {
    msg <- sprintf("furthest neighbor %.1f km > %.1f km for %s", max(dists, na.rm = TRUE), dist_max_km, field)
    decl_failures <<- rbind(decl_failures,
                            data.frame(unique_seg = seg_id, query_lon = lon, query_lat = lat, message = msg,
                                       stringsAsFactors = FALSE))
    return(NA_real_)
  }
  fit <- lm(as.formula(paste(field, "~ grid_lon + grid_lat")), data = nbrs)
  predict(fit, newdata = data.frame(grid_lon = lon, grid_lat = lat))
}


interp_spacetime_safe <- function(lon, lat, query_time, decl_sub, seg_id, field,k = 4, dist_max_km = 500, time_tol = 31*60) {
  
  dt_secs <- abs(as.numeric(difftime(decl_sub$cell_time, query_time, units = "secs")))
  i_min <- which.min(dt_secs)
  min_dt <- dt_secs[i_min]
    if (min_dt > time_tol) {
    msg <- sprintf("nearest time Δt=%.0f min > %.0f min for %s", min_dt/60, time_tol/60, field)
    decl_failures <<- rbind(decl_failures,
                            data.frame(unique_seg = seg_id, query_lon = lon, query_lat = lat, message = msg,
                                       stringsAsFactors = FALSE))
    return(NA_real_)
  }
  
  t_near <- decl_sub$cell_time[i_min]
  slice <- decl_sub[decl_sub$cell_time == t_near, ]
  interp_space_safe(lon, lat, slice, seg_id, field, k, dist_max_km)
}

# ------------------------------------------------------------------------------
# Apparent dip heading solver (magnetoclinic)
# ------------------------------------------------------------------------------
solve_heading_for_apparent_dip <- function(A0, D_local, I_local, 
                                           H_center,
                                           span      = 90,
                                           grid_deg  = 0.2,
                                           tol_root  = 1e-7)
{
  f <- function(H){
    M_r <- ((H - D_local) %% 360) * pi/180
    I_r <-  I_local * pi/180
    asin( sin(I_r) / sqrt( sin(M_r)^2 * cos(I_r)^2 + sin(I_r)^2 ) )*180/pi - A0
  }
  
  lo <- H_center - span
  hi <- H_center + span
  cand <- numeric(0)
  root1 <- tryCatch(uniroot(f, c(lo, hi), extendInt="yes", tol=tol_root)$root,
                    error=function(e) NA_real_)
  if (!is.na(root1)) cand <- c(cand, root1)
  
  ang  <- seq(lo, hi, by=grid_deg)
  vals <- sapply(ang, f)
  idx  <- which(diff(sign(vals)) != 0)
  if (length(idx)){
    extra <- vapply(idx, function(i){
      tryCatch(uniroot(f, c(ang[i], ang[i+1]), tol=tol_root)$root,
               error=function(e) NA_real_)},
      numeric(1))
    cand <- c(cand, extra[!is.na(extra)])
  }
  
  if (length(cand)){
    cand <- cand %% 360
    dist <- abs(((cand - H_center + 180) %% 360) - 180)
    return(cand[which.min(dist)])
  }
  
  optimise(function(h) abs(f(h)), c(lo, hi))$minimum %% 360
}

# ------------------------------------------------------------------------------
# Steppers
# ------------------------------------------------------------------------------
# 6a) Magnetoclinic
step_magnetoclinic <- function(seg, beta0, times, step_km, decl_static, decl_dynamic) {
  lon <- seg$first_long
  lat <- seg$first_lat
    decl0 <- if (   is.na(seg$MaxKp) | seg$MaxKp <= 5     ) {   
    interp_space_safe(lon, lat, decl_static, seg$unique_seg, "D")
  } else {
    interp_spacetime_safe(lon, lat, times[1], decl_dynamic, seg$unique_seg, "D")
  }

  incl0 <- if (is.na(seg$MaxKp) | seg$MaxKp <= 5 ) {
    interp_space_safe(lon, lat, decl_static, seg$unique_seg, "I")
  } else {
    interp_spacetime_safe(lon, lat, times[1], decl_dynamic, seg$unique_seg, "I")
  }
  if (is.na(decl0) || is.na(incl0)) return(NULL)

  init_mag <- (beta0 - decl0 + 360) %% 360
  M0_r <- init_mag * pi/180
  I0_r <- incl0 * pi/180
  A0 <- asin(sin(I0_r) / sqrt(sin(M0_r)^2 * cos(I0_r)^2 + sin(I0_r)^2)) * 180/pi
  
  n <- length(step_km) + 1
  rows <- vector("list", n)
  rows[[1]] <- list(lon = lon, lat = lat, time = times[1])
  H_prev <- beta0
  
  for (i in seq_along(step_km)) {
    if (is.na(seg$MaxKp) | seg$MaxKp <= 5 ) {
      dloc <- interp_space_safe(lon, lat, decl_static, seg$unique_seg, "D")
      iloc <- interp_space_safe(lon, lat, decl_static, seg$unique_seg, "I")
    } else {
      ref_time <- if ( !is.na(seg$first_high_kp_time) && seg$first_high_kp_time <= times[i]) times[i] else times[1]
      dloc <- interp_spacetime_safe(lon, lat, ref_time, decl_dynamic, seg$unique_seg, "D")
      iloc <- interp_spacetime_safe(lon, lat, ref_time, decl_dynamic, seg$unique_seg, "I")
    }
    if (is.na(dloc) || is.na(iloc)) return(NULL)
    
    true_h <- solve_heading_for_apparent_dip(A0, dloc, iloc, H_prev)
    mv <- geosphere::destPointRhumb(c(lon, lat), true_h, step_km[i] * 1000)
    lon <- mv[1]
    lat <- mv[2]
    rows[[i + 1]] <- list(lon = lon, lat = lat, time = times[i + 1])
    H_prev <- true_h
  }
  
  trk <- data.table::rbindlist(rows)
  setnames(trk, c("lon", "lat", "time"))
  trk
}


# 6b) Geomagnetic Loxodrome
step_Geomag_loxodrome <- function(seg, beta0, times, step_km, decl_static, decl_dynamic) {
  lon <- seg$first_long; lat <- seg$first_lat
  decl0 <- if (is.na(seg$MaxKp) | seg$MaxKp <= 5) {
    interp_space_safe(lon, lat, decl_static, seg$unique_seg, "D")
  } else {
    interp_spacetime_safe(lon, lat, times[1], decl_dynamic, seg$unique_seg, "D")
  }
  
  if (is.na(decl0)) return(NULL)
  base_h <- (beta0 - decl0 + 360) %% 360
  n <- length(step_km) + 1
  rows <- vector("list", n)
  rows[[1]] <- list(lon=lon, lat=lat, time=times[1])
  
  for (i in seq_along(step_km)) {
      if (i == 1) {
      true_h <- beta0     
      } else {    
          if (is.na(seg$MaxKp) | seg$MaxKp <= 5 ) {
        dloc <- interp_space_safe(lon, lat, decl_static, seg$unique_seg, "D")
        if (is.na(dloc)) return(NULL)
        } else {
        ref_time <- if (seg$first_high_kp_time > times[i])  times[1] else times[i]
        dloc <- interp_spacetime_safe(lon, lat, ref_time,decl_dynamic, seg$unique_seg, "D")
        if (is.na(dloc)) return(NULL)
      }
      true_h <- (base_h + dloc) %% 360
      }
    mv  <- geosphere::destPointRhumb(c(lon,lat), true_h, step_km[i] * 1000)
    lon <- mv[1]; lat <- mv[2]
    rows[[i+1]] <- list(lon=lon, lat=lat, time=times[i+1])
  }
  trk <- data.table::rbindlist(rows)
  setnames(trk, c("lon","lat","time"))
  trk
}

# ------------------------------------------------------------------------------
# Simulation wrapper
# ------------------------------------------------------------------------------
simulate_path <- function(seg, heading_deg, heading_type,speed_mode, route_type){
  
  pts   <- filter(points_df, unique_seg == seg$unique_seg)
  times <- pts$timestamp
  dt_h  <- as.numeric(difftime(times[-1], times[-nrow(pts)], units = "hours"))
  
  step_km <- switch(speed_mode,
                    const = (seg$cumulative_distance / seg$duration) * dt_h,
                    dyn   = pmax(head(pts$speed_kmh, -1) * dt_h, 0))
  
  grid_static <- filter(GeoMag_data,unique_seg == seg$unique_seg, is.na(cell_time))
  grid_dyn    <- filter(GeoMag_data,unique_seg == seg$unique_seg,!is.na(cell_time))
  
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

# ------------------------------------------------------------------------------
# Heading set (use established initial mean heading)
# ------------------------------------------------------------------------------
get_heading_set <- function(seg) {
  c(ini100km = seg$init_mean_heading)
}

# ------------------------------------------------------------------------------
# Run simulations (routes × headings × speed modes) with failure filtering
# ------------------------------------------------------------------------------
route_types <- c("magnetoclinic","geomagnetic")
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
    if (any(decl_failures$unique_seg == seg$unique_seg)) {
      return(NULL)
    }
    sim_routes
  }
)

write.csv(filter(all_routes,route_type=="magnetoclinic"), "Magnetoclinic_routes.csv", row.names=FALSE)
write.csv(filter(all_routes,route_type=="geomagnetic"),   "Geomagnetic_routes.csv",   row.names=FALSE)

write.csv(decl_failures, "declination_failures_static.csv", row.names=FALSE)
