# ==============================================================================
# Title: Global wind-optimal route
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
# Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi | am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
# 1) Read segment table with endpoints + hourly resampled points
# 2) Build a fixed-width corridor around the great-circle (GC) line
# 3) Use a time-expanded DP to minimise travel time with ERA5 winds
# - Waits are inserted only where the empirical track did not move ≥ 5 km
# - Transitions constrained to a small lateral spread per slice
# - Ground speed = Va + wind-along-track (hourly, bilinear in space & time)
# 4) Backtrack optimal path, then expand in time to include explicit wait hours
# 5) Export raw nodes (1 per slice) + expanded hourly series
# Inputs:
# - segments_resample.csv
# - points_resampled.csv
# - ERA5_wind_<YEAR>.nc (longitude, latitude, valid_time, u100, v100)
# Outputs:
# - wind_optimal_routes_raw.csv (one node per slice)
# - wind_optimal_routes.csv (hourly-expanded with waits)
# Usage:
# Rscript R/09_5_Route_simulation_WO_DP.R
# Notes:
# - Wind data: ERA5 hourly single levels (u100, v100)
# - annual NetCDF files can be downloaded from the Copernicus Climate Data Store (CDS), DOI: 10.24381/cds.adbb2d47.

# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
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

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------
segments <- read.csv("segments_resample.csv", stringsAsFactors = FALSE) %>%
  mutate(
    first_timestamp = as.POSIXct(first_timestamp,format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"),
    unique_seg      = as.character(unique_seg),
    across(c(first_lat, first_long, last_lat, last_long), as.numeric)) 

segments <- segments %>%
  mutate(wind_year = lubridate::year(first_timestamp))  

points_df <- read.csv("points_resampled.csv", stringsAsFactors = FALSE) %>%
  mutate(
    timestamp  = as.POSIXct(timestamp,format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"),
    unique_seg = as.character(unique_seg),
    lat        = as.numeric(location.lat),
    lon        = as.numeric(location.long)
  ) %>%  arrange(unique_seg, timestamp)

# ------------------------------------------------------------------------------
# Corridor parameters (fixed width based on GC length; constant lateral step)
# ------------------------------------------------------------------------------
choose_grid_params <- function(L_seg) {
  pad_width_km    <- L_seg / 4
  lateral_step_km <- 10
  half_steps      <- pad_width_km / lateral_step_km
  M_lateral       <- 2 * ceiling(half_steps) + 1
  spread          <- floor(M_lateral / 4)
  list(
    pad_width_km = pad_width_km,
    M_lateral    = M_lateral,
    spread       = spread
  )
}

# ------------------------------------------------------------------------------
# ERA5 cache (one open handle per year); interpolation helpers
# ------------------------------------------------------------------------------
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

find_corners <- function(lon_pt, lat_pt) {
  ix_raw <- floor((lon_pt - all_lons[1]) / cell_d) + 1
  ix     <- pmin(pmax(ix_raw, 1), length(all_lons) - 1)
  iy_raw <- floor((all_lats[1] - lat_pt) / cell_d) + 1
  iy     <- pmin(pmax(iy_raw, 1), length(all_lats) - 1)
  list(ix = ix, iy = iy)
}

wind_uv <- function(lon, lat, t_posix) {
  cr <- find_corners(lon, lat)
  fx <- (lon - all_lons[cr$ix]) / cell_d
  fy <- (all_lats[cr$iy] - lat) / cell_d
  k  <- findInterval(t_posix, valid_time_posix)
  k  <- pmax(pmin(k, length(valid_time_posix) - 1), 1)
  k2 <- k + 1
  dt_secs <- as.numeric(
    difftime(valid_time_posix[k2], valid_time_posix[k], units = "secs")
  )
  a <- if (dt_secs == 0) 0 else
    as.numeric(difftime(t_posix, valid_time_posix[k], units = "secs")) / dt_secs
  a <- min(max(a, 0), 1)
  
  read_slice <- function(var, idx) {
    ncvar_get(nc_data, var,
              start = c(cr$ix, cr$iy, idx),
              count = c(2,      2,      1))
  }
  U <- (1 - a) * read_slice("u100", k)  +  a * read_slice("u100", k2)
  V <- (1 - a) * read_slice("v100", k)  +  a * read_slice("v100", k2)
  
  u <- (1 - fx)*(1 - fy)*U[1,1] + fx*(1 - fy)*U[2,1] +
    (1 - fx)*fy      *U[1,2] + fx*fy       *U[2,2]
  v <- (1 - fx)*(1 - fy)*V[1,1] + fx*(1 - fy)*V[2,1] +
    (1 - fx)*fy      *V[1,2] + fx*fy       *V[2,2]
  
  c(u = u, v = v)
}

ground_time_h <- function(lonA, latA, lonB, latB, u, v) {
  d_m <- distHaversine(c(lonA, latA), c(lonB, latB))
  if (d_m < 1) return(0)
  brg        <- bearing(c(lonA, latA), c(lonB, latB)) * pi/180
  wind_along <- u * sin(brg) + v * cos(brg)
  gs_ms      <- Va_ms + wind_along
  if (gs_ms < 0.1) return(NA_real_)    # NA instead of Inf
  (d_m / gs_ms) / 3600
}

make_slice_vertices <- function(lon_c, lat_c, pad_km, M, main_brg_deg) {
  perp <- (main_brg_deg + 90) %% 360
  tibble(j_index = 1:M) %>%
    mutate(
      dist_km = seq(-pad_km, pad_km, length.out = M),
      tmp     = geosphere::destPoint(c(lon_c, lat_c), perp, dist_km * 1000),
      lon     = tmp[,1], lat = tmp[,2]
    ) %>%
    select(-dist_km, -tmp)
}

# ------------------------------------------------------------------------------
# DP solver
# ------------------------------------------------------------------------------
compute_dp_route <- function(seg_row, vertices_by_slice,
                             wait_hours,
                             Va_ms,
                             M_lateral, spread) {
  
  L          <- length(vertices_by_slice)
  centre_idx <- (M_lateral + 1) %/% 2
  
  cost   <- lapply(vertices_by_slice, function(v) rep(Inf, nrow(v)))
  arr_t  <- lapply(vertices_by_slice,
                   function(v) rep(as.POSIXct(NA, tz="UTC"), nrow(v)))
  prev_j <- lapply(vertices_by_slice, function(v) rep(NA_integer_, nrow(v)))
  
  cost[[1]][1]  <- 0
  arr_t[[1]][1] <- seg_row$first_timestamp
  
  for (i in 2:L) {
    src      <- vertices_by_slice[[i-1]]
    trg      <- vertices_by_slice[[i]]
    eff_src  <- if (nrow(src)>1) src$j_index else rep(centre_idx, nrow(src))
    eff_trg  <- if (nrow(trg)>1) trg$j_index else rep(centre_idx, nrow(trg))
    
    for (j in seq_len(nrow(trg))) {
      best   <- Inf; best_t <- NA; best_k <- NA
      valid_k <- which(
        (nrow(src)==1 | abs(eff_src-eff_trg[j]) <= spread) &
          is.finite(cost[[i-1]])
      )
      for (k in valid_k) {
        t0         <- arr_t[[i-1]][k]
        wait_h     <- wait_hours[i-1]
        t_depart   <- t0 + wait_h * 3600
        base_cost  <- cost[[i-1]][k] + wait_h
        
        uv0 <- wind_uv(src$lon[k], src$lat[k], t_depart)
        h0  <- ground_time_h(src$lon[k], src$lat[k],
                             trg$lon[j], trg$lat[j],
                             uv0["u"], uv0["v"])
        if (!is.finite(h0)) next
        
        
        lon_mid <- (src$lon[k] + trg$lon[j]) / 2
        lat_mid <- (src$lat[k] + trg$lat[j]) / 2
        t_mid   <- t_depart + h0*3600/2
        uv_mid  <- wind_uv(lon_mid, lat_mid, t_mid)
        
        h    <- ground_time_h(src$lon[k], src$lat[k],
                              trg$lon[j], trg$lat[j],
                              uv_mid["u"], uv_mid["v"])
        if (!is.finite(h)) next
        
        cand_cost <- base_cost + h
        t_arr     <- t_depart + h*3600
        if (cand_cost < best) {
          best   <- cand_cost
          best_t <- t_arr
          best_k <- k
        }
      }
      
      cost[[i]][j]   <- best
      arr_t[[i]][j]  <- best_t
      prev_j[[i]][j] <- best_k
    }
  }
  
  if (is.infinite(cost[[L]][1])) stop("DP failed for segment ", seg_row$unique_seg)
  
  # Preallocated backtrack (O(L) instead of c())
  layers <- integer(L)
  idxs   <- integer(L)
  times  <- as.POSIXct(rep(NA, L), tz="UTC")
  cur_i  <- L; cur_j <- 1; pos <- L
  while (pos >= 1) {
    layers[pos] <- cur_i - 1
    idxs[pos]   <- cur_j
    times[pos]  <- arr_t[[cur_i]][cur_j]
    if (cur_i == 1) break
    tmp_j   <- prev_j[[cur_i]][cur_j]
    cur_i   <- cur_i - 1
    cur_j   <- tmp_j
    pos     <- pos - 1
  }
  ord <- seq_len(L)
  
  tibble(
    layer = layers[ord],
    lon   = map2_dbl(layer, idxs[ord], ~ vertices_by_slice[[.x+1]]$lon[.y]),
    lat   = map2_dbl(layer, idxs[ord], ~ vertices_by_slice[[.x+1]]$lat[.y]),
    time  = times[ord]
  )
}

# ------------------------------------------------------------------------------
# Prepare segment (slices + waits)
# ------------------------------------------------------------------------------
prepare_segment <- function(seg, points_df, D_min_km = 5) {
  pts   <- filter(points_df, unique_seg == seg$unique_seg)
  d_km  <- distHaversine(
    cbind(pts$lon[-1], pts$lat[-1]),
    cbind(pts$lon[-nrow(pts)], pts$lat[-nrow(pts)])
  ) / 1000
  
    L_seg <- distHaversine(
    c(seg$first_long, seg$first_lat),
    c(seg$last_long,  seg$last_lat)
  ) / 1000
  
  gc_pts <- gcIntermediate(
    c(seg$first_long, seg$first_lat),
    c(seg$last_long,  seg$last_lat),
    n             = ceiling(L_seg) + 1,
    addStartEnd   = TRUE,
    breakAtDateLine = FALSE
  )
  
  idx_gc     <- 1
  centres    <- list(gc_pts[1,])
  wait_hours <- 0
  pending    <- 0
  
  for (step in d_km) {
    if (step < D_min_km) {
      pending <- pending + 1
      next
    }
    wait_hours[length(wait_hours)] <- wait_hours[length(wait_hours)] + pending
    pending <- 0
    
    rem_km <- nrow(gc_pts) - idx_gc
    adv    <- min(floor(step), rem_km)
    idx_gc <- idx_gc + adv
    
    centres      <- append(centres, list(gc_pts[idx_gc,]))
    wait_hours   <- c(wait_hours, 0)
    if (idx_gc == nrow(gc_pts)) break
  }
  wait_hours <- as.numeric(wait_hours)
  
  grid    <- choose_grid_params(L_seg)
  main_brg <- bearing(
    c(seg$first_long, seg$first_lat),
    c(seg$last_long,  seg$last_lat)
  )
  
  vertices_by_slice <- map2(
    map_dbl(centres, 1), map_dbl(centres, 2),
    make_slice_vertices,
    pad_km       = grid$pad_width_km,
    M            = grid$M_lateral,
    main_brg_deg = main_brg
  )
  
  vertices_by_slice[[1]]             <- tibble(
    j_index = 1,
    lon     = seg$first_long,
    lat     = seg$first_lat
  )
  vertices_by_slice[[length(centres)]] <- tibble(
    j_index = 1,
    lon     = seg$last_long,
    lat     = seg$last_lat
  )
  
  list(
    wait_hours        = wait_hours,
    vertices_by_slice = vertices_by_slice,
    grid              = grid
  )
}

# ------------------------------------------------------------------------------
# Simulate one wind-optimal route
# ------------------------------------------------------------------------------
simulate_wind_route <- function(seg, prep) {
  dp_raw <- compute_dp_route(
    seg_row           = seg,
    vertices_by_slice = prep$vertices_by_slice,
    wait_hours        = prep$wait_hours,
    Va_ms             = Va_ms,
    M_lateral         = prep$grid$M_lateral,
    spread            = max(prep$grid$spread, 4)
  ) %>%
    mutate(
      unique_seg   = seg$unique_seg,
      heading_init = bearing(c(seg$first_long, seg$first_lat),
                             c(seg$last_long,  seg$last_lat)),
      route_type   = "wind_optimal",
      speed_mode   = "dyn",
      heading_type = "wind_optimal",
      id_year      = seg$id_year,
      season       = seg$season,
      segment_id   = seg$segment_id
    )
  
  dp_2 <- dp_raw %>%
    mutate(wait = prep$wait_hours[layer + 1L]) %>%
    mutate(wait = replace(wait, row_number() == n(), 0L)) %>%
    mutate(times = map2(time, wait, ~ .x + hours(0:.y))) %>%
    unnest(times) %>%
    select(-time, -wait) %>%
    rename(time = times)
  
  list(raw = dp_raw, expanded = dp_2)
}

# ------------------------------------------------------------------------------
# Run all segments & save results
# ------------------------------------------------------------------------------
n <- nrow(segments)
raw_list      <- vector("list", n)
expanded_list <- vector("list", n)
n_total <- length(split(segments, segments$unique_seg))

i <- 1
for (seg in split(segments, segments$unique_seg)) {
  
  segments_left <- n_total - i
  cat(sprintf("[Iteration %3d/%3d] → %3d segments left\n",
              i, n_total, segments_left))
  
    meta             <- get_era5_for_year(seg$wind_year)
  nc_data          <- meta$nc            # <- visible to wind_uv()
  all_lons         <- meta$lons
  all_lats         <- meta$lats
  valid_time_posix <- meta$vt_posix

  prep <- prepare_segment(seg, points_df, D_min_km)
  out  <- try(simulate_wind_route(seg, prep), silent = TRUE)
  if (inherits(out, "try-error")) {
    message("Segment ", seg$unique_seg, " failed: ", out)
    next
  }
  raw_list[[i]]      <- out$raw
  expanded_list[[i]] <- out$expanded
  i <- i + 1
}

raw_tbl      <- bind_rows(raw_list)
expanded_tbl <- bind_rows(expanded_list)

write.csv(raw_tbl,      "wind_optimal_routes_raw.csv",      row.names = FALSE)
write.csv(expanded_tbl, "wind_optimal_routes.csv", row.names = FALSE)
