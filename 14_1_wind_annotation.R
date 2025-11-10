# ==============================================================================
# Title: ERA5 wind annotation + segment-level wind metrics
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   Part A — ERA5 wind annotation
#     1) Open yearly ERA5 NetCDF files with a small in-memory cache
#     2) For each GPS fix, do bilinear interpolation in space + linear in time
#     3) Append u100/v100 as u100_interp/v100_interp to the points table
#     4) Save annotated points ordered by unique_seg and timestamp
#
#   Part B — Segment wind metrics summary
#     5) Flag flight vs stationary fixes (≥5 km displacement)
#     6) Compute per-fix wind_speed, wind_support, crosswind
#     7) Summarise mean metrics per unique_seg
#     8) Save segment summary and augmented point-level CSVs
#
# Inputs:
#   - points_resampled.csv
#   - ERA5_wind_<YYYY>.nc  (longitude, latitude, valid_time, u100, v100)
# Outputs:
#   - points_resampled_annotated.csv
#   - segments_wind_summary.csv
#   - points_wind_metrics.csv
# Usage:
#   Rscript R/14_1_wind_annotation.R
# Notes:
#   - Grid step assumed 0.25° (ERA5-HRES single levels); timestamps are UTC
#   - Ref: Hersbach et al. 2020, QJRMS (ERA5 reanalysis)
# ==============================================================================

# ==============================================================================
# Part A — ERA5 wind annotation
# ==============================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------
pkgs <- c("dplyr", "ncdf4", "lubridate", "progress", "purrr")
invisible(lapply(
  pkgs,
  \(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
))
lapply(pkgs, library, character.only = TRUE)

# ------------------------------------------------------------------------------
# Set parameters
# ------------------------------------------------------------------------------
cell_d      <- 0.25                 # ERA5 grid step (°)
vars_needed <- c("u100", "v100")
outfile     <- "points_resampled_annotated.csv"

# ------------------------------------------------------------------------------
# Build NetCDF cache
# ------------------------------------------------------------------------------
nc_cache <- new.env(parent = emptyenv())

get_era5_for_year <- function(yr) {
  key <- as.character(yr)
  if (!exists(key, envir = nc_cache, inherits = FALSE)) {
    f <- sprintf("ERA5_wind_%d.nc", yr)
    if (!file.exists(f))
      stop("ERA‑5 file ", f, " not found")
    nc <- nc_open(f)
    meta <- list(
      nc       = nc,
      lons     = ncvar_get(nc, "longitude"),
      lats     = ncvar_get(nc, "latitude"),
      vt_posix = as.POSIXct(ncvar_get(nc, "valid_time"),
                            origin = "1970-01-01", tz = "UTC")
    )
    assign(key, meta, envir = nc_cache)
  }
  get(key, envir = nc_cache, inherits = FALSE)
}

on.exit({
  invisible(lapply(as.list(nc_cache),
                   \(m) try(nc_close(m$nc), silent = TRUE)))
}, add = TRUE)

# ------------------------------------------------------------------------------
# Define interpolation (bilinear space + linear time)
# ------------------------------------------------------------------------------
bilinear_uv <- function(lon_pt, lat_pt, t_pt, meta) {
  lons <- meta$lons; lats <- meta$lats; vt <- meta$vt_posix; nc <- meta$nc
  ix_raw <- floor((lon_pt - lons[1]) / cell_d) + 1
  iy_raw <- floor((lats[1] - lat_pt) / cell_d) + 1
  ix <- pmin(pmax(ix_raw, 1), length(lons) - 1)
  iy <- pmin(pmax(iy_raw, 1), length(lats) - 1)
  fx <- (lon_pt - lons[ix])   / cell_d
  fy <- (lats[iy]  - lat_pt)  / cell_d
  k1 <- findInterval(t_pt, vt)
  k1 <- pmin(pmax(k1, 1), length(vt) - 1)
  k2 <- k1 + 1
  dt_tot <- as.numeric(difftime(vt[k2], vt[k1], units = "secs"))
  a <- if (dt_tot == 0) 0 else
    as.numeric(difftime(t_pt, vt[k1], units = "secs")) / dt_tot
  a <- min(max(a, 0), 1)
  slice <- function(var, k) {
    ncvar_get(nc, var,
              start = c(ix, iy, k),
              count = c(2,  2,  1))
  }
  
  U <- (1 - a) * slice("u100", k1) + a * slice("u100", k2)
  V <- (1 - a) * slice("v100", k1) + a * slice("v100", k2)
  
  u <- (1 - fx)*(1 - fy)*U[1,1] + fx*(1 - fy)*U[2,1] +
    (1 - fx)*fy      *U[1,2] + fx*fy       *U[2,2]
  v <- (1 - fx)*(1 - fy)*V[1,1] + fx*(1 - fy)*V[2,1] +
    (1 - fx)*fy      *V[1,2] + fx*fy       *V[2,2]
  c(u100 = u, v100 = v)
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
pts <- read.csv("points_resampled.csv") %>%
  mutate(
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"),
    location.lat       = as.numeric(location.lat),
    location.long       = as.numeric(location.long),
    year      = lubridate::year(timestamp)
  ) %>%  arrange(timestamp)

pts$u100_interp <- NA_real_
pts$v100_interp <- NA_real_

# ------------------------------------------------------------------------------
# Annotate points
# ------------------------------------------------------------------------------
pb <- progress::progress_bar$new(
  format = "  Annotating [:bar] :percent | :current/:total | eta :eta",
  total = nrow(pts), clear = FALSE, width = 60)

for (i in seq_len(nrow(pts))) {
  r    <- pts[i, ]
  meta <- get_era5_for_year(r$year)
  uv   <- bilinear_uv(r$location.long, r$location.lat, r$timestamp, meta)
  pts$u100_interp[i] <- uv["u100"]
  pts$v100_interp[i] <- uv["v100"]
  pb$tick()
}

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
pts<-pts%>%  arrange(unique_seg, timestamp)

write.csv(pts, outfile, row.names = FALSE)
cat("✔  Finished. Results saved to", outfile, "\n")


# ==============================================================================
# Part B — Segment-level wind summaries & point wind metrics
# ==============================================================================

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------
required_pkgs <- c("dplyr", "readr", "geosphere")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# ------------------------------------------------------------------------------
# Set parameters
# ------------------------------------------------------------------------------
infile      <- "points_resampled_annotated.csv"
outfile_seg <- "segments_wind_summary.csv"
outfile_pts <- "points_wind_metrics.csv"
D_min_m     <- 5000  # minimum displacement (m) to count as flight

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
pts <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    timestamp     = as.POSIXct(timestamp, tz = "UTC"),
    location.lat  = as.numeric(location.lat),
    location.long = as.numeric(location.long),
    u100_interp   = as.numeric(u100_interp),
    v100_interp   = as.numeric(v100_interp),
    CalHeading    = as.numeric(CalHeading)
  ) %>%  arrange(unique_seg, timestamp)

# ------------------------------------------------------------------------------
# Compute displacement & flight flag
# ------------------------------------------------------------------------------
pts <- pts %>%
  group_by(unique_seg) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  mutate(
    disp_m    = distHaversine(
      cbind(lag(location.long), lag(location.lat)),
      cbind(location.long,      location.lat)
    ),
    in_flight = !is.na(disp_m) & disp_m >= D_min_m
  ) %>%  ungroup()

# ------------------------------------------------------------------------------
# Compute wind metrics per fix
# ------------------------------------------------------------------------------
pts <- pts %>%
  mutate(
    wind_speed   = sqrt(u100_interp^2 + v100_interp^2),  
    mov_u        = sin(CalHeading * pi/180),             
    mov_v        = cos(CalHeading * pi/180),            
    wind_support = mov_u * u100_interp + mov_v * v100_interp,
    crosswind    = abs(mov_u * v100_interp - mov_v * u100_interp)
  )

# ------------------------------------------------------------------------------
# Summarise metrics per segment
# ------------------------------------------------------------------------------
segment_summary <- pts %>%
  filter(is.na(disp_m) | in_flight) %>%
  group_by(unique_seg) %>%
  summarise(
    mean_wind_speed   = mean(wind_speed,   na.rm = TRUE),
    mean_wind_support = mean(wind_support, na.rm = TRUE),
    mean_crosswind    = mean(crosswind,    na.rm = TRUE),
    .groups = "drop")

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
write_csv(segment_summary, outfile_seg)
write_csv(pts, outfile_pts)
message("✔ Segment summary written to ", outfile_seg)
message("✔ Point-level metrics written to ", outfile_pts)
