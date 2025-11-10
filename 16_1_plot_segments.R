# ==============================================================================
# Title: Segment maps with routes, DEM, day/night & storms 
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read segments, hourly points, per-route CSVs (GC/WO/GL/ML/MC/SC/LW)
#   2) Join per-point Kp, classify day/night (–6° alt) and storm hours (Kp>5)
#   3) Draw DEM-shaded basemap + observed track + seven candidate routes
#   4) Highlight the per-segment best_route and annotate meta info
#   5) Export one PNG per segment to significant_patterns/all_segments/
# Inputs:
#   - segments_final.csv
#   - points_resampled.csv
#   - points_kpinfo.csv
#   - Geomagnetic_routes.csv, Geographic_routes.csv, Magnetoclinic_routes.csv,
#     SunCompass_classic_routes.csv, GreatCircle_routes.csv,
#     wind_optimal_routes.csv, local_wind_routes.csv
# Outputs:
#   - significant_patterns/all_segments/<unique_seg>.png
# Usage:
#   Rscript R/16_1_plot_segments.R

# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(geosphere)
  library(maps);      library(mapdata);   library(suncalc)
  library(raster);    library(lubridate); library(elevatr)
  library(sf);        library(readr);     library(dplyr)
})

# ------------------------------------------------------------------------------
# User files (paths)
# ------------------------------------------------------------------------------
# (Removed reliance on filtered_segments_*.csv files)
seg_file <- "segments_final.csv"
pts_file <- "points_resampled.csv"

files_routes <- list(
  gm = "Geomagnetic_routes.csv",
  gg = "Geographic_routes.csv",
  mc = "Magnetoclinic_routes.csv",
  sc = "SunCompass_classic_routes.csv",
  gc = "GreatCircle_routes.csv",
  wo = "wind_optimal_routes.csv",
  tw = "local_wind_routes.csv")

# ------------------------------------------------------------------------------
# Load static data
# ------------------------------------------------------------------------------
segments <- read_csv(seg_file, show_col_types = FALSE) %>%
  mutate(unique_seg      = as.character(unique_seg),
         first_timestamp = as.POSIXct(first_timestamp, "%Y-%m-%d %H:%M:%S", tz = "UTC"))

points_df <- read_csv(pts_file, show_col_types = FALSE) %>%
  mutate(unique_seg = as.character(unique_seg),
         time       = as.POSIXct(timestamp, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
         lon        = location.long,
         lat        = location.lat)

points_df_kp <- read_csv("points_kpinfo.csv", show_col_types = FALSE) %>%
  mutate(unique_seg = as.character(unique_seg),
         time       = as.POSIXct(timestamp, "%Y-%m-%d %H:%M:%S", tz = "UTC"))

points_df <- points_df %>%
  left_join(
    dplyr::select(points_df_kp, unique_seg, timestamp, Kp),
    by = c("unique_seg", "timestamp")
  )

routes <- lapply(files_routes, read_csv, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# Sun / daylight helper
# ------------------------------------------------------------------------------
SUN_THRESH <- -6
is_daylight <- function(lat, lon, time) {
  df  <- data.frame(date = time, lat = lat, lon = lon)
  alt <- suncalc::getSunlightPosition(data = df, keep = "altitude")$altitude * 180/pi
  alt >= SUN_THRESH}

# ------------------------------------------------------------------------------
# Output folder (ALL segments)
# ------------------------------------------------------------------------------
category_tag <- "all_segments"
main_out <- file.path("significant_patterns", category_tag)
dir.create(main_out, showWarnings = FALSE, recursive = TRUE)
sel_segs <- segments %>%
  dplyr::select(
    unique_seg, first_timestamp,
    first_long, first_lat, last_long, last_lat,
    bird_id, season, segment_type,
    cumulative_distance, great_circle_distance, duration, Linearity,
    best_route
  )

cat("• Total segments to plot:", nrow(sel_segs), "\n")

# ------------------------------------------------------------------------------
# Loop over segments → draw map, routes, track, and export PNG
# ------------------------------------------------------------------------------
grand_total <- 0
seg_count <- 0

for (seg_id in unique(sel_segs$unique_seg)) {
  seg_count <- seg_count + 1
  cat(sprintf("   [%3d/%3d]  %-25s\r",
              seg_count, nrow(sel_segs), seg_id))
  flush.console()
  
  seg_rows <- dplyr::filter(sel_segs, unique_seg == seg_id)
  seg_meta <- seg_rows[1, ]
  
  trk <- dplyr::filter(points_df, unique_seg == seg_id)
  if (nrow(trk) == 0) next
  
  gm <- dplyr::filter(routes$gm, unique_seg == seg_id)
  gg <- dplyr::filter(routes$gg, unique_seg == seg_id)
  mc <- dplyr::filter(routes$mc, unique_seg == seg_id)
  sc <- dplyr::filter(routes$sc, unique_seg == seg_id)
  gc <- dplyr::filter(routes$gc, unique_seg == seg_id)
  wo <- dplyr::filter(routes$wo, unique_seg == seg_id)
  tw <- dplyr::filter(routes$tw, unique_seg == seg_id)
  
  trk$day   <- is_daylight(trk$lat, trk$lon, trk$time)
  trk$storm <- !is.na(trk$Kp) & trk$Kp > 5
  day_no_storm   <-  trk$day & !trk$storm
  night_no_storm <- !trk$day & !trk$storm
  storm_idx      <-  trk$storm
  
  # info lines (stats only)
  meta_lines <- c(
    sprintf("Cum. Dist. : %.0f km", seg_meta$cumulative_distance),
    sprintf("Geo. Dist.  : %.0f km", seg_meta$great_circle_distance),
    sprintf("Duration    : % .1f h",  seg_meta$duration),
    sprintf("Linearity    :  %.2f",    seg_meta$Linearity)
  )
  
  f_png <- file.path(main_out, paste0(seg_id, ".png"))
  png(f_png, width = 4200, height = 3100, res = 500, type = "cairo-png")
  
  xr   <- range(c(trk$lon, gm$lon, gg$lon, mc$lon, sc$lon, gc$lon, wo$lon))
  yr   <- range(c(trk$lat, gm$lat, gg$lat, mc$lat, sc$lat, gc$lat, wo$lat))
  xpad <- diff(xr) * 0.3
  ypad <- diff(yr) * 0.65
  
  par(mar = c(2.5, 2.5, 1, 6),
      mgp = c(1.4, 0.35, 0),
      cex.axis = 0.55,
      cex.lab = 0.6,
      cex.main = .70)
  
  plot(NA, NA, asp = 1, type = "n",
       xlim = xr + c(-xpad, xpad),
       ylim = yr + c(-ypad, ypad),
       xaxs = "i", yaxs = "i",
       main = sprintf("%s. %d, %s, %s",
                      seg_meta$bird_id,
                      year(trk$timestamp[1]),
                      seg_meta$season,
                      seg_meta$segment_type),
       xlab = "Longitude (°E)", ylab = "Latitude (°N)")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "#dfeffc", border = NA)
  maps::map("worldHires", add = TRUE, fill = TRUE, col = "#FFFFFF", border = NA)
  
  # DEM shading
  bb <- rbind(
    c(xr[1]-xpad, yr[1]-ypad),
    c(xr[2]+xpad, yr[1]-ypad),
    c(xr[2]+xpad, yr[2]+ypad),
    c(xr[1]-xpad, yr[2]+ypad),
    c(xr[1]-xpad, yr[1]-ypad)
  )
  bb_poly <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(sf::st_polygon(list(bb)), crs = 4326)
  )
  
  try({
    dem <- elevatr::get_elev_raster(bb_poly, z = 5, clip = "locations")
    dem[dem <= 0] <- NA
    
    elev_breaks <- c(0,150,300,450,600,750,900)
    elev_cols   <- c("#FFFFFF",
                     colorRampPalette(c("#f0e6d2","#A0522D"))(length(elev_breaks)-3),
                     "#5C3317")
    
    dem_cl <- dem
    dem_cl[dem_cl < elev_breaks[1]] <- elev_breaks[1]
    dem_cl[dem_cl > elev_breaks[length(elev_breaks)]] <- elev_breaks[length(elev_breaks)]
    
    plot(dem_cl, breaks = elev_breaks,
         col    = c(elev_cols[1], adjustcolor(elev_cols[-1], alpha.f = .3)),
         legend = FALSE, add = TRUE)
    
    plot(dem_cl, legend.only = TRUE, breaks = elev_breaks, col = elev_cols,
         smallplot = c(0.9,0.92,0.55,0.85), legend.width = .3,
         axis.args = list(at = elev_breaks, labels = elev_breaks, cex.axis = .5),
         legend.args = list(text = "Elevation (m)", side = 3, line = .5, cex = .6))
  }, silent = TRUE)
  
  maps::map("worldHires", add = TRUE, fill = FALSE,
            col = adjustcolor("grey80", .6), lwd = .5)
  
  # Routes
  lines(gg$lon, gg$lat, col = "darkblue",    lwd = 1)
  lines(gm$lon, gm$lat, col = "firebrick",   lwd = 1)
  lines(mc$lon, mc$lat, col = "magenta",     lwd = 1)
  lines(sc$lon, sc$lat, col = "#FFD700",     lwd = 1)
  lines(tw$lon, tw$lat, col = "#00FF00",     lwd = 1)
  lines(gc$lon, gc$lat, col = "tomato",      lty = 2, lwd = 1)
  lines(wo$lon, wo$lat, col = "chartreuse4", lty = 2, lwd = 1)
  
  route_dfs  <- list(geographic=gg, geomagnetic=gm, magnetoclinic=mc, sun=sc,
                     local_wind=tw, GC=gc, wind_optimal=wo)
  route_cols <- c(geographic="darkblue", geomagnetic="firebrick", magnetoclinic="magenta",
                  sun="#FFD700", local_wind="springgreen", GC="tomato", wind_optimal="chartreuse4")
  route_lty  <- c(geographic=1, geomagnetic=1, magnetoclinic=1, sun=1, local_wind=1, GC=2, wind_optimal=2)
  route_lwd  <- c(geographic=1.1, geomagnetic=1.1, magnetoclinic=1, sun=1, local_wind=1, GC=1, wind_optimal=1)
  
  bn <- seg_rows$best_route[1]
  lines((df <- route_dfs[[bn]])$lon, df$lat,
        col = route_cols[bn], lty = route_lty[bn], lwd = route_lwd[bn] + 1)
  
  lines(trk$lon, trk$lat, col = "black", lwd = 1)
  points(seg_meta$first_long, seg_meta$first_lat, pch = 21, bg = "orange2", col = "orange2", cex = 1.7)
  points(seg_meta$last_long,  seg_meta$last_lat,  pch = 24, bg = "orange2", col = "orange2", cex = 1.7)
  
  points(trk$lon[day_no_storm],   trk$lat[day_no_storm],   pch = 1,  col = "black", cex = 0.65, lwd = 0.5)
  points(trk$lon[night_no_storm], trk$lat[night_no_storm], pch = 21, bg  = adjustcolor("black", 0.3),
         col = "black", cex = 0.65, lwd = 0.5)
  points(trk$lon[storm_idx],      trk$lat[storm_idx],      pch = 8,  col = "red",   cex = 0.6,  lwd = 0.4)
  
  usr <- par("usr")
  w   <- diff(usr[1:2]) * 0.16
  h   <- diff(usr[3:4]) * 0.13
  x0  <- usr[2] + diff(usr[1:2]) * 0.01
  x1  <- x0 + w
  y0  <- usr[3] + diff(usr[3:4]) * 0.32
  y1  <- y0 + h
  
  rect(x0, y0, x1, y1, col = adjustcolor("white", 0.85), border = "gray90", xpd = TRUE)
  text(x = x0 + 0.01 * diff(usr[1:2]),
       y = y1 - (0:3) * (h / 4) - 0.21,
       labels = meta_lines, adj = c(0, 1), cex = 0.6, xpd = TRUE)
  
  dev.off()
  grand_total <- grand_total + 1
}

cat("\n=== ALL DONE —", grand_total, "PNG files written ===\n")
