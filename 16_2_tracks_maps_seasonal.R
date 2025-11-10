# ==============================================================================
# Title: Seasonal track maps (Spring vs Autumn) 
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read hourly-resampled, filtered points and build LINESTRING per segment
#   2) Extract start/end fixes per segment (for endpoint markers)
#   3) Fetch Esri World Imagery tiles per season-specific extent
#   4) Plot tracks for SPRING (left) and AUTUMN (right) with clean axes
#   5) Save a two-panel PNG figure for journal-ready background map
# Inputs:
#   - points_resampled.csv  (must contain: unique_seg, season, lon, lat, timestamp)
# Outputs:
#   - Tracks_Spring_vs_Autumn.png
# Usage:
#   Rscript R/16_2_tracks_maps_seasonal.R

# ==============================================================================

rm(list=ls())

# ------------------------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------------------------
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(maptiles)
library(ggspatial)
library(patchwork)
library(units)

# ------------------------------------------------------------------------------
# 1) Read & prepare data
# ------------------------------------------------------------------------------
tracks <- read.csv("points_resampled.csv", stringsAsFactors = FALSE) |>
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"))

tracks_sf <- st_as_sf(
  tracks,
  coords = c("lon", "lat"),
  crs    = 4326,
  remove = FALSE
)

lines_sf <- tracks_sf |>
  arrange(unique_seg, timestamp) |>
  group_by(unique_seg, season) |>
  summarise(
    geometry = st_cast(st_combine(geometry), "LINESTRING"),
    .groups  = "drop"
  )

endpoints_sf <- tracks_sf |>
  group_by(unique_seg, season) |>
  arrange(timestamp) |>
  slice(c(1, n())) |>
  ungroup()

# ------------------------------------------------------------------------------
# 2) Mapping function (one season)
# ------------------------------------------------------------------------------
make_season_map <- function(season_name,
                            buffer_deg   = 2,
                            zoom         = 5,
                            asp_ratio    = 2/3) {
  
  seg_lines <- filter(lines_sf,     season == season_name)
  ep_pts    <- filter(endpoints_sf, season == season_name)
  
  bb_exp <- st_bbox(seg_lines) + c(-buffer_deg, -1.0, buffer_deg, 1.0)
  bb_sfc <- st_as_sfc(bb_exp)
  
  tiles_3857 <- get_tiles(
    bb_sfc,
    provider = "Esri.WorldImagery",
    crop     = TRUE,
    zoom     = zoom
  )
  
  tiles_ll <- terra::project(tiles_3857, "EPSG:4326")
  col_line <- if (season_name == "spring") "#440154FF" else "#ff2b6d"
  ggplot() +
    annotation_spatial(tiles_ll, alpha = 0.6) +
    geom_sf(data = seg_lines, colour = col_line, linewidth = 0.3) +
    geom_sf(data = ep_pts, shape = 21, fill = col_line,
            colour = "white", size = 2.4, stroke = 0.3) +
    coord_sf(
      crs = 4326,
      xlim = bb_exp[c("xmin", "xmax")],
      ylim = bb_exp[c("ymin", "ymax")],
      expand = FALSE
    ) +
    scale_x_continuous(expand = c(0,0), breaks = NULL, labels = NULL, name = NULL) +
    scale_y_continuous(expand = c(0,0), breaks = NULL, labels = NULL, name = NULL) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid        = element_blank(),
      legend.position   = "none",
      aspect.ratio      = asp_ratio,
      plot.margin       = margin(0, 0, 0, 0),
      axis.text         = element_blank(),
      axis.title        = element_blank(),
      axis.ticks        = element_blank(),
      axis.ticks.length = unit(0, "pt")
    )
}

# ------------------------------------------------------------------------------
# 3) Build panels & save figure
# ------------------------------------------------------------------------------
map_autumn <- make_season_map("autumn")
map_spring <- make_season_map("spring")

final_plot <- map_spring + map_autumn + plot_layout(ncol = 2, widths = c(0.5, 0.5))
print(final_plot)

ggsave(
  filename = "Tracks_Spring_vs_Autumn.png",
  plot     = final_plot,
  width    = 15, height = 5,
  dpi      = 600, units = "in"
)
