# ==============================================================================
# Title: Stacked Sankey (MGD → DTW → DIR) by season (SPRING | AUTUMN)
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
# Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi | am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
# 1) For each metric (MGD, DTW, DIR) and season (spring, autumn), find the
# winner route per segment within efficiency (GC, WO) and compass sets
# 2) Build two-column Sankey for each metric (SPRING | AUTUMN)
# 3) Stack three panels with headers and export final figure
# Inputs:
# - per_segment_compass_metrics_spring.csv
# - per_segment_compass_metrics_autumn.csv
# - per_segment_efficiency_metrics_spring.csv
# - per_segment_efficiency_metrics_autumn.csv
# Outputs:
# - sankey_plots/two_col_{med|dtw|dir}.png
# - sankey_plots/Fig7_sankey_stack.png
# Usage:
# Rscript R/11_sankey_stack.R
# Notes:
# - MGD = median geodesic distance; DTW = spatiotemporal; DIR = directional.

# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(networkD3)
library(htmlwidgets)
library(htmltools)
library(webshot)
library(magick)
# Ensure PhantomJS is available for webshot
if (!webshot::is_phantomjs_installed()) webshot::install_phantomjs()

# ------------------------------------------------------------------------------
# Route sets, labels, and style flags
# ------------------------------------------------------------------------------
route_types_compass    <- c("geographic","geomagnetic","magnetoclinic","sun","local_wind")
route_types_efficiency <- c("GC","wind_optimal")
route_label_map <- c(
  GC            = "GC",
  wind_optimal  = "WO",
  geographic    = "GL",
  geomagnetic   = "ML",
  magnetoclinic = "MC",
  sun           = "SC",
  local_wind    = "LW")

USE_SPATIOTEMPORAL_NO_HYPHEN <- TRUE     # set FALSE if journal prefers “spatio-temporal”
ST_SPATIOTEMPORAL <- if (USE_SPATIOTEMPORAL_NO_HYPHEN) "spatiotemporal" else "spatio-temporal"
EN_DASH <- " — "

# Methods order: spatial → spatiotemporal → directional
metrics_df <- tibble::tribble(
  ~metric, ~smaller_is_better, ~panel_title,                                  ~panel_letter,
  "med",   TRUE,                paste0("MGD", EN_DASH, "spatial similarity"),  "(a)",
  "dtw",   TRUE,                paste0("DTW", EN_DASH, ST_SPATIOTEMPORAL, " similarity"),     "(b)",
  "dir",   FALSE,               paste0("DIR", EN_DASH, "directional similarity"),         "(c)"
)
seasons <- c("spring","autumn")  # left = spring, right = autumn

# ------------------------------------------------------------------------------
# Output folder & global styling
# ------------------------------------------------------------------------------
out_dir <- "sankey_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

FONT_SIZE        <- 14   
NODE_WIDTH       <- 10
NODE_PADDING     <- 12
SVG_WIDTH        <- 700  
SVG_HEIGHT       <- 320
SNAPSHOT_ZOOM    <- 2    
SPACER_W         <- 80   
STRIP_HEIGHT     <- 40   
PANEL_TITLE_SIZE <- 26  
GLOBAL_HDR_H     <- 50   
GLOBAL_HDR_SIZE  <- 28   
GLOBAL_HDR_COL   <- "#333333"

# ------------------------------------------------------------------------------
# Helper: build ONE season’s Sankey → PNG
# ------------------------------------------------------------------------------
build_png <- function(season, metric, sib, outfile) {
  compass_segs <- read_csv(sprintf("per_segment_compass_metrics_%s.csv", season),show_col_types = FALSE)
  eff_segs     <- read_csv(sprintf("per_segment_efficiency_metrics_%s.csv", season),show_col_types = FALSE)
  all_eval     <- bind_rows(compass_segs, eff_segs)
  
# best efficiency for each segment
best_eff <- all_eval %>%filter(route %in% route_types_efficiency) %>%
    group_by(unique_seg) %>%
    { if (sib) slice_min(., .data[[metric]], with_ties = FALSE)
      else     slice_max(., .data[[metric]], with_ties = FALSE) } %>%
    ungroup() %>% transmute(unique_seg, best_eff = route)
  
# best compass for each segment
best_comp <- all_eval %>%filter(route %in% route_types_compass) %>%
    group_by(unique_seg) %>%
    { if (sib) slice_min(., .data[[metric]], with_ties = FALSE)
      else     slice_max(., .data[[metric]], with_ties = FALSE) } %>%
    ungroup() %>%transmute(unique_seg, best_comp = route)
  
sankey_df  <- inner_join(best_eff, best_comp, by = "unique_seg")
total_segs <- dplyr::n_distinct(sankey_df$unique_seg)
  
# links & nodes
links <- sankey_df %>%count(best_eff, best_comp, name = "value") %>%
    rename(source = best_eff, target = best_comp)
  
# Node order: efficiency first, then compass
node_order <- c("GC","wind_optimal","geographic","geomagnetic","magnetoclinic","local_wind","sun")
nodes <- data.frame(name = node_order, stringsAsFactors = FALSE)
  
links <- links %>%
    mutate(
      source_id = match(source, nodes$name) - 1,
      target_id = match(target, nodes$name) - 1)
  
# Node labels with percentages
eff_lab  <- sankey_df %>% count(best_eff)  %>% rename(name = best_eff,  count = n)
comp_lab <- sankey_df %>% count(best_comp) %>% rename(name = best_comp, count = n)
node_labels <- bind_rows(eff_lab, comp_lab) %>%
    mutate(
      pct     = 100 * count / total_segs,
      display = dplyr::recode(name, !!!route_label_map),
      label   = paste0(display, " (", sprintf("%.1f%%", pct), ")")
    ) %>% select(name, label)
  nodes <- dplyr::left_join(nodes, node_labels, by = "name")
  
# palette
palette <- c(
    GC            = "#2C7BB6",
    wind_optimal  = "#FEB24C",
    geographic    = "#009E73",
    geomagnetic   = "#56B4E9",
    magnetoclinic = "#F0E442",
    local_wind    = "#E69F00",
    sun           = "#CC79A7")

links$linkGroup <- links$source
nodes$nodeGroup <- nodes$name
colScale <- paste0(
    "d3.scaleOrdinal()",
    ".domain([\"", paste(names(palette), collapse = "\",\""), "\"])",
    ".range([\"",  paste(palette,       collapse = "\",\""), "\"])")
  
# Sankey widget
sankey <- sankeyNetwork(
    Links       = links, Nodes       = nodes,
    Source      = "source_id", Target = "target_id",
    Value       = "value",    NodeID  = "label",
    NodeGroup   = "nodeGroup", LinkGroup = "linkGroup",
    fontFamily  = "Arial",
    fontSize    = FONT_SIZE,
    nodeWidth   = NODE_WIDTH,
    nodePadding = NODE_PADDING,
    sinksRight  = FALSE,
    width       = SVG_WIDTH,
    height      = SVG_HEIGHT,
    colourScale = JS(colScale),
    units       = "segments"
  ) %>% onRender("
      function(el,x){
        var gap=8, svg=d3.select(el).select('svg');
        function fix(){
          svg.selectAll('.node').each(function(d){
            var t=d3.select(this).select('text');
            if(d.targetLinks.length===0){
              t.attr('x',-gap).attr('text-anchor','end');
            } else {
              t.attr('x',d.dx+gap).attr('text-anchor','start');
            }
          });
        }
        fix(); setTimeout(fix,0);
      }")
  
# snapshot to PNG
tmp_html <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(sankey, tmp_html, selfcontained = TRUE)
webshot(tmp_html, outfile,
          selector = "svg",
          vwidth   = SVG_WIDTH,
          vheight  = SVG_HEIGHT,
          zoom     = SNAPSHOT_ZOOM)
}

# ------------------------------------------------------------------------------
# Build metric-specific two-column images (SPRING | AUTUMN)
# ------------------------------------------------------------------------------
for (i in seq_len(nrow(metrics_df))) {
  metric <- metrics_df$metric[i]
  sib    <- metrics_df$smaller_is_better[i]
  message("⟹ Building panels for metric: ", metric)
  
  tmpS <- tempfile(fileext = ".png")
  tmpA <- tempfile(fileext = ".png")
  
  build_png("spring", metric, sib,  tmpS)
  build_png("autumn", metric, sib,  tmpA)
  
  imgS <- image_trim(image_read(tmpS))
  imgA <- image_trim(image_read(tmpA))
  
  spacer <- image_blank(
    width  = SPACER_W,
    height = max(image_info(imgS)$height, image_info(imgA)$height),
    color  = "white"
  )
  
  two_col <- image_append(c(imgS, spacer, imgA), stack = FALSE)
  two_col <- image_trim(two_col)
  
  out_png <- file.path(out_dir, paste0("two_col_", metric, ".png"))
  image_write(two_col, out_png)
  
  unlink(c(tmpA, tmpS))
  message("   ✔ saved: ", out_png)
}

# ------------------------------------------------------------------------------
# Stack the three panels with a single global (SPRING | AUTUMN) header
# ------------------------------------------------------------------------------
annotate_panel <- function(img, label, title) {
img <- image_annotate(img, label, gravity = "northwest",
                        size = 36, weight = 700, location = "+12+8")
strip <- image_blank(width = image_info(img)$width, height = STRIP_HEIGHT, color = "white")
strip <- image_annotate(strip, title, gravity = "center",size = PANEL_TITLE_SIZE, weight = 600)
  image_append(c(strip, img), stack = TRUE)
}

# Read the three per-metric images (Methods order)
p_spatial <- image_read(file.path(out_dir, "two_col_med.png"))
p_time    <- image_read(file.path(out_dir, "two_col_dtw.png"))
p_dir     <- image_read(file.path(out_dir, "two_col_dir.png"))

target_w <- min(image_info(p_spatial)$width, image_info(p_time)$width, image_info(p_dir)$width)
p_spatial <- image_resize(p_spatial, paste0(target_w, "x"))
p_time    <- image_resize(p_time,    paste0(target_w, "x"))
p_dir     <- image_resize(p_dir,     paste0(target_w, "x"))

p_spatial <- annotate_panel(p_spatial, metrics_df$panel_letter[metrics_df$metric=="med"],
                            metrics_df$panel_title[metrics_df$metric=="med"])
p_time    <- annotate_panel(p_time,    metrics_df$panel_letter[metrics_df$metric=="dtw"],
                            metrics_df$panel_title[metrics_df$metric=="dtw"])
p_dir     <- annotate_panel(p_dir,     metrics_df$panel_letter[metrics_df$metric=="dir"],
                            metrics_df$panel_title[metrics_df$metric=="dir"])
panel_gap <- 40  # adjust this number to control spacing (e.g., 40–80 px)
spacer_v  <- image_blank(width = target_w, height = panel_gap, color = "white")
fig_body <- image_append(c(p_spatial, spacer_v, p_time, spacer_v, p_dir), stack = TRUE)
fig_body <- image_trim(fig_body)
W <- image_info(fig_body)$width
header <- image_blank(width = W, height = GLOBAL_HDR_H, color = "white")
left_w  <- (W - SPACER_W) / 2
x_left  <- round(left_w/2)
x_right <- round(left_w + SPACER_W + left_w/2)

header <- image_annotate(header, "Spring",
                         gravity  = "northwest",
                         location = paste0("+", x_left  - 40, "+", 12),
                         size = GLOBAL_HDR_SIZE, weight = 500, color = GLOBAL_HDR_COL)

header <- image_annotate(header, "Autumn",
                         gravity  = "northwest",
                         location = paste0("+", x_right - 48, "+", 12),
                         size = GLOBAL_HDR_SIZE, weight = 500, color = GLOBAL_HDR_COL)
fig7 <- image_append(c(header, fig_body), stack = TRUE)
fig7 <- image_border(fig7, color = "white", geometry = "0x8")
fig7 <- image_resize(fig7, "1800x")
final_path <- file.path(out_dir, "sankey_stack.png")
image_write(fig7, final_path)
message("✔ Final stacked figure written: ", final_path)
