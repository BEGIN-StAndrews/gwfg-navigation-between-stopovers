# ==============================================================================
# Title: Triple correlation heatmaps (MGD, DTW, DIR) with panel titles & centered legend
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
# Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi | am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
# 1) Read per-segment evaluation metrics (med, dtw, dir)
# 2) Build three correlation heatmaps across compass routes for each metric
# 3) Stack panels (a) MGD, (b) DTW, (c) DIR with a single centered legend
# 4) Save PNG and PDF outputs sized for journal figures
# Inputs:
# - evaluation_metrics.csv (columns: unique_seg, route, med, dtw, dir)
# Outputs:
# - corr_outputs/corr_triple.png
# - corr_outputs/corr_triple.pdf
# Usage:
# Rscript R/13_correlation_heatmaps.R

# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(Hmisc) 
  library(ggplot2)
  library(viridisLite)
  library(grid) 
  library(cowplot) 
  library(scales) 
})

# ------------------------------------------------------------------------------
# I/O paths
# ------------------------------------------------------------------------------
in_file <- "evaluation_metrics.csv"
out_dir <- "corr_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# Load & filter data (compass routes only)
# ------------------------------------------------------------------------------
all_eval <- read_csv(in_file, show_col_types = FALSE)
route_types_compass <- c("geographic","geomagnetic","magnetoclinic","sun","local_wind")
label_map  <- c(geographic="GL", geomagnetic="ML", magnetoclinic="MC", sun="SC", local_wind="LW")
route_levels <- names(label_map)

compass_eval <- all_eval %>%filter(route %in% route_types_compass) %>%
  mutate(route = factor(route, levels = route_levels))

wide_eval <- compass_eval %>%
  pivot_longer(cols = c(med, dtw, dir), names_to = "metric", values_to = "value") %>%
  unite("metric_route", metric, route) %>%
  pivot_wider(names_from = metric_route, values_from = value)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
corr_mats <- function(prefix){
  cols <- grep(paste0("^", prefix, "_"), names(wide_eval), value = TRUE)
  cols <- cols[order(match(sub(paste0("^", prefix, "_"), "", cols), route_levels))]
  mat  <- as.matrix(wide_eval[cols])
  rc   <- rcorr(mat, type = "pearson")
  labs <- label_map[sub(paste0("^", prefix, "_"), "", colnames(mat))]
  rownames(rc$r) <- rownames(rc$P) <- colnames(rc$r) <- colnames(rc$P) <- labs
  list(R = rc$r, P = rc$P)
}

label_color_for_r <- function(r_vals, palette_fn = viridis){
  cols <- palette_fn(256)[pmax(1, pmin(256, round(r_vals*255)+1))]
  rgb <- t(col2rgb(cols))/255
  lum <- 0.2126*rgb[,1] + 0.7152*rgb[,2] + 0.0722*rgb[,3]
  ifelse(lum < 0.5, "white", "black")
}

corr_plot <- function(prefix, limits=c(0,1), use_sqrt=TRUE, palette_fn=viridis,
                      show_legend = TRUE){
  cp <- corr_mats(prefix); R <- cp$R; P <- cp$P
  df <- as.data.frame(as.table(R)); names(df) <- c("Var1","Var2","r")
  pdf <- as.data.frame(as.table(P)); names(pdf) <- c("Var1","Var2","p")
  df <- dplyr::left_join(df, pdf, by=c("Var1","Var2")) %>%
    dplyr::mutate(i = as.integer(factor(Var1, levels = colnames(R))),
                  j = as.integer(factor(Var2, levels = colnames(R)))) %>%
    dplyr::filter(j >= i)
  
  df$lab_col <- label_color_for_r(scales::rescale(df$r, to = c(0,1), from = limits),
                                  palette_fn = palette_fn)
  
  pal <- palette_fn(256)
  
  gg <- ggplot(df, aes(Var2, Var1, fill = r)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", r), color = I(lab_col)), size = 3.2) +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_family = "Arial") +
    theme(
      plot.title  = element_blank(),             # titles added later
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      legend.margin      = margin(0,0,0,0),
      legend.box.margin  = margin(0,0,0,0),
      legend.spacing     = unit(0, "pt"),
      legend.key.height  = unit(1.2, "cm"),
      legend.key.width   = unit(0.5, "cm"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 2, 4, 4)          # narrow right margin -> legend close
    )
  
  gg <- gg + scale_fill_gradientn(colors = pal, limits = limits,
                                  oob = scales::squish, name = "r",
                                  trans = if (use_sqrt) "sqrt" else "identity",
                                  breaks = seq(0,1,0.25))
  
  if (!show_legend) gg <- gg + theme(legend.position = "none")
  gg
}

p_med <- corr_plot("med", use_sqrt=TRUE,  show_legend = FALSE) +
  labs(title = "(a) MGD — spatial similarity") +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))

p_dtw <- corr_plot("dtw", use_sqrt=TRUE,  show_legend = TRUE) +
  labs(title = "(b) DTW — spatiotemporal") +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))

p_dir <- corr_plot("dir", use_sqrt=TRUE,  show_legend = FALSE) +
  labs(title = "(c) DIR — directional") +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))

legend_g <- cowplot::get_legend(
  p_dtw + theme(
    legend.position    = "right",
    legend.margin      = margin(0,0,0,0),
    legend.box.margin  = margin(0,0,0,0),
    legend.spacing     = unit(0, "pt")
  ))

p_dtw_nl <- p_dtw + theme(legend.position = "none")

gap <- 0.08
left_stack <- cowplot::plot_grid(
  p_med,
  NULL,
  p_dtw_nl,
  NULL,
  p_dir,
  ncol = 1, rel_heights = c(1, gap, 1, gap, 1), align = "v")

right_col <- cowplot::plot_grid(
  NULL,
  legend_g,
  NULL,
  ncol = 1,
  rel_heights = c(1 + gap/2, 1, 1 + gap/2))

final_plot <- cowplot::plot_grid(
  left_stack, right_col,
  ncol = 2, rel_widths = c(1, 0.10), align = "h")

ggsave(file.path(out_dir, "corr_triple.png"), final_plot, width = 7.2, height = 9.0, dpi = 300)
ggsave(file.path(out_dir, "corr_triple.pdf"), final_plot, width = 7.2, height = 9.0, device = cairo_pdf)

message("✔ Saved: ", normalizePath(out_dir))
