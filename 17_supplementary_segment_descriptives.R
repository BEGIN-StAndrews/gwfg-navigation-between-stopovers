# ==============================================================================
# Title: Supplementary segment summaries and figures (Table + 2 Figures)
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Script No.: 17
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
#
# Summary:
#   This script reproduces three supplementary outputs in one place:
#   1) Build filtered segment dataset (common segments across route files), classify
#      departure light, and generate seasonal summary table (Table-style CSV)
#   2) Plot distribution of number of flight segments per seasonal migration (PNG/PDF)
#   3) Plot individual variation (per-individual medians) by season for duration,
#      cumulative distance, and straightness (PNG/PDF)
#
# Inputs (working directory):
#   - segments_resample_filtered.csv
#   - points_resampled_filtered.csv
#   - Geomagnetic_routes.csv
#   - Geographic_routes.csv
#   - Magnetoclinic_routes.csv
#   - SunCompass_classic_routes.csv
#   - local_wind_routes.csv
#   - wind_optimal_routes_raw.csv
#
# Outputs (single folder):
#   - 17_supplementary_segment_descriptives/segments_final.csv
#   - 17_supplementary_segment_descriptives/seasonal_summary_table.csv
#   - 17_supplementary_segment_descriptives/segments_per_migration_discrete_by_season.png
#   - 17_supplementary_segment_descriptives/segments_per_migration_discrete_by_season.pdf
#   - 17_supplementary_segment_descriptives/FigS_individual_variation_duration_distance_straightness.png
#   - 17_supplementary_segment_descriptives/FigS_individual_variation_duration_distance_straightness.pdf
#
# Usage:
#   Rscript R/17_supplementary_segment_descriptives.R
# ==============================================================================

options(width = 120)

# ------------------------------------------------------------------------------
# Load libraries (auto-install if missing)
# ------------------------------------------------------------------------------
req <- c(
  "dplyr", "readr", "geosphere", "RANN", "data.table", "lubridate", "suncalc",
  "purrr", "tidyverse", "scales", "ggplot2", "tidyr", "patchwork"
)
inst <- setdiff(req, rownames(installed.packages()))
if (length(inst)) install.packages(inst, dependencies = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# Shared output folder (all outputs from all three blocks)
# ------------------------------------------------------------------------------
outdir <- file.path(getwd(), "17_supplementary_segment_descriptives")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# BLOCK 1 — Build filtered segment dataset + classify departure light + seasonal
#           summary table (table-style output)
# ==============================================================================

# -----------------------------
# STEP 1.1: Compute metrics from segments_resample_filtered.csv
# -----------------------------
segments <- read.csv("segments_resample_filtered.csv", stringsAsFactors = FALSE) %>%
  mutate(
    first_timestamp       = as.POSIXct(first_timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    first_lat             = as.numeric(first_lat),
    first_long            = as.numeric(first_long),
    last_lat              = as.numeric(last_lat),
    last_long             = as.numeric(last_long),
    duration              = as.numeric(duration),
    cumulative_distance   = as.numeric(cumulative_distance),
    great_circle_distance = as.numeric(great_circle_distance),
    MaxKp                 = as.numeric(MaxKp),
    unique_seg            = as.character(unique_seg),
    straightness          = great_circle_distance / cumulative_distance,
    geomag_storm          = ifelse(MaxKp > 5, 1, 0)
  )

# -----------------------------
# STEP 1.2: Get intersection of unique_seg across all required files
# -----------------------------
common_segs <- segments$unique_seg

for (f in c(
  "points_resampled_filtered.csv",
  "Geomagnetic_routes.csv",
  "Geographic_routes.csv",
  "Magnetoclinic_routes.csv",
  "SunCompass_classic_routes.csv",
  "local_wind_routes.csv",
  "wind_optimal_routes_raw.csv"
)) {
  segs <- read_csv(f, show_col_types = FALSE)$unique_seg
  common_segs <- intersect(common_segs, segs)
}

# -----------------------------
# STEP 1.3: Filter and save final segment dataset
# -----------------------------
segments_filtered <- filter(segments, unique_seg %in% common_segs)  # expected n = 1524
segments_final_file <- file.path(outdir, "segments_final.csv")
write.csv(segments_filtered, segments_final_file, row.names = FALSE)

# -----------------------------
# STEP 1.4: Classify light at departure
# -----------------------------
classify_start_light <- function(lat, lon, t0) {
  st <- getSunlightTimes(
    date = as_date(t0), lat = lat, lon = lon,
    keep = c("sunrise", "sunset", "solarNoon", "nadir")
  )
  alt_now <- getSunlightPosition(t0, lat, lon)$altitude
  
  if (!is.na(st$sunrise) && !is.na(st$sunset)) {
    if (alt_now < -6 * pi / 180) return("night") else return("day")
  } else {
    alt_noon <- getSunlightPosition(st$solarNoon, lat, lon)$altitude
    if (alt_noon > 0) return("polar_day") else return("polar_night")
  }
}

segments_filtered <- segments_filtered %>%
  rowwise() %>%
  mutate(light_start = classify_start_light(first_lat, first_long, first_timestamp)) %>%
  ungroup() %>%
  mutate(night_departure = ifelse(light_start %in% c("night", "polar_night"), 1, 0))

# -----------------------------
# STEP 1.5: Summary stats function
# -----------------------------
summary_stats <- function(df) {
  data.frame(
    n_individuals = n_distinct(df$individual.local.identifier),
    n_segments = nrow(df),
    
    duration = sprintf(
      "%.1f ± %.1f (%.1f–%.1f)",
      mean(df$duration), sd(df$duration), min(df$duration), max(df$duration)
    ),
    
    cumulative_distance = sprintf(
      "%.1f ± %.1f (%.1f–%.1f)",
      mean(df$cumulative_distance), sd(df$cumulative_distance),
      min(df$cumulative_distance), max(df$cumulative_distance)
    ),
    
    great_circle_distance = sprintf(
      "%.1f ± %.1f (%.1f–%.1f)",
      mean(df$great_circle_distance), sd(df$great_circle_distance),
      min(df$great_circle_distance), max(df$great_circle_distance)
    ),
    
    straightness = sprintf(
      "%.2f ± %.2f (%.2f–%.2f)",
      mean(df$straightness), sd(df$straightness),
      min(df$straightness), max(df$straightness)
    ),
    
    ini_heading = sprintf(
      "%.2f ± %.2f (%.2f–%.2f)",
      mean(df$init_mean_heading), sd(df$init_mean_heading),
      min(df$init_mean_heading), max(df$init_mean_heading)
    ),
    
    pct_geomag_storm = percent(mean(df$geomag_storm, na.rm = TRUE), accuracy = 0.1),
    pct_nocturnal_departure = percent(mean(df$night_departure, na.rm = TRUE), accuracy = 0.1)
  )
}

# -----------------------------
# STEP 1.6: Summarise per season and save table
# -----------------------------
spring_summary <- summary_stats(filter(segments_filtered, season == "spring"))
autumn_summary <- summary_stats(filter(segments_filtered, season == "autumn"))

seasonal_summary <- bind_rows(
  spring = spring_summary,
  autumn = autumn_summary,
  .id = "season"
)

print(seasonal_summary)

seasonal_summary_file <- file.path(outdir, "seasonal_summary_table.csv")
write.csv(seasonal_summary, seasonal_summary_file, row.names = FALSE)

# ==============================================================================
# BLOCK 2 — Figure: number of flight segments per seasonal migration
# ==============================================================================

# -----------------------------
# STEP 2.1: Inputs / outputs
# -----------------------------
infile_segments_final <- segments_final_file

# -----------------------------
# STEP 2.2: Read + checks
# -----------------------------
segments2 <- read_csv(infile_segments_final, show_col_types = FALSE)

req_cols <- c("unique_seg", "season", "first_timestamp", "individual.local.identifier")
missing <- setdiff(req_cols, names(segments2))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

segments2 <- segments2 %>%
  mutate(
    first_timestamp = as.POSIXct(first_timestamp, tz = "UTC"),
    season = tolower(season),
    mig_year = lubridate::year(first_timestamp)
  ) %>%
  filter(
    !is.na(first_timestamp),
    season %in% c("spring", "autumn"),
    !is.na(.data[["individual.local.identifier"]])
  )

# -----------------------------
# STEP 2.3: Define "migration" = bird × season × year
# -----------------------------
per_migration <- segments2 %>%
  group_by(individual.local.identifier, season, mig_year) %>%
  summarise(n_segments = n(), .groups = "drop")

# -----------------------------
# STEP 2.4: Plot (discrete frequency bars; faceted by season)
# -----------------------------
p_segments_per_migration <- ggplot(per_migration, aes(x = factor(n_segments), fill = season)) +
  geom_bar(width = 0.85, colour = "grey25", linewidth = 0.25) +
  facet_wrap(~ season, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    spring = "#2C7FB8",
    autumn = "#D95F0E"
  )) +
  labs(
    x = "Number of flight segments per seasonal migration",
    y = "Count of migrations"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    panel.spacing = unit(0.8, "lines"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.6),
    legend.position = "none"
  )

ggsave(
  file.path(outdir, "segments_per_migration_discrete_by_season.png"),
  p_segments_per_migration, width = 8.8, height = 6.2, dpi = 350
)
ggsave(
  file.path(outdir, "segments_per_migration_discrete_by_season.pdf"),
  p_segments_per_migration, width = 8.8, height = 6.2
)

print(p_segments_per_migration)

# ==============================================================================
# BLOCK 3 — Figure: individual variation (per-individual medians) by season
#           Panels: duration, cumulative distance, straightness
# ==============================================================================

# -----------------------------
# STEP 3.1: Inputs / outputs
# -----------------------------
in_file_individual <- segments_final_file

# -----------------------------
# STEP 3.2: Read + minimal checks
# -----------------------------
segments3 <- read_csv(in_file_individual, show_col_types = FALSE)

req_cols_ind <- c(
  "individual.local.identifier", "season",
  "duration", "cumulative_distance", "straightness"
)
missing_ind <- setdiff(req_cols_ind, names(segments3))
if (length(missing_ind) > 0) {
  stop("Missing required columns in ", in_file_individual, ":\n- ", paste(missing_ind, collapse = "\n- "))
}

segments3 <- segments3 %>%
  filter(season %in% c("spring", "autumn")) %>%
  mutate(season = factor(season, levels = c("autumn", "spring")))

# -----------------------------
# STEP 3.3: Per-individual summaries (median across segments)
# -----------------------------
ind_season <- segments3 %>%
  group_by(individual.local.identifier, season) %>%
  summarise(
    n_segments = n(),
    duration_med = median(duration, na.rm = TRUE),
    dist_med     = median(cumulative_distance, na.rm = TRUE),
    straight_med = median(straightness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(duration_med), is.finite(dist_med), is.finite(straight_med))

# Long format for plotting
ind_long <- ind_season %>%
  dplyr::select(individual.local.identifier, season, duration_med, dist_med, straight_med) %>%
  pivot_longer(
    cols = c(duration_med, dist_med, straight_med),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("duration_med", "dist_med", "straight_med"),
      labels = c(
        "A:  Duration (h)",
        "B:  Cumulative distance (km)",
        "C:  Straightness index"
      )
    )
  )

# -----------------------------
# STEP 3.4: Plot helper + figure
# -----------------------------
base_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    plot.margin = margin(6, 8, 6, 8),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88"),
    panel.grid.minor = element_blank()
  )

p_individual_variation <- ggplot(ind_long, aes(x = season, y = value, fill = season)) +
  geom_violin(trim = TRUE, alpha = 0.55, colour = "grey25", linewidth = 0.25) +
  geom_point(
    aes(colour = season),
    position = position_jitter(width = 0.12, height = 0),
    size = 1.4, alpha = 0.65
  ) +
  stat_summary(fun = median, geom = "point", size = 2.3, colour = "grey15") +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    autumn = "#D95F0E",
    spring = "#2C7FB8"
  )) +
  scale_colour_manual(values = c(
    autumn = "#B84A05",
    spring = "#1F5D86"
  )) +
  labs(x = "Season", y = NULL) +
  base_theme

# -----------------------------
# STEP 3.5: Save figure
# -----------------------------
ggsave(
  file.path(outdir, "FigS_individual_variation_duration_distance_straightness.png"),
  p_individual_variation, width = 11.5, height = 4.2, dpi = 400
)
ggsave(
  file.path(outdir, "FigS_individual_variation_duration_distance_straightness.pdf"),
  p_individual_variation, width = 11.5, height = 4.2
)

print(p_individual_variation)

# ------------------------------------------------------------------------------
# Final messages
# ------------------------------------------------------------------------------
message("Done. All outputs written to: ", normalizePath(outdir))
message("Key files:")
message(" - ", normalizePath(segments_final_file))
message(" - ", normalizePath(seasonal_summary_file))
message(" - ", normalizePath(file.path(outdir, "segments_per_migration_discrete_by_season.png")))
message(" - ", normalizePath(file.path(outdir, "FigS_individual_variation_duration_distance_straightness.png")))