# ==============================================================================
# Title: Route similarity measures
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
# Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi | am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
# 1) Load segment metadata, actual hourly points, actual equal-interval points
# 2) Load simulated routes (5 compass, 2 efficiency)
# 3) Compute 3 metrics per segment × route:
# • med = median geodesic distance (km)
# • dtw = DTW cost / warping-path length (km per step)
# • dir = mean cos(Δbearing) with internal stationarity masking
# 4) Summarise winners (per segment) within compass and efficiency sets
# 5) Save per-segment tables (all + seasonal splits) and set-level summaries
# Inputs:
# - segments_final.csv
# - points_resampled.csv 
# - actual_equal_interval.csv 
# - Geographic_routes.csv
# - Geomagnetic_routes.csv
# - Magnetoclinic_routes.csv
# - SunCompass_classic_routes.csv
# - local_wind_routes.csv
# - GreatCircle_equal_interval.csv
# - wind_optimal_equal_interval.csv
# Outputs:
# - evaluation_metrics.csv
# - summary_compass_routes.csv
# - summary_efficiency_routes.csv
# - per_segment_compass_metrics.csv
# - per_segment_compass_metrics_autumn.csv
# - per_segment_compass_metrics_spring.csv
# - per_segment_efficiency_metrics.csv
# - per_segment_efficiency_metrics_autumn.csv
# - per_segment_efficiency_metrics_spring.csv
# Usage:
# Rscript R/10_2_Route_similarity_measures.R

# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(geosphere)
library(dtw)
library(readr)
library(dplyr)
library(purrr)
library(rlang)
library(tidyr)

# ------------------------------------------------------------------------------
# Parameters & route sets
# ------------------------------------------------------------------------------
D_min_km <- 5
route_types_compass <- c("geographic","geomagnetic","magnetoclinic","sun","local_wind")
route_types_efficiency <- c("GC","wind_optimal")
route_types_all <- c(route_types_compass, route_types_efficiency)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
paths <- list(
  segments = "segments_final.csv",
  hourly_actual = "points_resampled.csv",
  eq_actual = "actual_equal_interval.csv",
  geomagnetic = "Geomagnetic_routes.csv",
  geographic = "Geographic_routes.csv",
  magnetoclinic = "Magnetoclinic_routes.csv",
  sun = "SunCompass_classic_routes.csv",
  local_wind = "local_wind_routes.csv",
  GC = "GreatCircle_equal_interval.csv",
  wind_optimal = "wind_optimal_equal_interval.csv")

# ------------------------------------------------------------------------------
# Helper functions
# - compute_med_geo_distance(): median geodesic error (km)
# - compute_dtw_position_geo(): DTW cost normalized by path length (km/step)
# - compute_dir_idx(): directional similarity with stationarity masking
# ------------------------------------------------------------------------------
compute_med_geo_distance <- function(real_df, sim_df) {
    n <- min(nrow(real_df), nrow(sim_df))
    pts_r <- real_df[1:n, c("lon","lat")]
    pts_s <- sim_df[1:n, c("lon","lat")]
    d     <- distHaversine(pts_r, pts_s)
    median(d, na.rm=TRUE)/1000
  }
  
compute_dtw_position_geo <- function(real_df, sim_df) {
    mat_r <- as.matrix(real_df[, c("lon","lat")])
    mat_s <- as.matrix(sim_df[, c("lon","lat")])
    C <- geosphere::distm(mat_r, mat_s, fun = geosphere::distHaversine)
    alignment <- dtw::dtw(
      x             = C,                  # your precomputed cost-matrix
      y             = NULL,               # tells dtw() C is already the local costs
      step.pattern  = dtw::symmetric1,    # your chosen step pattern
      keep.internals = FALSE,             # no need to retain the full cost-matrix
      distance.only = FALSE               # ensure backtracking & index1/index2 are computed
    )
    
    (alignment$distance / length(alignment$index1)) / 1000
    
  }

compute_dir_idx <- function(real_df, sim_df) {
    n <- min(nrow(real_df), nrow(sim_df))
    if (n < 2) return(NA_real_)
    rd <- real_df[1:n, ]
    sd <- sim_df[1:n, ]
    d_m   <- distHaversine(rd[-1, c("lon","lat")], rd[-n, c("lon","lat")]) / 1000
    is_stat <- d_m <= D_min_km
    br <- bearing(rd[-n, c("lon","lat")], rd[-1, c("lon","lat")])
    bs <- bearing(sd[-n, c("lon","lat")], sd[-1, c("lon","lat")])
    cosd <- cos((br - bs) * pi/180)
    cosd[is_stat] <- 0
    mean(cosd, na.rm = TRUE)
  }
  
# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------
segments_df   <- read.csv(paths$segments,  stringsAsFactors=FALSE)
actual_hourly <- read.csv(paths$hourly_actual,  stringsAsFactors=FALSE) %>%arrange(unique_seg, timestamp)
actual_equal  <- read.csv(paths$eq_actual,    stringsAsFactors=FALSE)%>%arrange(unique_seg, idx)
  
route_names <- c("geomagnetic","geographic","magnetoclinic","sun","local_wind","GC","wind_optimal")
  
route_dfs <- setNames(
    lapply(route_names, function(rt) {
      read_csv(paths[[rt]], show_col_types = FALSE)
    }),
    route_names)
  
# ------------------------------------------------------------------------------
# Filter to segments available across all compass routes
# ------------------------------------------------------------------------------
available_segs <- unique(actual_hourly$unique_seg)
for (rt in route_types_compass) {
  available_segs <- intersect(available_segs, unique(route_dfs[[rt]]$unique_seg))
}

actual_hourly <- filter(actual_hourly, unique_seg %in% available_segs)
actual_equal <- filter(actual_equal, unique_seg %in% available_segs)
segments_df <- filter(segments_df, unique_seg %in% available_segs)
route_dfs <- lapply(route_dfs, function(df) filter(df, unique_seg %in% available_segs))

# ------------------------------------------------------------------------------
# Compute metrics (per segment × route)
# ------------------------------------------------------------------------------
results <- list()
  cnt <- 1L
  
  for (seg in available_segs) {
    for (rt in route_types_all) {
      # choose real reference
      real_df <- if (rt %in% route_types_compass) actual_hourly else actual_equal
      real_sub <- filter(real_df, unique_seg==seg)
      sim_sub  <- filter(route_dfs[[rt]], unique_seg==seg) 
      
      results[[cnt]] <- tibble(
        unique_seg = seg,
        route      = rt,
        med        = compute_med_geo_distance(real_sub, sim_sub),
        dtw        = compute_dtw_position_geo(real_sub, sim_sub),
        dir        = compute_dir_idx(real_sub, sim_sub)
      )
      cnt <- cnt + 1L
    }
  }
  
all_eval <- bind_rows(results)
write.csv(all_eval, "evaluation_metrics.csv", row.names = FALSE)
  
  
  
# ------------------------------------------------------------------------------
# Set-level summaries: winner proportion per metric
# ------------------------------------------------------------------------------ 
n_segs   <- n_distinct(all_eval$unique_seg)
  
metrics_df <- tribble(
    ~metric, ~smaller_is_better,
    "med",    TRUE,
    "dtw",    TRUE,
    "dir",   FALSE)
  
summarize_set <- function(eval_df, route_set) {
    map_dfr(metrics_df$metric, function(m) {
      sib <- metrics_df$smaller_is_better[metrics_df$metric==m]
      sub <- filter(eval_df, route %in% route_set)
      best <- if (sib) {
        sub %>% group_by(unique_seg) %>% slice_min(.data[[m]], with_ties=FALSE)
      } else {
        sub %>% group_by(unique_seg) %>% slice_max(.data[[m]], with_ties=FALSE)
      }
      best %>% ungroup() %>%
        count(route, name="wins") %>%
        mutate(metric=m, proportion=wins/n_segs) %>%
        dplyr::select(metric, route, proportion)
    }) %>%
      pivot_wider(names_from=route, values_from=proportion, values_fill=0)
  }
  
summary_compass    <- summarize_set(all_eval, route_types_compass)
summary_efficiency <- summarize_set(all_eval, route_types_efficiency)
  
write.csv(summary_compass,    "summary_compass_routes.csv", row.names = FALSE)
write.csv(summary_efficiency, "summary_efficiency_routes.csv", row.names = FALSE)
  
# ------------------------------------------------------------------------------
# Save per-segment metrics + metadata (all & seasonal splits)
# ------------------------------------------------------------------------------
seg_meta <- segments_df %>%
  transmute(
    unique_seg,
    bird_id               = individual.local.identifier,
    season,
    cumulative_distance,
    great_circle_distance,
    Linearity             = great_circle_distance / cumulative_distance,
    MaxKp,
    duration
  )

# Compass per-segment
compass_perseg <- all_eval %>%
  filter(route %in% route_types_compass) %>%      # keep only the 5 compass routes
  left_join(seg_meta, by = "unique_seg") %>%      # attach metadata
  arrange(unique_seg, route)
write.csv(compass_perseg, "per_segment_compass_metrics.csv", row.names = FALSE)

compass_perseg_autumn <- compass_perseg %>%filter(season == "autumn") 
compass_perseg_spring <- compass_perseg %>%filter(season == "spring") 

write.csv(compass_perseg_autumn, "per_segment_compass_metrics_autumn.csv", row.names = FALSE)
write.csv(compass_perseg_spring, "per_segment_compass_metrics_spring.csv", row.names = FALSE)

# Efficiency per-segment
eff_perseg <- all_eval %>%
  filter(route %in% route_types_efficiency) %>%   # GC + wind_optimal
  left_join(seg_meta, by = "unique_seg") %>%      # attach metadata
  arrange(unique_seg, route)

write.csv(eff_perseg, "per_segment_efficiency_metrics.csv", row.names = FALSE)

eff_perseg_autumn <- eff_perseg %>% filter(season == "autumn") 
eff_perseg_spring <- eff_perseg %>% filter(season == "spring") 

write.csv(eff_perseg_autumn, "per_segment_efficiency_metrics_autumn.csv", row.names = FALSE)
write.csv(eff_perseg_spring, "per_segment_efficiency_metrics_spring.csv", row.names = FALSE)
