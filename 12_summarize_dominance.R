# ==============================================================================
# Title: Segment-level dominance & runner-up summaries
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
# Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi | am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
# 1) Build per-segment dominance table by metric and level (compass/efficiency)
# 2) Attach segment metadata and export seasonal splits
# 3) Produce runner-up tables (compass–spring–DIR; efficiency–autumn–DIR)
# 4) Produce detailed runner-up breakdown for compass–autumn–DTW
# Inputs:
# - evaluation_metrics.csv (unique_seg × route → med, dtw, dir)
# - segments_final.csv (segment metadata)
# - points_resampled_filtered.csv (per-point flags for stopovers)
# Outputs:
# - dominance_by_segment.csv
# - dominance_by_segment_autumn.csv
# - dominance_by_segment_spring.csv
# - compass_spring_DIR_runnerup.csv (summary table)
# - efficiency_autumn_DIR_runnerup.csv (summary table)
# - compass_autumn_DTW_runnerup_breakdown.csv
# Usage:
# Rscript R/12_summarize_dominance.R

# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(rlang)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
paths <- list(
  eval_metrics = "evaluation_metrics.csv",
  seg_table = "segments_final.csv",
  hourly_points = "points_resampled_filtered.csv")

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------
all_eval <- read_csv("evaluation_metrics.csv", show_col_types = FALSE)

seg_meta1 <- read_csv("segments_final.csv", show_col_types = FALSE) %>%
  transmute(
    unique_seg,
    bird_id               = individual.local.identifier,
    season,
    cumulative_distance,
    great_circle_distance,
    Linearity             = great_circle_distance / cumulative_distance,
    MaxKp,
    duration )

seg_meta2 <- read_csv("points_resampled_filtered.csv", show_col_types = FALSE) %>%
  distinct(
    unique_seg,
    segment_type,
    contain_short_stop,
    previous_stopover_duration,
    next_stopover_duration)

seg_meta <- seg_meta1 %>% left_join(seg_meta2, by = "unique_seg")

# ------------------------------------------------------------------------------
# Route sets (levels) & metric meta
# ------------------------------------------------------------------------------
route_types_compass    <- c("geographic","geomagnetic","magnetoclinic","sun","local_wind")
route_types_efficiency <- c("GC","wind_optimal")

levels_list <- list(
  compass    = route_types_compass,
  efficiency = route_types_efficiency
)

metrics_all <- tribble(
  ~metric, ~smaller_is_better,
  "med",    TRUE,
  "dtw",    TRUE,
  "dir",   FALSE
)

# ------------------------------------------------------------------------------
# Compute dominance (best & runner-up) per unique_seg × metric × level
# ------------------------------------------------------------------------------
dominance_df <- imap_dfr(levels_list, function(route_set, lvl) {
  all_eval %>%
    filter(route %in% route_set) %>%
    pivot_longer(
      cols      = c(med, dtw, dir),
      names_to  = "metric",
      values_to = "score"
    ) %>%
    left_join(metrics_all, by = "metric") %>%
    group_by(unique_seg, metric) %>%
    mutate(
      rank = if_else(
        smaller_is_better,
        rank(score, ties.method = "first"),
        rank(-score, ties.method = "first")
      )
    ) %>%
    filter(rank <= 2) %>%
    summarise(
      level               = lvl,
      best_route          = route[rank == 1],
      second_route        = route[rank == 2],
      v1                  = score[rank == 1],
      v2                  = score[rank == 2],
      smaller_is_better   = first(smaller_is_better),
      .groups             = "drop"
    ) %>%
    mutate(
      raw_diff  = if_else(smaller_is_better, v2 - v1, v1 - v2),
      # imp_pct   = if_else(v2 == 0, NA_real_, 100 * raw_diff / v2),
      
      imp_pct   = if_else(v2 == 0, NA_real_, 100 * raw_diff / abs(v2)),
      imp_perkm = raw_diff / seg_meta1$cumulative_distance[match(unique_seg, seg_meta1$unique_seg)]
    )
})

# ------------------------------------------------------------------------------
# Attach metadata & export dominance tables (all + seasonal splits)
# ------------------------------------------------------------------------------
dominance_with_meta <- dominance_df %>%
  left_join(seg_meta, by = "unique_seg") %>%
  dplyr::select(
    unique_seg, bird_id, season, level, metric,
    best_route, second_route,
    v1, v2,
    cumulative_distance, great_circle_distance, Linearity, MaxKp, duration,
    segment_type, contain_short_stop, previous_stopover_duration, next_stopover_duration,
    raw_diff, imp_pct, imp_perkm
  ) %>% arrange(unique_seg, metric, level)

## 8) Write out ---------------------------------------------------------------
write_csv(dominance_with_meta, "dominance_by_segment.csv")
message("✔︎ Written dominance_by_segment.csv (8 rows per segment: 2 levels × 4 metrics).")

dominance_with_meta_autumn <- dominance_with_meta %>% filter(season == "autumn") 
dominance_with_meta_spring <- dominance_with_meta %>% filter(season == "spring") 

write_csv(dominance_with_meta_autumn, "dominance_by_segment_autumn.csv")
write_csv(dominance_with_meta_spring, "dominance_by_segment_spring.csv")

# ------------------------------------------------------------------------------
# Runner-up summary: COMPASS × SPRING × DIR
# ------------------------------------------------------------------------------
dominance_df <- read_csv("dominance_by_segment.csv", show_col_types = FALSE)

final_table <- dominance_df %>% filter(level == "compass", season == "spring", metric == "dir") %>%
group_by(best_route, second_route) %>%
  summarise(
    runnerup_n    = n(),
    mean_imp_top  = mean(imp_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
left_join(
    dominance_df %>%
      filter(level == "compass", season == "spring", metric == "dir") %>%
      group_by(best_route) %>%
      summarise(
        n_segments    = n(),
        mean_imp_all  = mean(imp_pct, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "best_route"
  ) %>%  group_by(best_route) %>%
  slice_max(runnerup_n, n = 1, with_ties = FALSE) %>%  summarise(
    `Best-Fitting Strategy`                = best_route,
    `Top Second-Best`                      = second_route,
    `Top Second-Best Frequency (%)`        = round(100 * runnerup_n / n_segments, 1),
    `Mean Improvement Over Top Second-Best` = round(mean_imp_top, 1),
    `Mean Improvement Over All Second-Best` = round(mean_imp_all, 1),
    .groups = "drop"
  )
print(final_table)

# ------------------------------------------------------------------------------
# Runner-up summary: EFFICIENCY × AUTUMN × DIR
# ------------------------------------------------------------------------------
dominance_df <- read_csv("dominance_by_segment.csv", show_col_types = FALSE)

final_table_efficiency <- dominance_df %>%filter(level == "efficiency", season == "autumn", metric == "dir") %>%
  
group_by(best_route, second_route) %>%  summarise(
    runnerup_n    = n(),
    mean_imp_top  = mean(imp_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
    left_join(
    dominance_df %>%
      filter(level == "efficiency", season == "autumn", metric == "dir") %>%
      group_by(best_route) %>%
      summarise(
        n_segments    = n(),
        mean_imp_all  = mean(imp_pct, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "best_route"
  ) %>%  group_by(best_route) %>%
  slice_max(runnerup_n, n = 1, with_ties = FALSE) %>%  summarise(
    `Best-Fitting Strategy`                 = best_route,
    `Top Second-Best`                       = second_route,
    `Top Second-Best Frequency (%)`         = round(100 * runnerup_n / n_segments, 1),
    `Mean Improvement Over Top Second-Best` = round(mean_imp_top, 1),
    `Mean Improvement Over All Second-Best` = round(mean_imp_all, 1),
    .groups = "drop"
  )
print(final_table_efficiency)

# ------------------------------------------------------------------------------
# Detailed runner-up breakdown: COMPASS × AUTUMN × DTW (per runner-up)
# ------------------------------------------------------------------------------
dominance_df <- read_csv("dominance_by_segment.csv", show_col_types = FALSE)

df <- dominance_df %>%filter(level == "compass", season == "autumn", metric == "dtw")
total_autumn_compass_dtw <- nrow(df)
compass_routes <- c("geographic", "geomagnetic", "magnetoclinic", "sun", "local_wind")

result <- map_dfr(compass_routes, function(main_strategy) {
  df_main <- df %>% filter(best_route == main_strategy)
  n_segments <- nrow(df_main)
  win_pct <- 100 * n_segments / total_autumn_compass_dtw
  
  runnerup_stats <- map_dfr(setdiff(compass_routes, main_strategy), function(runner_up) {
    df_pair <- df_main %>% filter(second_route == runner_up)
    runnerup_n <- nrow(df_pair)
    runnerup_pct <- ifelse(n_segments > 0, 100 * runnerup_n / n_segments, NA)
    mean_imp <- mean(df_pair$imp_pct, na.rm = TRUE)
    
    tibble(
      strategy       = main_strategy,
      n_segments     = n_segments,
      win_pct        = round(win_pct, 1),
      runner_up      = runner_up,
      runnerup_n     = runnerup_n,
      runnerup_pct   = round(runnerup_pct, 1),
      mean_imp_pct   = round(mean_imp, 1)
    )
  })
  
  return(runnerup_stats)
})
print(result, n = Inf)
