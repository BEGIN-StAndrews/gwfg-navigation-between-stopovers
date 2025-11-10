# ==============================================================================
# Title: Repeatability of compass-route choice (one-vs-rest, DTW reference)
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read segments_all_meta_compass_dtw.csv and keep compass routes only
#   2) Filter to birds with ≥ 2 segments overall (report exclusions)
#   3) For each compass strategy, fit rptR binary (one-vs-rest) model with logit link
#   4) Extract LINK-scale ICC, bootstrap/empirical CI, and permutation p-values
#   5) Save a tidy summary CSV and an ICC errorbar plot
# Inputs:
#   - segments_all_meta_compass_dtw.csv
# Outputs:
#   - compass_rpt/compass_repeatability_summary.csv
#   - compass_rpt/compass_repeatability_ICC.png
#   - compass_rpt/compass_repeatability_ICC.pdf
# Usage:
#   Rscript R/15_3_repeatability_compass.R
# ==============================================================================

rm(list = ls()); set.seed(42)

# ------------------------------------------------------------------------------
# Load libraries (auto-install if missing)
# ------------------------------------------------------------------------------
pkgs <- c("dplyr","readr","rptR","lme4","ggplot2","tidyr","purrr","tibble")
inst <- setdiff(pkgs, rownames(installed.packages()))
if (length(inst)) install.packages(inst, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# I/O
# ------------------------------------------------------------------------------
INFILE <- "segments_all_meta_compass_dtw.csv"
OUTDIR <- "compass_rpt"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# ------------------------------------------------------------------------------
# Canonical order & labels
# ------------------------------------------------------------------------------
route_types_compass <- c("geographic","geomagnetic","magnetoclinic","sun","local_wind")
label_map <- c(
  geographic    = "GL",
  geomagnetic   = "ML",
  magnetoclinic = "MC",
  sun           = "SC",
  local_wind    = "LW")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
segments <- read_csv(INFILE, show_col_types = FALSE)

req <- c("individual.local.identifier","best_route")
stopifnot(all(req %in% names(segments)))
stopifnot(all(segments$best_route %in% route_types_compass))

segments <- segments %>%
  filter(best_route %in% route_types_compass) %>%
  mutate(
    id         = factor(individual.local.identifier),
    best_route = factor(best_route, levels = route_types_compass)
  )

# ------------------------------------------------------------------------------
# Keep birds with ≥ 2 segments overall (report exclusions)
# ------------------------------------------------------------------------------
total_segments_before_filter <- nrow(segments)
segments_ge2 <- segments %>%
  group_by(id) %>% filter(n() >= 2) %>% ungroup()
excluded_segments <- total_segments_before_filter - nrow(segments_ge2)
excluded_pct <- ifelse(total_segments_before_filter > 0,
                       100 * excluded_segments / total_segments_before_filter, NA_real_)

cat("Birds (≥2 segments overall):", n_distinct(segments_ge2$id),
    "| Segments:", nrow(segments_ge2), "\n")

cat("\n========== EXCLUSIONS DUE TO <2 SEGMENTS PER INDIVIDUAL ==========\n")
cat("Excluded segments:", excluded_segments,
    sprintf("(%.1f%% of all segments before filtering)\n", excluded_pct))

# ------------------------------------------------------------------------------
# Helpers to extract LINK-scale ICC, CI, perm-p
# ------------------------------------------------------------------------------
.pick_row <- function(x, target = "Link") {
  if (is.null(x)) return(NA_real_)
  if (is.matrix(x) || is.data.frame(x)) {
    rn <- rownames(x)
    if (!is.null(rn) && any(rn == target)) return(as.numeric(x[which(rn == target)[1], ]))
    return(as.numeric(x[1, ]))
  }
  if (is.list(x)) {
    if (!is.null(x[[target]])) return(as.numeric(x[[target]]))
    return(as.numeric(x[[1]]))
  }
  return(as.numeric(x))
}
.pick_scalar <- function(x, target = "Link") {
  v <- .pick_row(x, target)
  if (length(v)) return(v[1])
  NA_real_
}
.get_perm_p_link <- function(obj) {
  smry <- summary(obj)
  if (!is.null(smry$P_permut)) {
    p <- .pick_scalar(smry$P_permut, "Link")
    if (!is.na(p)) return(p)
  }
  # fallback: parse printed output
  txt <- capture.output(print(smry))
  idx <- grep("^Repeatability estimation overview", txt)
  if (length(idx)) {
    block <- txt[(idx+1):min(length(txt), idx+20)]
    link_line_idx <- grep("^\\s*Link\\b", block)
    if (length(link_line_idx)) {
      line <- block[link_line_idx[1]]
      toks <- strsplit(line, "\\s+")[[1]]
      nums <- suppressWarnings(as.numeric(toks))
      nums <- nums[!is.na(nums)]
      if (length(nums)) return(tail(nums, 1))
    }
  }
  NA_real_
}
extract_link_icc <- function(obj) {
  if (is.null(obj)) return(c(R = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P_perm = NA_real_))
  R_link   <- .pick_scalar(obj$R, "Link")
  CI_link  <- .pick_row(obj$CI_boot, "Link")
  if (all(is.na(CI_link))) CI_link <- .pick_row(obj$CI_emp, "Link")
  CI_low  <- ifelse(length(CI_link) >= 1, CI_link[1], NA_real_)
  CI_high <- ifelse(length(CI_link) >= 2, CI_link[2], NA_real_)
  P_perm   <- .get_perm_p_link(obj)
  c(R = R_link, CI_low = CI_low, CI_high = CI_high, P_perm = P_perm)
}

# ------------------------------------------------------------------------------
# Fit function (single model, no season)
# ------------------------------------------------------------------------------
fit_repeatability <- function(data, focal_route, B = 1000, P = 1000) {
  dat <- data %>%
    mutate(y = as.integer(best_route == focal_route)) %>%
    select(id, y) %>%
    tidyr::drop_na()
  
  prop1 <- mean(dat$y)
  if (prop1 %in% c(0, 1)) {
    return(list(
      focal = focal_route, prop1 = prop1,
      model = NULL,
      note = "Degenerate (all 0 or all 1); no ICC estimable."
    ))
  }
  
  ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  
  model <- tryCatch(
    rpt(y ~ (1 | id),
        grname   = "id",
        data     = dat,
        datatype = "Binary",
        link     = "logit",
        nboot    = B,
        npermut  = P,
        control  = ctrl,
        parallel = FALSE),
    error = function(e) NULL
  )
  
  list(focal = focal_route, prop1 = prop1,
       model = model, note = NA_character_)
}

# ------------------------------------------------------------------------------
# Run per-strategy
# ------------------------------------------------------------------------------
results <- purrr::map(route_types_compass, ~fit_repeatability(segments_ge2, .x))
if (!exists("results") || !exists("segments_ge2") || !exists("label_map")) {
  stop("❌ Required objects missing: make sure 'results', 'segments_ge2', and 'label_map' exist in memory.\n",
       "Re-run the ICC modelling section first.")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------------------
# Build tidy summary table
# ------------------------------------------------------------------------------
tab <- tibble::tibble(
  route      = purrr::map_chr(results, "focal"),
  label      = unname(label_map[route]),
  prevalence = purrr::map_dbl(results, "prop1"),
  R          = purrr::map_dbl(results, ~extract_link_icc(.x$model)["R"]),
  CI_low     = purrr::map_dbl(results, ~extract_link_icc(.x$model)["CI_low"]),
  CI_high    = purrr::map_dbl(results, ~extract_link_icc(.x$model)["CI_high"]),
  P_perm     = purrr::map_dbl(results, ~extract_link_icc(.x$model)["P_perm"]),
  note       = purrr::map_chr(results, ~.x$note %||% "")
)

counts <- segments_ge2 %>%
  group_by(best_route) %>%
  summarise(n_segs = dplyr::n(), .groups = "drop") %>%
  rename(route = best_route) %>%
  mutate(route = as.character(route))

tab <- tab %>%
  left_join(counts, by = "route") %>%
  mutate(n_birds = n_distinct(segments_ge2$id)) %>%
  arrange(match(route, route_types_compass))

# ------------------------------------------------------------------------------
# Save CSV
# ------------------------------------------------------------------------------
write_csv(tab, file.path(OUTDIR, "compass_repeatability_summary.csv"))
print(tab)

# ------------------------------------------------------------------------------
# Plot ICCs (single model)
# ------------------------------------------------------------------------------
plot_df <- tab %>%
  mutate(
    strategy = factor(label, levels = unname(label_map[route_types_compass])),
    R = R, lo = CI_low, hi = CI_high)

hi_max <- suppressWarnings(max(plot_df$hi, na.rm = TRUE))
hi_max <- ifelse(is.finite(hi_max), hi_max, 0.001)

p <- ggplot(plot_df, aes(x = strategy, y = R)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "#bbbbbb") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.20, linewidth = 0.6) +
  geom_point(size = 2.8) +
  coord_cartesian(ylim = c(0, max(0.001, hi_max))) +
  labs(
    x = "Compass option (one-vs-rest)",
    y = "Repeatability (ICC, latent/logit)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUTDIR, "compass_repeatability_ICC.png"), p,
       width = 150, height = 95, units = "mm", dpi = 300)
ggsave(file.path(OUTDIR, "compass_repeatability_ICC.pdf"), p,
       width = 150, height = 95, units = "mm")

cat("✔ Table and plots saved to:", normalizePath(OUTDIR), "\n")
