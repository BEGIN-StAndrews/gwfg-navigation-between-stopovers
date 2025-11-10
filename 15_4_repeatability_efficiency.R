# ==============================================================================
# Title: Repeatability of efficiency-route choice (GC vs wind-optimal, DTW ref)
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read segments_all_meta_efficiency_dtw.csv and code WO=1, GC=0
#   2) Keep birds with ≥ 2 segments; print sample sizes and exclusions
#   3) Fit rptR binary ICC on logit link with random intercept (bird ID)
#   4) Print detailed rptR summary and LINK-scale ICC one-liner (with CI & P_perm)
# Inputs:
#   - segments_all_meta_efficiency_dtw.csv
# Outputs:
#   - (console) sample sizes, exclusions, prevalence, balance, ICC summaries
# Usage:
#   Rscript R/15_4_repeatability_efficiency.R

# ==============================================================================

rm(list = ls()); set.seed(42)

# ------------------------------------------------------------------------------
# Load libraries (auto-install if missing)
# ------------------------------------------------------------------------------
pkgs <- c("rptR","dplyr","lme4")
inst <- setdiff(pkgs, rownames(installed.packages()))
if (length(inst)) install.packages(inst, dependencies = TRUE)
suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# Load & prepare data
# ------------------------------------------------------------------------------
segments <- read.csv("segments_all_meta_efficiency_dtw.csv", stringsAsFactors = FALSE)

stopifnot(all(c("best_route","individual.local.identifier") %in% names(segments)))
stopifnot(all(segments$best_route %in% c("GC","wind_optimal")))

# Coding: wind_optimal = 1 (success), GC = 0 (failure)
segments <- segments %>%
  mutate(
    best_route_binary = ifelse(best_route == "wind_optimal", 1L, 0L),
    individual.local.identifier = factor(individual.local.identifier)
  ) %>%  filter(!is.na(best_route_binary), !is.na(individual.local.identifier))

# ------------------------------------------------------------------------------
# Keep birds with ≥ 2 segments overall (report exclusions)
# ------------------------------------------------------------------------------
total_segments_before_filter <- nrow(segments)
segments_ge2 <- segments %>%
  group_by(individual.local.identifier) %>% filter(n() >= 2) %>% ungroup()
excluded_segments <- total_segments_before_filter - nrow(segments_ge2)
excluded_pct <- ifelse(total_segments_before_filter > 0,
                       100 * excluded_segments / total_segments_before_filter, NA_real_)

# ------------------------------------------------------------------------------
# Descriptive prints
# ------------------------------------------------------------------------------
cat("\n========== SAMPLE SIZES ==========\n")
n_birds <- n_distinct(segments_ge2$individual.local.identifier)
n_segments <- nrow(segments_ge2)
seg_per_bird <- as.integer(table(segments_ge2$individual.local.identifier))
cat("Birds (≥2 segments):", n_birds,
    "\nSegments:", n_segments,
    "\nMedian segments/bird:", median(seg_per_bird),
    "\nMean segments/bird:", round(mean(seg_per_bird), 1), "\n")

cat("\n========== EXCLUSIONS DUE TO <2 SEGMENTS PER INDIVIDUAL ==========\n")
cat("Excluded segments:", excluded_segments, 
    sprintf("(%.1f%% of all segments before filtering)\n", excluded_pct))

cat("\n========== PREVALENCE ==========\n")
cat("Overall Pr(wind_optimal):", round(mean(segments_ge2$best_route_binary), 3), "\n")

cat("\n========== BALANCE / SEPARATION ==========\n")
balance <- segments_ge2 %>% group_by(individual.local.identifier) %>%
  summarise(n = n(),
            n_windopt = sum(best_route_binary),
            n_gc = n - n_windopt,
            prop_wind_opt = mean(best_route_binary),
            .groups = "drop")
cat("Birds all GC (all 0):", sum(balance$prop_wind_opt == 0),
    " | Birds all wind-optimal (all 1):", sum(balance$prop_wind_opt == 1), "\n")
print(head(balance, 10))

# ------------------------------------------------------------------------------
# Model controls
# ------------------------------------------------------------------------------
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
B <- 1000; P <- 1000

# ------------------------------------------------------------------------------
# Fit ICC (rptR) model
# ------------------------------------------------------------------------------
cat("\n========== FITTING ICC MODEL ==========\n")

rpt_model <- rpt(
  best_route_binary ~ (1 | individual.local.identifier),
  grname   = "individual.local.identifier",
  data     = segments_ge2,
  datatype = "Binary",
  link     = "logit",
  nboot    = B,
  npermut  = P,
  control  = ctrl,
  parallel = FALSE
)

# ------------------------------------------------------------------------------
# Detailed summaries
# ------------------------------------------------------------------------------
cat("\n========== DETAILED MODEL SUMMARY ==========\n")
print(summary(rpt_model))

# --- Robust helpers for LINK-scale extraction with permutation-p fallback ----
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
  return(NA_real_)
}

# Fallback: parse permutation p from printed summary if summary(obj)$P_permut is NULL
.get_perm_p_link <- function(obj) {
  smry <- summary(obj)
  if (!is.null(smry$P_permut)) {
    p <- .pick_scalar(smry$P_permut, "Link")
    if (!is.na(p)) return(p)
  }
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

print_scales <- function(obj, label) {
  R_link   <- .pick_scalar(obj$R, "Link")
  CI_link  <- .pick_row(obj$CI_boot, "Link"); if (all(is.na(CI_link))) CI_link <- .pick_row(obj$CI_emp, "Link")
  P_perm   <- .get_perm_p_link(obj)
  
  cat(sprintf("%s [Link] R=%.3f CI=[%.3f, %.3f] P_perm=%.2f\n",
              label,
              ifelse(is.na(R_link), NA_real_, R_link),
              ifelse(length(CI_link)>=1, CI_link[1], NA_real_),
              ifelse(length(CI_link)>=2, CI_link[2], NA_real_),
              ifelse(is.na(P_perm), NA_real_, P_perm)))
}

one_liner <- function(obj, label) {
  R_link   <- .pick_scalar(obj$R, "Link")
  CI_link  <- .pick_row(obj$CI_boot, "Link"); if (all(is.na(CI_link))) CI_link <- .pick_row(obj$CI_emp, "Link")
  CI_low   <- ifelse(length(CI_link)>=1, CI_link[1], NA_real_)
  CI_high  <- ifelse(length(CI_link)>=2, CI_link[2], NA_real_)
  P_perm   <- .get_perm_p_link(obj)
  
  vc <- tryCatch(as.data.frame(lme4::VarCorr(obj$mod)), error = function(e) NULL)
  id_var <- if (!is.null(vc)) vc$vcov[vc$grp == "individual.local.identifier"] else NA_real_
  sing  <- tryCatch(lme4::isSingular(obj$mod, tol = 1e-5), error = function(e) NA)
  
  cat(sprintf(
    "%s: ICC (Link/latent) R = %.3f (95%% CI: %.3f–%.3f), P_perm = %.2f; Var(ID)=%.5f%s\n",
    label,
    ifelse(is.na(R_link), NA_real_, R_link),
    ifelse(is.na(CI_low),  NA_real_, CI_low),
    ifelse(is.na(CI_high), NA_real_, CI_high),
    ifelse(is.na(P_perm),  NA_real_, P_perm),
    ifelse(is.na(id_var),  NA_real_, id_var),
    if (isTRUE(sing)) " ; singular fit (Var(ID) ≈ 0)" else ""
  ))
}

# ------------------------------------------------------------------------------
# Output (LINK-scale one-liner + full scales)
# ------------------------------------------------------------------------------
cat("\n========== ONE-LINE SUMMARY (LINK scale) ==========\n")
one_liner(rpt_model,  "Model")

cat("\n========== FULL SCALE PRINT ==========\n")
print_scales(rpt_model,  "Model")

cat("\n[link='logit', nboot=", B, ", npermut=", P,
    ", optimizer='bobyqa'. Outcome: wind_optimal=1, GC=0]\n", sep = "")
