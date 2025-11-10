# ==============================================================================
# Title: Binary logistic GLM (DTW reference) with diagnostics
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read segments_all_meta_dtw.csv; set factors and collapse light_start to day/night
#   2) Fit full binomial GLM (GC vs wind_optimal) + stepwise AIC to final additive model
#   3) Compute diagnostics: McFadden R2, Tjur R2, AUC, Hosmer–Lemeshow, VIF
#   4) Report coefficients, odds ratios with 95% profile-likelihood CI, and AIC path
#   5) Write text report to GLM_output/GLM_binomial_report.txt
# Inputs:
#   - segments_all_meta_dtw.csv
# Outputs:
#   - GLM_output/GLM_binomial_report.txt
# Usage:
#   Rscript R/15_1_glm_binary.R
# ==============================================================================
options(width = 120)

# ------------------------------------------------------------------------------
# Load libraries (auto-install if missing)
# ------------------------------------------------------------------------------
req <- c("dplyr","readr","forcats","car","pROC","ResourceSelection","broom","tibble")
inst <- setdiff(req, rownames(installed.packages()))
if (length(inst)) install.packages(inst, dependencies = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# Paths & files
# ------------------------------------------------------------------------------
seg_fn <- "segments_all_meta_dtw.csv"
outdir <- "GLM_output"; dir.create(outdir, showWarnings = FALSE)
logf   <- file.path(outdir, "GLM_binomial_report.txt")

# ------------------------------------------------------------------------------
# Load data and set factors
# ------------------------------------------------------------------------------
seg_fn <- "segments_all_meta_dtw.csv"
outdir <- "GLM_output"; dir.create(outdir, showWarnings = FALSE)
logf   <- file.path(outdir, "GLM_binomial_report.txt")

segments <- readr::read_csv(seg_fn, show_col_types = FALSE) |>
  dplyr::mutate(
    best_route = factor(best_route, levels = c("GC","wind_optimal")), # GC=0, WO=1
    season     = forcats::fct_relevel(factor(season), "spring","autumn"),
    contain_short_stop = forcats::fct_relevel(factor(contain_short_stop), "FALSE","TRUE"),
    geomag_storm       = factor(geomag_storm),
    light_start        = factor(light_start),
    segment_type       = factor(segment_type))

# ------------------------------------------------------------------------------
# Collapse polar_day → day (two levels: day/night)
# ------------------------------------------------------------------------------
segments$light_start <- forcats::fct_collapse(
  segments$light_start,
  day   = c("day","polar_day"),
  night = "night")

# ------------------------------------------------------------------------------
# Model specification + stepwise AIC to final additive model
# ------------------------------------------------------------------------------
full_form <- best_route ~ season + mean_wind_support + mean_crosswind +
  geomag_storm + light_start + great_circle_distance + duration + segment_type +
  contain_short_stop + previous_stopover_duration + next_stopover_duration

mod_full <- glm(full_form, family = binomial, data = segments)
model    <- step(mod_full, trace = FALSE)  # final additive model

# ------------------------------------------------------------------------------
# Align NA handling & extract fitted components
# ------------------------------------------------------------------------------
mf   <- model.frame(model)                               
p_hat <- fitted(model)                                   
obs   <- as.numeric(mf$best_route == "wind_optimal")    

# ------------------------------------------------------------------------------
# Diagnostics: pseudo-R2, AUC, HL test, VIF
# ------------------------------------------------------------------------------
mf_r2   <- 1 - model$deviance / model$null.deviance
tjur_r2 <- mean(p_hat[obs == 1]) - mean(p_hat[obs == 0])

roc_obj <- pROC::roc(response = obs, predictor = p_hat, quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))

hl <- ResourceSelection::hoslem.test(x = obs, y = p_hat, g = 10)

vif_tbl <- tryCatch({
  v <- car::vif(model)
  if (is.matrix(v)) {
    as.data.frame(v) |>
      tibble::rownames_to_column("term") |>
      dplyr::rename(VIF = dplyr::last_col())
  } else {
    tibble::tibble(term = names(v), VIF = as.numeric(v))
  }
}, error = function(e) tibble::tibble())

# ------------------------------------------------------------------------------
# Odds ratios (native per-unit) with 95% profile-likelihood CI
# ------------------------------------------------------------------------------
ci_prof <- suppressMessages(confint(model))                 # log-odds scale
coefs   <- broom::tidy(model)                               # estimate/SE/z/p

row_match <- match(coefs$term, rownames(ci_prof))
or_tbl <- coefs |>
  dplyr::mutate(
    OR      = exp(estimate),
    CI_low  = exp(ci_prof[,1][row_match]),
    CI_high = exp(ci_prof[,2][row_match])
  ) |>  dplyr::select(term, estimate, std.error, statistic, p.value, OR, CI_low, CI_high)

# ------------------------------------------------------------------------------
# Write text report
# ------------------------------------------------------------------------------
if (file.exists(logf)) file.remove(logf)
sink(logf)

cat("Final additive model — ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
cat("Formula:\n")
print(formula(model)); cat("\n")

cat("Residual deviance:", round(model$deviance, 1), "on", df.residual(model), "df\n")
cat("AIC:", round(AIC(model), 1), "\n")
cat("McFadden R2:", round(mf_r2, 3), "\n")
cat("Tjur R2:", round(tjur_r2, 3), "\n")
cat("AUC:", round(auc_val, 3), "\n")
cat(sprintf("Hosmer–Lemeshow: HL X2=%.3f, df=%d, p=%.3f\n\n", hl$statistic, hl$parameter, hl$p.value))

cat("Coefficient summary (log-odds):\n")
print(summary(model)$coefficients); cat("\n")

cat("Odds ratios (native per-unit) with 95% profile-likelihood CI:\n")
print(or_tbl, row.names = FALSE, digits = 4); cat("\n")

if (nrow(vif_tbl)) {
  cat("VIF (final model):\n")
  print(vif_tbl, row.names = FALSE, digits = 4); cat("\n")
}

cat("Stepwise AIC path:\n")
print(model$anova); cat("\n")

sink()

message("Wrote one report: ", normalizePath(logf))
