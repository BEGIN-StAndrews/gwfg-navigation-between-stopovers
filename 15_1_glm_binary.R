# ==============================================================================
# Title: Binary logistic GLM (DTW reference) with cluster-robust inference
#        and GLMM sensitivity check
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
#
# Summary:
#   1) Read segments_all_meta_dtw.csv; set factors and collapse light_start to day/night
#   2) Fit full binomial GLM (GC vs wind_optimal) + stepwise AIC to final additive model
#   3) Compute GLM diagnostics (McFadden R2, Tjur R2, AUC, Hosmer–Lemeshow, VIF)
#   4) Estimate cluster-robust coefficient table (goose-ID clustered SEs; main inference)
#   5) Fit GLMM random-intercept sensitivity model (same fixed effects as final GLM)
#   6) Report GLMM singularity / random-intercept variance and fixed effects
#   7) Write one text report and compact outputs to BINOMIAL_output/
#
# Inputs:
#   - segments_all_meta_dtw.csv
#
# Outputs:
#   - BINOMIAL_output/binomial_glm_cluster_robust_report.txt
#   - BINOMIAL_output/binomial_glm_cluster_robust_coefficients.csv
#   - BINOMIAL_output/binomial_glmm_sensitivity_coefficients.csv
#   - BINOMIAL_output/model_binomial_glm_final.rds
#   - BINOMIAL_output/model_binomial_glmm_sensitivity.rds
#
# Usage:
#   Rscript R/15_1_glm_binary_cluster_robust.R
# ==============================================================================

options(width = 120)
set.seed(1)

# ------------------------------------------------------------------------------
# Load libraries 
# ------------------------------------------------------------------------------
req <- c(
  "dplyr", "readr", "forcats", "tibble", "broom",
  "car", "pROC", "ResourceSelection",
  "sandwich", "lmtest", "lme4"
)
inst <- setdiff(req, rownames(installed.packages()))
if (length(inst)) install.packages(inst, dependencies = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# Paths & files
# ------------------------------------------------------------------------------
seg_fn <- "segments_all_meta_dtw.csv"

outdir <- "BINOMIAL_output"
dir.create(outdir, showWarnings = FALSE)

logf      <- file.path(outdir, "binomial_glm_cluster_robust_report.txt")
tab_rob   <- file.path(outdir, "binomial_glm_cluster_robust_coefficients.csv")
tab_glmm  <- file.path(outdir, "binomial_glmm_sensitivity_coefficients.csv")
rds_glm   <- file.path(outdir, "model_binomial_glm_final.rds")
rds_glmm  <- file.path(outdir, "model_binomial_glmm_sensitivity.rds")

# ------------------------------------------------------------------------------
# Load data and set factors
# ------------------------------------------------------------------------------
segments <- readr::read_csv(seg_fn, show_col_types = FALSE) |>
  dplyr::mutate(
    best_route = factor(best_route, levels = c("GC", "wind_optimal")), # GC=0, WO=1
    season = forcats::fct_relevel(factor(season), "spring", "autumn"),
    contain_short_stop = forcats::fct_relevel(factor(contain_short_stop), "FALSE", "TRUE"),
    geomag_storm = factor(geomag_storm),
    light_start = factor(light_start),
    segment_type = factor(segment_type)
  )

# ------------------------------------------------------------------------------
# ID column
# ------------------------------------------------------------------------------
id_col <- "individual.local.identifier"

if (!id_col %in% names(segments)) {
  stop(
    "Required ID column not found: ", id_col, "\n",
    "Available columns include: ",
    paste(head(names(segments), 40), collapse = ", "),
    if (length(names(segments)) > 40) " ... (truncated)" else ""
  )
}
segments[[id_col]] <- factor(segments[[id_col]])

# ------------------------------------------------------------------------------
# Collapse polar_day -> day (two levels: day/night)
# ------------------------------------------------------------------------------
segments$light_start <- forcats::fct_collapse(
  segments$light_start,
  day   = c("day", "polar_day"),
  night = "night"
)

# ------------------------------------------------------------------------------
# Model specification + stepwise AIC to final additive model
# ------------------------------------------------------------------------------
full_form <- best_route ~ season + mean_wind_support + mean_crosswind +
  geomag_storm + light_start + great_circle_distance + duration + segment_type +
  contain_short_stop + previous_stopover_duration + next_stopover_duration

mod_full  <- glm(full_form, family = binomial, data = segments)
glm_final <- step(mod_full, trace = FALSE)

saveRDS(glm_final, rds_glm)

# ------------------------------------------------------------------------------
# Align NA handling and extract fitted components (exact rows used by final GLM)
# ------------------------------------------------------------------------------
mf_glm <- model.frame(glm_final)
p_hat  <- fitted(glm_final)
obs    <- as.numeric(mf_glm$best_route == "wind_optimal")

# Row indices in original data used by glm_final (for robust clustering / GLMM)
mf_rows <- as.integer(attr(mf_glm, "row.names"))
segments_used <- segments[mf_rows, , drop = FALSE]

# ------------------------------------------------------------------------------
# Diagnostics: pseudo-R2, AUC, HL test, VIF
# ------------------------------------------------------------------------------
mf_r2   <- 1 - glm_final$deviance / glm_final$null.deviance
tjur_r2 <- mean(p_hat[obs == 1]) - mean(p_hat[obs == 0])

roc_obj <- pROC::roc(response = obs, predictor = p_hat, quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))

hl <- ResourceSelection::hoslem.test(x = obs, y = p_hat, g = 10)

vif_tbl <- tryCatch({
  v <- car::vif(glm_final)
  if (is.matrix(v)) {
    as.data.frame(v) |>
      tibble::rownames_to_column("term") |>
      dplyr::rename(VIF = dplyr::last_col())
  } else {
    tibble::tibble(term = names(v), VIF = as.numeric(v))
  }
}, error = function(e) tibble::tibble())

# ------------------------------------------------------------------------------
# GLM with goose-ID cluster-robust SEs (main inference)
# ------------------------------------------------------------------------------
cluster_used <- segments_used[[id_col]]

# HC0 retained to reproduce revision checks exactly.
V_cl <- sandwich::vcovCL(glm_final, cluster = cluster_used, type = "HC0")
ct   <- lmtest::coeftest(glm_final, vcov. = V_cl)

tidy_rob <- tibble::tibble(
  term      = rownames(ct),
  estimate  = ct[, 1],
  std.error = ct[, 2],
  statistic = ct[, 3],
  p.value   = ct[, 4]
) |>
  dplyr::mutate(
    model = "GLM_cluster_robust",
    OR = exp(estimate)
  ) |>
  dplyr::select(model, term, estimate, std.error, statistic, p.value, OR)

readr::write_csv(tidy_rob, tab_rob)

# ------------------------------------------------------------------------------
# GLMM sensitivity check: same fixed effects as final GLM + goose-ID random intercept
# ------------------------------------------------------------------------------
rand_term <- as.formula(paste0(". ~ . + (1|`", id_col, "`)"))
glmm_form <- update(formula(glm_final), rand_term)

glmm_final <- lme4::glmer(
  glmm_form,
  family = binomial,
  data = segments_used,
  control = lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

saveRDS(glmm_final, rds_glmm)

glmm_sum <- summary(glmm_final)
fe <- glmm_sum$coefficients

tidy_glmm <- tibble::tibble(
  term      = rownames(fe),
  estimate  = fe[, "Estimate"],
  std.error = fe[, "Std. Error"],
  statistic = fe[, "z value"],
  p.value   = fe[, "Pr(>|z|)"]
) |>
  dplyr::mutate(
    model = "GLMM_random_intercept",
    OR = exp(estimate)
  ) |>
  dplyr::select(model, term, estimate, std.error, statistic, p.value, OR)

readr::write_csv(tidy_glmm, tab_glmm)

# Random intercept diagnostics
re_var_tbl <- as.data.frame(lme4::VarCorr(glmm_final))
is_sing <- lme4::isSingular(glmm_final, tol = 1e-5)

# ------------------------------------------------------------------------------
# Write one consolidated text report (robust-centred)
# ------------------------------------------------------------------------------
if (file.exists(logf)) file.remove(logf)
sink(logf)

cat("Binary logistic GLM (DTW reference): cluster-robust inference and GLMM sensitivity check\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("Input file:\n")
cat("  ", seg_fn, "\n\n", sep = "")

cat("Data rows:\n")
cat("  Rows in input data:", nrow(segments), "\n")
cat("  Rows used by final GLM (after NA handling):", nrow(mf_glm), "\n")
cat("  Clustering variable:", id_col, "\n\n")

cat("Final GLM formula (stepwise AIC):\n")
print(formula(glm_final)); cat("\n")

cat("GLMM sensitivity formula (same fixed effects + random intercept):\n")
print(formula(glmm_final)); cat("\n")

cat("GLM diagnostics (rows used by final model):\n")
cat("  Residual deviance:", round(glm_final$deviance, 1), "on", df.residual(glm_final), "df\n")
cat("  AIC:", round(AIC(glm_final), 1), "\n")
cat("  McFadden R2:", round(mf_r2, 3), "\n")
cat("  Tjur R2:", round(tjur_r2, 3), "\n")
cat("  AUC:", round(auc_val, 3), "\n")
cat(sprintf("  Hosmer–Lemeshow: HL X2=%.3f, df=%d, p=%.3f\n\n",
            hl$statistic, hl$parameter, hl$p.value))

if (nrow(vif_tbl)) {
  cat("VIF (final GLM):\n")
  print(vif_tbl, row.names = FALSE, digits = 4); cat("\n")
}

cat("Coefficient summary (GLM cluster-robust SEs; main inference, log-odds):\n")
print(tidy_rob, row.names = FALSE, digits = 4); cat("\n")

cat("GLMM random-intercept sensitivity diagnostics:\n")
cat("  isSingular:", is_sing, "\n")
cat("  Random-effects variance components:\n")
print(re_var_tbl, row.names = FALSE, digits = 6); cat("\n")

cat("GLMM fixed effects (Wald z tests; sensitivity check):\n")
print(tidy_glmm, row.names = FALSE, digits = 4); cat("\n")

cat("Stepwise AIC path (final GLM):\n")
print(glm_final$anova); cat("\n")

sink()

message("Wrote outputs to: ", normalizePath(outdir))