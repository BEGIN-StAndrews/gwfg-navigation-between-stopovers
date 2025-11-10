# ==============================================================================
# Title: Multinomial GLM (compass routes, DTW reference) + diagnostics report
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read segments_all_meta_compass_dtw.csv; set factor baselines & collapse light_start
#   2) Fit multinomial model (nnet::multinom) and stepwise AIC (MASS::stepAIC)
#   3) Compute metrics: McFadden R2, multiclass AUC, accuracy, per-class recall
#   4) Produce VIF (dummy-response), confusion matrix, and OR with 95% Wald CI
#   5) Write one text report to MNL_outputs/MNL_compass_report.txt
# Inputs:
#   - segments_all_meta_compass_dtw.csv
# Outputs:
#   - MNL_outputs/MNL_compass_report.txt
# Usage:
#   Rscript R/15_2_glm_multinomial_compass.R
# ==============================================================================

rm(list = ls()); options(width = 120)
set.seed(123)

# ------------------------------------------------------------------------------
# Load libraries (auto-install if missing)
# ------------------------------------------------------------------------------
req <- c("dplyr","readr","forcats","car","nnet","MASS","pROC","tibble","tidyr","stringr")
inst <- setdiff(req, rownames(installed.packages()))
if (length(inst)) install.packages(inst, dependencies = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# Paths & files
# ------------------------------------------------------------------------------
seg_fn <- "segments_all_meta_compass_dtw.csv"
outdir <- "MNL_outputs"; dir.create(outdir, showWarnings = FALSE)
logf   <- file.path(outdir, "MNL_compass_report.txt")

# ------------------------------------------------------------------------------
# Load data and factor baselines
# ------------------------------------------------------------------------------
segments <- readr::read_csv(seg_fn, show_col_types = FALSE) |>
  dplyr::mutate(
    best_route = forcats::fct_relevel(
      factor(best_route),
      "geographic","geomagnetic","magnetoclinic","sun","local_wind"
    ),
    season             = forcats::fct_relevel(factor(season), "spring","autumn"),
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
# Predictors & model formula
# ------------------------------------------------------------------------------
num_cols  <- c("mean_wind_support","mean_crosswind","great_circle_distance",
               "duration","previous_stopover_duration","next_stopover_duration")
rhs_terms <- c("season", num_cols, "geomag_storm","light_start","segment_type","contain_short_stop")
full_form <- as.formula(paste("best_route ~", paste(rhs_terms, collapse = " + ")))

# ------------------------------------------------------------------------------
# NA handling transparency
# ------------------------------------------------------------------------------
mf   <- model.frame(full_form, data = segments, na.action = na.omit)
n_in <- nrow(segments); n_eff <- nrow(mf)

# ------------------------------------------------------------------------------
# VIF (dummy response, outcome-agnostic)
# ------------------------------------------------------------------------------
vif_rhs <- paste(rhs_terms, collapse = " + ")
vif_df  <- segments; vif_df$.__y <- rnorm(nrow(segments))
vif_lm  <- lm(as.formula(paste(".__y ~", vif_rhs)), data = vif_df)
v       <- car::vif(vif_lm)

vif_tbl <- if (is.matrix(v)) {
  tibble::tibble(term = rownames(v),
                 GVIF = v[, "GVIF"],
                 Df   = v[, "Df"],
                 GVIF_adj = v[, "GVIF"]^(1/(2*v[, "Df"])))
} else {
  tibble::tibble(term = names(v), GVIF = as.numeric(v), Df = 1, GVIF_adj = as.numeric(v))
}

# ------------------------------------------------------------------------------
# Fit model + stepwise AIC
# ------------------------------------------------------------------------------
mod_full <- nnet::multinom(full_form, data = segments, trace = FALSE)
model    <- MASS::stepAIC(mod_full, trace = FALSE)

# ------------------------------------------------------------------------------
# Metrics (use model’s own frame to align truth/probabilities)
# ------------------------------------------------------------------------------
# McFadden pseudo-R2
ll_null <- as.numeric(logLik(update(model, . ~ 1)))
ll_mod  <- as.numeric(logLik(model))
mf_r2   <- 1 - ll_mod/ll_null

mf_mod     <- model.frame(model)
truth      <- mf_mod$best_route
pred_class <- predict(model, newdata = mf_mod, type = "class")
prob       <- as.data.frame(predict(model, newdata = mf_mod, type = "probs"))
prob       <- prob[, levels(truth), drop = FALSE]

cm   <- table(Truth = truth, Pred = pred_class)
acc  <- sum(diag(cm)) / sum(cm)
rec  <- diag(prop.table(cm, 1)); rec[!is.finite(rec)] <- NA
rec_df <- tibble::tibble(class = names(rec), recall = as.numeric(rec))

mroc    <- pROC::multiclass.roc(truth, as.matrix(prob))
auc_val <- as.numeric(pROC::auc(mroc))

# ------------------------------------------------------------------------------
# Odds ratios, Wald CI, p-values (keyed join; no misalignment)
# ------------------------------------------------------------------------------
s      <- summary(model)
co_mat <- s$coefficients
se_mat <- s$standard.errors
stopifnot(!is.null(dim(co_mat)), !is.null(dim(se_mat)), identical(dim(co_mat), dim(se_mat)))

coef_long <- as.data.frame(co_mat) |>
  tibble::rownames_to_column("y.level") |>
  tidyr::pivot_longer(-y.level, names_to = "term", values_to = "coef")

se_long <- as.data.frame(se_mat) |>
  tibble::rownames_to_column("y.level") |>
  tidyr::pivot_longer(-y.level, names_to = "term", values_to = "se")

wald <- dplyr::left_join(coef_long, se_long, by = c("y.level","term")) |>
  dplyr::mutate(
    z       = coef / se,
    p.value = 2 * (1 - pnorm(abs(z))),
    OR      = exp(coef),
    CI_low  = exp(coef - 1.96 * se),
    CI_high = exp(coef + 1.96 * se)
  ) |>
  dplyr::select(y.level, term, OR, CI_low, CI_high, p.value)

p_eps <- 1e-4
or_print <- wald |>
  dplyr::mutate(
    OR      = round(OR, 3),
    CI_low  = round(CI_low, 3),
    CI_high = round(CI_high, 3),
    p_str   = ifelse(is.na(p.value), "NA",
                     format.pval(p.value, eps = p_eps, digits = 3))
  ) |>
  dplyr::select(y.level, term, OR, CI_low, CI_high, p_str) |>
  dplyr::rename(`p.value` = p_str)

# ------------------------------------------------------------------------------
# Write ONE neat text report
# ------------------------------------------------------------------------------
if (file.exists(logf)) file.remove(logf)
con <- file(logf, open = "wt")
wl <- function(...) writeLines(paste0(...), con)

wl("Final multinomial (compass) model — ", format(Sys.time(), "%Y-%m-%d %H:%M"))
wl("Baseline class (reference): geographic")
wl(sprintf("Effective sample size: %d", n_eff))
wl("")
wl("Formula:"); wl(deparse(formula(model))); wl("")
wl("Factor codings (reference levels in parentheses):")
wl("  season: spring (ref), autumn")
wl("  light_start: day (ref), night")
lvl_seg <- paste(levels(segments$segment_type), collapse = " | ")
wl(sprintf("  segment_type: %s  (ref = %s)", lvl_seg, levels(segments$segment_type)[1]))
wl("  contain_short_stop: FALSE (ref), TRUE")
lvl_gs <- paste(levels(segments$geomag_storm), collapse = " | ")
wl(sprintf("  geomag_storm: %s  (ref = %s)", lvl_gs, levels(segments$geomag_storm)[1]))
wl("")
wl(sprintf("Rows in: %d | Rows used after NA omission: %d", n_in, n_eff))
wl(strrep("-", 78))
wl(sprintf("Residual deviance: %.1f", model$deviance))
wl(sprintf("AIC: %.1f", AIC(model)))
wl(sprintf("McFadden R2: %.3f", mf_r2))
wl(sprintf("Multiclass AUC (Hand–Till): %.3f", auc_val))
wl(sprintf("Overall accuracy (micro-recall): %.1f%%", 100*acc))
wl(sprintf("Macro-recall: %.3f", mean(rec, na.rm = TRUE)))
wl("")
wl("Per-class recall (row-wise):")
capture.output(print(rec_df, row.names = FALSE, digits = 3), file = con)
wl(strrep("-", 78))
wl("Confusion matrix (rows = truth, cols = pred):")
capture.output(print(cm), file = con)
wl(strrep("-", 78))
wl("VIF (dummy-response, outcome-agnostic):")
capture.output(print(vif_tbl, row.names = FALSE, digits = 3), file = con)
wl(strrep("-", 78))
if (!is.null(model$anova)) {
  wl("Stepwise AIC path:")
  capture.output(print(model$anova), file = con)
  wl(strrep("-", 78))
}
wl("Odds ratios (vs geographic), 95% Wald CI, p-values (multinomial):")
capture.output(print(or_print, row.names = FALSE, max = 1e6), file = con)
close(con)

message("Done. Report at: ", normalizePath(logf))
