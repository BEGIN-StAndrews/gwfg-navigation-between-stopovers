# ==============================================================================
# Title: Multinomial GLM (compass routes, DTW reference): cluster-robust inference
# Manuscript: Navigation Between Stopovers by Greater White-Fronted Geese:
#             Comparing Compass Mechanisms and Efficiency Benchmarks
# Project: Segment-level navigation in migratory geese (2014–2024)
# Author: Ali Moayedi  |  am636@st-andrews.ac.uk
# License: MIT (see LICENSE)
# Summary:
#   1) Read segments_all_meta_compass_dtw.csv; set factor baselines & collapse light_start
#   2) Fit multinomial model (nnet::multinom) and stepwise AIC (MASS::stepAIC)
#   3) Compute model fit metrics (McFadden R2, multiclass AUC, overall accuracy, macro-recall)
#   4) Produce VIF table (dummy-response, outcome-agnostic) and stepwise AIC path
#   5) Compute goose-ID clustered (sandwich) robust SEs manually from score contributions
#   6) Export class-wise OR tables with 95% Wald CI and clustered p-values
#   7) Write one text report + CSV outputs for Supplementary Material generation
# Inputs:
#   - segments_all_meta_compass_dtw.csv
# Outputs:
#   - MNL_outputs_cluster_robust/MNL_compass_cluster_robust_report.txt
#   - MNL_outputs_cluster_robust/MNL_compass_cluster_robust_OR_table.csv
#   - MNL_outputs_cluster_robust/MNL_compass_VIF.csv
#   - MNL_outputs_cluster_robust/MNL_compass_stepAIC_path.csv   (if available)
#   - MNL_outputs_cluster_robust/MNL_compass_metrics.csv
# Usage:
#   Rscript R/15_2_glm_multinomial_compass_cluster_robust.R
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
outdir <- "MNL_outputs_cluster_robust"
dir.create(outdir, showWarnings = FALSE)

logf        <- file.path(outdir, "MNL_compass_cluster_robust_report.txt")
tab_or      <- file.path(outdir, "MNL_compass_cluster_robust_OR_table.csv")
tab_vif     <- file.path(outdir, "MNL_compass_VIF.csv")
tab_step    <- file.path(outdir, "MNL_compass_stepAIC_path.csv")
tab_metrics <- file.path(outdir, "MNL_compass_metrics.csv")

# ------------------------------------------------------------------------------
# Clustering variable
# ------------------------------------------------------------------------------
id_col <- "individual.local.identifier"

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
    segment_type       = factor(segment_type)
  )

if (!id_col %in% names(segments)) {
  stop(
    "ID column not found: ", id_col, "\n",
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
  day   = c("day","polar_day"),
  night = "night"
)

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
mf_full <- model.frame(full_form, data = segments, na.action = na.omit)
n_in <- nrow(segments)
n_eff_full <- nrow(mf_full)

# ------------------------------------------------------------------------------
# VIF 
# ------------------------------------------------------------------------------
vif_rhs <- paste(rhs_terms, collapse = " + ")
vif_df  <- segments
vif_df$.__y <- rnorm(nrow(segments))
vif_lm  <- lm(as.formula(paste(".__y ~", vif_rhs)), data = vif_df)
v       <- car::vif(vif_lm)

vif_tbl <- if (is.matrix(v)) {
  tibble::tibble(
    term = rownames(v),
    GVIF = v[, "GVIF"],
    Df   = v[, "Df"],
    GVIF_adj = v[, "GVIF"]^(1 / (2 * v[, "Df"]))
  )
} else {
  tibble::tibble(term = names(v), GVIF = as.numeric(v), Df = 1, GVIF_adj = as.numeric(v))
}
readr::write_csv(vif_tbl, tab_vif)

# ------------------------------------------------------------------------------
# Fit multinomial model + stepwise AIC
# ------------------------------------------------------------------------------
mod_full <- nnet::multinom(full_form, data = segments, trace = FALSE)
model    <- MASS::stepAIC(mod_full, trace = FALSE)

# Stepwise AIC path (if available)
if (!is.null(model$anova)) {
  step_tbl <- as.data.frame(model$anova) |>
    tibble::rownames_to_column("step")
  readr::write_csv(step_tbl, tab_step)
}

# ------------------------------------------------------------------------------
# Final model frame 
# ------------------------------------------------------------------------------
mf_mod <- model.frame(model)
used_rows <- as.integer(rownames(mf_mod))
segments_used <- segments[used_rows, , drop = FALSE]
cluster_used  <- segments_used[[id_col]]

n_eff_final <- nrow(mf_mod)
n_clusters  <- dplyr::n_distinct(cluster_used)

# ------------------------------------------------------------------------------
# Metrics 
# ------------------------------------------------------------------------------
truth      <- mf_mod$best_route
pred_class <- predict(model, newdata = mf_mod, type = "class")
prob       <- as.data.frame(predict(model, newdata = mf_mod, type = "probs"))
prob       <- prob[, levels(truth), drop = FALSE]

acc <- mean(pred_class == truth)

# Per-class recall 
cm  <- table(Truth = truth, Pred = pred_class)
rec <- diag(prop.table(cm, 1))
rec[!is.finite(rec)] <- NA

# McFadden pseudo-R2
ll_null <- as.numeric(logLik(update(model, . ~ 1)))
ll_mod  <- as.numeric(logLik(model))
mf_r2   <- 1 - ll_mod / ll_null

# Multiclass AUC 
mroc    <- pROC::multiclass.roc(truth, as.matrix(prob))
auc_val <- as.numeric(pROC::auc(mroc))

metrics_tbl <- tibble::tibble(
  n_in = n_in,
  n_eff_full_candidate = n_eff_full,
  n_eff_final = n_eff_final,
  geese_eff = n_clusters,
  baseline = levels(truth)[1],
  residual_deviance = model$deviance,
  AIC = AIC(model),
  McFadden_R2 = mf_r2,
  multiclass_AUC = auc_val,
  accuracy = acc,
  macro_recall = mean(rec, na.rm = TRUE)
)
readr::write_csv(metrics_tbl, tab_metrics)

# ------------------------------------------------------------------------------
# Extract coefficients in vcov order
# ------------------------------------------------------------------------------
coef_mat <- coef(model)   # rows = non-baseline classes, cols = terms (incl intercept)
if (is.vector(coef_mat)) coef_mat <- matrix(coef_mat, nrow = 1)

y_levs     <- rownames(coef_mat)
coef_terms <- colnames(coef_mat)

beta_vec <- as.vector(t(coef_mat))
names(beta_vec) <- as.vector(t(outer(y_levs, coef_terms, paste, sep = ":")))

V_naive <- vcov(model)

# Align covariance to coefficient vector names if needed
if (!identical(names(beta_vec), colnames(V_naive))) {
  V_naive <- V_naive[names(beta_vec), names(beta_vec), drop = FALSE]
}

# ------------------------------------------------------------------------------
# Cluster-robust sandwich covariance 
# ------------------------------------------------------------------------------

X <- model.matrix(delete.response(terms(model)), data = mf_mod)  
p_all <- as.matrix(prob)                                         
levs_all <- levels(truth)
base <- levs_all[1]
nonbase <- levs_all[-1]

q <- length(beta_vec)
S <- matrix(0, nrow = nrow(X), ncol = q)
colnames(S) <- names(beta_vec)

for (k in nonbase) {
  y_ind <- as.numeric(truth == k)
  pk    <- p_all[, k]
  diff  <- y_ind - pk
  
  for (t in seq_len(ncol(X))) {
    term_name <- colnames(X)[t]
    col_name  <- paste0(k, ":", term_name)
    if (col_name %in% colnames(S)) {
      S[, col_name] <- X[, t] * diff
    }
  }
}

u_by_c <- rowsum(S, group = as.character(cluster_used), reorder = FALSE)

meat <- matrix(0, nrow = q, ncol = q)
for (g in seq_len(nrow(u_by_c))) {
  u <- as.numeric(u_by_c[g, ])
  meat <- meat + tcrossprod(u)
}

bread <- V_naive
V_rob <- bread %*% meat %*% bread

# Numerical guard: enforce symmetry
V_rob <- (V_rob + t(V_rob)) / 2

se_rob <- sqrt(pmax(diag(V_rob), 0))
names(se_rob) <- names(beta_vec)

se_summary <- tibble::tibble(
  n_parameters = length(se_rob),
  n_clusters_used = n_clusters,
  min_cluster_robust_SE = min(se_rob, na.rm = TRUE),
  max_cluster_robust_SE = max(se_rob, na.rm = TRUE)
)

# ------------------------------------------------------------------------------
# Build cluster-robust OR table (Wald z tests + 95% Wald CI)
# ------------------------------------------------------------------------------
term_split <- stringr::str_split_fixed(names(beta_vec), ":", 2)
y_level <- term_split[, 1]
term    <- term_split[, 2]

z_rob <- beta_vec / se_rob
p_rob <- 2 * (1 - pnorm(abs(z_rob)))

or_tbl <- tibble::tibble(
  y.level = y_level,
  term = term,
  beta = as.numeric(beta_vec),
  SE = as.numeric(se_rob),
  z = as.numeric(z_rob),
  p.value = as.numeric(p_rob),
  OR = exp(beta_vec),
  CI_low = exp(beta_vec - 1.96 * se_rob),
  CI_high = exp(beta_vec + 1.96 * se_rob)
) |>
  dplyr::arrange(y.level, term)

readr::write_csv(or_tbl, tab_or)

p_eps <- 1e-4
or_print <- or_tbl |>
  dplyr::mutate(
    OR      = round(OR, 3),
    CI_low  = round(CI_low, 3),
    CI_high = round(CI_high, 3),
    p.value = ifelse(is.na(p.value), "NA",
                     format.pval(p.value, eps = p_eps, digits = 3))
  ) |>
  dplyr::select(y.level, term, OR, CI_low, CI_high, p.value)

# ------------------------------------------------------------------------------
# Write ONE neat text report 
# ------------------------------------------------------------------------------
if (file.exists(logf)) file.remove(logf)
con <- file(logf, open = "wt")
wl <- function(...) writeLines(paste0(...), con)

wl("Multinomial GLM (compass routes, DTW reference): cluster-robust inference")
wl("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"))
wl("")
wl("Input file:")
wl("  ", seg_fn)
wl("")
wl("Data rows:")
wl("  Rows in input data: ", n_in)
wl("  Rows used after NA omission (full candidate formula): ", n_eff_full)
wl("  Rows used by final model (after stepwise selection / model frame): ", n_eff_final)
wl("  Clustering variable: ", id_col)
wl("  Distinct geese in final model rows: ", n_clusters)
wl("")
wl("Baseline class (reference): ", base)
wl("")
wl("Final model formula (stepwise AIC):")
wl(deparse(formula(model)))
wl("")
wl("Factor codings (reference levels in parentheses):")
wl("  season: spring (ref), autumn")
wl("  light_start: day (ref), night")
lvl_seg <- paste(levels(segments$segment_type), collapse = " | ")
wl(sprintf("  segment_type: %s  (ref = %s)", lvl_seg, levels(segments$segment_type)[1]))
wl("  contain_short_stop: FALSE (ref), TRUE")
lvl_gs <- paste(levels(segments$geomag_storm), collapse = " | ")
wl(sprintf("  geomag_storm: %s  (ref = %s)", lvl_gs, levels(segments$geomag_storm)[1]))
wl(strrep("-", 78))
wl("Model fit metrics:")
wl(sprintf("  Residual deviance: %.1f", model$deviance))
wl(sprintf("  AIC: %.1f", AIC(model)))
wl(sprintf("  McFadden R2: %.3f", mf_r2))
wl(sprintf("  Multiclass AUC (Hand-Till): %.3f", auc_val))
wl(sprintf("  Overall accuracy (micro-recall): %.1f%%", 100 * acc))
wl(sprintf("  Macro-recall: %.3f", mean(rec, na.rm = TRUE)))
wl(strrep("-", 78))
wl("VIF (dummy-response, outcome-agnostic):")
capture.output(print(vif_tbl, row.names = FALSE, digits = 3, n = Inf), file = con)
wl(strrep("-", 78))
if (!is.null(model$anova)) {
  wl("Stepwise AIC path:")
  capture.output(print(model$anova), file = con)
  wl(strrep("-", 78))
}
wl("Cluster-robust covariance method (multinomial):")
wl("  - manual sandwich covariance from per-observation score contributions")
wl("  - clustered by goose ID")
wl("  - Wald z tests and 95% Wald confidence intervals use cluster-robust SEs")
wl("")
wl("SE summary (cluster-robust):")
capture.output(print(se_summary, row.names = FALSE, digits = 6), file = con)
wl(strrep("-", 78))
wl("Odds ratios (vs geographic), 95% Wald CI, p-values (cluster-robust multinomial):")
capture.output(print(or_print, row.names = FALSE, n = Inf, max = 1e6), file = con)
wl(strrep("-", 78))
wl("CSV outputs:")
wl("  - ", normalizePath(tab_or))
wl("  - ", normalizePath(tab_vif))
if (file.exists(tab_step)) wl("  - ", normalizePath(tab_step))
wl("  - ", normalizePath(tab_metrics))

close(con)

message("Done. Report at: ", normalizePath(logf))