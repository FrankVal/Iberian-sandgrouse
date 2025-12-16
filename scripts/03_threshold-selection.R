# ============================================================
# Script: 03_threshold-selection.R
# Purpose: Derive a probability cut-off for binary SDM maps using three
#          threshold criteria discussed in the paper:
#          (1) Sensitivity = Specificity
#          (2) Predicted prevalence = Observed prevalence
#          (3) Mean predicted probability
#
# Inputs :
#   - output/spatial_sdm/models/model_spatial_repeat_<label>.rds  (spatialRF rf_repeat)
#   - output/spatial_sdm/objects/model_table_<label>.rds          (contains 'pa')
# Optional:
#   - output/spatial_sdm/objects/coordinates_<label>.rds          (not needed here)
#
# Outputs:
#   - output/thresholds/thresholds_by_repetition_<label>.csv
#   - output/thresholds/threshold_summary_<label>.txt
#   - output/thresholds/sessionInfo_<label>.txt
#
# Notes:
#   - Uses PresenceAbsence::optimal.thresholds with methods:
#       "Sens=Spec", "PredPrev=Obs", "MeanProb"
#   - For each repetition, we take the median of the three method-specific
#     thresholds; then we take the median across repetitions as the final cut-off.
# ============================================================

# ---- Packages ----
library(PresenceAbsence)

# ---- Project-root check ----
if (!dir.exists("output") || !dir.exists("scripts")) {
  stop("Run this script from the project root (the folder containing /output and /scripts).")
}

# ---- User settings ----
dataset_label <- "sdm"  # keep species-agnostic; must match the label used in 02_spatial-SDMs.R

model_path <- file.path("output", "spatial_sdm", "models",
                        paste0("model_spatial_repeat_", dataset_label, ".rds"))
table_path <- file.path("output", "spatial_sdm", "objects",
                        paste0("model_table_", dataset_label, ".rds"))

out_dir <- file.path("output", "thresholds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load inputs ----
if (!file.exists(model_path)) stop("Model not found: ", model_path)
if (!file.exists(table_path)) stop("Model table not found: ", table_path)

spatial_model_repeated <- readRDS(model_path)
model_table <- readRDS(table_path)

if (!("pa" %in% names(model_table))) stop("Column 'pa' not found in model_table.")
observed_pa <- as.integer(as.character(model_table$pa))  # expects 0/1 or factor convertible to 0/1

if (!is.list(spatial_model_repeated$predictions)) {
  stop("Unexpected model object structure: model$predictions is not a list.")
}

# In spatialRF rf_repeat, predictions are typically stored per repetition; we handle flexibly.
predictions_obj <- spatial_model_repeated$predictions

# ---- Helper: extract per-repetition prediction vector ----
extract_prediction_vector <- function(predictions, repetition_index) {
  # Most common: a data.frame with columns pred.values.per.repetition.repetition_1, ...
  if (is.data.frame(predictions)) {
    col_name <- paste0("pred.values.per.repetition.repetition_", repetition_index)
    if (!(col_name %in% names(predictions))) {
      stop("Prediction column not found: ", col_name)
    }
    return(predictions[[col_name]])
  }
  
  # Alternative: list of numeric vectors
  if (is.list(predictions) && length(predictions) >= repetition_index) {
    vec <- predictions[[repetition_index]]
    if (!is.numeric(vec)) stop("Predictions for repetition ", repetition_index, " are not numeric.")
    return(vec)
  }
  
  stop("Could not extract predictions for repetition ", repetition_index, ".")
}

# ---- Helper: compute a single repetition threshold ----
compute_repetition_threshold <- function(observed, predicted) {
  # PresenceAbsence expects: Id, Observed (0/1), Predicted (probability)
  threshold_input <- data.frame(
    Id = seq_along(observed),
    Observed = observed,
    Predicted = predicted
  )
  
  threshold_table <- PresenceAbsence::optimal.thresholds(
    threshold_input,
    opt.method = c("Sens=Spec", "PredPrev=Obs", "MeanProb")
  )
  
  # Some versions return a data.frame with method rows and a "predicted" / "Predicted" column
  pred_col <- intersect(c("Predicted", "predicted"), names(threshold_table))
  if (length(pred_col) != 1) {
    stop("Could not find predicted threshold column in optimal.thresholds output.")
  }
  
  method_thresholds <- threshold_table[[pred_col]]
  
  # One threshold per method; take median across methods for this repetition
  stats::median(method_thresholds, na.rm = TRUE)
}

# ---- Determine number of repetitions ----
# Prefer an explicit slot if present; otherwise infer from columns/list length.
n_reps <- NA_integer_

if (is.data.frame(predictions_obj)) {
  rep_cols <- grep("^pred\\.values\\.per\\.repetition\\.repetition_\\d+$", names(predictions_obj), value = TRUE)
  n_reps <- length(rep_cols)
} else if (is.list(predictions_obj)) {
  n_reps <- length(predictions_obj)
}

if (is.na(n_reps) || n_reps < 1) stop("Could not infer number of repetitions from model predictions.")

# ---- Compute thresholds across repetitions ----
set.seed(123)  # not strictly necessary, but keeps any downstream randomness stable

threshold_by_repetition <- data.frame(
  repetition = seq_len(n_reps),
  threshold = NA_real_
)

for (r in seq_len(n_reps)) {
  predicted_prob <- extract_prediction_vector(predictions_obj, r)
  
  if (length(predicted_prob) != length(observed_pa)) {
    stop("Length mismatch at repetition ", r, ": predicted=", length(predicted_prob),
         " observed=", length(observed_pa))
  }
  
  threshold_by_repetition$threshold[r] <- compute_repetition_threshold(observed_pa, predicted_prob)
}

# Final cut-off = median across repetitions
final_cutoff <- stats::median(threshold_by_repetition$threshold, na.rm = TRUE)

# ---- Save outputs ----
out_csv <- file.path(out_dir, paste0("thresholds_by_repetition_", dataset_label, ".csv"))
write.csv(threshold_by_repetition, out_csv, row.names = FALSE)

out_txt <- file.path(out_dir, paste0("threshold_summary_", dataset_label, ".txt"))
summary_lines <- c(
  paste0("Dataset label: ", dataset_label),
  paste0("Repetitions: ", n_reps),
  paste0("Threshold methods: Sens=Spec; PredPrev=Obs; MeanProb"),
  paste0("Per-repetition threshold = median across the 3 methods"),
  paste0("Final cut-off = median across repetitions"),
  "",
  paste0("Final cut-off (median across repetitions): ", signif(final_cutoff, 6)),
  paste0("Per-repetition threshold summary:"),
  paste0("  mean  = ", signif(mean(threshold_by_repetition$threshold, na.rm = TRUE), 6)),
  paste0("  sd    = ", signif(stats::sd(threshold_by_repetition$threshold, na.rm = TRUE), 6)),
  paste0("  min   = ", signif(min(threshold_by_repetition$threshold, na.rm = TRUE), 6)),
  paste0("  max   = ", signif(max(threshold_by_repetition$threshold, na.rm = TRUE), 6))
)
writeLines(summary_lines, out_txt)

writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, paste0("sessionInfo_", dataset_label, ".txt")))

message("Done. Final cut-off = ", signif(final_cutoff, 6),
        "\nSaved:\n- ", out_csv, "\n- ", out_txt)
