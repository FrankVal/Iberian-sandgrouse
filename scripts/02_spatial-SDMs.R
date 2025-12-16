# ============================================================
# Script: 02_spatial-SDMs.R
# Purpose: Spatial SDM workflow using spatialRF (Moran diagnostics, tuning,
#          spatial eigenvector filtering, repeated fitting, evaluation)
#
# Input  : CSV table with:
#          - response column: pa  (presence/absence; will be coerced to factor)
#          - coordinate columns: x, y  (required for distance matrix)
#          - predictor columns (numeric and/or factors)
#
# Place the CSV in: data/sdm/presence_absence.csv
#
# Columns excluded from predictors (case-insensitive):
#   date, species, x, y
#
# Documented predictor sets used previously (as comments only; framework is agnostic):
#   Columns retained for the P_alchata dataset were:
#     pa, Bio1, Bio2, GPP, LST, Shee_Ext_M, HV, SC, Pig_Ext_M
#
#   Columns retained for the P_orientalis dataset were:
#     pa, GPP, Bio1, Bio2, Bio12, Shee_Ext_M, MSAVI2, OCD, Farm_Ext_M, LST, BD, HV
#
# Output : output/spatial_sdm/...
#          figs/spatial_sdm/...
#
# Run    : from project root (folder containing /data and /scripts)
# ============================================================

# ---- Packages ----
library(spatialRF)
library(ggplot2)
library(viridis)

# ---- Project-root check ----
if (!dir.exists("data") || !dir.exists("scripts")) {
  stop("Run this script from the project root (the folder containing /data and /scripts).")
}

# ---- User settings ----
dataset_label <- "sdm"  # used only to label outputs (keep species-agnostic)

input_csv_path <- file.path("data", "sdm", "presence_absence.csv")  # edit if needed

distance_thresholds <- c(0, 1000, 3000, 7000, 10000)  # same units as x/y
random_seed <- 1234

# Hyperparameter grids for tuning
tuning_grid_mtry          <- c(2, 3, 5)
tuning_grid_min_node_size <- c(5, 10, 20, 40)
tuning_grid_num_trees     <- c(500, 1000, 2000)

# Final hyperparameters for spatial models (edit after inspecting tuning)
final_mtry          <- 5
final_min_node_size <- 5
final_num_trees     <- 2000

repetitions <- 30
training_fraction <- 0.75

n_cores <- max(1, parallel::detectCores() - 1)

# ---- Output folders ----
out_root   <- file.path("output", "spatial_sdm")
out_models <- file.path(out_root, "models")
out_objs   <- file.path(out_root, "objects")
out_tables <- file.path(out_root, "tables")
fig_dir    <- file.path("figs", "spatial_sdm")

dir.create(out_root,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_models, recursive = TRUE, showWarnings = FALSE)
dir.create(out_objs,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,    recursive = TRUE, showWarnings = FALSE)

# ---- Helper: consistent plot styling (portable; avoids system fonts) ----
apply_plot_style <- function(p) {
  p +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

# ---- Load data ----
if (!file.exists(input_csv_path)) stop("Input CSV not found: ", input_csv_path)
occurrence_table <- read.csv(input_csv_path, stringsAsFactors = FALSE)

# ---- Validate essential columns ----
if (!("pa" %in% names(occurrence_table))) stop("Response column 'pa' not found in the CSV.")

name_lower <- tolower(names(occurrence_table))
if (!all(c("x", "y") %in% name_lower)) {
  stop("Coordinate columns x and y are required (case-insensitive) to compute the distance matrix.")
}

x_col <- names(occurrence_table)[match("x", name_lower)]
y_col <- names(occurrence_table)[match("y", name_lower)]

# ---- Prepare response and coordinates ----
occurrence_table$pa <- as.factor(occurrence_table$pa)

coordinates <- data.frame(
  x = occurrence_table[[x_col]],
  y = occurrence_table[[y_col]]
)

# ---- Exclude non-predictor columns ----
excluded_from_predictors <- c("date", "species", "x", "y")
columns_to_exclude <- names(occurrence_table)[tolower(names(occurrence_table)) %in% excluded_from_predictors]

# Model table retains response + predictors (but not x/y/date/species)
model_table <- occurrence_table[, setdiff(names(occurrence_table), columns_to_exclude), drop = FALSE]

# ---- Drop incomplete rows BEFORE distance matrix (to preserve alignment) ----
# We need complete cases across response + predictors + coordinates
complete_rows <- stats::complete.cases(model_table) & stats::complete.cases(coordinates)

if (!all(complete_rows)) {
  warning("Dropping rows with missing values in response/predictors/coordinates before computing distances.")
  model_table <- model_table[complete_rows, , drop = FALSE]
  coordinates <- coordinates[complete_rows, , drop = FALSE]
}

# ---- Compute distance matrix (Euclidean; units follow x/y) ----
distance_matrix <- as.matrix(stats::dist(coordinates, diag = TRUE, upper = TRUE))
saveRDS(distance_matrix, file = file.path(out_objs, paste0("distance_matrix_", dataset_label, ".rds")))
saveRDS(model_table,     file = file.path(out_objs, paste0("model_table_", dataset_label, ".rds")))
saveRDS(coordinates,     file = file.path(out_objs, paste0("coordinates_", dataset_label, ".rds")))

# ---- Predictor names + mtry ----
dependent_var <- "pa"
predictor_names <- setdiff(names(model_table), dependent_var)
if (length(predictor_names) < 1) stop("No predictors available after exclusions.")

mtry_default <- max(1, floor(sqrt(length(predictor_names))))

# ---- Moran test on training data (exploratory) ----
moran_training_plot <- spatialRF::plot_training_df_moran(
  data = model_table,
  dependent.variable.name = dependent_var,
  predictor.variable.names = predictor_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  fill.color = viridis::viridis(100, direction = -1),
  point.color = "gray40"
)

moran_training_plot <- apply_plot_style(moran_training_plot) + labs(color = "Feature")

ggsave(
  filename = file.path(fig_dir, paste0("moran_training_", dataset_label, ".tiff")),
  plot     = moran_training_plot,
  width    = 22, height = 22, units = "cm", dpi = 600
)

# ---- Binary check (useful guardrail) ----
binary_check <- spatialRF::is_binary(data = model_table, dependent.variable.name = dependent_var)
writeLines(
  paste0("is_binary(pa) = ", binary_check),
  con = file.path(out_tables, paste0("binary_check_", dataset_label, ".txt"))
)

# ---- Fit non-spatial RF (baseline) ----
set.seed(random_seed)
non_spatial_model <- spatialRF::rf(
  data = model_table,
  xy = coordinates,
  dependent.variable.name = dependent_var,
  predictor.variable.names = predictor_names,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  seed = random_seed,
  verbose = FALSE,
  ranger.arguments = list(
    mtry = mtry_default,
    num.trees = final_num_trees
  )
)

saveRDS(non_spatial_model, file = file.path(out_models, paste0("model_non_spatial_", dataset_label, ".rds")))

# Residual diagnostics (baseline)
baseline_residuals_plot <- spatialRF::plot_residuals_diagnostics(non_spatial_model, verbose = FALSE)
baseline_residuals_plot <- apply_plot_style(baseline_residuals_plot)

ggsave(
  filename = file.path(fig_dir, paste0("residuals_non_spatial_", dataset_label, ".tiff")),
  plot     = baseline_residuals_plot,
  width    = 22, height = 22, units = "cm", dpi = 600
)

# ---- Hyperparameter tuning (optional but recommended) ----
set.seed(1)
tuning_results <- spatialRF::rf_tuning(
  model = non_spatial_model,
  mtry = tuning_grid_mtry,
  min.node.size = tuning_grid_min_node_size,
  xy = coordinates,
  num.trees = tuning_grid_num_trees,
  repetitions = repetitions,
  training.fraction = training_fraction,
  seed = 1,
  verbose = TRUE,
  n.cores = n_cores,
  cluster = NULL
)

saveRDS(tuning_results, file = file.path(out_objs, paste0("tuning_results_", dataset_label, ".rds")))

tuning_plot <- spatialRF::plot_tuning(
  tuning_results,
  point.color = viridis::viridis(100, option = "F"),
  verbose = TRUE
)
tuning_plot <- apply_plot_style(tuning_plot)

ggsave(
  filename = file.path(fig_dir, paste0("tuning_hyperparameters_", dataset_label, ".tiff")),
  plot     = tuning_plot,
  width    = 22, height = 22, units = "cm", dpi = 600
)

# ---- Fit spatial RF (MEM filtering; defaults for eigenvector selection) ----
set.seed(random_seed)
spatial_model <- spatialRF::rf_spatial(
  model = non_spatial_model,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  xy = coordinates,
  method = "mem.moran.sequential",
  ranger.arguments = list(
    mtry = final_mtry,
    num.trees = final_num_trees,
    min.node.size = final_min_node_size
  ),
  verbose = FALSE,
  seed = random_seed
)

saveRDS(spatial_model, file = file.path(out_models, paste0("model_spatial_", dataset_label, ".rds")))

# Residual diagnostics (spatial)
spatial_residuals_plot <- spatialRF::plot_residuals_diagnostics(spatial_model, verbose = FALSE)
spatial_residuals_plot <- apply_plot_style(spatial_residuals_plot)

ggsave(
  filename = file.path(fig_dir, paste0("residuals_spatial_", dataset_label, ".tiff")),
  plot     = spatial_residuals_plot,
  width    = 22, height = 22, units = "cm", dpi = 600
)

# ---- Repeat spatial model (quantify variation across resamples/trees) ----
set.seed(random_seed)
spatial_model_repeated <- spatialRF::rf_repeat(
  model = spatial_model,
  xy = coordinates,
  repetitions = repetitions,
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds,
  seed = random_seed,
  ranger.arguments = list(
    mtry = final_mtry,
    num.trees = final_num_trees,
    min.node.size = final_min_node_size
  ),
  verbose = TRUE,
  n.cores = n_cores,
  cluster = NULL
)

saveRDS(spatial_model_repeated, file = file.path(out_models, paste0("model_spatial_repeat_", dataset_label, ".rds")))

# Residual diagnostics (repeated)
repeat_residuals_plot <- spatialRF::plot_residuals_diagnostics(spatial_model_repeated, verbose = FALSE)
repeat_residuals_plot <- apply_plot_style(repeat_residuals_plot)

ggsave(
  filename = file.path(fig_dir, paste0("residuals_spatial_repeat_", dataset_label, ".tiff")),
  plot     = repeat_residuals_plot,
  width    = 22, height = 22, units = "cm", dpi = 600
)

# ---- Moran residuals plot ----
moran_residuals_plot <- spatialRF::plot_moran(
  spatial_model_repeated,
  point.color = viridis::viridis(100, option = "F", direction = -1),
  line.color = "gray30",
  option = 1,
  ncol = 1,
  verbose = TRUE
)
moran_residuals_plot <- apply_plot_style(moran_residuals_plot) + labs(color = "Feature")

ggsave(
  filename = file.path(fig_dir, paste0("moran_residuals_", dataset_label, ".tiff")),
  plot     = moran_residuals_plot,
  width    = 22, height = 22, units = "cm", dpi = 600
)

# ---- Permutation importance plot (error increase) ----
permutation_importance_plot <- spatialRF::plot_importance(
  spatial_model_repeated,
  verbose = FALSE
)
permutation_importance_plot <- apply_plot_style(permutation_importance_plot) + labs(color = "Feature")

ggsave(
  filename = file.path(fig_dir, paste0("permutation_importance_", dataset_label, ".tiff")),
  plot     = permutation_importance_plot,
  width    = 22, height = 22, units = "cm", dpi = 600
)

# ---- Cross-validated importance (transferability) ----
# (This returns an object with importance summaries; you can plot it later if needed.)
cv_importance <- spatialRF::rf_importance(
  model = spatial_model,
  xy = coordinates,
  repetitions = repetitions,
  training.fraction = training_fraction,
  metric = c("auc"),
  fill.color = viridis::viridis(100, option = "F", direction = -1, alpha = 1, end = 0.9),
  line.color = "white",
  seed = 1,
  verbose = TRUE,
  n.cores = n_cores,
  cluster = NULL
)

saveRDS(cv_importance, file = file.path(out_objs, paste0("cv_importance_", dataset_label, ".rds")))

# ---- Model evaluation (AUC across repetitions) ----
evaluation_results <- spatialRF::rf_evaluate(
  model = spatial_model_repeated,
  xy = coordinates,
  repetitions = repetitions,
  training.fraction = training_fraction,
  metrics = c("auc"),
  seed = random_seed,
  verbose = TRUE,
  n.cores = n_cores,
  cluster = NULL
)

saveRDS(evaluation_results, file = file.path(out_objs, paste0("evaluation_", dataset_label, ".rds")))

# Plot evaluation
evaluation_plot <- spatialRF::plot_evaluation(
  evaluation_results,
  fill.color = viridis::viridis(3, option = "F", alpha = 0.8, direction = -1),
  line.color = "gray30",
  verbose = TRUE,
  notch = TRUE
)
evaluation_plot <- apply_plot_style(evaluation_plot)

ggsave(
  filename = file.path(fig_dir, paste0("accuracy_auc_", dataset_label, ".tiff")),
  plot     = evaluation_plot,
  width    = 22, height = 25, units = "cm", dpi = 600
)

# ---- MEM thresholds object (optional diagnostic) ----
mem_thresholds <- spatialRF::mem_multithreshold(
  distance.matrix = distance_matrix,
  distance.thresholds = distance_thresholds
)
saveRDS(mem_thresholds, file = file.path(out_objs, paste0("mem_thresholds_", dataset_label, ".rds")))

# ---- Record session info (reproducibility) ----
writeLines(capture.output(sessionInfo()), con = file.path(out_root, paste0("sessionInfo_", dataset_label, ".txt")))

message("Done. Outputs written under:\n- ", out_root, "\n- ", fig_dir)
