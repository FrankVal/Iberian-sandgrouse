# ============================================================
# Script: 01_screening-boruta_csv.R
# Purpose: Predictor screening for a sandgrouse-style occurrence framework
# Method : Boruta feature selection on a presence/absence response
#
# Input  : CSV table with response column 'pa' and predictor columns
#          Place the CSV in: data/processed/
#
# Exclude: date, Species, x, y columns (case-insensitive; not used as predictors)
#
# Output : output/screening/boruta_importance.csv
#          output/screening/boruta_selected_predictors.txt
#          output/screening/sessionInfo.txt
#          output/screening/correlation_matrix.csv
#          figs/screening_boruta_importance.tiff
#          figs/screening_correlation_matrix.tiff
#
# Run    : from project root (folder containing /data, /scripts, /R)
# ============================================================

# ---- Packages ----
library(Boruta)
library(ggplot2)
library(ggcorrplot)

# ---- Project-root check (portable relative paths) ----
if (!dir.exists("data") || !dir.exists("scripts")) {
  stop("Run this script from the project root (the folder containing /data and /scripts).")
}

# ---- Paths (relative) ----
input_csv_path <- file.path("data", "processed", "presence_absence.csv")  # edit if needed
screening_out_dir <- file.path("output", "screening")
figures_dir <- "figs"

dir.create(screening_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
if (!file.exists(input_csv_path)) {
  stop(
    "Input CSV not found: ", input_csv_path,
    "\nPut your file in data/processed/ and name it 'presence_absence.csv',",
    "\nor edit 'input_csv_path' to match your filename."
  )
}

occurrence_table <- read.csv(input_csv_path, stringsAsFactors = FALSE)

# ---- Validate response column ----
if (!("pa" %in% names(occurrence_table))) {
  stop("Response column 'pa' not found. Ensure the CSV contains a 'pa' column (0/1 or presence/absence).")
}

# ---- Drop non-predictor columns (case-insensitive) ----
excluded_columns <- c("date", "species", "x", "y")
columns_to_exclude <- names(occurrence_table)[tolower(names(occurrence_table)) %in% excluded_columns]

model_table <- occurrence_table[, setdiff(names(occurrence_table), columns_to_exclude), drop = FALSE]

# ---- Prepare response (classification) ----
model_table$pa <- as.factor(model_table$pa)

# ---- Handle missing values (screening step) ----
if (anyNA(model_table)) {
  warning("NAs detected: dropping rows with any NA for Boruta screening.")
  model_table <- na.omit(model_table)
}

# ---- Correlation plot (numeric predictors only) ----
numeric_predictors <- model_table[, setdiff(names(model_table), "pa"), drop = FALSE]
numeric_predictors <- numeric_predictors[, vapply(numeric_predictors, is.numeric, logical(1)), drop = FALSE]

# Remove zero-variance predictors (cor() can fail or create NA-heavy matrices)
if (ncol(numeric_predictors) >= 2) {
  non_constant <- vapply(numeric_predictors, function(x) stats::sd(x, na.rm = TRUE) > 0, logical(1))
  numeric_predictors <- numeric_predictors[, non_constant, drop = FALSE]
}

if (ncol(numeric_predictors) >= 2) {
  correlation_matrix <- stats::cor(numeric_predictors, use = "pairwise.complete.obs")
  
  correlation_plot <- ggcorrplot::ggcorrplot(
    correlation_matrix,
    type = "lower",
    lab = FALSE,
    show.diag = TRUE,
    hc.order = TRUE
  ) +
    theme_bw() +
    labs(x = NULL, y = NULL)
  
  ggsave(
    filename = file.path(figures_dir, "screening_correlation_matrix.tiff"),
    plot     = correlation_plot,
    width    = 30, height = 30, units = "cm", dpi = 600
  )
  
  write.csv(
    correlation_matrix,
    file = file.path(screening_out_dir, "correlation_matrix.csv"),
    row.names = TRUE
  )
} else {
  warning("Not enough numeric, non-constant predictors to compute a correlation matrix.")
}

# ---- Choose mtry from number of predictors ----
predictor_names <- setdiff(names(model_table), "pa")
if (length(predictor_names) < 1) stop("No predictors remain after excluding columns. Check the CSV headers.")

mtry_value <- max(1, floor(sqrt(length(predictor_names))))

# ---- Run Boruta ----
set.seed(123)
boruta_model <- Boruta(
  pa ~ .,
  data    = model_table,
  maxRuns = 2000,
  doTrace = 1,
  ntree   = 2000,
  mtry    = mtry_value
)

boruta_model_fixed <- TentativeRoughFix(boruta_model)

# ---- Selected predictors ----
selected_predictors <- getSelectedAttributes(boruta_model_fixed)
writeLines(selected_predictors, con = file.path(screening_out_dir, "boruta_selected_predictors.txt"))

# ---- Importance table ----
importance_stats <- attStats(boruta_model_fixed)
importance_stats$predictor <- rownames(importance_stats)

retained_importance <- importance_stats[
  importance_stats$decision != "Rejected" & importance_stats$predictor != "group",
  c("predictor", "meanImp", "decision")
]

retained_importance$importance_pct <-
  100 * retained_importance$meanImp / sum(retained_importance$meanImp, na.rm = TRUE)

retained_importance <- retained_importance[order(-retained_importance$importance_pct), ]

write.csv(
  retained_importance,
  file = file.path(screening_out_dir, "boruta_importance.csv"),
  row.names = FALSE
)

# ---- Importance plot ----
importance_plot <- ggplot(
  retained_importance,
  aes(x = reorder(predictor, importance_pct), y = importance_pct)
) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Importance (%)") +
  theme_bw()

ggsave(
  filename = file.path(figures_dir, "screening_boruta_importance.tiff"),
  plot     = importance_plot,
  width    = 30, height = 30, units = "cm", dpi = 600
)

# ---- Record session info (reproducibility) ----
writeLines(capture.output(sessionInfo()), con = file.path(screening_out_dir, "sessionInfo.txt"))

message("Done. Outputs written to:\n- ", screening_out_dir, "\n- ", figures_dir)
