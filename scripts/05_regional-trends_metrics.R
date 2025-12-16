# ============================================================
# Script: 05_regional-trends_metrics.R
# Purpose: Regional (NUTS-2) changes in multiple habitat metrics derived from
#          Spain-wide binary habitat maps (https://doi.org/10.6084/m9.figshare.30898223), comparing:
#            - Initial period: 2005–2007
#            - Ending period : 2019–2022
#
# Metrics (FRAGSTATS-style; class = habitat, value == 1):
#   - CA      (km²)   Total class area (habitat amount)
#   - NP      (count) Number of habitat patches (8-neighbour connectivity)
#   - AREA_MN (ha)    Mean patch area
#   - LPI     (%)     Largest patch index
#   - ED      (m/ha)  Edge density
#
# Input rasters (GeoTIFF; values 0/1, NA outside Spain), downloaded from Figshare:
#   Pterocles_alchata_<YEAR>_Mean_Thresholds.tif     (2005–2022)
#   Pterocles_orientalis_<YEAR>_Mean_Thresholds.tif  (2005–2022)
#
# Expected repository structure (BES/reproducible style):
#   data/binary_maps/      (input GeoTIFFs)
#   output/regional_trends/ (created by this script)
#
# Output (single table for all metrics):
#   output/regional_trends/metrics_NUTS2_summary_2005-2007_vs_2019-2022.csv
#   output/regional_trends/metrics_NUTS2_summary_2005-2007_vs_2019-2022.rds
#   (optional) output/regional_trends/metrics_NUTS2_summary_2005-2007_vs_2019-2022.xlsx
#   output/regional_trends/sessionInfo_regional_trends.txt
#
# Run: from project root (folder containing /data and /scripts)
# ============================================================

# ---- Packages ----
library(terra)
library(dplyr)
library(sf)
library(giscoR)
library(tidyr)

# ---- Project-root check ----
if (!dir.exists("data") || !dir.exists("scripts")) {
  stop("Run this script from the project root (the folder containing /data and /scripts).")
}

# ------------------------------------------------------------
# 0) Settings
# ------------------------------------------------------------
binary_raster_dir <- file.path("data", "binary_maps")
out_dir <- file.path("output", "regional_trends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raster_sets <- c("Pterocles_alchata", "Pterocles_orientalis")

periods <- list(
  list(label = "2005–2007", years = 2005:2007),
  list(label = "2019–2022", years = 2019:2022)
)

patch_connectivity <- 8  # 8-neighbour (queen); keep consistent across scripts

metric_units <- c(
  CA = "km²",
  NP = "count",
  AREA_MN = "ha",
  LPI = "%",
  ED = "m/ha"
)

# ------------------------------------------------------------
# 1) Helpers
# ------------------------------------------------------------
raster_path_for_year <- function(raster_set, year, raster_dir) {
  file.path(raster_dir, sprintf("%s_%d_Mean_Thresholds.tif", raster_set, year))
}

assert_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) stop("Missing raster(s):\n", paste(missing, collapse = "\n"))
}

assert_projected_metre_crs <- function(r) {
  if (is.lonlat(r)) {
    stop(
      "Rasters appear to be lon/lat. ED (perimeter) requires a projected CRS in metres.\n",
      "Please reproject rasters (or compute ED in a projected CRS) before running."
    )
  }
  invisible(TRUE)
}

assert_binary_values <- function(r) {
  v <- unique(values(r))
  v <- v[!is.na(v)]
  if (length(v) == 0) return(invisible(TRUE))
  if (!all(v %in% c(0, 1))) {
    stop(
      "Raster is not strictly binary (0/1). Example unique non-NA values: ",
      paste(utils::head(sort(v), 20), collapse = ", ")
    )
  }
  invisible(TRUE)
}

# Build a period-composite binary raster (0/1/NA) using:
#   median across years -> presence if median >= 0.5
make_period_binary <- function(raster_set, years, raster_dir) {
  paths <- vapply(years, function(y) raster_path_for_year(raster_set, y, raster_dir), character(1))
  assert_files_exist(paths)
  
  r_stack <- rast(paths)
  assert_binary_values(r_stack[[1]])
  
  r_med <- app(r_stack, median, na.rm = TRUE)
  
  # Keep NA outside Spain; inside Spain threshold to 0/1
  ifel(is.na(r_med), NA, ifel(r_med >= 0.5, 1, 0))
}

# Compute FRAGSTATS-style metrics for a given landscape raster (0/1/NA)
# Optionally restrict to a polygon (SpatVector of length 1).
compute_metrics <- function(binary_raster, region_poly = NULL, directions = 8) {
  
  r <- binary_raster
  if (!is.null(region_poly)) {
    r <- mask(crop(r, region_poly), region_poly)
  }
  
  # If no raster coverage in region (all NA), return zeros
  if (isTRUE(all(is.na(values(r))))) {
    return(list(CA = 0, NP = 0, AREA_MN = 0, LPI = 0, ED = 0))
  }
  
  assert_projected_metre_crs(r)
  
  # Geometry (assumes metres)
  res_xy <- res(r)
  cell_area_m2 <- res_xy[1] * res_xy[2]
  cell_area_ha <- cell_area_m2 / 10000
  
  # Landscape area A (ha): all non-NA cells (class 0 + class 1)
  na_cells <- terra::global(is.na(r), "sum", na.rm = TRUE)[1, 1]
  landscape_cells <- terra::ncell(r) - na_cells
  landscape_area_ha <- landscape_cells * cell_area_ha
  
  # CA (km²): sum of cell areas where class == 1
  cell_area_km2 <- cell_area_m2 / 1e6
  class1_cells <- terra::global(r == 1, "sum", na.rm = TRUE)[1, 1]
  ca_km2 <- class1_cells * cell_area_km2
  
  # Habitat-only raster for patch metrics
  habitat <- r
  habitat[habitat != 1] <- NA
  
  # If no habitat present, patch metrics are zero
  if (isTRUE(all(is.na(values(habitat))))) {
    return(list(CA = ca_km2, NP = 0, AREA_MN = 0, LPI = 0, ED = 0))
  }
  
  # Patches (class 1)
  patch_id <- patches(habitat, directions = directions)
  patch_freq <- terra::freq(patch_id, useNA = "no")
  
  if (is.null(patch_freq) || nrow(patch_freq) == 0) {
    return(list(CA = ca_km2, NP = 0, AREA_MN = 0, LPI = 0, ED = 0))
  }
  
  patch_area_ha <- patch_freq$count * cell_area_ha
  np <- nrow(patch_freq)
  area_mn <- mean(patch_area_ha)
  lpi <- if (landscape_area_ha > 0) 100 * max(patch_area_ha) / landscape_area_ha else 0
  
  # ED (m/ha): total habitat edge length (m) / landscape area (ha)
  patch_polygons <- as.polygons(patch_id, dissolve = TRUE, na.rm = TRUE)
  perim_m <- expanse(patch_polygons, unit = "m", what = "perimeter")
  total_edge_m <- sum(perim_m, na.rm = TRUE)
  ed <- if (landscape_area_ha > 0) total_edge_m / landscape_area_ha else 0
  
  list(CA = ca_km2, NP = np, AREA_MN = area_mn, LPI = lpi, ED = ed)
}

metrics_to_long <- function(metrics_list, raster_set, region_name, nuts_id, period_label) {
  tibble(
    RasterSet = raster_set,
    Region    = region_name,
    NUTS_ID   = nuts_id,
    Period    = period_label,
    CA        = metrics_list$CA,
    NP        = metrics_list$NP,
    AREA_MN   = metrics_list$AREA_MN,
    LPI       = metrics_list$LPI,
    ED        = metrics_list$ED
  ) %>%
    pivot_longer(cols = c(CA, NP, AREA_MN, LPI, ED),
                 names_to = "Metric",
                 values_to = "Value") %>%
    mutate(Units = unname(metric_units[Metric]))
}

# ------------------------------------------------------------
# 2) NUTS-2 polygons (Spain)
# ------------------------------------------------------------
nuts2_sf <- gisco_get_nuts(
  year       = "2021",
  nuts_level = 2,
  country    = "ES",
  resolution = "10"
)

# ------------------------------------------------------------
# 3) Main: compute metrics for each raster set, period, and region
# ------------------------------------------------------------
all_long <- lapply(raster_sets, function(raster_set) {
  
  # Period rasters (0/1/NA)
  period_rasters <- lapply(periods, function(p) {
    r <- make_period_binary(raster_set, p$years, binary_raster_dir)
    assert_projected_metre_crs(r)
    list(label = p$label, raster = r)
  })
  
  # Align NUTS-2 to raster CRS
  target_crs <- crs(period_rasters[[1]]$raster)
  nuts2_proj <- st_transform(nuts2_sf, target_crs)
  nuts2_vect <- vect(nuts2_proj)
  
  nuts2_key <- st_drop_geometry(nuts2_proj) %>%
    mutate(.row_id = row_number()) %>%
    transmute(
      .row_id,
      NUTS_ID,
      Region = NUTS_NAME
    )
  
  # --- National ("All (Spain)") rows ---
  all_rows <- lapply(period_rasters, function(pr) {
    m <- compute_metrics(pr$raster, region_poly = NULL, directions = patch_connectivity)
    metrics_to_long(m, raster_set, "All (Spain)", NA_character_, pr$label)
  }) %>% bind_rows()
  
  # --- Per-region rows ---
  region_rows <- lapply(seq_len(nrow(nuts2_key)), function(i) {
    region_name <- nuts2_key$Region[i]
    nuts_id     <- nuts2_key$NUTS_ID[i]
    region_poly <- nuts2_vect[i]
    
    bind_rows(lapply(period_rasters, function(pr) {
      m <- compute_metrics(pr$raster, region_poly = region_poly, directions = patch_connectivity)
      metrics_to_long(m, raster_set, region_name, nuts_id, pr$label)
    }))
  }) %>% bind_rows()
  
  bind_rows(all_rows, region_rows)
}) %>% bind_rows()

# ------------------------------------------------------------
# 4) Summarise: Initial vs Ending + differences
# ------------------------------------------------------------
metrics_summary <- all_long %>%
  pivot_wider(names_from = Period, values_from = Value) %>%
  rename(
    Initial = `2005–2007`,
    Ending  = `2019–2022`
  ) %>%
  mutate(
    Change = Ending - Initial,
    `Change (%)` = if_else(Initial > 0, 100 * Change / Initial, NA_real_)
  ) %>%
  select(RasterSet, Region, NUTS_ID, Metric, Units, Initial, Ending, Change, `Change (%)`) %>%
  arrange(RasterSet, Metric, Region)

# Tidy rounding for output (keeps counts as integers)
metrics_summary <- metrics_summary %>%
  mutate(
    Initial = if_else(Metric == "NP", round(Initial, 0), round(Initial, 3)),
    Ending  = if_else(Metric == "NP", round(Ending,  0), round(Ending,  3)),
    Change  = if_else(Metric == "NP", round(Change,  0), round(Change,  3)),
    `Change (%)` = round(`Change (%)`, 2)
  )

# ------------------------------------------------------------
# 5) Save outputs
# ------------------------------------------------------------
csv_path <- file.path(out_dir, "metrics_NUTS2_summary_2005-2007_vs_2019-2022.csv")
rds_path <- file.path(out_dir, "metrics_NUTS2_summary_2005-2007_vs_2019-2022.rds")

write.csv(metrics_summary, csv_path, row.names = FALSE)
saveRDS(metrics_summary, rds_path)

if (requireNamespace("writexl", quietly = TRUE)) {
  xlsx_path <- file.path(out_dir, "metrics_NUTS2_summary_2005-2007_vs_2019-2022.xlsx")
  writexl::write_xlsx(metrics_summary, xlsx_path)
}

writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_regional_trends.txt"))

metrics_summary
