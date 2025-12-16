# ============================================================
# Script: 04_national-trends_fragstats-metrics.R
# Purpose: Compute FRAGSTATS-style class/landscape metrics from binary habitat maps.
#          (CA, NP, AREA_MN, LPI, ED) for 2005–2022.
#
# Input rasters (downloaded from Figshare: https://doi.org/10.6084/m9.figshare.30898223 and stored locally):
#   Pterocles_alchata_<YEAR>_Mean_Thresholds.tif     (2005–2022)
#   Pterocles_orientalis_<YEAR>_Mean_Thresholds.tif  (2005–2022)
#
# Expected properties:
#   - Binary values: 0/1
#   - NA outside Spain (landscape extent is the non-NA area)
#   - Projected CRS with metres (perimeter in m; ED in m/ha)
#
# Expected storage (BES/reproducible style):
#   data/binary_maps/   (GeoTIFFs)
#
# Outputs:
#   output/national_trends/fragstats_metrics.csv
#   output/national_trends/fragstats_metrics.rds
# ============================================================

library(terra)
library(dplyr)

# ---- Project-root check ----
if (!dir.exists("data") || !dir.exists("scripts")) {
  stop("Run this script from the project root (the folder containing /data and /scripts).")
}

# ----------------------------------------------------------
# 0) User settings
# ----------------------------------------------------------
years <- 2005:2022

# Where you placed the GeoTIFFs downloaded from Figshare
binary_raster_dir <- file.path("data", "binary_maps")

# Raster prefixes present in the Figshare archive (file naming convention)
raster_prefixes <- c("Pterocles_alchata", "Pterocles_orientalis")

# Patch connectivity (FRAGSTATS can be 4 or 8; keep 8 to match prior code)
patch_directions <- 8

# ----------------------------------------------------------
# 1) Helpers
# ----------------------------------------------------------
read_binary_map <- function(year, raster_prefix, raster_dir) {
  raster_path <- file.path(raster_dir, sprintf("%s_%d_Mean_Thresholds.tif", raster_prefix, year))
  if (!file.exists(raster_path)) stop("Raster not found: ", raster_path)
  rast(raster_path)
}

assert_binary_raster <- function(r) {
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

# ----------------------------------------------------------
# 2) Core metrics (FRAGSTATS-style) for one raster
#    CA (ha), NP, AREA_MN (ha), LPI (%), ED (m/ha)
# ----------------------------------------------------------
compute_fragstats_metrics <- function(year, raster_prefix, raster_dir, directions = 8) {
  message("Processing: ", raster_prefix, " | year: ", year)
  
  r <- read_binary_map(year, raster_prefix, raster_dir)
  assert_binary_raster(r)
  
  # Geometry
  res_xy <- res(r)               # assumed metres
  cell_area_m2 <- res_xy[1] * res_xy[2]
  cell_area_ha <- cell_area_m2 / 10000
  
  # Landscape area A (ha): all non-NA cells inside Spain (class 0 + class 1)
  na_cells <- terra::global(is.na(r), "sum", na.rm = TRUE)[1, 1]
  landscape_cells <- terra::ncell(r) - na_cells
  landscape_area_ha <- landscape_cells * cell_area_ha
  
  # CA (Total Class Area) for class 1 (ha)
  class1_cells <- terra::global(r == 1, "sum", na.rm = TRUE)[1, 1]
  class_area_ha <- class1_cells * cell_area_ha
  
  # Habitat-only raster for patch detection
  habitat <- r
  habitat[habitat != 1] <- NA
  
  # If no habitat present, return zeros (CA already computed)
  if (isTRUE(all(is.na(values(habitat))))) {
    return(data.frame(
      RasterSet = raster_prefix,
      Year      = year,
      CA        = class_area_ha,
      NP        = 0,
      AREA_MN   = 0,
      LPI       = 0,
      ED        = 0
    ))
  }
  
  # Patches (class 1), 8-neighbour by default
  patch_id_raster <- patches(habitat, directions = directions)
  
  # Patch sizes (counts) without pulling all cell values
  patch_freq <- terra::freq(patch_id_raster, useNA = "no")
  if (is.null(patch_freq) || nrow(patch_freq) == 0) {
    return(data.frame(
      RasterSet = raster_prefix,
      Year      = year,
      CA        = class_area_ha,
      NP        = 0,
      AREA_MN   = 0,
      LPI       = 0,
      ED        = 0
    ))
  }
  
  patch_cell_counts <- patch_freq$count
  patch_area_ha <- patch_cell_counts * cell_area_ha
  
  number_of_patches <- length(patch_area_ha)          # NP
  mean_patch_area_ha <- mean(patch_area_ha)           # AREA_MN
  largest_patch_ha <- max(patch_area_ha)              # for LPI
  
  # LPI (%): (largest patch area / landscape area) * 100
  largest_patch_index <- if (landscape_area_ha > 0) (largest_patch_ha / landscape_area_ha) * 100 else 0
  
  # ED (m/ha): total perimeter of all habitat patches (m) / landscape area (ha)
  patch_polygons <- as.polygons(patch_id_raster, dissolve = TRUE, na.rm = TRUE)
  patch_perimeters_m <- expanse(patch_polygons, unit = "m", what = "perimeter")
  total_edge_m <- sum(patch_perimeters_m, na.rm = TRUE)
  
  edge_density_m_per_ha <- if (landscape_area_ha > 0) total_edge_m / landscape_area_ha else NA_real_
  
  data.frame(
    RasterSet = raster_prefix,
    Year      = year,
    CA        = class_area_ha,
    NP        = number_of_patches,
    AREA_MN   = mean_patch_area_ha,
    LPI       = largest_patch_index,
    ED        = edge_density_m_per_ha
  )
}

# ----------------------------------------------------------
# 3) Run for all years and both raster sets
# ----------------------------------------------------------
metrics_by_set <- lapply(raster_prefixes, function(prefix) {
  lapply(years, compute_fragstats_metrics,
         raster_prefix = prefix,
         raster_dir    = binary_raster_dir,
         directions    = patch_directions) |>
    bind_rows()
}) |> bind_rows()

# Optional rounding (keeps output tidy)
fragstats_metrics <- metrics_by_set %>%
  mutate(
    CA      = round(CA, 3),
    LPI     = round(LPI, 3),
    ED      = round(ED, 3),
    AREA_MN = round(AREA_MN, 3),
    NP      = as.integer(round(NP, 0))
  ) %>%
  arrange(RasterSet, Year)

# ----------------------------------------------------------
# 4) Save outputs
# ----------------------------------------------------------
out_dir <- file.path("output", "national_trends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  fragstats_metrics,
  file = file.path(out_dir, "fragstats_metrics.csv"),
  row.names = FALSE
)

saveRDS(
  fragstats_metrics,
  file = file.path(out_dir, "fragstats_metrics.rds")
)

fragstats_metrics
