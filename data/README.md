# Data folder

This folder contains the input datasets used in the paper and the associated metadata.

## Metadata files
- `metadata/codebook_sdm_occurrence.csv`: data dictionary for SDM presence/absence inputs.
- `metadata/codebook_demography_fpca_glmm.csv`: data dictionary for demographic status (persistence/extinction) used in FPCA-GLMM analyses.
- `metadata/sites_lookup.csv`: lookup table of unique 10×10 km UTM grid cells.

## Coordinate system
SDM point coordinates are stored as `x`/`y` in meters, projected CRS: WGS 84 / UTM zone 29N (EPSG:32629).
