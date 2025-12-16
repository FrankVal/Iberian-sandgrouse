setwd("E:/0Test")
# Create folders
dir.create("data", showWarnings = FALSE)
dir.create(file.path("data", "metadata"), recursive = TRUE, showWarnings = FALSE)


readme <- c(
  "# Data folder",
  "",
  "This folder contains the input datasets used in the paper and the associated metadata.",
  "",
  "## Metadata files",
  "- `metadata/codebook_sdm_occurrence.csv`: data dictionary for SDM presence/absence inputs.",
  "- `metadata/codebook_demography_fpca_glmm.csv`: data dictionary for demographic inputs used in FPCA-GLMM analyses.",
  "- `metadata/sites_lookup.csv`: stable site_id lookup (and optional coarse spatial descriptors).",
  "",
  "## Coordinate system",
  "SDM point coordinates are stored as `x`/`y` in meters, projected CRS: WGS 84 / UTM zone 29N (EPSG:32629).",
  "",
  "## Notes on sensitive locations",
  "If precise coordinates are sensitive, public releases may use coarsened coordinates or grid IDs. Document any anonymisation here."
)

writeLines(readme, "data/README.md")



codebook_sdm <- data.frame(
  field_name     = c("species","site_id","year","pa","x","y","epsg"),
  type           = c("character","character","integer","integer","numeric","numeric","integer"),
  units          = c("","", "year","", "meters","meters",""),
  allowed_values = c("P_alchata|P_orientalis","","","0|1","","","32629"),
  description    = c(
    "Species code",
    "Stable site identifier used across datasets",
    "Observation year (if applicable; otherwise document fixed-year dataset)",
    "Presence (1) / absence (0)",
    "UTM easting (WGS84 / UTM 29N)",
    "UTM northing (WGS84 / UTM 29N)",
    "EPSG code of coordinate reference system"
  ),
  required       = c("yes","recommended","recommended","yes","yes","yes","recommended"),
  notes          = c(
    "If one file per species, species column can be omitted",
    "Prefer an existing unique ID; otherwise create and freeze it",
    "If no year column exists, document temporal meaning in README",
    "",
    "",
    "",
    "Your layer reports WGS 84 / UTM zone 29N"
  ),
  stringsAsFactors = FALSE
)

write.csv(codebook_sdm, "data/metadata/codebook_sdm_occurrence.csv", row.names = FALSE)






codebook_demo <- data.frame(
  field_name     = c("species","site_id","year","occupied","persistence","extinction","n_visits","effort_km","region"),
  type           = c("character","character","integer","integer","integer","integer","integer","numeric","character"),
  units          = c("","", "year","","","","visits","km",""),
  allowed_values = c("P_alchata|P_orientalis","","","0|1","0|1","0|1",">=0",">=0",""),
  description    = c(
    "Species code",
    "Stable site identifier used across years",
    "Year t",
    "Occupancy state in year t (if stored directly)",
    "1 if occupied at t and t+1 (if used as response)",
    "1 if occupied at t and empty at t+1 (if used as response)",
    "Number of visits/surveys (if available)",
    "Survey effort in km (if available)",
    "Region/stratum (if used as random effect)"
  ),
  required       = c("yes","yes","yes","conditional","optional","optional","optional","optional","optional"),
  notes          = c(
    "If one file per species, species column can be omitted",
    "",
    "",
    "If you store annual states, persistence/extinction can be derived in scripts",
    "Define the rule explicitly in scripts + README",
    "Define the rule explicitly in scripts + README",
    "",
    "",
    ""
  ),
  stringsAsFactors = FALSE
)

write.csv(codebook_demo, "data/metadata/codebook_demography_fpca_glmm.csv", row.names = FALSE)























library(sf)
library(dplyr)

# --- Inputs ---
bbs_path <- "E:/000JoaoPauloSilva/000Francois/Dados/0Extract/00Final_Filling/Filled_Daily_Year_Static_UTM.shp"
pts_path <- "E:/000JoaoPauloSilva/000Francois/000Pterocles_alchata/Dados/0Extract/00Final_Filling/Original/Filled_Daily_Year_Static_UTM.shp"

bbs <- read_sf(bbs_path) |> mutate(species = "BBS")  # orientalis
pts <- read_sf(pts_path) |> mutate(species = "PTS")  # alchata

occ <- bind_rows(bbs, pts)

# Ensure UTM 29N
if (is.na(st_crs(occ))) stop("CRS is missing in one of the layers.")
if (st_crs(occ)$epsg != 32629) occ <- st_transform(occ, 32629)

# 10x10 km cell (lower-left corner) from x/y
xy <- st_coordinates(occ)
e0 <- floor(xy[,1] / 10000) * 10000
n0 <- floor(xy[,2] / 10000) * 10000

sites <- tibble(
  site_id    = paste0("UTM29N_E", e0, "_N", n0),
  cell_e0    = e0,
  cell_n0    = n0,
  x_centroid = e0 + 5000,
  y_centroid = n0 + 5000,
  epsg       = 32629,
  species    = occ$species
) |>
  distinct() |>
  group_by(site_id, cell_e0, cell_n0, x_centroid, y_centroid, epsg) |>
  summarise(
    in_PTS = as.integer(any(species == "PTS")),
    in_BBS = as.integer(any(species == "BBS")),
    .groups = "drop"
  )

# Add lon/lat centroids (optional)
cent_sf <- st_as_sf(sites, coords = c("x_centroid","y_centroid"), crs = 32629)
cent_ll <- st_transform(cent_sf, 4326)
ll <- st_coordinates(cent_ll)
sites$lon_centroid <- ll[,1]
sites$lat_centroid <- ll[,2]

sites$cell_e0<-NULL
sites$cell_n0<-NULL

# Write to your repo
out_dir <- "E:/Iberian-sandgrouse/data/metadata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(sites, file.path(out_dir, "sites_lookup.csv"), row.names = FALSE)
message("Wrote: ", file.path(out_dir, "sites_lookup.csv"))


###############################
library(readr)
library(tibble)

codebook_path <- "E:/Iberian-sandgrouse/data/metadata/codebook_demography_fpca_glmm.csv"

codebook <- tibble(
  field_name = c("UTMCODE","Status","CA","AREA_MN","ED","LPI","NP","Year","RegSect"),
  type = c("character","character","numeric","numeric","numeric","numeric","integer","integer","character"),
  units = c("", "", "km^2", "ha", "m/ha", "%", "n", "year", ""),
  allowed_values = c("", "Persistent|Extinct", "", "", "", "", "", "", ""),
  description = c(
    "10×10 km UTM grid cell code (UTMCODE) used as spatial unit (Mougeot et al., 2021).",
    "Demographic status of the cell over the interval: Persistent or Extinct.",
    "Class area (CA) of the focal habitat class within the 10×10 km cell.",
    "Mean patch area (AREA_MN) of the focal habitat class within the 10×10 km cell.",
    "Edge density (ED) of the focal habitat class within the 10×10 km cell.",
    "Largest patch index (LPI) of the focal habitat class within the 10×10 km cell.",
    "Number of patches (NP) of the focal habitat class within the 10×10 km cell.",
    "Reference year associated with the landscape metrics/status.",
    "Regional sector identifier"
  ),
  required = c("yes","yes","yes","yes","yes","yes","yes","yes","recommended"),
  notes = c(
    "Primary key for joins across demography and landscape metrics.",
    "Define the exact rule for Persistent/Extinct in Methods.",
    "In the manuscript/figures this metric is labelled CA (km^2).",
    "Units follow the figure (ha).",
    "Units follow the figure (m/ha).",
    "Units follow the figure (%).",
    "Units follow the figure (count).",
    "",
    ""
  )
)

write_csv(codebook, codebook_path)
message("Updated: ", codebook_path)





















library(sf)

# 1) Read layer
P_alchata <- st_read(
  dsn   = "E:/000JoaoPauloSilva/000Francois/000Pterocles_alchata/Dados/0Extract/00Final_Filling/Original",
  layer = "Filled_Daily_Year_Static_UTM",
  quiet = TRUE
)

# 2) Drop unwanted columns (only if present)
drop_cols <- c(
  "long","lat","ID","timestamp2","Topographi","OBJECTID",
  "ed_500","ed_1000","pland_500","pland_1000","shdi_500","shdi_1000",
  "IMD","Slope"
)

P_alchata <- P_alchata[, setdiff(names(P_alchata), intersect(names(P_alchata), drop_cols))]

# 3) Rename livestock masked suffix: _E_1 -> _Ext_M
# (safe: does not change _E_10 etc.)
names(P_alchata) <- gsub("_E_1(?![0-9])", "_Ext_M", names(P_alchata), perl = TRUE)

# 4) Rename timestamp -> date (if exists)
if ("timestamp" %in% names(P_alchata)) {
  names(P_alchata)[names(P_alchata) == "timestamp"] <- "date"
}

if ("Dist_Roads" %in% names(P_alchata)) {
  names(P_alchata)[names(P_alchata) == "Dist_Roads"] <- "Dist_MR"
}


# 5) Add x/y as longitude/latitude from geometry (EPSG:4326)
xy_ll <- st_coordinates(st_transform(P_alchata, 4326))
P_alchata$x <- xy_ll[,1]  # longitude
P_alchata$y <- xy_ll[,2]  # latitude

# Check
str(P_alchata)
names(P_alchata)


library(sf)
library(readr)

# Ensure CRS is UTM 29N
if (is.na(st_crs(P_alchata))) stop("P_alchata has no CRS set.")
if (st_crs(P_alchata)$epsg != 32629) P_alchata <- st_transform(P_alchata, 32629)

# Add x/y in meters (UTM 29N)
xy_utm <- st_coordinates(P_alchata)
P_alchata$x <- xy_utm[,1]
P_alchata$y <- xy_utm[,2]

ncol(P_alchata)

P_alchata$Species<-"P_alchata"

# Export to CSV (drop geometry)
out_dir  <- "E:/Iberian-sandgrouse/data/sdm"
out_file <- file.path(out_dir, "P_alchata_presence_absence.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(st_drop_geometry(P_alchata), out_file)

message("Wrote: ", out_file)




















library(sf)

# 1) Read layer
P_orientalis <- st_read(
  dsn   = "E:/000JoaoPauloSilva/000Francois/Dados/0Extract/00Final_Filling",
  layer = "Filled_Daily_Year_Static_UTM",
  quiet = TRUE
)

# 2) Drop unwanted columns (only if present)
drop_cols <- c(
  "long","lat","ID","timestamp2","Topographi","OBJECTID",
  "ed_500","ed_1000","pland_500","pland_1000","shdi_500","shdi_1000",
  "IMD","Slope"
)

P_orientalis <- P_orientalis[, setdiff(names(P_orientalis), intersect(names(P_orientalis), drop_cols))]

# 3) Rename livestock masked suffix: _E_1 -> _Ext_M
# (safe: does not change _E_10 etc.)
names(P_orientalis) <- gsub("_E_1(?![0-9])", "_Ext_M", names(P_orientalis), perl = TRUE)

# 4) Rename timestamp -> date (if exists)
if ("timestamp" %in% names(P_orientalis)) {
  names(P_orientalis)[names(P_orientalis) == "timestamp"] <- "date"
}

if ("Dist_Roads" %in% names(P_orientalis)) {
  names(P_orientalis)[names(P_orientalis) == "Dist_Roads"] <- "Dist_MR"
}


# 5) Add x/y as longitude/latitude from geometry (EPSG:4326)
xy_ll <- st_coordinates(st_transform(P_orientalis, 4326))
P_orientalis$x <- xy_ll[,1]  # longitude
P_orientalis$y <- xy_ll[,2]  # latitude

# Check
str(P_orientalis)
names(P_orientalis)
ncol(P_orientalis)

library(sf)
library(readr)

# Ensure CRS is UTM 29N
if (is.na(st_crs(P_orientalis))) stop("P_orientalis has no CRS set.")
if (st_crs(P_orientalis)$epsg != 32629) P_orientalis <- st_transform(P_orientalis, 32629)

# Add x/y in meters (UTM 29N)
xy_utm <- st_coordinates(P_orientalis)
P_orientalis$x <- xy_utm[,1]
P_orientalis$y <- xy_utm[,2]

P_orientalis$Species<-"P_orientalis"
ncol(P_orientalis)

# Export to CSV (drop geometry)
out_dir  <- "E:/Iberian-sandgrouse/data/sdm"
out_file <- file.path(out_dir, "P_orientalis_presence_absence.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(st_drop_geometry(P_orientalis), out_file)

message("Wrote: ", out_file)






library(sf)
library(dplyr)
library(readr)

# ---- 0) Read & prep data ----
polys_metrics_all <- sf::st_read(
  "E:/000JoaoPauloSilva/000Francois/000Paper/Chapter_1/4StatusExplanation/Alchata/Extraction/Alchata_Status.shp",
  quiet = TRUE
)
polys_df <- sf::st_drop_geometry(polys_metrics_all)

# ---- 1) Keep only requested columns + rename ----
keep <- c("UTMCODE","Status","AREA_MN","ED","LPI","NP","Year","TA","ncPALCH")

polys_df <- polys_df %>%
  select(any_of(keep)) %>%
  rename(
    RegSect = ncPALCH,
    CA      = TA
  ) %>%
  mutate(Species = "P_alchata") %>%
  select(Species, everything())   # put Species first (optional)

# ---- 2) Write to repo folder ----
out_dir  <- "E:/Iberian-sandgrouse/data/demography_fpca_glmm"
out_file <- file.path(out_dir, "P_alchata_demography_long.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(polys_df, out_file)
message("Wrote: ", out_file)











library(sf)
library(dplyr)
library(readr)

# ---- 0) Read & prep data ----
polys_metrics_all <- sf::st_read(
  "E:/000JoaoPauloSilva/000Francois/000Paper/Chapter_1/4StatusExplanation/Alchata/Extraction/Alchata_Status.shp",
  quiet = TRUE
)
polys_df <- sf::st_drop_geometry(polys_metrics_all)

# ---- 1) Keep only requested columns + rename ----
keep <- c("UTMCODE","Status","AREA_MN","ED","LPI","NP","Year","TA","ncPALCH")

polys_df <- polys_df %>%
  select(any_of(keep)) %>%
  rename(
    RegSect = ncPALCH,
    CA      = TA
  ) %>%
  mutate(Species = "P_alchata") %>%
  select(Species, everything())   # put Species first (optional)

# ---- 2) Write to repo folder ----
out_dir  <- "E:/Iberian-sandgrouse/data/demography_fpca_glmm"
out_file <- file.path(out_dir, "P_alchata_demography_long.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(polys_df, out_file)
message("Wrote: ", out_file)






















library(sf)
library(dplyr)
library(readr)

# ---- 0) Read & prep data ----
polys_metrics_all <- sf::st_read(
  "E:/000JoaoPauloSilva/000Francois/000Paper/Chapter_1/4StatusExplanation/Orientalis/Extraction/Orientalis_Status.shp",
  quiet = TRUE
)
polys_df <- sf::st_drop_geometry(polys_metrics_all)

# ---- 1) Keep only requested columns + rename ----
keep <- c("UTMCODE","Status","AREA_MN","ED","LPI","NP","Year","TA","nclPORI")

polys_df <- polys_df %>%
  select(any_of(keep)) %>%
  rename(
    RegSect = nclPORI,
    CA      = TA
  ) %>%
  mutate(Species = "P_orientalis") %>%
  select(Species, everything())   # put Species first (optional)

# ---- 2) Write to repo folder ----
out_dir  <- "E:/Iberian-sandgrouse/data/demography_fpca_glmm"
out_file <- file.path(out_dir, "P_orientalis_demography_long.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(polys_df, out_file)
message("Wrote: ", out_file)











cd E:\Iberian-sandgrouse
git add data/demography_fpca_glmm/P_orientalis_demography_long.csv 
git commit -m "Add P. orientalis demography long table" 
git pull --rebase origin main
git push