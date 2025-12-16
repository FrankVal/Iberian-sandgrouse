# Iberian Sandgrouse (Spain, 2005–2022) — habitat dynamics & persistence modelling

This repository contains **data + reproducible R code** supporting our analyses of two declining Iberian sandgrouses:

🐤 **Focal species**
- **Pin‑tailed Sandgrouse** (*Pterocles alchata*) — PTS
- **Black‑bellied Sandgrouse** (*Pterocles orientalis*

The workflow links **species distribution modelling (SDMs)** with **annual habitat mapping (2005–2022)**, **landscape change metrics**, and **functional persistence/extinction modelling** using **FPCA–GLMM**.

---

## 🔗 Quick links

- 📦 **Dataset (Figshare)**: https://doi.org/10.6084/m9.figshare.30898223  
- 🛰️ **Google Earth Engine (GEE) code directory**: https://code.earthengine.google.com/?accept_repo=users/valeriofrank/CorticolIberia  
- 🐛 **Issues / questions**: please use GitHub Issues in this repository

---

## 🧭 Summary

Across Iberian agro-steppe landscapes, sandgrouse populations have declined alongside widespread land-use change and agricultural intensification. Here we combine SDMs with time-series habitat mapping and landscape metrics to quantify **where suitable habitat occurs**, **how it changed from 2005–2022**, and **how habitat dynamics relate to persistence/extinction patterns**.

At a high level, the pipeline:

1) fits SDMs using occurrence data and environmental predictors,  
2) produces annual suitability surfaces and binary habitat maps,  
3) derives landscape metrics (e.g., habitat amount and configuration) nationally and by region, and  
4) models persistence/extinction using FPCA–GLMM on demographic trajectories.

---

## ✨ Key features

- 🔎**SDMs** with Random Forest and spatial filtering (Moran eigenvector–style approaches when applicable)
- 🗺️ **Annual habitat maps** (2005–2022) for both species
- 📏 **Landscape trends** using FRAGSTATS-style metrics (national,regional,local)
- 🧬 **Demography + persistence modelling** via FRAGSTATS metrics and sandgrouse status using FPCA–GLMM

---

## 🗂️ Repository structure

```
.
├── data/
│   ├── demography_fpca_glmm/
│   │   ├── P_alchata_demography.csv
│   │   ├── P_orientalis_demography.csv
│   │   └── README.md
│   ├── metadata/
│   │   ├── codebook_demography_fpca_glmm.csv
│   │   ├── codebook_sdm_occurrence.csv
│   │   ├── sites_lookup.csv
│   │   └── README.md
│   └── sdm/
│       ├── P_alchata_presence_absence_*.csv
│       ├── P_orientalis_presence_absence_*.csv
│       └── README.md
├── scripts/
│   ├── 01_screening-boruta.R
│   ├── 02_spatial-SDMs.R
│   ├── 03_threshold-selection.R
│   ├── 04_national-trends_fragstats-*.R
│   ├── 05_regional-trends_metrics.R
│   └── 06_fpca-glmm.R
├── LICENSE
└── README.md
```

> Each `data/**/README.md` describes the files in that folder and expected formats.

---

## 🚀 Setup

### 1) Get the data
Download the dataset from Figshare and place the contents in the repository `data/` folder (keeping the same subfolder structure):

📦 https://doi.org/10.6084/m9.figshare.30898223

### 2) R environment
We recommend running with a recent R version (≥ 4.2). Install required packages before running the scripts.

Typical dependencies include (non-exhaustive):  
`terra`, `sf`, `dplyr`, `tidyr`, `ggplot2`, `readr`, `stringr`, `lubridate`,  
`ranger`, `spatialRF` (if used), `landscapemetrics`, `glmmTMB` (and/or similar), and FPCA utilities.

If your workflow uses **renv**, run:
```r
renv::restore()
```

---

## 🛰️ Google Earth Engine (optional but recommended)

Some inputs (e.g., remote-sensing predictors and/or mapped surfaces) can be generated in Google Earth Engine.

👉 Open the shared GEE repository here:  
https://code.earthengine.google.com/?accept_repo=users/valeriofrank/CorticolIberia

---

## 📟  Reproduce the analysis

Run scripts in order from the project root:

1. **`scripts/01_screening-boruta.R`**  
   Variable screening / selection and exploratory checks.

2. **`scripts/02_spatial-SDMs.R`**  
   Fit SDMs for each species (including spatial components when enabled) and generate predictions.

3. **`scripts/03_threshold-selection.R`**  
   Convert suitability to binary habitat (threshold optimisation) and prepare annual habitat rasters.

4. **`scripts/04_national-trends_fragstats-*.R`**  
   Compute national trajectories of landscape metrics (habitat amount + configuration).

5. **`scripts/05_regional-trends_metrics.R`**  
   Compute regional trajectories and summaries (e.g., by administrative regions / reporting units).

6. **`scripts/06_fpca-glmm.R`**  
   FPCA of demographic curves and GLMM modelling of persistence/extinction responses.

> **Tip:** Scripts are designed to be run sequentially. If you re-run only later scripts, ensure the expected outputs from earlier steps exist.

---

## ✅ Outputs

Most scripts write results (tables/figures/maps) to project output folders (created automatically or specified within each script).
If you prefer a standard layout, we recommend creating:

```
outputs/
├── sdm/
├── habitat_maps/
├── landscape_trends/
└── fpca_glmm/
```

---

## 🔖 Citation

If you use this code or dataset, please cite:

- **Dataset (Figshare)**: https://doi.org/10.6084/m9.figshare.30898223  
- **Manuscript**: Valerio, F. *et al.* (in preparation). *[title to be updated]*

A `CITATION.cff` file can be added once the final reference (journal / DOI) is available.

---

## ✉️ Contact

Francesco Valerio — fvalerio@cibio.up.pt

---

## 🧾 License

This repository is released under the terms in `LICENSE` (see file).
