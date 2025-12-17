![R](https://img.shields.io/badge/R-%E2%89%A5%204.2-blue)
[![Google Earth Engine](https://img.shields.io/badge/Google%20Earth%20Engine-GEE-4285F4?logo=googleearth&logoColor=white)](https://code.earthengine.google.com/?accept_repo=users/valeriofrank/CorticolIberia)
![License: MIT](https://img.shields.io/badge/License-MIT-green)

![Study area](https://img.shields.io/badge/Area-Continental%20Spain-informational)
![Period](https://img.shields.io/badge/Period-2005%E2%80%932022-informational)
![Resolution](https://img.shields.io/badge/Grain-250%20m-informational)


# Sandgrouse analyses — habitat dynamics & persistence modelling

This repository contains **data + reproducible R code** supporting our analyses of two declining Iberian sandgrouses:

🐤 **Focal species**
- **Pin‑tailed Sandgrouse** (*Pterocles alchata*) — PTS
- **Black‑bellied Sandgrouse** (*Pterocles orientalis*) — BBS

The workflow links **species distribution modelling (SDMs)** with **annual habitat mapping**, **landscape class-level metrics**, and **functional persistence/extinction modelling** using **FPCA–GLMM**.


---

## Summary

Across Iberian agro-steppe landscapes, sandgrouse populations have declined alongside widespread land-use change and agricultural intensification. Here we combine SDMs with time-series habitat mapping and landscape metrics to quantify **where suitable habitat occurs**, **how it changed from 2005–2022**, and **how habitat dynamics relate to persistence/extinction patterns**.

At a high level, the pipeline:

1) fits SDMs using occurrence data and environmental predictors,  
2) produces annual habitat maps,  
3) derives landscape metrics (e.g., habitat amount and configuration) and  
4) models persistence/extinction on landscape trajectories.

---

## Key features

- 🔎 **SDMs** with Random Forests and spatial filtering (Moran eigenvector)
- 🗺️ **Annual habitat maps** (2005–2022) for both species
- 〽️ **Landscape trends** using FRAGSTATS-style metrics (national,regional,local)
- 🚩 **Persistence modelling** via local class-level metrics and sandgrouse status using FPCA–GLMM

---

## Repository structure

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
│       ├── P_alchata_presence_absence.csv
│       ├── P_orientalis_presence_absence.csv
│       └── README.md
├── scripts/
│   ├── 01_screening-boruta.R
│   ├── 02_spatial-SDMs.R
│   ├── 03_threshold-selection.R
│   ├── 04_national-trends_fragstats.R
│   ├── 05_regional-trends_metrics.R
│   └── 06_fpca-glmm.R
├── LICENSE
└── README.md
```

> Each `data/**/README.md` describes the files in that folder and expected formats.

---

## Setup

### 1) Get the data 🗂️ 
Download the dataset from Figshare and place the contents in the repository `data/` folder (keeping the same subfolder structure):

https://doi.org/10.6084/m9.figshare.30898223

### 2) R environment 💻 
We recommend running with a recent R version (≥ 4.2). Install required packages before running the scripts.

Typical dependencies include:  
`terra`, `sf`, `dplyr`, `tidyr`, `ggplot2`, `readr`, `stringr`, `lubridate`,  
`ranger`, `spatialRF` (if used), `landscapemetrics`, `glmmTMB` (and/or similar), and FPCA utilities.

---

## Google Earth Engine (optional but recommended) 🛰️ 

Some inputs (e.g., remote-sensing predictors) can be generated in Google Earth Engine.

Open the shared GEE repository here:  
https://code.earthengine.google.com/?accept_repo=users/valeriofrank/CorticolIberia

---

## Reproduce the analysis 📟 

Run scripts in order from the project root:

1. **`scripts/01_screening-boruta.R`**  
   Variable screening / selection and exploratory checks.

2. **`scripts/02_spatial-SDMs.R`**  
   Fit SDMs for each species (including spatial components).

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

## Outputs ✅ 

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

## Citation 🔖 

If you use this code or dataset, please cite:

- **Dataset (Figshare)**: https://doi.org/10.6084/m9.figshare.30898223  
- **Manuscript**: Valerio, F. *et al.* (in preparation). *[title to be updated]*

A `CITATION.cff` file can be added once the final reference (journal / DOI) is available.

---

## Contact ✉️ 

Francesco Valerio — fvalerio@cibio.up.pt

---

## License

🧾 This repository is released under the terms in `LICENSE` (see file).
