# Iberian-sandgrouse

Reproducible **data + R workflow** for the Iberian sandgrouse manuscript:  
(1) spatial habitat suitability modelling, (2) multi-scale habitat trends (2005–2022), and (3) time-varying landscape effects on persistence (FPCA-GLMM).

**Focal species**
- Pin-tailed Sandgrouse (*Pterocles alchata*) — **PTS**
- Black-bellied Sandgrouse (*Pterocles orientalis*) — **BBS**

---

## Overview

This repository accompanies the manuscript:

**Valerio, F.** et al. *Habitat composition and configuration influence persistence of two declining sandgrouse in Spanish agro-steppes.* (in review / preprint pending)

Core workflow:
1. Fit **spatially explicit SDMs** (random forest + Moran eigenvectors) from census occurrence data and environmental predictors.
2. Reconstruct **yearly habitat maps (2005–2022)** and quantify FRAGSTATS-style class metrics at national/regional/local scales.
3. Model **persistence vs extinction** using **FPCA of metric trajectories (2005–2019)** in a **binomial mixed model** framework.

> Tip: start in `scripts/` and run the numbered scripts in order.

---

## Repository structure

├── data/
│ ├── demography_fpca_glmm/
│ │ ├── P_alchata_demography.csv
│ │ ├── P_orientalis_demography.csv
│ │ └── README.md
│ ├── metadata/
│ │ ├── README.md
│ │ ├── codebook_demography_fpca_glmm.csv
│ │ ├── codebook_sdm_occurrence.csv
│ │ └── sites_lookup.csv
│ └── sdm/
│ ├── P_alchata_presence_absence.csv
│ ├── P_orientalis_presence_absence.csv
│ └── README.md
│
├── scripts/
│ ├── 01_screening-boruta.R
│ ├── 02_spatial-SDMs.R
│ ├── 03_threshold-selection.R
│ ├── 04_national-trends_fragstats-metrics.R
│ ├── 05_regional-trends_metrics.R
│ └── 06_fpca-glmm.R
│
├── LICENSE
└── README.md
