# Iberian-sandgrouse

![R](https://img.shields.io/badge/R-%E2%89%A54.2-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

Reproducible **data + R workflow** for an Iberian sandgrouse manuscript, covering:

1. **Spatial habitat suitability modelling** (SDMs)
2. **Habitat trends** from **FRAGSTATS-style class metrics** (2005–2022)
3. **Time-varying habitat effects on persistence/extinction** using **FPCA–GLMM** (2005–2019)

**Focal species**
- Pin-tailed Sandgrouse (*Pterocles alchata*) — **PTS**
- Black-bellied Sandgrouse (*Pterocles orientalis*) — **BBS**

> The repository is organised so analyses can be re-run from the project root using **relative paths**.

---

## Citation

Please cite the associated manuscript:

**Valerio, F.** et al. (Year). *Title*. Journal / preprint DOI.

A `CITATION.cff` file can be added once the reference/DOI is final.

---

## Requirements

- **R ≥ 4.2**
- Packages used across scripts (see each script header for exact needs):  
  `terra`, `sf`, `dplyr`, `tidyr`, `purrr`, `ggplot2`, `glmmTMB`, `pROC`, `performance`

Recommended:
- Run from the repository root (RStudio Project or set working directory to repo root).
- Keep large/derived outputs out of version control (see `.gitignore`).

---

## How to run (recommended order)

From the repository root:

```bash
Rscript scripts/01_screening-boruta.R
Rscript scripts/02_spatial-SDMs.R
Rscript scripts/03_threshold-selection.R
Rscript scripts/04_national-trends_fragstats-metrics.R
Rscript scripts/05_regional-trends_metrics.R
Rscript scripts/06_fpca-glmm.R
