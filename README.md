# Iberian-sandgrouse

![R](https://img.shields.io/badge/R-%E2%89%A54.2-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

# Iberian-sandgrouse

Reproducible **R + Google Earth Engine (GEE)** workflow supporting the manuscript:

**Valerio, F.**, Mougeot, F., Marques, A. T., Godinho, S., Mendes, T., Gameiro, J., Guise, I., Pacheco, C., Ferraz, G., Benítez‑López, A., & Silva, J. P. *Habitat composition and configuration influence persistence of two declining sandgrouse in Spanish agro‑steppes.* (in preparation).

The repository links:
1) **Spatially explicit SDMs** for Pin‑tailed Sandgrouse (*Pterocles alchata*) and Black‑bellied Sandgrouse (*P. orientalis*),
2) **Multi‑year habitat maps** (2005–2022) and FRAGSTATS‑style landscape metrics, and
3) **Time‑varying effects** of habitat change on local persistence/extinction (2005–2019).

---

## Reproducibility principles (BES-aligned)

- Run everything from the **project root** (relative paths).
- Keep **raw inputs**, **derived data**, and **outputs** separated.
- Document the workflow (this README + per-folder READMEs).
- Record software/package versions (recommended: `renv`).
- Archive code + analysis-ready data in a long-term repository for submission (e.g., Zenodo/Dryad/figshare), if/when allowed by data ownership and sensitivity.

---

## Repository structure

```
Iberian-sandgrouse/
├─ data/
│  ├─ raw/            # raw inputs (some may be excluded due to licensing)
│  ├─ derived/        # processed tables/rasters used by the analyses
│  └─ README.md       # data inventory + download instructions (recommended)
├─ scripts/           # numbered analysis scripts (run in order)
├─ outputs/
│  ├─ figures/
│  ├─ tables/
│  └─ logs/
├─ renv/              # optional: local package library (recommended)
├─ renv.lock          # optional: exact package versions (recommended)
├─ README.md
└─ LICENSE
```

> If your current repository uses slightly different folder names, keep the same logic (raw vs derived vs outputs) and update paths in scripts accordingly.

---

## What each script does

> **Run order matters**: scripts are numbered so outputs from earlier steps become inputs to later steps.

1. **`01_screening-boruta.R`**  
   Screens and reduces candidate predictors using Boruta feature selection and correlation filtering.

2. **`02_spatialRF_SDM.R`**  
   Fits spatially explicit Random Forest SDMs using `spatialRF` (Moran eigenvector maps / MEM) to address spatial autocorrelation; runs repeated train/test splits and exports suitability maps.

3. **`03_threshold_binary_habitat.R`**  
   Converts continuous suitability into habitat/non‑habitat for each year using robust threshold rules (e.g., sensitivity = specificity; prevalence-based thresholds), and applies a consensus threshold across runs.

4. **`04_landscape_metrics_fragstats.R`**  
   Computes annual class-level landscape metrics (CA, NP, AREA_MN, ED, LPI) from the binary habitat maps at national, regional and 10×10 km scales.

5. **`05_persistence_glmm_fpca.R`**  
   Derives demographic status (persistence vs extinction) at 10×10 km, summarises “Start/End/Change” metrics, and fits (a) univariate GLMMs and (b) FPCA‑GLMMs for time‑varying metric effects on persistence.

6. **`06_figures_tables.R`**  
   Reproduces manuscript figures/tables from saved intermediate objects (SDM diagnostics, habitat trends, FPCA‑GLMM coefficient curves, etc.).

---

## Requirements

### Software
- R (recommended ≥ 4.2)
- Optional: RStudio
- If using GEE extraction scripts: a Google account with Earth Engine enabled

### Recommended for fully reproducible runs
Use `renv` to freeze and restore exact package versions:

```r
install.packages("renv")
renv::init()      # first time (creates renv.lock)
renv::restore()   # on a fresh machine
```

---

## How to reproduce the analyses

1. **Clone the repository**
   ```bash
   git clone https://github.com/FrankVal/Iberian-sandgrouse.git
   cd Iberian-sandgrouse
   ```

2. **Restore packages (recommended)**
   ```r
   install.packages("renv")
   renv::restore()
   ```

3. **Prepare data**
   - Put raw inputs in `data/raw/` (or follow `data/README.md`).
   - If some datasets cannot be redistributed (e.g., restricted census data), document how to request/access them and how derived products were created.

4. **Run scripts in order**
   ```r
   source("scripts/01_screening-boruta.R")
   source("scripts/02_spatialRF_SDM.R")
   source("scripts/03_threshold_binary_habitat.R")
   source("scripts/04_landscape_metrics_fragstats.R")
   source("scripts/05_persistence_glmm_fpca.R")
   source("scripts/06_figures_tables.R")
   ```

Outputs (figures/tables/logs) are written to `outputs/` by default.

---

## Data and code availability

- **Code**: this GitHub repository (recommended: create a tagged release and archive it for a DOI at submission).
- **Data**: deposit analysis-ready datasets in a public archive at/near acceptance, or document third‑party/restricted data with access instructions.

---

## How to cite

Add the final manuscript citation here once the reference/DOI is available.

Recommended: add a `CITATION.cff` once you mint a DOI (Zenodo).

---

## License

MIT License (see `LICENSE`).

---

## Contact

Francesco Valerio — fvalerio@cibio.up.pt

