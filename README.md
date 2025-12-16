# Iberian sandgrouse — habitat dynamics & persistence modelling (Spain, 2005–2022)

![R](https://img.shields.io/badge/R-%E2%89%A5%204.2-blue) ![License: MIT](https://img.shields.io/badge/License-MIT-green)

🧭 **Purpose.** Reproducible **data + R workflow** supporting a manuscript on how **habitat composition and configuration** relate to **persistence/extinction** of two declining sandgrouse in Spanish agro‑steppes.

🐤 **Focal species**
- **Pin‑tailed Sandgrouse** (*Pterocles alchata*) — PTS
- **Black‑bellied Sandgrouse** (*Pterocles orientalis*) — BBS

---

## 💿 What this repository does

This pipeline links three components:

1) **Spatial SDMs (habitat suitability)** from presence/absence data and environmental predictors  
2) **Annual habitat maps (2005–2022)** + FRAGSTATS‑style **habitat amount / composition / configuration** metrics  
3) **Persistence/extinction modelling (2005–2019)** using FPCA–GLMM to estimate **time‑varying** landscape effects

> The repository is organised so analyses can be re‑run from the project root using relative paths.

### 🔎 Study design (at a glance)

- **Occurrences (SDMs):** presence/absence records aggregated to a common grid for modelling.
- **Persistence/extinction:** local status inferred by comparing an early vs a recent survey period.
- **Habitat trends:** yearly suitability/binary habitat maps summarised for **2005–2022**.

---

## 🚀 Quick start

Clone the repo and run the scripts in order:

```bash
git clone https://github.com/FrankVal/Iberian-sandgrouse.git
cd Iberian-sandgrouse

Rscript scripts/01_screening-boruta.R
Rscript scripts/02_spatial-SDMs.R
Rscript scripts/03_threshold-selection.R
Rscript scripts/04_national-trends_fragstats-metrics.R
Rscript scripts/05_regional-trends_metrics.R
Rscript scripts/06_fpca-glmm.R
```

⚠️ Some steps may rely on large rasters / third‑party layers that are not stored in GitHub. See `data/README.md` and the script headers for download / path instructions.

---

## 📁 Repository layout

```text
Iberian-sandgrouse/
├─ data/
│  ├─ sdm/                    # occurrence tables used for SDMs
│  ├─ demography_fpca_glmm/   # cell-level demographic status tables
│  ├─ metadata/               # codebooks + lookup tables
│  └─ README.md               # data inventory / provenance notes
├─ scripts/                   # numbered analysis scripts (run in order)
├─ README.md
├─ LICENSE
└─ .gitignore
```

---

## 📟 Scripts

🟦 **01 — screening & predictor filtering**  
`01_screening-boruta.R` screens candidate predictors (e.g., Boruta / correlation checks) to reduce redundancy before SDMs.

🟩 **02 — spatial SDMs**  
`02_spatial-SDMs.R` fits spatially explicit SDMs (Random Forest with spatial autocorrelation control) and exports suitability predictions.

🟨 **03 — thresholding & binary habitat**  
`03_threshold-selection.R` converts continuous suitability to yearly habitat / non‑habitat maps using an explicit thresholding approach.

🟧 **04 — national habitat trends (2005–2022)**  
`04_national-trends_fragstats-metrics.R` computes national‑level FRAGSTATS‑style class metrics from annual habitat maps.

🟥 **05 — regional trends**  
`05_regional-trends_metrics.R` repeats trend analyses at regional / sub‑national scales and exports summary tables.

🟪 **06 — persistence/extinction modelling**  
`06_fpca-glmm.R` links landscape dynamics to persistence/extinction (10×10 km cells) using FPCA‑GLMM and produces effect curves and model outputs.

---

## 🗃️ Data (what’s included here)

✅ This repository includes key **analysis‑ready tables** and **metadata**, with codebooks in `data/metadata/`.

- `data/sdm/` — presence/absence tables per species  
- `data/demography_fpca_glmm/` — demographic status tables used in persistence/extinction analyses  
- `data/metadata/` — codebooks + lookup tables used across scripts

If you add new datasets, please also update `data/README.md` (source, date, processing notes, and any access constraints).

---

## 🔓 Reproducibility notes

- Run from the **repository root** (avoid absolute paths).  
- Set / record random seeds where applicable (SDM resampling, RF fitting).  
- Consider freezing package versions with `renv` (optional but recommended):

```r
install.packages("renv")
renv::init()      # once, creates renv.lock
renv::restore()   # on a new machine
```

---

## 📌 Citation

Please cite the associated manuscript (details/DOI will be added when available).

Tip: once a DOI is minted (e.g., Zenodo release), add a `CITATION.cff` to enable “Cite this repository” on GitHub.

---

## 🪪 License

This project is released under the **MIT License** (see `LICENSE`).

---

## ✉️ Contact

Francesco Valerio — fvalerio@cibio.up.pt

