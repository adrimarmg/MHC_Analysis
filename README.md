# MHC_Analysis

This repository contains R scripts and data used for the analysis presented in Morales-Guerrero *et al.* (2026) on immune-related genomic signatures in the human MHC region.

The goal of this repository is to provide fully reproducible code for:
- Sensitivity analysis of historical predictors.
- Estimation of geographic distance from Africa for each population.
- Extraction of climatic variables across different historical time periods.
- Estimation of the area under the curve (AUC) for sociopolitical complexity trajectories.

---

## Repository structure
```text
MHC_Analysis/
├── README.md
├── scripts/
│   ├── 01_sensitivity_analysis.R
│   ├── 02_distance_from_africa.R          # (to be added)
│   ├── 03_climate_extraction.R            # (to be added)
│   └── 04_AUC_sociopolitical_complexity.R # (to be added)
└── data/
    ├── predictors.csv                     # predictors file for sensitivity analysis
    ├── response_variables.csv             # response variable (Tajima's D or heterozygosity)
    ├── coords.csv                         # population coordinates for distance-from-Africa
    ├── climate_rasters                    # (to be added)
    └── sociopolitical_timeseries.csv      # (to be added) complexity trajectories for AUC
---
```
## Scripts

### 1. `01_sensitivity_analysis.R`

**Purpose:**  
Performs the sensitivity analysis of the selected predictors by adding controlled random error to the variables and re-running the regression models.

**Inputs (from `data/`):**
- `predictors.csv` — matrix of predictors (e.g., climate PCs, ΔN–U, pathogen stress, α index, etc.)
- `response_variables.csv` — response variable (e.g., Tajima’s D ond heterozygosity for MHC)

---
