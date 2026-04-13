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
│   ├── 02_distance_from_africa.R
│   ├── 03_PCA_scores_climate_extraction.R
│   └── 04_AUC_sociopolitical_complexity.R # (to be added)
└── data/
    ├── predictors.csv                                  # predictors file for sensitivity analysis
    ├── response_variables.csv                          # response variable (Tajima's D or heterozygosity)
    ├── coords.csv                                      # population coordinates for distance-from-Africa
    ├── bio0000.txt, bio1000.txt, ..., bio20000.txt     # URLs for CHELSA-TraCE21k bioclimatic rasters (bio01–bio19) 
    │                                                   # sampled at 1,000-year intervals from 0 to 20,000 years BP
    └── sociopolitical_timeseries.csv                   # (to be added) complexity trajectories for AUC
---
```
## Scripts

### 1. `01_sensitivity_analysis.R`

**Purpose:**  
Performs the sensitivity analysis of the selected predictors by adding controlled random error to the variables and re-running the regression models.

**Inputs (from `data/`):**
- `predictors.csv` — matrix of predictors (e.g., climate PCs, ΔN–U, pathogen stress, α index, etc.)
- `response_variables.csv` — response variable (e.g., Tajima’s D and heterozygosity for MHC)

### 2. `02_distance_from_africa.R`

**Purpose:**  
Calculates geographic distances from African populations to all other populations by constraining routes through predefined waypoint locations that approximate major human migration paths.

**Inputs (from `data/`):**
- `coords.csv` — Geographic coordinates from each HGDP population

### 3. `03_PCA_scores_climate_extraction.R`

**Purpose:**  
Extracts bioclimatic variables (bio01–bio19) across time slices and performs principal component analysis (PCA) across all time periods.

**Inputs (from `data/`):**
- `bio0000.txt, bio1000.txt, ..., bio20000.txt` — URLs for CHELSA-TraCE21k rasters sampled every 1,000 years from 20,000 years BP to present

### 3. `04_AUC_sociopolitical_complexity.R`

**Purpose:**  
Calculate the area under the curve (AUC) for sociopolitical complexity trajectories of the 54 HGDP populations.

**Inputs (from `data/`):**
- `sociopolitical_timeseries.csv ` — trajectories for AUC

---
