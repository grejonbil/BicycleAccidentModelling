# Predicting bicycle accident locations from road infrastructure

Code and results for the article

> Biland, G., Fu, C., Brucks, W., Weibel, R. *Speed reduction alone may not be enough for bicycle safety: Predicting bicycle accident locations from road infrastructure using machine learning.* Geo-spatial Information Science (under review).

Police-reported bicycle accidents in the City of Zurich (2016–2025) are modelled against pseudo-absence points on the cycling network. Gradient Boosting (GB), Random Forest (RF) and XGBoost (XGB) classifiers are trained on 14 infrastructure and traffic features, evaluated with random and spatial block cross-validation and probabilistic scoring rules, interpreted with SHAP values and partial dependence, and applied to six infrastructure scenarios with bootstrap uncertainty estimates.

## Repository structure

```         
R/
  prepare_data.R         Downloads and processes all input data: accidents, infrastructure layers,
                         derived junctions, pseudo-absence points, modelling table, prediction grids
  features.R             Minimum distance intersection calculation (MDIC): features of any point set
  setup.R                Shared configuration: data loading, model specification, prediction grids,
                         scenario definitions (Table 1), file paths
  model_training.R       Hold-out evaluation (Table 3) and RMSE of the prediction grids (Table 2)
  model_evaluation.R     Validation, sensitivity, interpretation and scenario analyses
                         (Tables 4-6, S2-S7) and the data of the figures
  figures.Rmd            All figures (Figures 4-9, S1), one chunk per figure
  install_packages.R     Installs the required R packages (versions used)
data/
  complex_intersections_2023.gpkg   Complex intersections identified in the original study
  pseudo_absence_points_2023.gpkg   Pseudo-absence points of the original study
results/
  tables/                Results tables (CSV)
  figure_data/           Model outputs used by figures.Rmd (scenario predictions, SHAP values)
  figures/               Figures 4-9 and S1 (PNG)
  key_results.csv        Single values quoted in the text
docs/
  data.md                Data sources, variables and directory structure
```

## Data

The input data are open government data of the City and Canton of Zurich, the swissTLM3D landscape model and OpenStreetMap, plus one restricted dataset (the GVM-ZH 2018 traffic model, available on request from the Office of Transport of the Canton of Zurich). Accidents cover 2016–2025; the infrastructure layers represent the state of 2025/26 (bicycle network: 2022). `prepare_data.R` downloads and processes all open data; sources, versions, variables and the directory structure are documented in [`docs/data.md`](docs/data.md).

The data directory defaults to `Data/` next to `R/` and can be set with the environment variable `BAM_DATA_DIR`; results are written to `results/` (or `BAM_RESULTS_DIR`).

## Requirements

- R 4.3 (results were produced with R 4.3.3)
- System libraries for `sf` (GDAL, GEOS, PROJ, udunits) and GLPK (for `igraph`)
- R packages: `Rscript R/install_packages.R`. `caret` requires `xgboost` \< 2.0; the script installs xgboost 1.7.7.1.

## Workflow

Run all scripts from the `R/` directory:

``` bash
cd R
Rscript prepare_data.R          # 1. data preparation
Rscript model_training.R        # 2. Tables 2 and 3
Rscript model_evaluation.R      # 3. all other results (or selected stages, e.g. "cv metrics")
Rscript -e 'rmarkdown::render("figures.Rmd")'   # 4. figures
```

The figures can also be edited and drawn one by one by running the chunks of `figures.Rmd` in RStudio (working directory `R/`). They only read the outputs of steps 2 and 3 and the prepared data, so no model is re-fitted.

The number of cores used for cross-validation and bootstrapping can be set with `BAM_CORES` (default: all but one).

## Citation

If you use this code, please cite the article above.
