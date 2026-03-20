# GEO511 – Master Thesis Programming

This repository contains the code for modelling and predicting accident locations, along with the associated data processing pipeline for the master's thesis (GEO511).

## Contents

- `R/Location_Prediction_Clean.Rmd`  
  Central analysis notebook: data preparation, feature engineering, model training (GBM/RF/XGB), generation of random points, scenario predictions, and raster visualisation.
- `Data/`, `GIS/`, `JMP/`, `Literature/`  
  Data, GIS materials, supplementary analyses, and literature (stored locally; not versioned or only partially included).
- `random_grid.py`, `random_grid.ipynb`, `dissolve.ipynb`  
  Helper scripts/notebooks for supplementary processing steps.

## Requirements

The following R packages are required (minimum):

```r
sf, raster, rasterVis, ggplot2, dplyr, caret, gbm, randomForest, xgboost, future, future.apply, boot, gridExtra
```

> **Note:** Exact package versions are not pinned. If reproducibility is required, please add an `renv` or `packrat` setup.

## Data Structure (Expected)

The R Markdown file expects a defined folder structure, in particular:

- Input data located under `Data/` and/or `GIS/`
- Target output written to `Prediction_Grid/Random/` (created automatically if needed)

Path logic and file names are defined within the notebook. If your local directory structure differs, adjust the paths in `R/Location_Prediction_Clean.Rmd` accordingly.

## Workflow (Overview)

1. Load and clean data
2. Feature computation (e.g. distances to infrastructure objects)
3. Model training (GBM, RF, XGB)
4. Generation of random points within raster cells
5. Scenario predictions + rasterisation and plotting

## Raster and Point Parameters

The notebook iterates over combinations of raster cell sizes and numbers of modelled points:

| Parameter | Values |
|---|---|
| Raster cell sizes (metres) | `4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 70, 80` |
| Number of points | `15,000 · 30,000 · 45,000 · 100,000 · 200,000` |

Generated point files follow this naming convention:

```
point_network_random_<grid>m_<points>.geojson
```

**Example:** `point_network_random_8m_100000.geojson`

## Running the Analysis

Open `R/Location_Prediction_Clean.Rmd` in RStudio and execute the relevant chunks.  
Many steps are set to `eval=FALSE` to reduce runtime. Remove this flag deliberately when needed.

## Output

Output files and maps are configured within the notebook (e.g. written to `Prediction_Grid/Random/` and designated plot output folders).

## Notes

> ⚠️ Processing can be computationally intensive depending on the number of points and cell size selected. Ensure sufficient RAM and runtime, particularly when using 200,000 points. Parallel processing is handled via the `future` package.
