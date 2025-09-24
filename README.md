# Pleurotus pulmonarius — EB Prediction

This repository contains R scripts and data for substrate characterization and 
biological efficiency (EB%) modeling in *Pleurotus pulmonarius*.

## Contents
- **01_pipeline_EB.R** — Main pipeline:
  - Random Forest (Initial vs Consumed composition)
  - Decision Trees (pruned ~8 leaves, CV10 OOF)
  - Bootstrap stability + heatmaps
  - Bland–Altman comparison with new data
- **data/** — Example input data (if shared)
- **outputs/** — Generated results (PNG, CSV, TXT)

## Outputs (generated automatically)
Running `01_pipeline_EB.R` will generate:

- `RF_EB_INICIAL_metrics_split.csv | RF_EB_CONSUMO_metrics_split.csv`
- `RF_EB_INICIAL_metrics_cv10.csv  | RF_EB_CONSUMO_metrics_cv10.csv`
- `RF_EB_INICIAL_importance_bootstrap.csv | RF_EB_CONSUMO_importance_bootstrap.csv`
- `RF_EB_INICIAL_heatmap_top2.png  | RF_EB_CONSUMO_heatmap_top2.png`
- `DT_EB_INICIAL_8leaves.png       | DT_EB_CONSUMO_8leaves.png`
- `DT_EB_INICIAL_rules.txt         | DT_EB_CONSUMO_rules.txt`
- `Predicciones_EB_nuevos_datos.csv`
- `BA_RF_inicial.png | BA_RF_consumo.png | BA_DT_inicial.png | BA_DT_consumo.png`
- `BlandAltman_summary_nuevos_datos.csv`

## Requirements
- R >= 4.2
- Packages: `ranger`, `caret`, `dplyr`, `ggplot2`, `viridisLite`, 
  `readr`, `tidyr`, `rpart`, `rpart.plot`, `scales`

## Usage
Clone the repo and run the pipeline:

```r
source("01_pipeline_EB.R")

