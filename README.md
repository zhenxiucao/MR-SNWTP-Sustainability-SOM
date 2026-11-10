# MR-SNWTP Sustainability SOM
Research code for SOM-based classification of sustainability response patterns using ecological risk (ERI), grain production (GP) and economic output (GDP), together with spatial-block cluster-robust multinomial logistic regression for pathway association analysis in the MR-SNWTP recipient region.

## Scope

`SOM_classification.R` processes **three annual response-coefficient rasters** (ecological risk index, grain production, and GDP) for each of **2010, 2015, and 2020**. It trains a separate SOM for each year, groups the SOM codebook vectors using K-means, exports cluster maps, and produces diagnostics for cluster-number selection and initialization stability.

This repository provides code for SOM-based sustainability pathway classification and spatial-block cluster-robust multinomial logistic regression. It does not reproduce the upstream response-coefficient rasters, InVEST simulations, constraint-envelope analysis, hotspot analysis, enrichment analysis or landscape-fragmentation analysis. These steps were conducted separately.

## Software requirements

The SOM classification analysis was implemented in R.

The following software environment is currently configured:

| Software / Package | Version |
|---|---|
| R | 4.5.2 |
| terra | 1.8.93 |
| kohonen | 3.0.13 |
| cluster | 2.1.8.1 |
| Operating system | Windows 11 x64 |

The SOM and K-means analyses use a fixed random seed of 2020.

The exact R session information can be recorded using `sessionInfo()`.

## Input files

Set `input_dir` in `SOM_classification.R` to the directory containing the following single-band GeoTIFF files:

| Year | ERI response coefficient | Grain-production response coefficient | GDP response coefficient |
| --- | --- | --- | --- |
| 2010 | `wyeri10.tif` | `wyls10.tif` | `wygdp10.tif` |
| 2015 | `wyeri15.tif` | `wyls15.tif` | `wygdp15.tif` |
| 2020 | `wyeri20.tif` | `wyls20.tif` | `wygdp20.tif` |

The coefficient rasters must have a documented geographic reference and compatible units/meaning across the three features. The code aligns predictor rasters to the first (ERI) raster where necessary, and excludes cells with missing or non-finite values. **Raster creation is not included in this script.** To independently reproduce the maps, users need access to these nine input rasters or an independently specified method for recreating them.

### Optional reference rasters for checking the manuscript's pathway labels

The script looks for these previously generated raw K-means maps in `reference_dir`:

- `SOM_k6_2010.tif`
- `SOM_k6_2015.tif`
- `SOM_k6_2020.tif`

If these files are absent, the script can still produce raw cluster maps and validation diagnostics, but it **will not export the manuscript's paper-type maps**. Where present, reference rasters must match the input geometry and the rerun must agree with them exactly before paper-type export continues.

## Analysis settings

The settings below are specified in the submitted script:

| Parameter | Value |
| --- | --- |
| Analysis years | 2010, 2015, 2020 |
| SOM grid | 8 × 8, hexagonal |
| SOM training iterations (`rlen`) | 200 |
| SOM learning rate (`alpha`) | 0.05 to 0.01 |
| Random seed | 2020 (reinitialized for each year) |
| Maximum number of SOM training cells | 50,000 per year |
| Predictor scaling | Per-year Z-score, based on the sampled SOM training cells |
| Main K-means solution | K = 6, `nstart = 30` |
| Candidate cluster counts | K = 2–10, `nstart = 100` |
| Initialization stability | 20 repeated K = 6 K-means fits on each fixed SOM codebook |

Candidate K diagnostics include silhouette width, within-cluster sum of squares, cluster sizes and response profiles. The script **does not establish that K = 6 maximizes silhouette width in every year**. Stability checks measure K-means initialization variability on a fixed trained SOM; they do not evaluate stability across repeated SOM training runs.

## Running the analysis

1. Install R and the required packages.
2. Download this repository and make the nine coefficient rasters available locally.
3. Edit `input_dir`, `reference_dir`, and `output_dir` near the top of `SOM_classification.R`. The paths in the originally supplied script are specific to the authors' computer and must be changed for other environments.
4. Run `SOM_classification.R` in R or RStudio.
5. Inspect the exported raster audit, cluster response profiles, cluster-number diagnostics, and validation summaries.

## Main outputs

For each year, the script writes results to an output subdirectory named after that year, including:

- `SOM_K6_Raw_<year>.tif`: raw K-means cluster IDs (not standardized to manuscript Type 1–6 labels).
- `SOM_K6_Model_<year>.rds`: trained SOM, K-means model, feature scaler, and seed.
- `Original_K6_Response_Profiles_<year>.csv`: cluster counts and feature response summaries.
- `Candidate_K_<year>.csv` and `Candidate_K_Diagnostics_<year>.png`: K = 2–10 diagnostics.
- `Kmeans_Initialization_Stability_<year>.csv`: repeat-fit K-means stability.
- `Raster_Audit_<year>.csv`: actual raster resolution, dimensions, and valid cell count.
- `Validation_Summary_<year>.csv`: selected solution and comparison with any original reference map.

Combined summary files, `R_Session_Info.txt`, and `REPRODUCIBILITY_NOTES.txt` are saved to the configured output directory.

## Important: manuscript pathway crosswalk is not yet confirmed

Raw K-means class numbers are arbitrary and are **not directly comparable among years**. The script creates a `pathway_type_crosswalk.csv` containing a **provisional** translation from original cluster numbers to the manuscript's six sustainability pathway types, and initializes `semantic_mapping_confirmed <- FALSE`.

Before releasing verified Fig. 4a paper-type maps, the authors must check the crosswalk against the original cluster response profiles, manuscript pathway definitions, and downstream Fig. 4 analyses. Once verified, mark its `review_status` values as `CONFIRMED`, set `semantic_mapping_confirmed <- TRUE`, and rerun with the original reference rasters available. The provisional mapping **must not be presented as independently validated**.

## Spatial resolution

This script uses the **native geometry of the ERI input raster** and records its actual resolution in `Raster_Audit_<year>.csv`. It does not force a 1 km classification grid. The input raster used in the original workflow has been described as approximately 2.516 km, but the definitive resolution must be confirmed from the audit output and reported consistently with the manuscript.

## Data and reuse

The repository's code and the data required to reproduce its outputs have different availability requirements. Do not redistribute third-party or restricted data without authorization. The source code is released under the MIT License (see LICENSE). The input datasets required to reproduce the analyses are subject to their respective data-access and redistribution conditions. Users must obtain the necessary inputs from their original sources or authorized providers.

## Citation

The software archive (v1.0.0) is available on Zenodo:
https://doi.org/10.5281/zenodo.23883896
This DOI identifies the archived v1.0.0 release. The latest development version is available from the GitHub repository. When citing the software used for a specific analysis, the corresponding archived version should be referenced..
