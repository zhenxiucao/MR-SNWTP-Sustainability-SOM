# MR-SNWTP-Sustainability-SOM
R code for SOM-based sustainability pathway classification in the MR-SNWTP recipient region
MR-SNWTP Sustainability SOM
R code for self-organizing map (SOM) clustering and K-means classification of ecological, agricultural, and economic sustainability response patterns in the recipient region of the Middle Route of China's South-to-North Water Transfer Project (MR-SNWTP).
Scope
SOM_classification.R processes three annual response-coefficient rasters (ecological risk index, grain production, and GDP) for each of 2010, 2015, and 2020. It trains a separate SOM for each year, groups the SOM codebook vectors using K-means, exports cluster maps, and produces diagnostics for cluster-number selection and initialization stability.

This repository does not reproduce the upstream regression coefficients, InVEST simulations, constraint-envelope analysis, multinomial logistic regression, hotspot analysis, enrichment analysis, or landscape-fragmentation analysis. Those steps were conducted separately.
Software requirements
R (record the exact version used for the published analysis)
R packages: terra, kohonen, and cluster
Base/recommended R package: stats
The script writes R_Session_Info.txt after a successful execution. Retain this file with the archived analysis materials to document the actual R and package versions.
Input files
Set input_dir in SOM_classification.R to the directory containing the following single-band GeoTIFF files:
Year	ERI response coefficient	Grain-production response coefficient	GDP response coefficient
2010	wyeri10.tif	wyls10.tif	wygdp10.tif
2015	wyeri15.tif	wyls15.tif	wygdp15.tif
2020	wyeri20.tif	wyls20.tif	wygdp20.tif


The coefficient rasters must have a documented geographic reference and compatible units/meaning across the three features. The code aligns predictor rasters to the first (ERI) raster where necessary, and excludes cells with missing or non-finite values. Raster creation is not included in this script. To independently reproduce the maps, users need access to these nine input rasters or an independently specified method for recreating them.
Optional reference rasters for checking the manuscript's pathway labels
The script looks for these previously generated raw K-means maps in reference_dir:
SOM_k6_2010.tif
SOM_k6_2015.tif
SOM_k6_2020.tif
If these files are absent, the script can still produce raw cluster maps and validation diagnostics, but it will not export the manuscript's paper-type maps. Where present, reference rasters must match the input geometry and the rerun must agree with them exactly before paper-type export continues.
Analysis settings
The settings below are specified in the submitted script:
Parameter	Value
Analysis years	2010, 2015, 2020
SOM grid	8 × 8, hexagonal
SOM training iterations (rlen)	200
SOM learning rate (alpha)	0.05 to 0.01
Random seed	2020 (reinitialized for each year)
Maximum number of SOM training cells	50,000 per year
Predictor scaling	Per-year Z-score, based on the sampled SOM training cells
Main K-means solution	K = 6, nstart = 30
Candidate cluster counts	K = 2–10, nstart = 100
Initialization stability	20 repeated K = 6 K-means fits on each fixed SOM codebook


Candidate K diagnostics include silhouette width, within-cluster sum of squares, cluster sizes and response profiles. The script does not establish that K = 6 maximizes silhouette width in every year. Stability checks measure K-means initialization variability on a fixed trained SOM; they do not evaluate stability across repeated SOM training runs.
Running the analysis
1. Install R and the required packages.
2. Download this repository and make the nine coefficient rasters available locally.
3. Edit input_dir, reference_dir, and output_dir near the top of SOM_classification.R. The paths in the originally supplied script are specific to the authors' computer and must be changed for other environments.
4. Run SOM_classification.R in R or RStudio.
5. Inspect the exported raster audit, cluster response profiles, cluster-number diagnostics, and validation summaries.
Main outputs
For each year, the script writes results to an output subdirectory named after that year, including:
SOM_K6_Raw_<year>.tif: raw K-means cluster IDs (not standardized to manuscript Type 1–6 labels).
SOM_K6_Model_<year>.rds: trained SOM, K-means model, feature scaler, and seed.
Original_K6_Response_Profiles_<year>.csv: cluster counts and feature response summaries.
Candidate_K_<year>.csv and Candidate_K_Diagnostics_<year>.png: K = 2–10 diagnostics.
Kmeans_Initialization_Stability_<year>.csv: repeat-fit K-means stability.
Raster_Audit_<year>.csv: actual raster resolution, dimensions, and valid cell count.
Validation_Summary_<year>.csv: selected solution and comparison with any original reference map.
Combined summary files, R_Session_Info.txt, and REPRODUCIBILITY_NOTES.txt are saved to the configured output directory.
Important: manuscript pathway crosswalk is not yet confirmed
Raw K-means class numbers are arbitrary and are not directly comparable among years. The script creates a pathway_type_crosswalk.csv containing a provisional translation from original cluster numbers to the manuscript's six sustainability pathway types, and initializes semantic_mapping_confirmed <- FALSE.

Before releasing verified Fig. 4a paper-type maps, the authors must check the crosswalk against the original cluster response profiles, manuscript pathway definitions, and downstream Fig. 4 analyses. Once verified, mark its review_status values as CONFIRMED, set semantic_mapping_confirmed <- TRUE, and rerun with the original reference rasters available. The provisional mapping must not be presented as independently validated.
Spatial resolution
This script uses the native geometry of the ERI input raster and records its actual resolution in Raster_Audit_<year>.csv. It does not force a 1 km classification grid. The input raster used in the original workflow has been described as approximately 2.516 km, but the definitive resolution must be confirmed from the audit output and reported consistently with the manuscript.
Data and reuse
The repository's code and the data required to reproduce its outputs have different availability requirements. Do not redistribute third-party or restricted data without authorization. Add a public or reviewer-accessible input dataset archive and an explicit software license, with the consent of all relevant authors, before claiming complete independent reproducibility or permissive reuse.
Citation
A persistent Zenodo DOI and a citation to the associated manuscript can be added once the authors have confirmed the code version, pathway mappings, and public release. Do not cite a DOI before the archive exists.
