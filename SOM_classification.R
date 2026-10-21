# SOM pathway classification and cluster-number validation
# Input: annual ERI, grain-production, and GDP response-coefficient GeoTIFFs
# Years: 2010, 2015, 2020
# This script reproduces the original year-specific SOM + K-means analysis.
# It does not construct the upstream regression-coefficient rasters.

required_packages <- c("terra", "kohonen", "cluster")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(terra)
  library(kohonen)
  library(cluster)
})

# User configuration
input_dir <- "C:/Users/DELL/Desktop/SOM/01data"
reference_dir <- "C:/Users/DELL/Desktop/SOM/SOM_outputs_separate"
output_dir <- "C:/Users/DELL/Desktop/SOM/SOM_outputs_publication"
years <- c(2010L, 2015L, 2020L)

random_seed <- 2020L
max_training_pixels <- 50000L
som_grid_x <- 8L
som_grid_y <- 8L
som_grid_topology <- "hexagonal"
som_training_length <- 200L
som_learning_rate <- c(0.05, 0.01)
selected_k <- 6L
candidate_k <- 2:10
original_kmeans_nstart <- 30L
validation_kmeans_nstart <- 100L
kmeans_stability_repeats <- 20L

# The mapping in pathway_type_crosswalk.csv is a scientific interpretation.
# It must be checked against the authors' original pathway definitions and
# downstream Fig. 4 analysis before a publication-type map is considered final.
semantic_mapping_confirmed <- FALSE
crosswalk_path <- file.path(output_dir, "pathway_type_crosswalk.csv")

if (!selected_k %in% candidate_k) stop("selected_k must occur in candidate_k")
if (!length(years)) stop("No years specified")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pathway_names <- c(
  "Agricultural gains and ecological improvement",
  "Weak coupling",
  "Economic and agricultural gains with ecological pressure",
  "Economic losses and agricultural gains",
  "Economic gains and agricultural losses",
  "Multidimensional sustainability synergy"
)

# Provisional map from the previously reported scientific pathway definitions
# to the original year-specific raw K-means cluster IDs.
provisional_mapping <- data.frame(
  year = rep(years, each = selected_k),
  paper_type = rep(seq_len(selected_k), times = length(years)),
  original_cluster = c(
    6L, 4L, 5L, 1L, 2L, 3L,
    5L, 6L, 2L, 1L, 3L, 4L,
    5L, 4L, 2L, 1L, 6L, 3L
  ),
  pathway_name = rep(pathway_names, times = length(years)),
  review_status = "PROVISIONAL_VERIFY_WITH_MANUSCRIPT",
  stringsAsFactors = FALSE
)
if (!file.exists(crosswalk_path)) {
  write.csv(provisional_mapping, crosswalk_path, row.names = FALSE)
}

mapping <- read.csv(crosswalk_path, stringsAsFactors = FALSE)
needed_columns <- c("year", "paper_type", "original_cluster", "pathway_name")
if (!all(needed_columns %in% names(mapping))) {
  stop("The crosswalk CSV is missing required columns")
}
for (yr in years) {
  sub <- mapping[mapping$year == yr, , drop = FALSE]
  if (nrow(sub) != selected_k ||
      !setequal(sub$paper_type, seq_len(selected_k)) ||
      !setequal(sub$original_cluster, seq_len(selected_k))) {
    stop("Crosswalk must be a one-to-one six-class mapping for year ", yr)
  }
}
if (semantic_mapping_confirmed &&
    ("review_status" %in% names(mapping)) &&
    any(mapping$review_status != "CONFIRMED")) {
  stop("Set every crosswalk review_status to CONFIRMED before publishing type maps")
}

input_paths <- function(year) {
  suffix <- sprintf("%02d", year %% 100L)
  c(
    b_ERI = file.path(input_dir, paste0("wyeri", suffix, ".tif")),
    b_GP = file.path(input_dir, paste0("wyls", suffix, ".tif")),
    b_GDP = file.path(input_dir, paste0("wygdp", suffix, ".tif"))
  )
}

read_predictors <- function(year) {
  paths <- input_paths(year)
  if (!all(file.exists(paths))) {
    stop("Missing input raster(s) for ", year, ": ",
         paste(paths[!file.exists(paths)], collapse = "; "))
  }
  rr <- lapply(paths, terra::rast)
  if (any(vapply(rr, terra::nlyr, integer(1)) != 1L)) {
    stop("Each input must have exactly one layer")
  }
  for (j in 2:3) {
    if (!terra::same.crs(rr[[1]], rr[[j]])) {
      message("Projecting predictor ", names(paths)[j], " in ", year)
      rr[[j]] <- terra::project(rr[[j]], rr[[1]], method = "bilinear")
    } else if (!terra::compareGeom(rr[[1]], rr[[j]], stopOnError = FALSE)) {
      message("Resampling predictor ", names(paths)[j], " in ", year)
      rr[[j]] <- terra::resample(rr[[j]], rr[[1]], method = "bilinear")
    }
  }
  s <- c(rr[[1]], rr[[2]], rr[[3]])
  if (!inherits(s, "SpatRaster") || terra::nlyr(s) != 3L) {
    stop("Failed to build the three-layer predictor raster")
  }
  names(s) <- names(paths)
  s
}

read_feature_matrix <- function(s) {
  x <- terra::values(s, mat = TRUE)
  if (!is.matrix(x)) x <- matrix(x, ncol = 3L)
  colnames(x) <- names(s)
  storage.mode(x) <- "double"
  x
}

valid_feature_rows <- function(x) {
  complete.cases(x) & rowSums(!is.finite(x)) == 0L
}

make_scaler <- function(x) {
  mu <- colMeans(x)
  sigma <- apply(x, 2L, stats::sd)
  if (any(!is.finite(mu)) || any(!is.finite(sigma)) || any(sigma <= 0)) {
    stop("Invalid scaling parameters")
  }
  list(mu = mu, sigma = sigma)
}

scale_features <- function(x, scaler) {
  sweep(sweep(as.matrix(x), 2L, scaler$mu, "-"), 2L, scaler$sigma, "/")
}

choose_two <- function(z) z * (z - 1) / 2

adjusted_rand_index <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]
  b <- b[ok]
  if (length(a) < 2L) return(NA_real_)
  tab <- table(a, b)
  n <- sum(tab)
  a_pairs <- sum(choose_two(rowSums(tab)))
  b_pairs <- sum(choose_two(colSums(tab)))
  both <- sum(choose_two(tab))
  expectation <- a_pairs * b_pairs / choose_two(n)
  maximum <- (a_pairs + b_pairs) / 2
  if (abs(maximum - expectation) < 1e-12) return(1)
  (both - expectation) / (maximum - expectation)
}

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(j) {
    tail <- permutations(x[-j])
    cbind(rep(x[j], nrow(tail)), tail)
  }))
}

align_cluster_labels <- function(old, recent, k) {
  tab <- table(factor(old, levels = seq_len(k)),
               factor(recent, levels = seq_len(k)))
  ps <- permutations(seq_len(k))
  matching <- apply(ps, 1L, function(z) sum(tab[cbind(z, seq_len(k))]))
  best <- as.integer(ps[which.max(matching), ])
  list(new_to_original = best,
       matched_fraction = max(matching) / length(old),
       ari = adjusted_rand_index(old, recent),
       contingency = tab)
}

save_class_raster <- function(s, valid_indices, labels, filename, layer_name) {
  v <- rep(NA_integer_, terra::ncell(s))
  v[valid_indices] <- as.integer(labels)
  out <- terra::setValues(terra::rast(s[[1L]]), v)
  names(out) <- layer_name
  terra::writeRaster(out, filename, overwrite = TRUE,
                     datatype = "INT1U", NAflag = 255)
  invisible(filename)
}

response_profiles <- function(x, labels, year, k, method) {
  rows <- vector("list", k)
  for (group in seq_len(k)) {
    subset <- x[which(labels == group), , drop = FALSE]
    n_group <- nrow(subset)
    rec <- data.frame(
      year = year, k = k, method = method, cluster = group,
      n_pixels = n_group,
      percentage = 100 * n_group / nrow(x)
    )
    for (feature in colnames(x)) {
      z <- subset[, feature]
      rec[[paste0(feature, "_mean")]] <- if (n_group) mean(z) else NA_real_
      rec[[paste0(feature, "_median")]] <- if (n_group) median(z) else NA_real_
      rec[[paste0(feature, "_sd")]] <- if (n_group > 1L) sd(z) else NA_real_
      rec[[paste0(feature, "_positive_fraction")]] <- if (n_group) mean(z > 0) else NA_real_
      rec[[paste0(feature, "_negative_fraction")]] <- if (n_group) mean(z < 0) else NA_real_
    }
    rows[[group]] <- rec
  }
  do.call(rbind, rows)
}

compare_reference <- function(year, s, valid_indices, new_labels, year_dir) {
  fn <- file.path(reference_dir, sprintf("SOM_k6_%d.tif", year))
  result <- list(
    year = year, reference_found = FALSE,
    ari = NA_real_, best_match_fraction = NA_real_,
    new_to_original = NULL
  )
  if (!file.exists(fn)) {
    message("Reference raster not found; label mapping skipped: ", fn)
    return(result)
  }
  old <- terra::rast(fn)
  if (!terra::compareGeom(old, s[[1L]], stopOnError = FALSE)) {
    stop("Reference raster geometry does not match the predictor raster for ", year)
  }
  z <- as.vector(terra::values(old, mat = TRUE)[, 1L])
  old_labels <- z[valid_indices]
  keep <- is.finite(old_labels) & old_labels %in% seq_len(selected_k)
  if (!all(keep)) {
    stop("Original reference has missing or unexpected class IDs for year ", year)
  }
  match <- align_cluster_labels(as.integer(old_labels), new_labels, selected_k)
  write.csv(as.data.frame(match$contingency),
            file.path(year_dir, sprintf("Rerun_vs_Original_Contingency_%d.csv", year)),
            row.names = FALSE)
  write.csv(data.frame(
    year = year,
    new_cluster = seq_len(selected_k),
    original_cluster = match$new_to_original
  ), file.path(year_dir, sprintf("Rerun_to_Original_Clusters_%d.csv", year)),
  row.names = FALSE)
  result$reference_found <- TRUE
  result$ari <- match$ari
  result$best_match_fraction <- match$matched_fraction
  result$new_to_original <- match$new_to_original
  result
}

run_year <- function(year) {
  message("--- Processing ", year, " ---")
  year_dir <- file.path(output_dir, as.character(year))
  dir.create(year_dir, showWarnings = FALSE, recursive = TRUE)
  s <- read_predictors(year)
  x <- read_feature_matrix(s)
  valid <- valid_feature_rows(x)
  indices <- which(valid)
  X <- x[valid, , drop = FALSE]
  n_valid <- nrow(X)
  if (n_valid < 100L) stop("Too few valid pixels in ", year)

  resolution <- terra::res(s)
  pixel_area <- if (!terra::is.lonlat(s)) prod(resolution) / 1e6 else NA_real_
  audit <- data.frame(
    year = year,
    input_resolution_x = resolution[1L],
    input_resolution_y = resolution[2L],
    crs = terra::crs(s, proj = TRUE),
    n_total_cells = terra::ncell(s),
    n_valid_cells = n_valid,
    pixel_area_km2 = pixel_area,
    valid_area_km2 = n_valid * pixel_area
  )
  write.csv(audit, file.path(year_dir, sprintf("Raster_Audit_%d.csv", year)),
            row.names = FALSE)
  paths <- input_paths(year)
  inventory <- data.frame(
    year = year,
    feature = names(paths),
    filename = basename(paths),
    file_md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
  write.csv(inventory, file.path(year_dir, sprintf("Input_MD5_%d.csv", year)),
            row.names = FALSE)

  # Keep this random-number sequence identical to the validated original workflow.
  set.seed(random_seed)
  n_train <- min(max_training_pixels, n_valid)
  training_idx <- sample.int(n_valid, n_train)
  training <- X[training_idx, , drop = FALSE]
  scaler <- make_scaler(training)
  training_scaled <- scale_features(training, scaler)
  X_scaled <- scale_features(X, scaler)
  som_grid <- kohonen::somgrid(
    xdim = som_grid_x, ydim = som_grid_y, topo = som_grid_topology
  )
  som_model <- kohonen::som(
    training_scaled,
    grid = som_grid,
    rlen = som_training_length,
    alpha = som_learning_rate,
    keep.data = TRUE
  )
  codebook <- as.matrix(som_model$codes[[1L]])
  colnames(codebook) <- colnames(X)
  unit <- as.integer(kohonen::map(som_model, newdata = X_scaled)$unit.classif)
  if (length(unit) != n_valid || anyNA(unit)) {
    stop("SOM unit assignment failed in year ", year)
  }

  # Original K=6 model, before any validation uses additional random numbers.
  original_k6 <- stats::kmeans(
    codebook, centers = selected_k, nstart = original_kmeans_nstart
  )
  original_labels <- original_k6$cluster[unit]
  saveRDS(list(
    year = year, som_model = som_model,
    kmeans_model = original_k6,
    scaler = scaler,
    features = colnames(X),
    seed = random_seed
  ), file.path(year_dir, sprintf("SOM_K6_Model_%d.rds", year)))
  save_class_raster(
    s, indices, original_labels,
    file.path(year_dir, sprintf("SOM_K6_Raw_%d.tif", year)), "raw_cluster"
  )
  write.csv(response_profiles(X, original_labels, year, selected_k, "original_k6"),
            file.path(year_dir, sprintf("Original_K6_Response_Profiles_%d.csv", year)),
            row.names = FALSE)

  comparison <- compare_reference(year, s, indices, original_labels, year_dir)
  if (comparison$reference_found) {
    old_ids <- comparison$new_to_original[original_labels]
    crosswalk <- mapping[mapping$year == year, , drop = FALSE]
    map_old_to_paper <- integer(selected_k)
    map_old_to_paper[crosswalk$original_cluster] <- crosswalk$paper_type
    paper_type <- map_old_to_paper[old_ids]
    if (!identical(sort(unique(paper_type)), as.integer(seq_len(selected_k)))) {
      stop("Pathway type mapping failed in ", year)
    }
    if (comparison$best_match_fraction < 1 - 1e-10) {
      stop("Reproduced cluster raster differs from reference in year ", year,
           ". Scientific type export stopped; review the model before publication")
    }
    map_name <- if (semantic_mapping_confirmed) {
      sprintf("SOM_Paper_Type_%d.tif", year)
    } else {
      sprintf("SOM_Provisional_Type_UNVERIFIED_%d.tif", year)
    }
    save_class_raster(s, indices, paper_type,
                      file.path(year_dir, map_name), "paper_type")
    write.csv(response_profiles(X, paper_type, year, selected_k, "paper_type_mapping"),
              file.path(year_dir, sprintf("Paper_Type_Response_Profiles_%d.csv", year)),
              row.names = FALSE)
  }

  distance <- stats::dist(codebook)
  candidate_models <- list()
  candidate_scores <- vector("list", length(candidate_k))
  for (i in seq_along(candidate_k)) {
    k <- candidate_k[i]
    set.seed(random_seed + year * 100L + k)
    km <- stats::kmeans(codebook, centers = k,
                        nstart = validation_kmeans_nstart, iter.max = 200L)
    candidate_models[[as.character(k)]] <- km
    pixel_labels <- km$cluster[unit]
    silhouette <- cluster::silhouette(km$cluster, distance)
    candidate_scores[[i]] <- data.frame(
      year = year, k = k,
      silhouette = mean(silhouette[, "sil_width"]),
      withinss = km$tot.withinss,
      between_total_ratio = km$betweenss / km$totss,
      min_units_per_cluster = min(tabulate(km$cluster, nbins = k)),
      min_pixels_per_cluster = min(tabulate(pixel_labels, nbins = k)),
      min_pixel_share_pct = 100 * min(tabulate(pixel_labels, nbins = k)) / n_valid
    )
    if (k %in% c(5L, 6L)) {
      write.csv(response_profiles(X, pixel_labels, year, k, "candidate_kmeans"),
                file.path(year_dir, sprintf("K%d_Response_Profiles_%d.csv", k, year)),
                row.names = FALSE)
    }
  }
  scores <- do.call(rbind, candidate_scores)
  write.csv(scores, file.path(year_dir, sprintf("Candidate_K_%d.csv", year)),
            row.names = FALSE)
  if (5L %in% candidate_k && 6L %in% candidate_k) {
    labels5 <- candidate_models[["5"]]$cluster[unit]
    labels6 <- candidate_models[["6"]]$cluster[unit]
    overlap <- as.data.frame(table(K5 = labels5, K6 = labels6))
    names(overlap)[3L] <- "n_pixels"
    overlap$year <- year
    write.csv(overlap, file.path(year_dir, sprintf("K5_vs_K6_Overlap_%d.csv", year)),
              row.names = FALSE)
  }

  stability <- do.call(rbind, lapply(seq_len(kmeans_stability_repeats), function(i) {
    set.seed(random_seed + 10000L + year * 100L + i)
    trial <- stats::kmeans(codebook, centers = selected_k,
                           nstart = 10L, iter.max = 200L)
    data.frame(
      year = year, run = i,
      ari_codebook = adjusted_rand_index(original_k6$cluster, trial$cluster),
      ari_pixels = adjusted_rand_index(original_labels, trial$cluster[unit]),
      withinss = trial$tot.withinss
    )
  }))
  write.csv(stability,
            file.path(year_dir, sprintf("Kmeans_Initialization_Stability_%d.csv", year)),
            row.names = FALSE)

  png(file.path(year_dir, sprintf("Candidate_K_Diagnostics_%d.png", year)),
      width = 1400, height = 650, res = 150)
  oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  plot(scores$k, scores$silhouette, type = "b", pch = 19,
       xlab = "Number of clusters", ylab = "Mean silhouette width",
       main = paste(year, "SOM codebook"))
  abline(v = selected_k, lty = 2)
  plot(scores$k, scores$withinss, type = "b", pch = 19,
       xlab = "Number of clusters", ylab = "Within-cluster sum of squares")
  abline(v = selected_k, lty = 2)
  par(oldpar)
  dev.off()

  chosen <- scores[scores$k == selected_k, , drop = FALSE]
  summary <- data.frame(
    year = year,
    valid_pixels = n_valid,
    best_silhouette_k = scores$k[which.max(scores$silhouette)],
    selected_k = selected_k,
    selected_silhouette = chosen$silhouette,
    selected_withinss = chosen$withinss,
    mean_initialization_ari = mean(stability$ari_pixels),
    original_reference_found = comparison$reference_found,
    rerun_reference_ari = comparison$ari,
    rerun_reference_match = comparison$best_match_fraction
  )
  write.csv(summary,
            file.path(year_dir, sprintf("Validation_Summary_%d.csv", year)),
            row.names = FALSE)
  message("Completed ", year, ": ", n_valid, " valid cells, ",
          "K=6 silhouette=", round(chosen$silhouette, 4))
  list(summary = summary, scores = scores, audit = audit)
}

results <- lapply(years, run_year)
combined_summary <- do.call(rbind, lapply(results, function(z) z$summary))
combined_scores <- do.call(rbind, lapply(results, function(z) z$scores))
combined_audit <- do.call(rbind, lapply(results, function(z) z$audit))
write.csv(combined_summary, file.path(output_dir, "Validation_Summary_All_Years.csv"),
          row.names = FALSE)
write.csv(combined_scores, file.path(output_dir, "Candidate_K_All_Years.csv"),
          row.names = FALSE)
write.csv(combined_audit, file.path(output_dir, "Raster_Audit_All_Years.csv"),
          row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "R_Session_Info.txt"))
notes <- c(
  "Original year-specific SOM + K-means classification; three input coefficient rasters per year.",
  "The SOM is retrained separately for each year and uses a distinct annual Z-score scaler.",
  "The raw K-means cluster IDs are NOT directly comparable across years.",
  "The reported paper Type mapping is provisional unless authors verify crosswalk and downstream figure calculations.",
  paste0("semantic_mapping_confirmed: ", semantic_mapping_confirmed),
  "Candidate K comparisons assess sensitivity; K=6 is not the maximum-silhouette solution in all years.",
  "K-means initialization ARI is computed on a fixed trained SOM, not across full SOM retraining.",
  "Raster_Audit_All_Years.csv reports the actual SOM input grid size (around 2516 m), not 1 km.",
  "Upstream regression-coefficient raster creation, spatial fragmentation, MNL, enrichment and urban-rural analyses are not part of this script.",
  paste0("R version: ", R.version.string),
  paste0("terra: ", as.character(packageVersion("terra"))),
  paste0("kohonen: ", as.character(packageVersion("kohonen"))),
  paste0("cluster: ", as.character(packageVersion("cluster")))
)
writeLines(notes, file.path(output_dir, "REPRODUCIBILITY_NOTES.txt"))
message("Completed. Outputs saved to: ", normalizePath(output_dir, winslash = "/"))
