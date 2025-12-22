# Helper function to check and install missing packages
check_and_install_packages <- function(pkg_list) {
  new_pkgs <- setdiff(pkg_list, rownames(installed.packages()))
  if (length(new_pkgs) > 0) {
    message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
    install.packages(new_pkgs, repos = "https://cloud.r-project.org/")
  }
  failed <- pkg_list[!sapply(pkg_list, require, character.only = TRUE, quietly = TRUE)]
  if (length(failed) > 0) {
    stop("Failed to load required packages: ", paste(failed, collapse = ", "))
  }
}

# Helper function to source all required analysis functions from the current directory
source_required_functions <- function() {
  function_files <- c(
    "prepare_dataset.R",
    "perform_umap_projection.R",
    "perform_ward_clustering.R",
    "plot_umap_with_voronoi.R",
    "plot_misclassification_heatmap.R"
  )
  missing <- function_files[!file.exists(function_files)]
  if (length(missing) > 0) {
    stop("Required function files not found: ", paste(missing, collapse = ", "))
  }
  lapply(function_files, function(file) {
    source(file)
    message("Sourced function from ", file)
  })
}


#' Perform UMAP projection and misclassification analysis
#'
#' @description
#' A comprehensive analysis workflow that integrates data preparation,
#' UMAP projection, clustering, and visualization. This function handles
#' the entire process from raw data to final visualizations and outputs
#' both a Voronoi UMAP plot and a misclassification heatmap.
#'
#' @param data Data frame or matrix containing input data (samples as rows)
#' @param target Optional vector of target class values for classification (length must match nrow(data))
#' @param labels Optional vector of labels for data points (length must match nrow(data))
#' @param determine_cluster_number Logical, whether to automatically determine cluster number and membership. 
#' @param                         if TRUE, then 
#' @param voronoi_targets Logical, whether Voronoi cells are colored for prior targets, else for clusters
#' @param output_dir Character, directory for saving plots (default: "results")
#' @param file_prefix Character, prefix for output filenames (default: "umap_analysis")
#' @param file_format Character, "svg" or "png" (default: "svg")
#' @param label_points Logical, whether to display point labels in plots (default: TRUE)
#' @param row_font_size Numeric, font size for heatmap row labels (default: 6)
#' @param width Numeric, plot width in inches (default: 12)
#' @param height Numeric, plot height in inches (default: 9)
#' @param dpi Integer, resolution for PNG output (default: 300)
#' @param n_neighbors Integer, number of nearest neighbors used in UMAP (default: 15)
#'
#' @return A list containing all analysis components:
#'   \item{prepared_data}{Data frame with prepared input data}
#'   \item{umap_result}{UMAP projection results}
#'   \item{cluster_result}{Clustering results}
#'   \item{voronoi_plot}{UMAP Voronoi plot}
#'   \item{heatmap_result}{Misclassification heatmap results}
#'   \item{combined_plot}{Combined visualization of both plots}
#'   \item{misclassification_rate}{Numeric value indicating the proportion of misclassified samples}
#'   \item{misclassified_samples}{Data frame containing labels, expected classes, and assigned clusters for misclassified samples}
#'
#' @importFrom ggplot2 ggsave
#' @importFrom gridExtra grid.arrange

umap_ward_misclassification_analysis <- function(data,
                                                 target = NULL,
                                                 labels = NULL,
                                                 determine_cluster_number = FALSE,
                                                 voronoi_targets = TRUE,
                                                 output_dir = "results",
                                                 file_prefix = "umap_analysis",
                                                 file_format = "svg",
                                                 label_points = TRUE,
                                                 row_font_size = 6,
                                                 width = 12,
                                                 height = 9,
                                                 dpi = 300,
                                                 n_neighbors = 15) {

  # --- Input Validation ---
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data frame or matrix.")
  }
  n_samples <- nrow(data)

  if (!is.null(target) && length(target) != n_samples) {
    stop(sprintf("Length of 'target' (%d) must equal number of samples in data (%d).", length(target), n_samples))
  }

  if (!is.null(labels) && length(labels) != n_samples) {
    stop(sprintf("Length of 'labels' (%d) must equal number of samples in data (%d).", length(labels), n_samples))
  }

  if (!file_format %in% c("svg", "png")) {
    stop("file_format must be 'svg' or 'png'.")
  }

  if (!is.numeric(width) || !is.numeric(height) || width <= 0 || height <= 0) {
    stop("Both 'width' and 'height' must be positive numbers.")
  }

  if (!is.numeric(dpi) || dpi <= 0) {
    stop("'dpi' must be a positive integer.")
  }

  # --- Check and load required packages ---
  required_packages <- c("ggplot2", "tidyr", "scales", "deldir", "umap", "gridExtra")
  check_and_install_packages(required_packages)

  # --- Create output directory if it doesn't exist ---
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  } else if (!file.access(output_dir, 2) == 0) {
    stop("No write permission in output directory: ", output_dir)
  }

  # --- Source required functions ---
  source_required_functions()

  # --- Prepare data ---
  message("Preparing data...")
  prepared_data <- prepare_dataset(data, Target = target, Label = labels)

  # --- Perform UMAP projection ---
  message("Performing UMAP projection...")
  umap_result <- perform_umap_projection(data = prepared_data, seed = 42, scaleX = TRUE, n_neighbors = n_neighbors)

  # --- Perform clustering ---
  message("Performing clustering...")
  cluster_result <- perform_ward_clustering(
    projection_data = umap_result$Projected,
    target = umap_result$UniqueData$Target,
    determine_cluster_number = determine_cluster_number
  )

  # --- Combine data for visualization ---
  combined_data <- data.frame(
    umap_result$Projected,
    Target = umap_result$UniqueData$Target,
    Label = umap_result$UniqueData$Label,
    Cluster = cluster_result$clusters
  )

  # --- Create Voronoi plot ---
  message("Creating UMAP Voronoi plot...")
  voronoi_plot <- plot_umap_with_voronoi(
    umap_projection = combined_data,
    targets = combined_data$Target,
    clusters = if (voronoi_targets) NULL else combined_data$Cluster,
    labels = combined_data$Label,
    label_points = label_points
  )

  # --- Create misclassification heatmap ---
  message("Creating misclassification heatmap...")
  heatmap_result <- plot_misclassification_heatmap(
    data = combined_data,
    prior_class_col = "Target",
    cluster_col = "Cluster",
    label_col = "Label",
    row_font_size = row_font_size
  )

  # --- Identify misclassified samples (using 1:1 assignment, not re-mapping) ---
  message("Identifying misclassified samples based on provided cluster labels...")
  expected_class <- as.character(combined_data$Target)
  assigned_cluster <- as.character(combined_data$Cluster)
  misclassified_indices <- which(expected_class != assigned_cluster)
  misclassified_samples <- data.frame(
    Label = combined_data$Label[misclassified_indices],
    ExpectedClass = combined_data$Target[misclassified_indices],
    AssignedCluster = combined_data$Cluster[misclassified_indices]
  )

  # --- Create combined plot ---
  message("Creating combined visualization...")
  combined_plot <- gridExtra::grid.arrange(
    voronoi_plot,
    heatmap_result$plot,
    ncol = 2,
    widths = c(1.5, 1)
  )

  # --- Save plots (each, and combined) ---
  file_extension <- file_format
  voronoi_file <- file.path(output_dir, paste0(file_prefix, "_voronoi.", file_extension))
  heatmap_file <- file.path(output_dir, paste0(file_prefix, "_heatmap.", file_extension))
  combined_file <- file.path(output_dir, paste0(file_prefix, "_combined.", file_extension))

  # Check file overwrite and save (warn if file exists)
  for (f in c(voronoi_file, heatmap_file, combined_file)) {
    if (file.exists(f)) warning("Overwriting existing file: ", f)
  }

  message(paste("Saving Voronoi plot to", voronoi_file))
  ggplot2::ggsave(
    filename = voronoi_file,
    plot = voronoi_plot,
    width = width * 0.6,
    height = height,
    dpi = dpi
  )

  message(paste("Saving heatmap to", heatmap_file))
  ggplot2::ggsave(
    filename = heatmap_file,
    plot = heatmap_result$plot,
    width = width * 0.4,
    height = height,
    dpi = dpi
  )

  message(paste("Saving combined plot to", combined_file))
  # Saving grid object (grob), use explicit device
  if (file_format == "png") {
    png(combined_file, width = width, height = height, units = "in", res = dpi)
    grid::grid.draw(combined_plot)
    dev.off()
  } else {
    svg(combined_file, width = width, height = height)
    grid::grid.draw(combined_plot)
    dev.off()
  }

  # --- Write misclassified samples to file ---
  if (nrow(misclassified_samples) > 0) {
    misclassified_file <- file.path(output_dir, paste0(file_prefix, "_misclassified_samples.csv"))
    message(paste("Writing misclassified samples to", misclassified_file))
    write.csv(misclassified_samples, file = misclassified_file, row.names = FALSE)
  }

  # --- Return results ---
  message("Analysis complete!")
  return(list(
    prepared_data = prepared_data,
    umap_result = umap_result,
    cluster_result = cluster_result,
    voronoi_plot = voronoi_plot,
    heatmap_result = heatmap_result,
    combined_plot = combined_plot,
    misclassification_rate = heatmap_result$misclassification_rate,
    misclassified_samples = misclassified_samples
  ))
}

# Example usage
run_examples <- FALSE

if (run_examples) {
  # Load your data frame: each row is a sample, each column a feature (e.g., lipid species)
  lipid_profiles <- read.csv("lipid_profiles.csv")

  # Load or extract your class labels for each sample
  sample_metadata <- read.csv("sample_metadata.csv")
  sample_types <- sample_metadata$SampleType

  # Run the integrated UMAP projection and clustering/misclassification analysis
  set.seed(42)
  results <- umap_ward_misclassification_analysis(
    data = lipid_profiles[, !names(lipid_profiles) %in% c("SampleID")], # Features data
    target = sample_types, # Ground truth (prior classes)
    labels = sample_metadata$SampleID, # Optional: row labels for plots, if available
    output_dir = "qc_results" # Directory to save QC visualizations and results
  )

  # Output misclassification rate and list which samples were misclassified
  cat("Sample misclassification rate:",
      sprintf("%.2f%%", results$misclassification_rate * 100), "\n")

  if (!is.null(results$misclassified_samples) && nrow(results$misclassified_samples) > 0) {
    cat("Misclassified samples:\n")
    print(results$misclassified_samples)
  } else {
    cat("No misclassified samples.\n")
  }

  # Optionally: view UMAP plot object, if provided
  print(results$umap_plot)

  # Run the integrated UMAP projection and clustering/misclassification analysis on randomly permuted data
  SEED <- 42
  lipid_profiles_rndm <- mapply(function(col, j) {
    set.seed(SEED + j)
    perm <- sample(seq_along(col), length(col))
    col[perm]
  }, lipid_profiles, seq_len(ncol(lipid_profiles)), SIMPLIFY = FALSE)

  lipid_profiles_rndm <- as.data.frame(lipid_profiles_rndm, stringsAsFactors = FALSE)

  set.seed(42)
  results_rndm <- umap_ward_misclassification_analysis(
    data = lipid_profiles_rndm[, !names(lipid_profiles) %in% c("SampleID")], # Features data
    target = sample_types, # Ground truth (prior classes)
    labels = sample_metadata$SampleID, # Optional: row labels for plots, if available
    output_dir = "qc_results_rndm" # Directory to save QC visualizations and results
  )

  # Output misclassification rate and list which samples were misclassified
  cat("Sample misclassification rate in permuted data:",
      sprintf("%.2f%%", results_rndm$misclassification_rate * 100), "\n")

  if (!is.null(results_rndm$misclassified_samples) && nrow(results_rndm$misclassified_samples) > 0) {
    cat("Misclassified samples in permuted data:\n")
    print(results_rndm$misclassified_samples)
  } else {
    cat("No misclassified samples in permuted data.\n")
  }

  # Optionally: view UMAP plot object of permuted data, if provided
  print(results_rndm$umap_plot)
}
