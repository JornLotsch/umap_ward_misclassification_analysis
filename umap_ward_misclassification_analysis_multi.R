#' Perform UMAP projection and misclassification analysis for multiple target variables
#'
#' @description
#' Extends the single-target analysis to handle multiple target variables.
#' The UMAP projection is computed once, then clustering and visualization
#' are performed separately for each target variable. All heatmaps are combined
#' into a single multi-panel figure.
#'
#' @param data Data frame or matrix containing input data (samples as rows)
#' @param targets Data frame where each column is a different target variable.
#'                Column names will be used as target names.
#' @param labels Optional vector of labels for data points
#' @param determine_cluster_number Logical, whether to automatically determine cluster number
#' @param voronoi_targets Logical, whether Voronoi cells are colored for prior targets
#' @param output_dir Character, directory for saving plots
#' @param file_prefix Character, prefix for output filenames
#' @param file_format Character, "svg" or "png"
#' @param label_points Logical, whether to display point labels
#' @param row_font_size Numeric, font size for heatmap row labels
#' @param width Numeric, plot width in inches
#' @param height Numeric, plot height in inches
#' @param dpi Integer, resolution for PNG output
#' @param n_neighbors Integer, number of nearest neighbors for UMAP
#' @param heatmap_ncol Integer, number of columns in the combined heatmap layout.
#'   If NULL (default), it is chosen automatically to have more columns than rows.
#'   For example: 3 plots → 2 columns, 5 plots → 3 columns, etc.
#'
#' @return A list with:
#'   \item{umap_result}{UMAP projection (shared across all targets)}
#'   \item{target_results}{List of results for each target}
#'   \item{summary_stats}{Data frame with misclassification rates for all targets}
#'   \item{combined_heatmap_plot}{Combined multi-panel heatmap figure}
#'   \item{prepared_data}{Prepared input data}
#'
umap_ward_misclassification_analysis_multi <- function(data,
                                                       targets,
                                                       labels = NULL,
                                                       determine_cluster_number = FALSE,
                                                       voronoi_targets = TRUE,
                                                       output_dir = "results_multi",
                                                       file_prefix = "umap_analysis",
                                                       file_format = "svg",
                                                       label_points = TRUE,
                                                       row_font_size = 6,
                                                       width = 16,
                                                       height = 12,
                                                       dpi = 300,
                                                       n_neighbors = 15,
                                                       heatmap_ncol = NULL) {

  # --- Input Validation ---
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data frame or matrix.")
  }
  n_samples <- nrow(data)

  # Convert targets to data frame if it's a vector
  if (is.vector(targets) || is.factor(targets)) {
    targets <- data.frame(Target1 = targets)
  }

  if (!is.data.frame(targets)) {
    stop("'targets' must be a data frame, matrix, or vector.")
  }

  if (nrow(targets) != n_samples) {
    stop(sprintf("Number of rows in 'targets' (%d) must equal number of samples in data (%d).",
                 nrow(targets), n_samples))
  }

  # Use column names from targets data frame
  target_names <- colnames(targets)
  if (is.null(target_names)) {
    target_names <- paste0("Target", seq_len(ncol(targets)))
    colnames(targets) <- target_names
  }

  # --- Check and load required packages ---
  required_packages <- c("ggplot2", "tidyr", "scales", "deldir", "umap", "gridExtra")
  check_and_install_packages(required_packages)

  # --- Create output directory ---
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  }

  # --- Source required functions ---
  source_required_functions()

  # --- Prepare data (use first target for initial preparation) ---
  message("Preparing data...")
  prepared_data <- prepare_dataset(data, Target = targets[, 1], Label = labels)

  # --- Perform UMAP projection ONCE (independent of targets) ---
  message("Performing UMAP projection (shared across all targets)...")
  umap_result <- perform_umap_projection(data = prepared_data, seed = 42, scaleX = TRUE, n_neighbors = n_neighbors)

  # --- Store results for each target ---
  all_results <- list()
  all_heatmap_plots <- list()
  summary_stats <- data.frame(
    TargetName = character(),
    MisclassificationRate = numeric(),
    NumMisclassified = integer(),
    NumSamples = integer(),
    stringsAsFactors = FALSE
  )

  # --- Loop through each target variable ---
  for (i in seq_along(target_names)) {
    target_name <- target_names[i]
    current_target <- targets[, i]

    message(sprintf("\n=== Processing target: %s (%d/%d) ===", target_name, i, length(target_names)))

    # Create subdirectory for this target
    target_output_dir <- file.path(output_dir, target_name)
    if (!dir.exists(target_output_dir)) {
      dir.create(target_output_dir, recursive = TRUE)
    }

    # --- Perform clustering for this target ---
    message(sprintf("  Performing clustering for %s...", target_name))
    cluster_result <- perform_ward_clustering(
      projection_data = umap_result$Projected,
      target = current_target,
      determine_cluster_number = determine_cluster_number
    )

    # --- Combine data for visualization ---
    combined_data <- data.frame(
      umap_result$Projected,
      Target = current_target,
      Label = if (is.null(labels)) umap_result$UniqueData$Label else labels,
      Cluster = cluster_result$clusters
    )

    # --- Create Voronoi plot ---
    message(sprintf("  Creating Voronoi plot for %s...", target_name))
    voronoi_plot <- plot_umap_with_voronoi(
      umap_projection = combined_data,
      targets = combined_data$Target,
      clusters = if (voronoi_targets) NULL else combined_data$Cluster,
      labels = combined_data$Label,
      label_points = label_points)
    
    voronoi_plot <- voronoi_plot + labs(title = paste("UMAP Projection -", target_name))

    # --- Create misclassification heatmap ---
    message(sprintf("  Creating heatmap for %s...", target_name))
    heatmap_result <- plot_misclassification_heatmap(
      data = combined_data,
      prior_class_col = "Target",
      cluster_col = "Cluster",
      label_col = "Label",
      title = target_name,  # Use target name as title
      row_font_size = row_font_size
    )

    # Store the heatmap plot for combined figure
    all_heatmap_plots[[target_name]] <- heatmap_result$plot

    # --- Identify misclassified samples ---
    expected_class <- as.character(combined_data$Target)
    assigned_cluster <- as.character(combined_data$Cluster)
    misclassified_indices <- which(expected_class != assigned_cluster)
    misclassified_samples <- data.frame(
      Label = combined_data$Label[misclassified_indices],
      ExpectedClass = combined_data$Target[misclassified_indices],
      AssignedCluster = combined_data$Cluster[misclassified_indices]
    )

    # --- Create combined plot (Voronoi + Heatmap) ---
    message(sprintf("  Creating combined visualization for %s...", target_name))
    combined_plot <- gridExtra::grid.arrange(
      voronoi_plot,
      heatmap_result$plot,
      ncol = 2,
      widths = c(1.5, 1)
    )

    # --- Save individual plots ---
    file_extension <- file_format
    target_prefix <- paste0(file_prefix, "_", gsub("[^[:alnum:]_]", "_", target_name))
    voronoi_file <- file.path(target_output_dir, paste0(target_prefix, "_voronoi.", file_extension))
    heatmap_file <- file.path(target_output_dir, paste0(target_prefix, "_heatmap.", file_extension))
    combined_file <- file.path(target_output_dir, paste0(target_prefix, "_combined.", file_extension))

    message(sprintf("  Saving plots for %s...", target_name))
    ggplot2::ggsave(filename = voronoi_file, plot = voronoi_plot,
                    width = width * 0.6, height = height * 0.5, dpi = dpi)
    ggplot2::ggsave(filename = heatmap_file, plot = heatmap_result$plot,
                    width = width * 0.4, height = height * 0.5, dpi = dpi)

    if (file_format == "png") {
      png(combined_file, width = width * 0.5, height = height * 0.5, units = "in", res = dpi)
      grid::grid.draw(combined_plot)
      dev.off()
    } else {
      svg(combined_file, width = width * 0.5, height = height * 0.5)
      grid::grid.draw(combined_plot)
      dev.off()
    }

    # --- Save misclassified samples ---
    if (nrow(misclassified_samples) > 0) {
      misclassified_file <- file.path(target_output_dir, paste0(target_prefix, "_misclassified.csv"))
      write.csv(misclassified_samples, file = misclassified_file, row.names = FALSE)
    }

    # --- Store results ---
    all_results[[target_name]] <- list(
      target_name = target_name,
      cluster_result = cluster_result,
      voronoi_plot = voronoi_plot,
      heatmap_result = heatmap_result,
      combined_plot = combined_plot,
      combined_data = combined_data,
      misclassification_rate = heatmap_result$misclassification_rate,
      misclassified_samples = misclassified_samples
    )

    # --- Add to summary statistics ---
    summary_stats <- rbind(summary_stats, data.frame(
      TargetName = target_name,
      MisclassificationRate = heatmap_result$misclassification_rate,
      NumMisclassified = nrow(misclassified_samples),
      NumSamples = nrow(combined_data),
      stringsAsFactors = FALSE
    ))
  }

  # --- Create combined heatmap figure ---
  message("\n=== Creating combined heatmap figure ===")
  n_targets <- length(all_heatmap_plots)
  
  if (is.null(heatmap_ncol)) {
    heatmap_ncol <- ceiling(sqrt(n_targets))
    message(sprintf("Auto layout: %d plots → %d columns", n_targets, heatmap_ncol))
  } else {
    message(sprintf("Using specified layout: %d columns", heatmap_ncol))
  }
  
  # EXPLICIT nrow calculation + direct grid.arrange call (no do.call)
  nrow_layout <- ceiling(n_targets / heatmap_ncol)
  
  message(sprintf("Final layout: %dx%d grid (%d plots, %d empty cells)", 
                  nrow_layout, heatmap_ncol, n_targets, nrow_layout*heatmap_ncol - n_targets))
  
  combined_heatmap_plot <- gridExtra::grid.arrange(
    grobs = all_heatmap_plots,  # Use grobs= parameter
    ncol = heatmap_ncol,
    nrow = nrow_layout
  )
  
  
  # --- Save combined heatmap figure ---
  combined_heatmap_file <- file.path(output_dir, paste0(file_prefix, "_all_heatmaps.", file_format))
  message(sprintf("Saving combined heatmap figure to %s", combined_heatmap_file))

  if (file_format == "png") {
    png(combined_heatmap_file, width = width, height = height, units = "in", res = dpi)
    grid::grid.draw(combined_heatmap_plot)
    dev.off()
  } else {
    svg(combined_heatmap_file, width = width, height = height)
    grid::grid.draw(combined_heatmap_plot)
    dev.off()
  }

  # --- Save summary statistics ---
  summary_file <- file.path(output_dir, paste0(file_prefix, "_summary_all_targets.csv"))
  message(sprintf("Saving summary statistics to %s", summary_file))
  write.csv(summary_stats, file = summary_file, row.names = FALSE)

  # --- Print summary ---
  message("\n=== Analysis Complete ===")
  message("Summary of misclassification rates:")
  print(summary_stats)

  # --- Return results ---
  return(list(
    umap_result = umap_result,
    target_results = all_results,
    summary_stats = summary_stats,
    combined_heatmap_plot = combined_heatmap_plot,
    prepared_data = prepared_data
  ))
}

# Example usage
run_multi_examples <- FALSE

if (run_multi_examples) {
  # Load data
  lipid_profiles <- read.csv("lipid_profiles.csv")
  sample_metadata <- read.csv("sample_metadata.csv")

  # Extract sample labels (assuming there's a SampleID column)
  sample_labels <- sample_metadata$SampleID

  # Select target columns from metadata (excluding label columns)
  # Option 1: Select specific columns by name
  # target_cols <- c("SampleType", "BatchID", "TimePoint", "AgeGroup")
  # targets_df <- sample_metadata[, target_cols, drop = FALSE]

  # Option 2: Select all columns except SampleID (more flexible)
  targets_df <- sample_metadata[, !names(sample_metadata) %in% c("SampleID"), drop = FALSE]
  
  # Put some missings into two columns
  n <- nrow(targets_df)
  na_count <- max(1, round(0.10 * n))
  set.seed(123)
  na_indices <- sample(seq_len(n), na_count, replace = FALSE)
  targets_df[na_indices, "AgeGroup"] <- NA
  targets_df[na_indices, "SampleType"] <- NA
  
  # Run multi-target analysis
  set.seed(42)
  results_multi <- umap_ward_misclassification_analysis_multi(
    data = lipid_profiles,
    targets = targets_df,  # Multiple target variables from metadata
    labels = sample_labels,  # Separate labels parameter
    output_dir = "results_multi_target",
    file_prefix = "lipid_analysis",
    file_format = "svg",
    width = 16,
    height = 12,
    heatmap_ncol = NULL  # 2 columns in the combined heatmap layout
  )

  # View summary
  print(results_multi$summary_stats)

  # Access results for specific target (use actual column name from your data)
  if ("SampleType" %in% names(results_multi$target_results)) {
    sample_type_results <- results_multi$target_results[["SampleType"]]
    cat(sprintf("SampleType misclassification rate: %.2f%%\n",
                sample_type_results$misclassification_rate))
  }
}