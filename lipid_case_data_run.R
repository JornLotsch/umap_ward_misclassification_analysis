# =============================================================================
# UMAP-Ward Misclassification Analysis Script
# Lipid Profiles (PsA vs Controls)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. SETUP: Source helper functions and set working directory
# -----------------------------------------------------------------------------

# Helper function to source all required analysis functions from current directory
source_required_functions <- function() {
  function_files <- c("umap_ward_misclassification_analysis.R")
  missing <- function_files[!file.exists(function_files)]
  if (length(missing) > 0) {
    stop("Required function files not found: ", paste(missing, collapse = ", "))
  }
  lapply(function_files, function(file) {
    source(file)
    message("Sourced function from ", file)
  })
}

# Set working directory and source functions
setwd("/home/joern/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/projectionsPlots/umap_ward_misclassification_analysis/umap_ward_misclassification_analysis")
source_required_functions()

# Set seed
SEED <- 42

# -----------------------------------------------------------------------------
# 2. REAL DATA ANALYSIS
# -----------------------------------------------------------------------------

# Load lipid profiles (rows = samples, columns = lipid features)
lipid_profiles <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished/PsA_lipids_BC.csv", row.names = 1)

# Load sample metadata (class labels)
lipid_metadata <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished/PsA_classes.csv")
sample_types <- lipid_metadata$PsA

# Run UMAP + Ward clustering analysis (two variants)
set.seed(SEED)
results_lipids <- umap_ward_misclassification_analysis(
  data = lipid_profiles, # Feature matrix
  target = sample_types, # Ground truth classes
  labels = lipid_metadata$ID, # Sample labels for plots
  output_dir = "lipid_results",
  determine_cluster_number = FALSE
)

set.seed(SEED)
results_lipids_2 <- umap_ward_misclassification_analysis(
  data = lipid_profiles,
  target = sample_types,
  labels = lipid_metadata$ID,
  output_dir = "lipid_results",
  determine_cluster_number = FALSE,
  voronoi_targets = FALSE
)

# Report misclassification results
cat("PsA/controls misclassification rate:",
    sprintf("%.2f%%", results_lipids$misclassification_rate * 100), "\n")

if (!is.null(results_lipids$misclassified_samples) && nrow(results_lipids$misclassified_samples) > 0) {
  cat("Misclassified PsA/controls:\n")
  print(results_lipids$misclassified_samples)
} else {
  cat("No misclassified PsA/controls.\n")
}

print(results_lipids$umap_plot)

# Combine and save real data plots
lipids_combined_plot_1 <- cowplot::plot_grid(
  results_lipids$voronoi_plot + labs(title = "UMAP projection of data by study group", subtitle = ""),
  results_lipids_2$voronoi_plot + labs(title = "UMAP projection of data by cluster", subtitle = ""),
  results_lipids$heatmap_result$plot +
    theme(legend.position.inside = TRUE, legend.position = c(.5, .5), legend.direction = "vertical") +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
    labs(title = "Groups and misclassifieds"),
  labels = "AUTO", nrow = 1, rel_widths = c(2, 2, 1),
  align = "h", axis = "tb"
)

print(lipids_combined_plot_1)
ggsave(plot = lipids_combined_plot_1, filename = "lipids_combined_plot_1.svg", width = 18, height = 8)

# -----------------------------------------------------------------------------
# 3. PERMUTED DATA ANALYSIS (Null model)
# -----------------------------------------------------------------------------

# Create permuted lipid profiles (columns shuffled within each feature)
SEED <- 42
lipid_profiles_rndm <- mapply(function(col, j) {
  set.seed(SEED + j)
  perm <- sample(seq_along(col), length(col))
  col[perm]
}, lipid_profiles, seq_len(ncol(lipid_profiles)), SIMPLIFY = FALSE)

lipid_profiles_rndm <- as.data.frame(lipid_profiles_rndm, stringsAsFactors = FALSE)

# Run analysis on permuted data
set.seed(SEED)
results_lipids_rndm <- umap_ward_misclassification_analysis(
  data = lipid_profiles_rndm,
  target = sample_types,
  labels = lipid_metadata$ID,
  output_dir = "lipid_results_rndm",
  determine_cluster_number = FALSE,
  voronoi_targets = FALSE
)

set.seed(SEED)
results_lipids_rndm_2 <- umap_ward_misclassification_analysis(
  data = lipid_profiles_rndm,
  target = sample_types,
  labels = lipid_metadata$ID,
  output_dir = "lipid_results_rndm",
  determine_cluster_number = TRUE,
  voronoi_targets = FALSE
)

# Report permuted data results
cat("PsA/controls misclassification rate in permuted data:",
    sprintf("%.2f%%", results_lipids_rndm$misclassification_rate * 100), "\n")

if (!is.null(results_lipids_rndm$misclassified_samples) && nrow(results_lipids_rndm$misclassified_samples) > 0) {
  cat("Misclassified PsA/controls in permuted data:\n")
  print(results_lipids_rndm$misclassified_samples)
} else {
  cat("No misclassified PsA/controls in permuted data.\n")
}

print(results_lipids_rndm$umap_plot)

# -----------------------------------------------------------------------------
# 4. HEATMAP FOR PERMUTED DATA UMAP PROJECTION
# -----------------------------------------------------------------------------

library(ComplexHeatmap)
pal0 <- colorRampPalette(c("cornsilk", "cornsilk2", "cornsilk3", "cornsilk4"))

# Prepare heatmap matrix from UMAP projection
hm_mat <- as.data.frame(results_lipids_rndm$umap_result$Projected)
rownames(hm_mat) <- results_lipids_rndm$umap_result$UniqueData$Label

# Bottom margin annotation (empty space)
ha_bottom <- ComplexHeatmap::HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  which = "column",
  height = grid::unit(1, "cm")
)

# Create and capture heatmap
heat_map_lipids_rdm_umap <- grid::grid.grabExpr({
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      as.matrix(hm_mat),
      col = pal0(123),
      cluster_rows = TRUE,
      clustering_method_rows = "ward.D2",
      show_row_dend = TRUE,
      cluster_columns = FALSE,
      row_dend_width = unit(6, "cm"),
      column_dend_height = unit(4, "cm"),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 6),
      bottom_annotation = ha_bottom,
      show_heatmap_legend = FALSE,
      border = FALSE,
      rect_gp = grid::gpar(col = NA),
      column_title = "Dendrogram\n",
      column_title_gp = grid::gpar(fontsize = 16, fontface = "plain")
    )
  )
})

# Combine and save permuted data plots
lipids_combined_plot_2 <- cowplot::plot_grid(
  results_lipids_rndm$voronoi_plot + labs(title = "UMAP projection of permuted data by study group", subtitle = ""),
  results_lipids_rndm_2$voronoi_plot + labs(title = "UMAP projection of permuted data by clusters", subtitle = ""),
  results_lipids_rndm_2$heatmap_result$plot +
    theme(legend.position.inside = TRUE, legend.position = c(.5, .5),
          legend.direction = "vertical", legend.text = element_text(size = 4.5), plot.subtitle = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
    labs(title = "Misclassifieds"),
  heat_map_lipids_rdm_umap,
  labels = "AUTO", nrow = 1, rel_widths = c(3, 3, 1, 2),
  align = "h", axis = "tb"
)

print(lipids_combined_plot_2)
ggsave(plot = lipids_combined_plot_2, filename = "lipids_combined_plot_2.svg", width = 18, height = 8)


# =============================================================================
# 5. Supervised Classification Analysis 
# =============================================================================

# Load required libraries
required_pkgs <- c("caret", "randomForest", "pbmcapply", "caTools", "parallel", "dplyr", "ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# Source the main function
source("perform_supervised_classification.R")

# Analysis parameters
controls_only <- TRUE
RF_only_plotted <- TRUE
training_size <- 0.67
n_iterations <- 100
# permute_data <- TRUE

supervised_results <- lapply(c(FALSE, TRUE), function(permute_data) {

  # -----------------------------------------------------------------------------
  # 5.1. READ DATA
  # -----------------------------------------------------------------------------

  lipid_profiles <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished/PsA_lipids_BC.csv", row.names = 1)
  lipid_metadata <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished/PsA_classes.csv", row.names = 1)

  if (controls_only) {
    control_lab_dates <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/09Originale/Controls_Lab_Dates.csv", row.names = 2, stringsAsFactors = FALSE)

    control_lab_dates <- control_lab_dates %>%
      mutate(Sampling.date = as.Date(Sampling.date, format = "%m/%d/%Y"),
             month = format(Sampling.date, "%m"))

    control_ids <- rownames(control_lab_dates)
    lipid_profiles <- lipid_profiles[rownames(lipid_profiles) %in% control_ids,]
    lipid_metadata <- lipid_metadata[rownames(lipid_metadata) %in% control_ids,]

    control_lab_dates$ID <- rownames(control_lab_dates)
    lipid_metadata$ID <- rownames(lipid_metadata)
    lipid_metadata <- lipid_metadata %>%
      left_join(control_lab_dates[, c("ID", "month")], by = "ID") %>%
      select(-ID)
  }

  # Permute data
  if (permute_data) {
    # Create permuted lipid profiles (columns shuffled within each feature)
    SEED <- 42
    lipid_profiles_rndm <- mapply(function(col, j) {
      set.seed(SEED + j)
      perm <- sample(seq_along(col), length(col))
      col[perm]
    }, lipid_profiles, seq_len(ncol(lipid_profiles)), SIMPLIFY = FALSE)

    lipid_profiles_rndm <- as.data.frame(lipid_profiles_rndm, stringsAsFactors = FALSE)
    lipid_profiles_for_classification <- lipid_profiles_rndm
  } else lipid_profiles_for_classification <- lipid_profiles

  # Run analysis
  message("Starting supervised classification analysis...")
  results <- perform_supervised_classification(
    X = lipid_profiles_for_classification,
    Y = lipid_metadata,
    seed = SEED,
    n_iter = n_iterations,
    training_size = training_size,
    skip_tuning = if (permute_data) TRUE else FALSE
  )

  # Save results
  write.csv(results$pooled_summary, "Supervised_Analysis_Summary_AllTargets_pooledBA.csv", row.names = FALSE)
  write.csv(results$final_summary, "Supervised_Analysis_Summary_AllTargets_detailed.csv", row.names = FALSE)
  write.csv(results$class_distributions, "Class_Distributions.csv", row.names = FALSE)

  message("Analysis completed — results saved")

  # -----------------------------------------------------------------------------
  # 5.2. PLOT RESULTS
  # -----------------------------------------------------------------------------

  plot_df <- results$pooled_summary |>
    filter(!is.na(Median)) |>
    mutate(
      Target = factor(Target),
      Model = factor(Model, levels = c("RandomForest", "SVM"))
    ) |>
    group_by(Target) |>
    mutate(
      max_median = max(Median),
  # CI_lower corresponding to that max Median; adjust if needed
      ci_at_max = CI_lower[which.max(Median)]
    ) |>
    ungroup() |>
    mutate(
  # main sort: max_median; tie‑breaker: ci_at_max
      sort_key = max_median * 1000 + ci_at_max,
      Target_ordered = reorder(Target, sort_key)
    )

  # If only RF is plotted 
  if (RF_only_plotted) plot_df <- plot_df[plot_df$Model == "RandomForest",]

  p <- ggplot(plot_df, aes(x = Target_ordered, y = Median, color = Model)) +
    geom_hline(yintercept = 0.5, color = "salmon", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
    geom_rect(aes(xmin = as.numeric(Target_ordered) - 0.4, xmax = as.numeric(Target_ordered) + 0.4,
                  ymin = CI_lower, ymax = CI_upper, fill = Model),
              color = "black", alpha = 0.7, linewidth = 0.1) +
    facet_wrap(Model ~ .) +
    geom_segment(aes(x = as.numeric(Target_ordered) - 0.35, xend = as.numeric(Target_ordered) + 0.35,
                     y = Median, yend = Median),
                 linewidth = 1.4, alpha = 1) +
    scale_color_manual(values = c("RandomForest" = "cornsilk4", "SVM" = "cornsilk4")) +
    scale_fill_manual(values = c("RandomForest" = alpha("cornsilk3", 0.2),
                                 "SVM" = alpha("cornsilk2", 0.2))) +
    coord_flip() +
    labs(x = "Target", y = "Balanced accuracy", title = "Balanced accuracy (controls only)",
         color = "Model", fill = "Model") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 7), legend.position = "top",
          strip.background = element_rect(fill = "cornsilk"),
          strip.text = element_text(colour = "black"))

  print(p)

  ggsave(paste0("Supervised_Analysis_BA_Plot_permuted", permute_data, ".svg"), p, width = if (RF_only_plotted) 6 else 12, height = 8, dpi = 300)

  message("Plot saved as Supervised_Analysis_BA_Plot.svg")

  return(list(
    results = results,
    plot = p))
})
message("Supervised analyses completed.")

