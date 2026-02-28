#===============================================================================
# UMAP-Ward Misclassification and Supervised Classification Analysis
# Lipid Profiles: VALIDATION DATASET
# Reproducible Analysis Pipeline
#
# This script performs:
# 1. Unsupervised clustering (UMAP + Ward) with/without artificial bias
# 2. Supervised classification (RF/SVM) with/without permutation and bias
#===============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP: Libraries, Seeds, and Global Parameters
# -----------------------------------------------------------------------------

# Required libraries (installed automatically if missing)
required_packages <- c("ggplot2", "cowplot", "ComplexHeatmap", "grid", "umap", "deldir", "gridExtra")
for(pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)
  }
}

# Set global seed for reproducibility
SEED <- 42
set.seed(SEED)

# Global bias parameters (percentage increase for PLATE 1 samples)
BIAS_AMOUNT_PERCENT_UNSUPERVISED <- 95
BIAS_AMOUNT_PERCENT_SUPERVISED <- 95

# Set working directory
setwd("/home/joern/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/projectionsPlots/umap_ward_misclassification_analysis/umap_ward_misclassification_analysis/")

# -----------------------------------------------------------------------------
# CONFIGURATION: File paths and parameters
# -----------------------------------------------------------------------------
data_dir <- "/home/joern/Aktuell/RheumaMetabolomicsFFM/09Originale/ScienceDirect_files_18Feb2026_07-11-32.148/1-s2.0-S0003267015002962-mmc2/"
lipid2_file <- file.path(data_dir, "2groups.csv")

# -----------------------------------------------------------------------------
# 1. UTILITY FUNCTIONS
# -----------------------------------------------------------------------------
source_required_functions <- function() {
  function_files <- c(
    "umap_ward_misclassification_analysis.R",
    "umap_ward_misclassification_analysis_multi.R",
    "perform_supervised_classification.R",
    "plot_classification_stability_heatmap.R"
  )
  missing_files <- function_files[!file.exists(function_files)]
  if(length(missing_files) > 0) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
  }
  invisible(lapply(function_files, function(f) {
    source(f)
    message("✓ Sourced: ", f)
  }))
}

# Load required analysis functions
source_required_functions()

# -----------------------------------------------------------------------------
# 2. DATA LOADING AND PREPROCESSING
# -----------------------------------------------------------------------------
message("Loading lipidomics data...")
lipid2_profiles_raw <- read.csv(lipid2_file, row.names = 1, check.names = FALSE)

# Rotate matrix and separate metadata

lipid2_profiles <- as.data.frame(t(lipid2_profiles_raw[3:nrow(lipid2_profiles_raw),]))
names(lipid2_profiles) <- paste0("Lipid_", 1:ncol(lipid2_profiles))
lipid2_profiles <- sapply(lipid2_profiles, as.numeric)
lipid2_metadata <- as.data.frame(t(lipid2_profiles_raw[1:2,]))
rownames(lipid2_profiles) <- rownames(lipid2_metadata)


message(sprintf("✓ Loaded %d samples × %d lipid features", 
                nrow(lipid2_profiles), ncol(lipid2_profiles)))

# -----------------------------------------------------------------------------
# 3. UNSUPERVISED ANALYSIS: UMAP-WARD CLUSTERING (VALIDATION DATA)
# -----------------------------------------------------------------------------
# This section performs unsupervised clustering analysis with two conditions:
# - Original data (no bias)
# - Biased data (PLATE 1 samples increased by BIAS_AMOUNT_PERCENT_UNSUPERVISED)
# -----------------------------------------------------------------------------

message("\n=== UMAP-Ward Analysis: Validation Data ===")

unsupervised_results_val <- lapply(c(FALSE, TRUE), function(introduce_bias) {

  # ---------------------------------------------------------------------------
  # 3.1. DATA PREPARATION: Introduce bias (if enabled)
  # ---------------------------------------------------------------------------
  message(sprintf("  Processing unsupervised analysis (introduce_bias=%s)...", introduce_bias))

  if (introduce_bias) {
    rows_plate1 <- lipid2_metadata$Plate == "PLATE 1"
    lipid2_profiles_bias <- lipid2_profiles
    lipid2_profiles_bias[rows_plate1, ] <- lipid2_profiles_bias[rows_plate1, ] * (1 + BIAS_AMOUNT_PERCENT_UNSUPERVISED / 100)
    message(sprintf("    ✓ %s%% bias introduced for PLATE 1 samples", BIAS_AMOUNT_PERCENT_UNSUPERVISED))
  } else {
    lipid2_profiles_bias <- lipid2_profiles
    message("    ✓ Using original data (no bias)")
  }

  # ---------------------------------------------------------------------------
  # 3.2. RUN UMAP-WARD CLUSTERING ANALYSIS
  # ---------------------------------------------------------------------------
  # Performs UMAP dimensionality reduction followed by Ward clustering
  # Voronoi tessellation colored by target groups (Treatment, Plate)
  results_real_val_1 <- umap_ward_misclassification_analysis_multi(
    data = lipid2_profiles_bias,
    targets = lipid2_metadata,
    labels = lipid2_metadata$ID,
    output_dir = ifelse(introduce_bias, "lipid2_results_bias", "lipid2_results_real"),
    determine_cluster_number = FALSE,
    n_neighbors = 2
  )

  # ---------------------------------------------------------------------------
  # 3.3. REPORT MISCLASSIFICATION RATES
  # ---------------------------------------------------------------------------
  cat(sprintf("\n%s - Treatment misclassification rate: %.1f%%\n",
              ifelse(introduce_bias, "BIASED DATA", "ORIGINAL DATA"),
              results_real_val_1$target_results$Treatment$misclassification_rate * 100))

  # ---------------------------------------------------------------------------
  # 3.4. CREATE PUBLICATION-QUALITY FIGURE PANELS
  # ---------------------------------------------------------------------------
  # Create title plot
  title_unsupervised_string <- ifelse(introduce_bias,
                                      sprintf("Validation data: Unsupervised results - %s%% bias introduced", BIAS_AMOUNT_PERCENT_UNSUPERVISED),
                                      "Validation data: Unsupervised results - Original data")
  p_title_unsupervised <- ggplot() +
    labs(title = title_unsupervised_string) +
    theme_void() +
    theme(
      plot.title = element_text(face = "plain", size = 16, color = "#222222",
                                hjust = 0),
      plot.background = element_rect(fill = "white", color = NA)
    )

  # Create 4-panel row: Treatment (UMAP + heatmap) + Plate (UMAP + heatmap)
  panel_row <- cowplot::plot_grid(
    results_real_val_1$target_results$Treatment$voronoi_plot +
      labs(title = "UMAP: Treatment groups"),
    results_real_val_1$target_results$Treatment$heatmap_result$plot +
      theme(legend.position.inside = TRUE, legend.position = c(0.5, 0.5),
            legend.direction = "vertical") +
      guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
      labs(title = "Misclassifications - Treatment"),
    results_real_val_1$target_results$Plate$voronoi_plot +
      labs(title = "UMAP: Plates"),
    results_real_val_1$target_results$Plate$heatmap_result$plot +
      theme(legend.position.inside = TRUE, legend.position = c(0.5, 0.5),
            legend.direction = "vertical") +
      guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
      labs(title = "Misclassifications - Plate"),
    labels = c("A", "B", "C", "D"), nrow = 1, rel_widths = c(2, 1, 2, 1),
    align = "h", axis = "tb"
  )

  # Combine title and panels vertically
  fig1a_val <- cowplot::plot_grid(
    p_title_unsupervised,
    panel_row,
    ncol = 1, rel_heights = c(0.05, 1)
  )

  return(fig1a_val)
})

# ---------------------------------------------------------------------------
# 3.5. SAVE UNSUPERVISED ANALYSIS FIGURES
# ---------------------------------------------------------------------------
ggsave(unsupervised_results_val[[1]],
       filename = "Figure1_lipids_validation_analysis_unsupervised_original.svg",
       width = 22, height = 9, dpi = 300)
message("✓ Saved: Figure1_lipids_validation_analysis_unsupervised_original.svg")

ggsave(unsupervised_results_val[[2]],
       filename = "Figure1_lipids_validation_analysis_unsupervised_bias.svg",
       width = 22, height = 9, dpi = 300)
message("✓ Saved: Figure1_lipids_validation_analysis_unsupervised_bias.svg")

# -----------------------------------------------------------------------------
# 4. SUPERVISED CLASSIFICATION ANALYSIS
# -----------------------------------------------------------------------------
# This section performs supervised classification with four conditions:
# - Real data, original (no permutation, no bias)
# - Real data, biased (no permutation, PLATE 1 increased)
# - Permuted data, original (permuted features, no bias) - NULL MODEL
# - Permuted data, biased (permuted features, PLATE 1 increased) - NULL MODEL
# -----------------------------------------------------------------------------

message("\n=== Supervised Classification Analysis ===")

# Analysis parameters
training_size <- 0.67
n_iterations <- 100
RF_only_plotted <- TRUE  # Set to FALSE to include SVM results

supervised_results_val <- lapply(c(FALSE, TRUE), function(permute_data) {
  lapply(c(FALSE, TRUE), function(introduce_bias) {

    # -------------------------------------------------------------------------
    # 4.1. DATA PREPARATION
    # -------------------------------------------------------------------------
    message(sprintf("\n  Processing supervised analysis (permute=%s, introduce_bias=%s)...",
                    permute_data, introduce_bias))

    # Step 1: Introduce bias (if enabled)
    if (introduce_bias) {
      rows_plate1 <- lipid2_metadata$Plate == "PLATE 1"
      lipid2_profiles_bias <- lipid2_profiles
      lipid2_profiles_bias[rows_plate1, ] <- lipid2_profiles_bias[rows_plate1, ] * (1 + BIAS_AMOUNT_PERCENT_SUPERVISED / 100)
      message(sprintf("    ✓ %s%% bias introduced for PLATE 1 samples", BIAS_AMOUNT_PERCENT_SUPERVISED))
    } else {
      lipid2_profiles_bias <- lipid2_profiles
      message("    ✓ Using original data (no bias)")
    }

    # Step 2: Permute features (null model)
    if (permute_data) {
      set.seed(SEED)
      lipid2_profiles_perm <- mapply(function(col, j) {
        set.seed(SEED + j)
        col[sample(length(col))]
      }, as.data.frame(lipid2_profiles_bias), seq_along(lipid2_profiles_bias), SIMPLIFY = FALSE)
      lipid2_profiles_for_classification <- as.data.frame(lipid2_profiles_perm, stringsAsFactors = FALSE)
      message("    ✓ Features permuted (null model)")
    } else {
      lipid2_profiles_for_classification <- lipid2_profiles_bias
      message("    ✓ Features not permuted (real data)")
    }

    # -------------------------------------------------------------------------
    # 4.2. RUN CLASSIFICATION ANALYSIS
    # -------------------------------------------------------------------------
    message(sprintf("  Running %d iterations of supervised classification...", n_iterations))
    results_val <- perform_supervised_classification(
      X = lipid2_profiles_for_classification,
      Y = lipid2_metadata,
      seed = SEED,
      n_iter = n_iterations,
      training_size = training_size,
      skip_tuning = permute_data  # Skip hyperparameter tuning for permuted data
    )

    # -------------------------------------------------------------------------
    # 4.3. CREATE CLASSIFICATION STABILITY HEATMAPS
    # -------------------------------------------------------------------------
    # Generate file prefix based on conditions
    prefix <- paste0(ifelse(permute_data, "permuted", "real"), "_",
                     ifelse(introduce_bias, "bias", "original"))
    target_heatmaps <- list()

    for(target_name in names(results_val$rf_assignments)) {
      message(sprintf("  Creating heatmap for target: %s", target_name))

      # Random Forest heatmap (always generated)
      rf_heatmap <- plot_classification_stability_heatmap(
        class_assignments = results_val$rf_assignments[[target_name]],
        true_classes = lipid2_metadata[[target_name]],
        title = paste(target_name, "(RF)"),
        row_font_size = 6
      )
      target_heatmaps[[paste0(target_name, "_RF")]] <- rf_heatmap$plot

      # Save individual RF heatmap
      ggsave(sprintf("Stability_Heatmap_%s_%s_RF_val.svg", prefix, target_name),
             rf_heatmap$plot, width = 10, height = 8, dpi = 300)

      # SVM heatmap (only if RF_only_plotted = FALSE)
      if(!RF_only_plotted && target_name %in% names(results_val$svm_assignments)) {
        svm_heatmap <- plot_classification_stability_heatmap(
          class_assignments = results_val$svm_assignments[[target_name]],
          true_classes = lipid2_metadata[[target_name]],
          title = paste(target_name, "(SVM)"),
          row_font_size = 6
        )
        target_heatmaps[[paste0(target_name, "_SVM")]] <- svm_heatmap$plot

        # Save individual SVM heatmap
        ggsave(sprintf("Stability_Heatmap_%s_%s_SVM_val.svg", prefix, target_name),
               svm_heatmap$plot, width = 10, height = 8, dpi = 300)
      }
    }

    # -------------------------------------------------------------------------
    # 4.4. COMBINED HEATMAP FIGURE
    # -------------------------------------------------------------------------
    n_plots <- length(target_heatmaps)
    ncol_plot <- ifelse(RF_only_plotted, 3, 4)

    combined_heatmap_val <- cowplot::plot_grid(
      plotlist = target_heatmaps,
      ncol = ncol_plot,
      nrow = ceiling(n_plots / ncol_plot),
      label_size = 14,
      align = "hv"
    )

    print(combined_heatmap_val)

    ggsave(sprintf("Figure3_Supervised_Heatmaps_val_%s.svg", prefix),
           combined_heatmap_val,
           width = ifelse(RF_only_plotted, 20, 25),
           height = 20, dpi = 300)
    message(sprintf("  ✓ Saved: Figure3_Supervised_Heatmaps_val_%s.svg", prefix))

    # -------------------------------------------------------------------------
    # 4.5. CREATE BALANCED ACCURACY PERFORMANCE PLOT
    # -------------------------------------------------------------------------
    # Prepare data for plotting
    plot_df <- results_val$pooled_summary |>
      filter(!is.na(Median)) |>
      mutate(
        Target = factor(Target),
        Model = factor(Model, levels = c("RandomForest", "SVM"))
      ) |>
      group_by(Target) |>
      mutate(n_label = paste0("n=", first(N_complete))) |>
      ungroup() |>
      group_by(Target) |>
      mutate(
        max_median = max(Median),
        ci_at_max = CI_lower[which.max(Median)]
      ) |>
      ungroup() |>
      mutate(
        sort_key = max_median * 1000 + ci_at_max,
        Target_ordered = reorder(Target, sort_key)
      )

    # Filter to RF only if specified
    if (RF_only_plotted) {
      plot_df <- plot_df[plot_df$Model == "RandomForest", ]
    }

    # Create balanced accuracy plot
    p <- ggplot(plot_df, aes(x = Target_ordered, y = Median)) +
      geom_hline(yintercept = 0.5, color = "salmon", linetype = "dashed",
                 linewidth = 0.8, alpha = 0.7) +
      geom_text(aes(y = -0.05, label = n_label), hjust = 1, size = 3,
                fontface = "bold", color = "black") +
      geom_rect(aes(xmin = as.numeric(Target_ordered) - 0.4,
                    xmax = as.numeric(Target_ordered) + 0.4,
                    ymin = CI_lower, ymax = CI_upper),
                fill = "cornsilk3", color = "black", alpha = 0.7, linewidth = 0.1) +
      facet_wrap(Model ~ ., scales = "free_y") +
      geom_segment(aes(x = as.numeric(Target_ordered) - 0.35,
                       xend = as.numeric(Target_ordered) + 0.35,
                       y = Median, yend = Median), linewidth = 1.4, alpha = 1) +
      coord_flip(ylim = c(-0.12, 1.05)) +
      labs(
        x = "Metadata target",
        y = "Balanced accuracy",
        title = sprintf("Classification Performance (%s, %s, n=%d)",
                        ifelse(permute_data, "permuted", "real"),
                        ifelse(introduce_bias, "bias", "original"),
                        nrow(lipid2_profiles_for_classification))
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 7),
        legend.position = "top",
        strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black"),
        plot.margin = ggplot2::margin(l = 60, r = 10, t = 10, b = 10, unit = "pt")
      )

    print(p)
    ggsave(sprintf("Supervised_BA_Plot_%s.svg", prefix),
           plot = p, width = ifelse(RF_only_plotted, 6, 12), height = 8, dpi = 300)
    message(sprintf("  ✓ Saved: Supervised_BA_Plot_%s.svg", prefix))

    # -------------------------------------------------------------------------
    # 4.6. SAVE SUMMARY TABLES
    # -------------------------------------------------------------------------
    write.csv(results_val$pooled_summary,
              sprintf("Supervised_val_%s_pooled_summary.csv", prefix),
              row.names = FALSE)
    write.csv(results_val$final_summary,
              sprintf("Supervised_val_%s_detailed_summary.csv", prefix),
              row.names = FALSE)
    write.csv(results_val$class_distributions,
              sprintf("Class_Distributions_val_%s.csv", prefix),
              row.names = FALSE)
    message(sprintf("  ✓ Saved: Summary tables for %s condition", prefix))

    # Return results for this condition
    return(list(
      results = results_val,
      ba_plot = p,
      heatmaps = target_heatmaps,
      prefix = prefix,
      combined_heatmap = combined_heatmap_val
    ))
  })
})

# -----------------------------------------------------------------------------
# 5. CREATE COMBINED SUPERVISED ANALYSIS FIGURE (2×2 PANEL)
# -----------------------------------------------------------------------------
# Combines all 4 balanced accuracy plots into a single publication figure:
# - Top row: Real data (original vs. bias)
# - Bottom row: Permuted data (original vs. bias)
# -----------------------------------------------------------------------------

message("\n=== Creating combined supervised analysis figure ===")

# Create title plot
p_title_supervised <- ggplot() +
  labs(title = "Validation data: Supervised classification results") +
  theme_void() +
  theme(
    plot.title = element_text(face = "plain", size = 16, color = "#222222",
                              hjust = 0),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Create 2×2 panel grid
# Row 1: Real data (original, bias)
# Row 2: Permuted data (original, bias)
panel_row_supervised <- cowplot::plot_grid(
  supervised_results_val[[1]][[1]]$ba_plot,  # Real, original
  supervised_results_val[[1]][[2]]$ba_plot,  # Real, bias
  supervised_results_val[[2]][[1]]$ba_plot,  # Permuted, original
  supervised_results_val[[2]][[2]]$ba_plot,  # Permuted, bias
  labels = c("A", "B", "C", "D"),
  nrow = 2, ncol = 2,
  align = "h", axis = "tb"
)

# Combine title and panels vertically
fig1a_val_supervised <- cowplot::plot_grid(
  p_title_supervised,
  panel_row_supervised,
  ncol = 1, rel_heights = c(0.05, 1)
)

# Display and save combined figure
print(fig1a_val_supervised)
ggsave(fig1a_val_supervised,
       filename = "Figure1_lipids_validation_analysis_supervised.svg",
       width = 12, height = 9, dpi = 300)
message("✓ Saved: Figure1_lipids_validation_analysis_supervised.svg")

# -----------------------------------------------------------------------------
# 6. FINAL SUMMARY
# -----------------------------------------------------------------------------
message("\n" %+% strrep("=", 80))
message("ANALYSIS COMPLETE")
message(strrep("=", 80))

message("\n✓ Unsupervised analysis complete (2 conditions: original/bias)")
message("  - Figure1_lipids_validation_analysis_unsupervised_original.svg")
message("  - Figure1_lipids_validation_analysis_unsupervised_bias.svg")

message("\n✓ Supervised classification analysis complete (4 conditions: real/permuted × original/bias)")
message("\nGenerated files per condition:")
for(outer_res in supervised_results_val) {
  for(inner_res in outer_res) {
    message(sprintf("\n  [%s]", toupper(inner_res$prefix)))
    message(sprintf("    - Supervised_val_%s_pooled_summary.csv", inner_res$prefix))
    message(sprintf("    - Supervised_val_%s_detailed_summary.csv", inner_res$prefix))
    message(sprintf("    - Class_Distributions_val_%s.csv", inner_res$prefix))
    message(sprintf("    - Supervised_BA_Plot_%s.svg", inner_res$prefix))
    message(sprintf("    - Figure3_Supervised_Heatmaps_val_%s.svg", inner_res$prefix))
  }
}

message("\n✓ Combined supervised figure:")
message("  - Figure1_lipids_validation_analysis_supervised.svg")

message("\n" %+% strrep("=", 80))
message("All analyses completed successfully!")
message(strrep("=", 80) %+% "\n")
