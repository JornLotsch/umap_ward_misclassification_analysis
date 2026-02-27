#===============================================================================
# UMAP-Ward Misclassification and Supervised Classification Analysis
# Lipid Profiles: Psoriatic Arthritis (PsA) vs Controls
# Reproducible Analysis Pipeline
#===============================================================================

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

setwd("/home/joern/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/projectionsPlots/umap_ward_misclassification_analysis/umap_ward_misclassification_analysis/")

# -----------------------------------------------------------------------------
# CONFIGURATION: File paths and parameters
# -----------------------------------------------------------------------------
data_dir <- "/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished"
lipid_file <- file.path(data_dir, "PsA_lipids_BC.csv")
metadata_file <- file.path(data_dir, "PsA_classes.csv")

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
lipid_profiles <- read.csv(lipid_file, row.names = 1, check.names = FALSE)
lipid_metadata <- read.csv(metadata_file, row.names = 1)

# Extract processing month from sampling dates
lipid_metadata$Month <- format(
  as.Date(lipid_metadata$Sampling_date, format = "%m/%d/%Y"), "%B"
)

message(sprintf("✓ Loaded %d samples × %d lipid features", 
                nrow(lipid_profiles), ncol(lipid_profiles)))

# -----------------------------------------------------------------------------
# 3. UMAP-WARD CLUSTERING ANALYSIS (REAL DATA)
# -----------------------------------------------------------------------------
message("\n=== UMAP-Ward Analysis: Real Data ===")

# Analysis 1: Voronoi colored by target groups
results_real_1 <- umap_ward_misclassification_analysis_multi(
  data = lipid_profiles,
  targets = lipid_metadata,
  labels = lipid_metadata$ID,
  output_dir = "lipid_results_real",
  determine_cluster_number = FALSE
)

# Analysis 2: Voronoi colored by clusters
results_real_2 <- umap_ward_misclassification_analysis_multi(
  data = lipid_profiles,
  targets = lipid_metadata,
  labels = lipid_metadata$ID,
  output_dir = "lipid_results_real",
  determine_cluster_number = FALSE,
  voronoi_targets = FALSE
)

# Report primary results (PsA classification)
cat(sprintf("\nPsA classification misclassification rate: %.1f%%\n", 
            results_real_1$target_results$PsA$misclassification_rate * 100))

# Create publication-quality figure panels
fig1a_psa <- cowplot::plot_grid(
  results_real_1$target_results$PsA$voronoi_plot + 
    labs(title = "UMAP: Study groups"),
  results_real_2$target_results$PsA$voronoi_plot + 
    labs(title = "UMAP: Clusters"),
  results_real_1$target_results$PsA$heatmap_result$plot + 
    theme(legend.position.inside = TRUE, legend.position = c(0.5, 0.5), 
          legend.direction = "vertical") +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
    labs(title = "Misclassifications"),
  labels = c("A", "B", "C"), nrow = 1, rel_widths = c(2, 2, 1),
  align = "h", axis = "tb"
)

fig1b_month <- cowplot::plot_grid(
  results_real_1$target_results$Month$voronoi_plot + 
    labs(title = "UMAP: Processing month"),
  results_real_2$target_results$Month$voronoi_plot + 
    labs(title = "UMAP: Clusters"),
  results_real_1$target_results$Month$heatmap_result$plot + 
    theme(legend.position.inside = TRUE, legend.position = c(0.5, 0.5), 
          legend.direction = "vertical") +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
    labs(title = "Misclassifications"),
  labels = c("D", "E", "F"), nrow = 1, rel_widths = c(2, 2, 1),
  align = "h", axis = "tb"
)

# Save figures
ggsave(fig1a_psa, filename = "Figure1_PsA_analysis.svg", width = 18, height = 8, dpi = 300)
ggsave(fig1b_month, filename = "Figure1_Month_analysis.svg", width = 18, height = 8, dpi = 300)
ggsave(cowplot::ggdraw(results_real_1$combined_heatmap_plot), filename = "Figure2_all_metadata_heatmaps.svg", width = 25, height = 20, dpi = 300)

# -----------------------------------------------------------------------------
# 4. NULL MODEL: PERMUTED DATA ANALYSIS
# -----------------------------------------------------------------------------
message("\n=== Null Model: Permuted Data Analysis ===")

# Permute lipid features (destroy biological signal)
lipid_profiles_perm <- as.data.frame(lapply(lipid_profiles, function(col, j) {
  set.seed(SEED + j)
  col[sample(length(col))]
}, j = seq_along(lipid_profiles)))

# Single-target permuted analysis
results_perm_1 <- umap_ward_misclassification_analysis(
  data = lipid_profiles_perm,
  target = lipid_metadata$PsA,
  labels = lipid_metadata$ID,
  output_dir = "lipid_results_perm",
  determine_cluster_number = FALSE,
  voronoi_targets = FALSE
)

results_perm_2 <- umap_ward_misclassification_analysis(
  data = lipid_profiles_perm,
  target = lipid_metadata$PsA,
  labels = lipid_metadata$ID,
  output_dir = "lipid_results_perm",
  determine_cluster_number = TRUE,
  voronoi_targets = FALSE
)

cat(sprintf("Permuted data misclassification rate: %.1f%%\n", 
            results_perm_1$misclassification_rate * 100))

# Permuted UMAP heatmap
hm_mat <- as.data.frame(results_perm_1$umap_result$Projected)
rownames(hm_mat) <- results_perm_1$umap_result$UniqueData$Label

pal0 <- colorRampPalette(c("cornsilk", "cornsilk2", "cornsilk3", "cornsilk4"))
ha_bottom <- ComplexHeatmap::HeatmapAnnotation(empty = anno_empty(border = FALSE),
                                               which = "column", height = unit(1, "cm"))

perm_umap_heatmap <- grid.grabExpr({
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
    as.matrix(hm_mat), col = pal0(123),
    cluster_rows = TRUE, clustering_method_rows = "ward.D2",
    show_row_dend = TRUE, cluster_columns = FALSE,
    row_dend_width = unit(6, "cm"), column_dend_height = unit(4, "cm"),
    row_names_side = "left", row_names_gp = gpar(fontsize = 6),
    bottom_annotation = ha_bottom, show_heatmap_legend = FALSE,
    border = FALSE, rect_gp = gpar(col = NA),
    column_title = "UMAP Dendrogram"
  ))
})

# Combined permuted data figure
fig2 <- cowplot::plot_grid(
  results_perm_1$voronoi_plot + labs(title = "Permuted: Study groups"),
  results_perm_2$voronoi_plot + labs(title = "Permuted: Clusters"),
  results_perm_2$heatmap_result$plot + 
    theme(legend.position.inside = TRUE, legend.position = c(0.5, 0.5),
          legend.direction = "vertical", legend.text = element_text(size = 4.5)) +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
    labs(title = "Misclassifications"),
  perm_umap_heatmap,
  labels = c("A", "B", "C", "D"), nrow = 1, 
  rel_widths = c(3, 3, 1, 2), align = "h", axis = "tb"
)

ggsave(fig2, filename = "Figure2_Permuted_analysis.svg", width = 22, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 5. SUPERVISED CLASSIFICATION ANALYSIS
# -----------------------------------------------------------------------------
message("\n=== Supervised Classification Analysis ===")

# Analysis parameters
SEED <- 42
training_size <- 0.67
n_iterations <- 100
RF_only_plotted <- TRUE

supervised_results <- lapply(c(FALSE, TRUE), function(permute_data) {
  lapply(c(FALSE, TRUE), function(controls_only) {
    
    # -----------------------------------------------------------------------------
    # 5.1. DATA PREPARATION
    # -----------------------------------------------------------------------------
    message(sprintf("  Preparing data (permute=%s, controls_only=%s)...", permute_data, controls_only))
    
    # Load data
    lipid_profiles <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished/PsA_lipids_BC.csv", row.names = 1)
    lipid_metadata <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished/PsA_classes.csv", row.names = 1)
    
    # Extract processing month
    lipid_metadata$Month <- format(
      as.Date(lipid_metadata$Sampling_date, format = "%m/%d/%Y"), "%B"
    )
    
    # Filter to controls only (if enabled)
    if (controls_only) {
      control_lab_dates <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/09Originale/Controls_Lab_Dates.csv", row.names = 2, stringsAsFactors = FALSE)
      control_ids <- rownames(control_lab_dates)
      lipid_profiles <- lipid_profiles[rownames(lipid_profiles) %in% control_ids, ]
      lipid_metadata <- lipid_metadata[rownames(lipid_metadata) %in% control_ids, ]
      message(sprintf("    ✓ Filtered to %d control samples", nrow(lipid_profiles)))
    }
    
    # Permute features (null model)
    if (permute_data) {
      set.seed(SEED)
      lipid_profiles_perm <- mapply(function(col, j) {
        set.seed(SEED + j)
        col[sample(length(col))]
      }, as.data.frame(lipid_profiles), seq_along(lipid_profiles), SIMPLIFY = FALSE)
      lipid_profiles_for_classification <- as.data.frame(lipid_profiles_perm, stringsAsFactors = FALSE)
    } else {
      lipid_profiles_for_classification <- lipid_profiles
    }
    
    # -----------------------------------------------------------------------------
    # 5.2. RUN CLASSIFICATION ANALYSIS
    # -----------------------------------------------------------------------------
    message(sprintf("  Running %d iterations...", n_iterations))
    results <- perform_supervised_classification(
      X = lipid_profiles_for_classification,
      Y = lipid_metadata,
      seed = SEED,
      n_iter = n_iterations,
      training_size = training_size,
      skip_tuning = permute_data
    )
    
    # -----------------------------------------------------------------------------
    # 5.2.1. CREATE STABILITY HEATMAPS (RF_ONLY + optional SVM)
    # -----------------------------------------------------------------------------
    prefix <- paste0(ifelse(permute_data, "permuted", "real"), "_", ifelse(controls_only, "controls", "all"))
    target_heatmaps <- list()
    
    for(target_name in names(results$rf_assignments)) {
      message(sprintf("  Processing target: %s", target_name))
      
      # **RF plot - ALWAYS**
      rf_heatmap <- plot_classification_stability_heatmap(
        class_assignments = results$rf_assignments[[target_name]],
        true_classes = lipid_metadata[[target_name]],
        title = paste(target_name, "(RF)"),
        row_font_size = 6
      )
      target_heatmaps[[paste0(target_name, "_RF")]] <- rf_heatmap$plot
      
      # Save individual RF
      ggsave(sprintf("Stability_Heatmap_%s_%s_RF.svg", prefix, target_name), 
             rf_heatmap$plot, width = 10, height = 8, dpi = 300)
      
      # **SVM plot - ONLY if RF_only_plotted = FALSE**
      if(!RF_only_plotted && target_name %in% names(results$svm_assignments)) {
        svm_heatmap <- plot_classification_stability_heatmap(
          class_assignments = results$svm_assignments[[target_name]],
          true_classes = lipid_metadata[[target_name]],
          title = paste(target_name, "(SVM)"),
          row_font_size = 6
        )
        target_heatmaps[[paste0(target_name, "_SVM")]] <- svm_heatmap$plot
        
        # Save individual SVM
        ggsave(sprintf("Stability_Heatmap_%s_%s_SVM.svg", prefix, target_name), 
               svm_heatmap$plot, width = 10, height = 8, dpi = 300)
      }
    }
    
    # -----------------------------------------------------------------------------
    # 5.2.2. COMBINED HEATMAP FIGURE
    # -----------------------------------------------------------------------------
    n_plots <- length(target_heatmaps)
    ncol_plot <- ifelse(RF_only_plotted, 3, 4)
    
    combined_heatmap <- cowplot::plot_grid(
      plotlist = target_heatmaps,
      ncol = ncol_plot,
      nrow = ceiling(n_plots / ncol_plot),
      label_size = 14,
      align = "hv"
    )
    
    print(combined_heatmap)
    
    ggsave(sprintf("Figure3_Supervised_Heatmaps_%s.svg", prefix), 
           combined_heatmap, width = ifelse(RF_only_plotted, 20, 25), height = 20, dpi = 300)
    message(sprintf("    ✓ Combined heatmap saved: Figure3_Supervised_Heatmaps_%s.svg", prefix))
    
    # -----------------------------------------------------------------------------
    # 5.3. CREATE PERFORMANCE PLOT
    # -----------------------------------------------------------------------------
    plot_df <- results$pooled_summary |>
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
    
    if (RF_only_plotted) {
      plot_df <- plot_df[plot_df$Model == "RandomForest", ]
    }
    
    p <- ggplot(plot_df, aes(x = Target_ordered, y = Median)) +
      geom_hline(yintercept = 0.5, color = "salmon", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
      geom_text(aes(y = -0.05, label = n_label), hjust = 1, size = 3, fontface = "bold", color = "black") +
      geom_rect(aes(xmin = as.numeric(Target_ordered) - 0.4, xmax = as.numeric(Target_ordered) + 0.4,
                    ymin = CI_lower, ymax = CI_upper), 
                fill = "cornsilk3", color = "black", alpha = 0.7, linewidth = 0.1) +
      facet_wrap(Model ~ ., scales = "free_y") +
      geom_segment(aes(x = as.numeric(Target_ordered) - 0.35, xend = as.numeric(Target_ordered) + 0.35,
                       y = Median, yend = Median), linewidth = 1.4, alpha = 1) +
      coord_flip(ylim = c(-0.12, 1.05)) +
      labs(
        x = "Metadata target", 
        y = "Balanced accuracy",
        title = sprintf("Classification Performance (%s, %s, n=%d)", 
                        ifelse(permute_data, "permuted", "real"), 
                        ifelse(controls_only, "controls", "all cases"), 
                        nrow(lipid_profiles_for_classification))
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
    
    # Save results
    write.csv(results$pooled_summary, sprintf("Supervised_%s_pooled_summary.csv", prefix), row.names = FALSE)
    write.csv(results$final_summary, sprintf("Supervised_%s_detailed_summary.csv", prefix), row.names = FALSE)
    write.csv(results$class_distributions, sprintf("Class_Distributions_%s.csv", prefix), row.names = FALSE)
    
    message(sprintf("    ✓ Results saved (%s)", prefix))
    
    return(list(results = results, ba_plot = p, heatmaps = target_heatmaps, prefix = prefix, combined_heatmap = combined_heatmap))
  })
})

# Final summary
message("\n✓ Supervised classification analysis complete (4 combinations)")
message("Generated files:")
for(res in supervised_results) {
  for(inner_res in res) {
    message(sprintf("- Supervised_%s_pooled_summary.csv", inner_res$prefix))
    message(sprintf("- Supervised_BA_Plot_%s.svg", inner_res$prefix))
    message(sprintf("- Figure3_Supervised_Heatmaps_%s.svg", inner_res$prefix))
  }
}
