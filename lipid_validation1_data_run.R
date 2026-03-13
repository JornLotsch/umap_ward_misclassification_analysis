#===============================================================================
# UMAP-Ward Misclassification and Supervised Classification Analysis
# Lipid Profiles: VALIDATION DATASET (Batch Effects)
# Reproducible Analysis Pipeline - March 2026
#===============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP: Libraries, Seeds, and Global Parameters
# -----------------------------------------------------------------------------
required_packages <- c("ggplot2", "cowplot", "ComplexHeatmap", "grid", "umap", 
                       "deldir", "gridExtra", "readxl", "dplyr", "tools")
for(pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)
  }
}

SEED <- 42
set.seed(SEED)

# Set working directory
setwd("/home/joern/Aktuell/ProjectionsBiomed/08AnalyseProgramme/R/projectionsPlots/umap_ward_misclassification_analysis/umap_ward_misclassification_analysis/")

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
source_required_functions()

# -----------------------------------------------------------------------------
# 2. DATA LOADING AND BATCH EFFECTS SIMULATION
# -----------------------------------------------------------------------------
message("Loading lipidomics data and simulating batch effects...")

data_dir <- "/home/joern/Dokumente/BiomarkerLipide/90RawData/"
lipid_val1_file <- file.path(data_dir, "BiomarkerLipide.xlsx")

lipid_val1_file_raw <- readxl::read_excel(lipid_val1_file, skip = 1)
lipids_val1 <- lipid_val1_file_raw[lipid_val1_file_raw$Group == "Healthy_S",]
lipids_val1_profiles <- as.data.frame(apply(lipid_val1_file_raw[,16:ncol(lipid_val1_file_raw)], 2, as.numeric))

# Remove variables with NaNs
nan_vars <- apply(lipids_val1_profiles, 2, function(x) sum(is.na(x)))
lipids_val1_profiles_complete <- lipids_val1_profiles[,-which(nan_vars > 0)]
rownames(lipids_val1_profiles_complete) <- paste0("sample_", 1:nrow(lipids_val1_profiles_complete))

message(sprintf("✓ Loaded %d samples × %d lipid features", 
                nrow(lipids_val1_profiles_complete), ncol(lipids_val1_profiles_complete)))

# CREATE SHARED BATCH METADATA (2 batches, 50/50 split for ML identification)
batch_metadata <- data.frame(
  sample_id = rownames(lipids_val1_profiles_complete),
  batch_day = sample(1:2, nrow(lipids_val1_profiles_complete), replace = TRUE, prob = c(0.6, 0.4))
)
batch_metadata$run_order <- sample(1:nrow(batch_metadata), nrow(batch_metadata))

lipid_cols <- names(lipids_val1_profiles_complete)

# Function to add batch effects
add_batch_effects <- function(data, magnitude) {
  result <- data
  for(col in lipid_cols) {
    day_effect <- (batch_metadata$batch_day - 1.5) * magnitude * 0.7
    order_effect <- scale(batch_metadata$run_order) * magnitude * 0.3
    total_effect <- day_effect + order_effect
    result[[col]] <- data[[col]] * exp(total_effect)
  }
  return(result)
}

# Generate 5 datasets: Original + 4 batch levels
lipids_val1_profiles_complete_batch_weak <- add_batch_effects(lipids_val1_profiles_complete, 0.10)    # 10% CV
lipids_val1_profiles_complete_batch_medium <- add_batch_effects(lipids_val1_profiles_complete, 0.20)  # 20% CV
lipids_val1_profiles_complete_batch_strong <- add_batch_effects(lipids_val1_profiles_complete, 0.30)  # 30% CV
lipids_val1_profiles_complete_batch_extreme <- add_batch_effects(lipids_val1_profiles_complete, 1.00) # 100% CV

# -----------------------------------------------------------------------------
# 3. UNSUPERVISED ANALYSIS: UMAP-WARD CLUSTERING
# -----------------------------------------------------------------------------
message("\n=== UMAP-Ward Analysis: Batch Effect Detection ===")

datasets <- c("lipids_val1_profiles_complete", "lipids_val1_profiles_complete_batch_weak", 
              "lipids_val1_profiles_complete_batch_medium", "lipids_val1_profiles_complete_batch_strong",
              "lipids_val1_profiles_complete_batch_extreme")

batch_levels <- c("Original", "Weak (10% CV)", "Medium (20% CV)", "Strong (30% CV)", "Extreme (100% CV)")

unsupervised_results_val <- lapply(seq_along(datasets), function(i) {
  actual_dataset <- datasets[i]
  batch_level <- batch_levels[i]
  
  message(sprintf("  Processing %s dataset", batch_level))
  
  results_batch_val <- umap_ward_misclassification_analysis_multi(
    data = get(actual_dataset),
    targets = subset(batch_metadata, select = "batch_day"),
    labels = batch_metadata$sample_id,
    output_dir = "lipids_val1",
    determine_cluster_number = FALSE, 
    label_points = FALSE
  )
  
  cat(sprintf("\n%s - Batch misclassification rate: %.1f%%\n",
              batch_level, results_batch_val$target_results$batch_day$misclassification_rate * 100))
  
  # Create figure panels
  title_str <- sprintf("Unsupervised batch detection - %s", batch_level)
  p_title <- ggplot() + labs(title = title_str) + theme_void() +
    theme(plot.title = element_text(size = 16, hjust = 0), plot.background = element_rect(fill = "white"))
  
  panel_row <- cowplot::plot_grid(
    results_batch_val$target_results$batch_day$voronoi_plot + labs(title = "UMAP: Batch structure"),
    results_batch_val$target_results$batch_day$heatmap_result$plot + 
      theme(legend.position = "top", legend.direction = "vertical") +
      guides(fill = guide_legend(nrow = 2)) + labs(title = "Misclassifications"),
    nrow = 1, rel_widths = c(2, 1), align = "h"
  )
  
  fig <- cowplot::plot_grid(p_title, panel_row, ncol = 1, rel_heights = c(0.1, 1))
  ggsave(sprintf("Unsupervised_val1_Batch_Detection_%s.svg", gsub(" \\(.*\\)", "", batch_level)), 
         fig, width = 12, height = 6, dpi = 300)
  
  return(fig)
})

# Combined unsupervised figure (2x3 layout for 5 panels)
summary_unsupervised <- cowplot::plot_grid(plotlist = unsupervised_results_val, 
                                           ncol = 3, nrow = 2, labels = "AUTO")
ggsave("Unsupervised_val1_Batch_Detection_Summary.svg", summary_unsupervised, width = 36, height = 24, dpi = 300)
ggsave("Unsupervised_val1_Batch_Detection_Summary.png", summary_unsupervised, width = 36, height = 24, dpi = 300)

# -----------------------------------------------------------------------------
# 4. SUPERVISED CLASSIFICATION ANALYSIS
# -----------------------------------------------------------------------------
message("\n=== Supervised Batch Classification Analysis ===")

training_size <- 0.67
n_iterations <- 100
RF_only_plotted <- TRUE

supervised_results_val <- lapply(seq_along(datasets), function(i) {
  actual_dataset <- datasets[i]
  batch_level <- gsub("_", "-", tolower(gsub("lipids_val1_profiles_complete_batch_", "batch-", datasets[i])))
  if(batch_level == "lipids-val1-profiles-complete") batch_level <- "original"
  
  message(sprintf("\n  Processing supervised analysis (%s)...", batch_level))
  
  results_val <- perform_supervised_classification(
    X = get(actual_dataset),
    Y = batch_metadata["batch_day"],
    seed = SEED, n_iter = n_iterations, training_size = training_size, skip_tuning = FALSE
  )
  
  # Heatmap
  rf_heatmap <- plot_classification_stability_heatmap(
    results_val$rf_assignments$batch_day, batch_metadata$batch_day, "Batch Day (RF)", row_font_size = 6
  )
  ggsave(sprintf("Stability_Heatmap_val1_%s_batch_day_RF.svg", batch_level), rf_heatmap$plot, 
         width = 10, height = 8, dpi = 300)
  
  # BA plot
  plot_df <- results_val$pooled_summary |> filter(!is.na(Median)) |> 
    mutate(Target = factor(Target), Model = factor(Model, levels = c("RandomForest", "SVM"))) |>
    group_by(Target) |> mutate(n_label = paste0("n=", first(N_complete))) |> ungroup()
  
  if(RF_only_plotted) plot_df <- plot_df[plot_df$Model == "RandomForest", ]
  
  p <- ggplot(plot_df, aes(x = Target, y = Median)) +
    geom_hline(yintercept = 0.5, color = "salmon", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
    geom_text(aes(y = -0.05, label = n_label), hjust = 1, size = 3, fontface = "bold") +
    geom_rect(aes(xmin = as.numeric(Target)-0.4, xmax = as.numeric(Target)+0.4, ymin = CI_lower, ymax = CI_upper),
              fill = "cornsilk3", alpha = 0.7, color = "black", linewidth = 0.1) +
    geom_segment(aes(x = as.numeric(Target)-0.35, xend = as.numeric(Target)+0.35, y = Median, yend = Median), 
                 linewidth = 1.4) + coord_flip(ylim = c(-0.12, 1.05)) +
    labs(x = "Batch target", y = "Balanced accuracy", title = sprintf("Batch Classification (%s)", batch_level)) +
    theme_bw() + theme(axis.text.y = element_text(size = 7), legend.position = "top",
                       strip.background = element_rect(fill = "cornsilk"))
  
  ggsave(sprintf("Supervised_val1_BA_Plot_%s.svg", batch_level), p, width = 8, height = 8, dpi = 300)
  
  # Save tables
  write.csv(results_val$pooled_summary, sprintf("Supervised_val1_%s_pooled_summary.csv", batch_level), row.names = FALSE)
  write.csv(results_val$final_summary, sprintf("Supervised_val1_%s_detailed_summary.csv", batch_level), row.names = FALSE)
  
  return(list(ba_plot = p, prefix = batch_level))
})

# -----------------------------------------------------------------------------
# 5. COMBINED FIGURES
# -----------------------------------------------------------------------------
# Supervised combined (horizontal row)
p_title_sup <- ggplot() + labs(title = "Supervised batch classification") + theme_void() + 
  theme(plot.title = element_text(size = 16, hjust = 0), plot.background = element_rect(fill = "white"))
ba_plots <- lapply(supervised_results_val, function(x) x$ba_plot)
panel_sup <- cowplot::plot_grid(plotlist = ba_plots, labels = c("A", "B", "C", "D", "E"), nrow = 1)
fig_sup_combined <- cowplot::plot_grid(p_title_sup, panel_sup, ncol = 1, rel_heights = c(0.05, 1))
fig_sup_combined
ggsave("Figure_Supervised_val1_Batch_Classification.svg", fig_sup_combined, width = 20, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 6. FINAL SUMMARY
# -----------------------------------------------------------------------------
message("\n", strrep("=", 80))
message("BATCH EFFECTS VALIDATION COMPLETE")
message(strrep("=", 80))

message("\n✓ Unsupervised analysis (5 levels):")
message("  Unsupervised_val1_Batch_Detection_{Original,Weak,Medium,Strong,Extreme}.svg")
message("  Unsupervised_val1_Batch_Detection_Summary.svg")

message("\n✓ Supervised analysis (5 levels):")
message("  Supervised_val1_{original,batch-weak,batch-medium,batch-strong,batch-extreme}_*.svg/csv")

message("\n✓ Combined figures:")
message("  Unsupervised_val1_Batch_Detection_Summary.svg")
message("  Figure_Supervised_val1_Batch_Classification.svg")

message("\n", strrep("=", 80), "\nAll analyses completed successfully!\n", strrep("=", 80))
