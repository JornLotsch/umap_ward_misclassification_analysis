#' create_sample_files.R
#'
#' Creates sample lipidomics data with multiple metadata variables for testing UMAP/clustering analysis.
#' Generates both feature data (lipid profiles) and comprehensive sample metadata with multiple target variables.
#' Data is designed to demonstrate both single-target and multi-target misclassification analysis.
#'
#' @description
#' Produces two CSV files:
#'   - `lipid_profiles.csv`: 100 samples × 8 lipid features with class-specific separation
#'   - `sample_metadata.csv`: Sample IDs + 5 metadata variables (SampleType, BatchID, TimePoint, AgeGroup, Gender)
#'
#' Features are constructed to:
#' - Show clear SampleType separation (Lipid1/2 vs Lipid5/6)
#' - Include independent metadata variables for multi-target testing
#' - Add controlled noise for realistic misclassification scenarios
#'
#' @return No return value. Creates two CSV files in current directory.
#'
#' @examples
#' source("create_sample_files.R")
#' # Then use:
#' lipid_data <- read.csv("lipid_profiles.csv")
#' metadata <- read.csv("sample_metadata.csv")
#'
#' # Single-target: umap_ward_misclassification_analysis(lipid_data, metadata$SampleType)
#' # Multi-target: umap_ward_misclassification_analysis_multi(lipid_data, metadata[,c("SampleType", "BatchID", "AgeGroup")])
#'
#' @author Adapted for comprehensive UMAP misclassification testing
#' @export
create_sample_files <- function() {
  
  # --- Set random seeds for reproducibility ---
  set.seed(123)
  
  # --- Define dataset parameters ---
  n_samples <- 100         # Total samples
  n_features <- 8          # Lipid features
  samples_per_class <- 50  # Balanced 2-class design
  
  message("=== Creating sample lipidomics data ===")
  message(sprintf("Dataset: %d samples × %d features", n_samples, n_features))
  
  # --- Generate sample identifiers ---
  sample_ids <- sprintf("S%03d", 1:n_samples)
  
  # --- Primary target variable (SampleType) ---
  # Drives main feature separation for single-target testing
  sample_types <- rep(c("ClassA", "ClassB"), each = samples_per_class)
  
  # --- Additional metadata variables (independent of SampleType) ---
  # Designed for multi-target analysis demonstrations
  
  # BatchID: Technical batches (4 levels, evenly distributed)
  batch_ids <- rep(c("Batch1", "Batch2", "Batch3", "Batch4"), each = 25)
  
  # TimePoint: Temporal variable (4 time points, sequential)
  time_points <- rep(c("T0", "T1", "T2", "T3"), times = 25)
  
  # AgeGroup: Biological variable, UNBALANCED (70% Young, 30% Senior), random
  set.seed(789)
  age_groups <- sample(c("Young", "Senior"), size = n_samples, 
                       replace = TRUE, prob = c(0.7, 0.3))
  
  # Gender: Biological variable, roughly balanced, random
  set.seed(456)
  gender <- sample(c("Male", "Female"), size = n_samples, replace = TRUE)
  
  # --- Generate lipid feature matrix ---
  # Base: multivariate normal noise
  feature_matrix <- matrix(rnorm(n_samples * n_features, mean = 0, sd = 1),
                           nrow = n_samples, ncol = n_features)
  
  # Add class-specific signals (SampleType separation)
  # ClassA: elevate Lipid1, Lipid2
  feature_matrix[1:samples_per_class, 1:2] <- 
    feature_matrix[1:samples_per_class, 1:2] + 2
  
  # ClassB: elevate Lipid5, Lipid6
  feature_matrix[(samples_per_class+1):n_samples, 5:6] <- 
    feature_matrix[(samples_per_class+1):n_samples, 5:6] + 2
  
  # Add controlled noise (simulates misclassification cases)
  n_errors <- 6
  set.seed(123)
  error_indices <- sample(1:n_samples, n_errors)
  feature_matrix[error_indices, ] <- 
    feature_matrix[error_indices, ] + rnorm(n_errors * n_features, mean = 0, sd = 2)
  
  # Name features
  colnames(feature_matrix) <- paste0("Lipid", 1:n_features)
  
  # --- Create and save lipid profiles ---
  lipid_profiles <- as.data.frame(feature_matrix)
  write.csv(lipid_profiles, "lipid_profiles.csv", row.names = FALSE)
  message("✓ Saved: lipid_profiles.csv (", n_samples, " samples × ", n_features, " lipids)")
  
  # --- Create comprehensive sample metadata ---
  sample_metadata <- data.frame(
    SampleID = sample_ids,
    SampleType = sample_types,
    BatchID = batch_ids,
    TimePoint = time_points,
    AgeGroup = age_groups,
    Gender = gender,
    stringsAsFactors = FALSE
  )
  write.csv(sample_metadata, "sample_metadata.csv", row.names = FALSE)
  message("✓ Saved: sample_metadata.csv (", n_samples, " samples, 6 variables)")
  
  # --- Print detailed summary statistics ---
  cat("\n=== DATASET SUMMARY ===\n")
  cat(sprintf("Features: %d samples × %d lipids\n", n_samples, n_features))
  cat("\nSampleType distribution (primary target):\n")
  print(table(sample_types))
  
  cat("\nMetadata distributions:\n")
  cat("  BatchID:", paste0("(", length(unique(batch_ids)), " batches)"), "\n")
  cat("  TimePoint:", paste0("(", length(unique(time_points)), " time points)"), "\n")
  cat("  AgeGroup:", sum(age_groups == "Young"), "Young,", 
      sum(age_groups == "Senior"), "Senior (70/30 unbalanced)\n")
  cat("  Gender:", sum(gender == "Male"), "Male,", 
      sum(gender == "Female"), "Female\n")
  
  cat("\n=== USAGE EXAMPLES ===\n")
  cat("Single-target analysis:\n")
  cat("  umap_ward_misclassification_analysis(lipid_profiles, sample_metadata$SampleType)\n")
  cat("\nMulti-target analysis:\n")
  cat("  targets <- sample_metadata[, c('SampleType', 'BatchID', 'TimePoint', 'AgeGroup')]\n")
  cat("  umap_ward_misclassification_analysis_multi(lipid_profiles, targets)\n")
  
  cat("\n=== FILES CREATED ===\n")
  cat("✓ lipid_profiles.csv\n")
  cat("✓ sample_metadata.csv\n")
  
  invisible(list(
    lipid_profiles = lipid_profiles,
    sample_metadata = sample_metadata
  ))
}

# --- Execute if script is run directly ---
if (sys.nframe() == 0) {
  create_sample_files()
}
