# example_umap_ward_misclassification_analysis.R

# Example usage of umap_ward_misclassification_analysis

# Load your data frame: each row is a sample, each column a feature (e.g., lipid species)
lipid_profiles <- read.csv("lipid_profiles.csv")

# Load or extract your class labels for each sample
sample_metadata <- read.csv("sample_metadata.csv")
sample_types <- sample_metadata$SampleType

# Run the integrated UMAP projection and clustering/misclassification analysis
set.seed(42)
results <- umap_ward_misclassification_analysis(
  data = lipid_profiles,                 # Features data
  target = sample_types,                 # Ground truth (prior classes)
  labels = sample_metadata$SampleID,     # Optional: row labels for plots, if available
  output_dir = "qc_results"              # Directory to save QC visualizations and results
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

# Optionally: view UMAP Voronoi plot object
print(results$voronoi_plot)