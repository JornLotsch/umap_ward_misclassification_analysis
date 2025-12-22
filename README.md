# UMAP Ward misclassification analysis
## When artifacts masquerade as discovery: A case study revealing hidden laboratory errors in lipidomics data with biologically apparently plausible results.

This R code provides a comprehensive workflow for performing UMAP projection, Ward clustering, and misclassification analysis on high-dimensional data. It's particularly useful for quality control and exploratory data analysis in omics studies.

## Installation

You can download this code to you local hard drive and run it from there. 

## Functions
### `check_and_install_packages(pkg_list)`
**Description**: Utility function that checks for missing R packages and automatically installs them from CRAN.
**Parameters**:
- `pkg_list`: Character vector of package names to check and install

**Returns**:
- No return value (NULL). Installs missing packages and loads all required packages
- Throws an error if any packages fail to load after installation

**Example**:
``` r
required_packages <- c("ggplot2", "umap", "cluster")
check_and_install_packages(required_packages)
```
### `source_required_functions()`
**Description**: Sources all required analysis functions from the current working directory.
**Parameters**: None
**Returns**:
- No return value (NULL). Sources external R scripts containing analysis functions
- Throws an error if any required function files are missing

**Required Files**:
- `prepare_dataset.R`
- `perform_umap_projection.R`
- `perform_ward_clustering.R`
- `plot_umap_with_voronoi.R`
- `plot_misclassification_heatmap.R`
- `perform_supervised_classification.R`

### `umap_ward_misclassification_analysis()`
**Description**: Main analysis function that performs a complete workflow including data preparation, UMAP projection, Ward clustering, visualization, and misclassification analysis.
#### Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `data` | data.frame/matrix | Required | Input data with samples as rows and features as columns |
| `target` | vector | NULL | Optional target class labels (length must match nrow(data)) |
| `labels` | vector | NULL | Optional sample labels for visualization (length must match nrow(data)) |
| `determine_cluster_number` | logical | FALSE | Whether to automatically determine cluster number and membership |
| `voronoi_targets` | logical | TRUE | Whether Voronoi cells are colored for prior targets (TRUE) or for clusters (FALSE) |
| `output_dir` | character | "results" | Directory for saving output files |
| `file_prefix` | character | "umap_analysis" | Prefix for output filenames |
| `file_format` | character | "svg" | Plot format: "svg" or "png" |
| `label_points` | logical | TRUE | Whether to display point labels in plots |
| `row_font_size` | numeric | 6 | Font size for heatmap row labels |
| `width` | numeric | 12 | Plot width in inches |
| `height` | numeric | 9 | Plot height in inches |
| `dpi` | integer | 300 | Resolution for PNG output |
| `n_neighbors` | integer | 15 | Number of nearest neighbors for UMAP |
#### Returns
A named list containing the following components:

| Component | Type | Description |
| --- | --- | --- |
| `prepared_data` | data.frame | Processed input data ready for analysis |
| `umap_result` | list | Complete UMAP projection results including coordinates and parameters |
| `cluster_result` | list | Ward clustering results including cluster assignments |
| `voronoi_plot` | ggplot2 | UMAP scatter plot with Voronoi tessellation overlay |
| `heatmap_result` | list | Misclassification heatmap and associated statistics |
| `combined_plot` | grob | Combined visualization of Voronoi plot and heatmap |
| `misclassification_rate` | numeric | Proportion of misclassified samples (0-1) |
| `misclassified_samples` | data.frame | Details of misclassified samples with expected vs assigned classes |
#### Output Files
The function generates several output files in the specified directory:
1. **`{file_prefix}_voronoi.{format}`**: UMAP plot with Voronoi tessellation
2. **`{file_prefix}_heatmap.{format}`**: Misclassification heatmap
3. **`{file_prefix}_combined.{format}`**: Combined visualization
4. **`{file_prefix}_misclassified_samples.csv`**: CSV file listing misclassified samples (if any)

## Usage Examples

### Basic Usage with Sample Data

``` r
# Load your data
lipid_profiles <- read.csv("lipid_profiles.csv")
sample_metadata <- read.csv("sample_metadata.csv")
sample_types <- sample_metadata$SampleType

# Run the analysis
results <- umap_ward_misclassification_analysis(
  data = lipid_profiles,
  target = sample_types,
  labels = sample_metadata$SampleID,
  output_dir = "qc_results",
  file_prefix = "lipid_analysis",
  file_format = "png",
  width = 14,
  height = 10
)

# Check misclassification rate
cat("Misclassification rate:", 
    sprintf("%.2f%%", results$misclassification_rate * 100), "\n")

# View misclassified samples
if (nrow(results$misclassified_samples) > 0) {
  print(results$misclassified_samples)
}

# Access individual components
umap_coordinates <- results$umap_result$Projected
cluster_assignments <- results$cluster_result$clusters
```

## Example Scripts

This repository includes several example R scripts demonstrating different use cases:

### Main Publication Example

**`lipid_case_data_run.R`**: Demonstrates UMAP-Ward clustering analysis on real lipidomics data (PsA patients vs. controls). This is the primary example associated with the publication (Lotsch et al., 2025) and shows both standard analysis and randomized permutation analyses for validation.

### Demonstration Examples

**`golfball_dataset_run.R`**: Illustrates the workflow using an artificial "golfball" dataset consisting of concentric spheres in 3D space. This example clearly demonstrates how UMAP can separate distinct geometric structures and how Ward clustering assigns points to clusters. The output includes a combined visualization (`golfballs_combined_plot.svg`) showing:
- **Left panel**: UMAP projection with Voronoi tessellation colored by cluster assignments
- **Right panel**: Misclassification heatmap comparing prior class assignments (spheres) with Ward-assigned clusters

![Golfball UMAP Analysis Results](golfballs_combined_plot.svg)

**`create_sample_files.R`**: Generates minimal synthetic lipidomics data for quick testing and learning purposes. Creates `lipid_profiles.csv` and `sample_metadata.csv` files with balanced class distributions and simulated class-specific feature variations.

**`example_run.R`**: A simplified introductory example showing the basic workflow with minimal configuration.
## Dependencies

### Core Analysis Packages
The main analysis functions automatically install and load the following R packages:
- `ggplot2`: Data visualization 
- `tidyr`: Data manipulation
- `scales`: Plot scaling functions
- `deldir`: Voronoi tessellation
- `umap`: UMAP dimensionality reduction
- `gridExtra`: Combining plots
- `ggrepel`: Text labeling in plots
- `clue`: Optimal cluster assignment (Hungarian algorithm)
- `NbClust`: Determining optimal number of clusters (optional, if `determine_cluster_number = TRUE`)

### Supervised Classification Packages
For supervised classification analysis (via `perform_supervised_classification.R`):
- `caret`: Machine learning framework
- `randomForest`: Random Forest implementation
- `e1071`: Support Vector Machines
- `caTools`: Utility functions
- `pbmcapply`: Parallel processing with progress bars
- `dplyr`: Data manipulation

### Example Script Dependencies
For running the example scripts, additional packages may be required:
- `ComplexHeatmap`: Advanced heatmap visualizations (used in lipid and golfball examples)
- `plot3D`: 3D plotting (golfball dataset generation)
- `car`: Companion to Applied Regression (golfball data preparation)
- `ggplotify`: Convert base plots to ggplot2 objects (golfball visualization) 

## Requirements
- R >= 3.6.0
- External function files (listed in `source_required_functions()`)
- Write permissions in the output directory

## Error Handling
The function includes comprehensive input validation:
- Checks data types and dimensions
- Validates parameter ranges and formats
- Ensures required files exist
- Verifies directory permissions
- Provides informative error messages for troubleshooting

## Additional Resources

### Supervised Classification Analysis

**`perform_supervised_classification.R`**: Provides supervised machine learning classification using Random Forest and Support Vector Machine (SVM) models. This function:
- Performs repeated random subsampling validation across multiple target variables
- Includes automated hyperparameter tuning for both RF and SVM models
- Computes class-balanced accuracy with 95% confidence intervals
- Handles small and imbalanced datasets robustly
- Supports parallel processing for efficient computation

This is useful for evaluating the predictive accuracy of features or comparing classification performance across different datasets and models.

### Alternative Projection Methods
While this package focuses on UMAP for dimensionality reduction, the accompanying file `projection_methods_snippet.R` provides implementations of multiple alternative projection methods that can be used for exploratory analysis, including:

**Unsupervised methods:**
- PCA (Principal Component Analysis)
- ICA (Independent Component Analysis)
- MDS (Multidimensional Scaling)
- t-SNE (t-Distributed Stochastic Neighbor Embedding)
- Isomap (Isometric Feature Mapping)
- LLE (Local Linear Embedding)
- NMF (Non-negative Matrix Factorization)

**Supervised methods:**
- LDA (Linear Discriminant Analysis)
- PLS-DA (Partial Least Squares Discriminant Analysis)
- PLS-LDA (PLS followed by LDA)
- IPCA (Independent Principal Component Analysis)

The snippet demonstrates the `perform_projection()` function with detailed parameter documentation and usage notes for each method. This can be useful for comparing different dimensionality reduction approaches on your data.

## License

GPL-3

## Citation

If you use this tool in your work, please cite:

Lotsch J, Kringel D, Hahnefeld L, Gurke R, Himmelspach A, Behrens F, and Geisslinger G. When artifacts align with biology: The case for investigative flexibility over standardized quality control in lipidomics datasets. *2025* (submitted).