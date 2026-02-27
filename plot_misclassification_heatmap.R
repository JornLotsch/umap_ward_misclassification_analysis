#' Create a heatmap showing cluster assignments and misclassifications
#'
#' @description
#' Generates a heatmap comparing prior class assignments with clustering results to
#' identify and visualize misclassified samples. Supports multiple clusters and uses
#' a colorblind-friendly palette.
#'
#' @param data Data frame containing prior classes and cluster assignments
#' @param prior_class_col Character, column name for prior classes (default: "Target")
#' @param cluster_col Character, column name for cluster assignments (default: "Cluster")
#' @param label_col Character, column name for row labels (default: row names)
#' @param title Character, plot title (default: "Cluster Assignments and Misclassifications")
#' @param row_font_size Numeric, font size for row labels (default: 6)
#'
#' @return A list containing:
#'   \item{plot}{ggplot2 object with the misclassification heatmap}
#'   \item{data}{Data frame with prior classes, clusters and misclassification indicators}
#'   \item{long_data}{Long-format data used for plotting}
#'   \item{misclassification_rate}{Numeric, percentage of misclassified samples}
#'
#' @importFrom ggplot2 ggplot geom_tile scale_fill_manual theme_minimal labs theme element_text element_rect element_blank
#' @importFrom tidyr pivot_longer
#' @importFrom scales alpha
plot_misclassification_heatmap <- function(
  data,
  prior_class_col = "Target",
  cluster_col = "Cluster",
  label_col = NULL,
  title = "Cluster Assignments and Misclassifications",
  row_font_size = 6
) {
  # Check required packages
  required_packages <- c("ggplot2", "tidyr", "scales")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed"))
    }
  }

  # Input validation
  if (!is.data.frame(data)) stop("Input must be a data frame")
  if (!(prior_class_col %in% colnames(data))) stop(paste("Column", prior_class_col, "not found in data"))
  if (!(cluster_col %in% colnames(data))) stop(paste("Column", cluster_col, "not found in data"))
  if (!is.null(label_col) && !(label_col %in% colnames(data))) stop(paste("Label column", label_col, "not found in data"))

  # Use row names if label column not specified
  row_label <- if (is.null(label_col)) {
    rn <- rownames(data)
    if (!is.null(rn) && all(rn != "")) rn else as.character(seq_len(nrow(data)))
  } else {
    as.character(data[[label_col]])
  }

  # Coerce prior_class and cluster to character for reliable comparison, do NOT mutate original data
  prior_class <- as.character(data[[prior_class_col]])
  cluster <- as.character(data[[cluster_col]])

  # Determine misclassification (robust to type mismatches)
  misclassified <- ifelse(prior_class == cluster, 0L, 1L)

  clean_data <- data.frame(
    row_label = row_label,
    prior_class = prior_class,
    cluster = cluster,
    misclassified = misclassified,
    stringsAsFactors = FALSE
  )

  # Prepare data for heatmap in long format
  prepare_long_data <- function(df) {
    df_plot <- data.frame(
      row = df$row_label,
      prior_class = df$prior_class,
      cluster = df$cluster,
      misclassified = ifelse(df$misclassified == 1L, "Misclassified", "Not misclassified"),
      stringsAsFactors = FALSE
    )
    # Gather to long format
    df_long <- tidyr::pivot_longer(
      df_plot,
      cols = c("prior_class", "cluster", "misclassified"),
      names_to = "column",
      values_to = "value"
    )
    # Reverse row order for heatmap display
    df_long$row <- factor(df_long$row, levels = rev(unique(df_long$row)))
    df_long$column <- factor(df_long$column, levels = c("prior_class", "cluster", "misclassified"),
                             labels = c("Prior Class", "Cluster", "Misclassified"))
    df_long
  }
  df_long <- prepare_long_data(clean_data)

  # Build colorblind-friendly palette (expand if needed) for all unique class/cluster values
  cluster_levels <- sort(unique(c(clean_data$prior_class, clean_data$cluster)))
  n_clusters <- length(cluster_levels)
  cb_palette_base <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
  if (n_clusters > length(cb_palette_base)) {
    more_cols <- grDevices::rainbow(n_clusters - length(cb_palette_base))
    cb_palette <- c(cb_palette_base, more_cols)
  } else {
    cb_palette <- cb_palette_base
  }
  col_vals <- setNames(cb_palette[seq_len(n_clusters)], cluster_levels)

  # Specify color mapping for legend values
  value_levels <- c(cluster_levels, "Not misclassified", "Misclassified")
  fill_colors <- c(col_vals, "Not misclassified" = "lightyellow2", "Misclassified" = "salmon")

  # Set the factor level order for fill
  df_long$value <- factor(df_long$value, levels = value_levels)

  # Legend labels mapping
  legend_labels <- c(setNames(paste("Class/cluster", cluster_levels), cluster_levels),
                     "Not misclassified" = "Not misclassified",
                     "Misclassified" = "Misclassified")

  # Calculate misclassification rate
  misclass_rate <- round(100 * mean(clean_data$misclassified == 1L), 1)
  subtitle <- paste0("Misclassification rate: ", misclass_rate, "%")

  plot <- ggplot2::ggplot(df_long, ggplot2::aes(x = column, y = row, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = fill_colors,
      na.value = "white",
      name = "Assignment",
      labels = legend_labels
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL, y = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = row_font_size),
      legend.position = "top",
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.5), color = NA),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  list(
    plot = plot,
    data = clean_data,
    long_data = df_long,
    misclassification_rate = misclass_rate
  )
}