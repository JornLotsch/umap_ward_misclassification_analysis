#' Classification stability heatmap (supervised/unsupervised counterpart)
#'
#' Visualizes classification stability across cross-validation iterations.
#' Companion to unsupervised `plot_misclassification_heatmap()`.
#'
#' Creates 3-column heatmap:
#' - Column 1: True/prior class labels (discrete colors)
#' - Column 2: Dominant assigned class (same colors, alpha-coded confidence)
#' - Column 3: Misclassification status (lightyellow2/salmon, alpha-coded error rate)
#'
#' Handles samples with no assignments (all-NA rows) by showing white tiles.
#'
#' @param class_assignments Data frame. Rows = samples, columns = CV iterations.
#'                         Contains class predictions from RF or SVM.
#' @param true_classes Character vector. True class labels for each sample.
#' @param label_col Character vector. Sample labels (defaults to rownames).
#' @param title Character. Plot title (auto-generates if NULL).
#' @param row_font_size Numeric. Row label font size (default: 6).
#'
#' @return List containing:
#'   \itemize{
#'     \item{plot}{Final ggplot2 object}
#'     \item{data}{Clean data frame with row_label, prior_class, assigned_class, misclassified_status}
#'     \item{n_complete}{Number of samples with assignments}
#'   }
#'
#' @import ggplot2 tidyr scales grid
#' @export
plot_classification_stability_heatmap <- function(
    class_assignments, true_classes, label_col = NULL,
    title = NULL, row_font_size = 6
) {
  
  # Load required libraries
  if (!require("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
  if (!require("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required")
  if (!require("scales", quietly = TRUE)) stop("Package 'scales' is required")
  if (!require("grid", quietly = TRUE)) stop("Package 'grid' is required")
  
  # Set default title
  if (is.null(title)) title <- "Classification Stability Analysis"
  
  # Ensure input alignment
  n_rows <- nrow(class_assignments)
  true_classes <- as.character(true_classes[1:n_rows])
  
  # Create row labels as ordered factor (preserves input row order)
  row_label <- factor(rownames(class_assignments), levels = rownames(class_assignments))
  
  # Convert assignments to matrix for processing
  rf_matrix <- as.matrix(class_assignments)
  
  # Identify unique true class levels
  true_class_levels <- unique(true_classes[!is.na(true_classes)])
  
  # Initialize result vectors (length matches input rows)
  assigned_class <- rep(NA_character_, n_rows)
  assigned_class_pct <- rep(0, n_rows)
  true_class_pct <- rep(0, n_rows)
  misclassified_status <- rep(NA_character_, n_rows)
  
  # Process each sample: compute dominant class and confidence
  for(i in 1:n_rows) {
    row_data <- rf_matrix[i, ]
    
    # Handle samples with no assignments
    if(all(is.na(row_data))) {
      assigned_class[i] <- NA_character_
      assigned_class_pct[i] <- 0
    } else {
      # Compute class frequencies across CV iterations
      class_table <- table(na.omit(row_data))
      total <- sum(class_table)
      pct_table <- prop.table(class_table) * 100
      
      # Ensure all true classes present (fill missing with 0%)
      for(cls in true_class_levels) {
        if(!cls %in% names(pct_table)) pct_table[cls] <- 0
      }
      
      # Dominant class = highest frequency
      max_class <- names(pct_table)[which.max(pct_table)]
      assigned_class[i] <- max_class
      assigned_class_pct[i] <- max(pct_table)
    }
    
    # Compute true class percentage and status for samples with assignments
    if(!is.na(assigned_class[i])) {
      true_class_pct[i] <- pct_table[true_classes[i]]
      misclassified_status[i] <- ifelse(assigned_class[i] == true_classes[i], 
                                        "Not misclassified", "Misclassified")
    }
  }
  
  # Compute misclassification percentages
  misclassified_pct <- 100 - true_class_pct
  
  # Calculate performance statistics
  complete_rows <- rowSums(!is.na(rf_matrix)) > 0
  n_complete <- sum(complete_rows)
  
  # Raw misclassification rate (iteration-level)
  raw_misclass_rate <- if(n_complete > 0) {
    valid_preds <- sum(!is.na(rf_matrix[complete_rows, ]))
    sum(rf_matrix[complete_rows, ] != true_classes[complete_rows], na.rm = TRUE) / valid_preds * 100
  } else 0
  
  # Dominant class misclassification rate (sample-level)
  dominant_misclass_rate <- mean(misclassified_status == "Misclassified", na.rm = TRUE) * 100
  
  # Subtitle with dual performance metrics
  subtitle <- sprintf("Dominant: %.1f%% | Raw: %.1f%% (%d/%d rows)", 
                      dominant_misclass_rate, raw_misclass_rate, n_complete, n_rows)
  
  # Create summary data frame
  clean_data <- data.frame(
    row_label = row_label,
    prior_class = true_classes,
    assigned_class = assigned_class,
    misclassified_status = misclassified_status,
    stringsAsFactors = FALSE
  )
  
  # Convert to long format for heatmap
  df_long <- tidyr::pivot_longer(
    clean_data,
    cols = c("prior_class", "assigned_class", "misclassified_status"),
    names_to = "column",
    values_to = "value"
  )
  
  # Set column factor levels
  df_long$row <- df_long$row_label
  df_long$column <- factor(df_long$column, 
                           levels = c("prior_class", "assigned_class", "misclassified_status"),
                           labels = c("Prior Class", "Assigned Class", "Misclassified"))
  
  # Colorblind-friendly palette matching sister function
  class_levels <- unique(c(true_classes, assigned_class[!is.na(assigned_class)]))
  n_classes <- length(class_levels)
  cb_palette_base <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                       "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  # Extend palette if needed
  if (n_classes > length(cb_palette_base)) {
    cb_palette <- c(cb_palette_base, rainbow(n_classes - length(cb_palette_base)))
  } else {
    cb_palette <- cb_palette_base[1:n_classes]
  }
  
  col_vals <- setNames(cb_palette, class_levels)
  
  # Define complete value levels and colors (matches sister function)
  value_levels <- c(class_levels, "Not misclassified", "Misclassified")
  fill_colors <- c(col_vals, 
                   "Not misclassified" = "lightyellow2", 
                   "Misclassified" = "salmon")
  
  df_long$value <- factor(df_long$value, levels = value_levels)
  
  # Create heatmap
  plot <- ggplot2::ggplot(df_long, aes(x = column, y = row, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = fill_colors,
      na.value = "white",
      name = "Assignment"
    ) +
    ggplot2::scale_y_discrete(limits = rev(levels(row_label))) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = row_font_size),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = grid::unit(0.4, "lines"),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.5), color = NA),
      plot.title = ggplot2::element_text(size = 14),
      plot.subtitle = ggplot2::element_text(size = 12),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  
  # Return plot and metadata
  list(
    plot = plot,
    data = clean_data,
    n_complete = n_complete
  )
}
