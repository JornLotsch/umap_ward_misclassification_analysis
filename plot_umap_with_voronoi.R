#' Plot UMAP projection with Voronoi cells colored by target
#'
#' @description
#' Creates a visualization of UMAP projection data with points and Voronoi cells
#' both colored by target values. Provides a clear visualization of data clusters
#' and their boundaries.
#'
#' @param umap_projection Data frame with projection coordinates or result from
#'                       perform_umap_projection()
#' @param targets Vector of target values for coloring points and also cells if cluster is NULL 
#' @param clusters Vector of cluster values for coloring cells if not NULL
#' @param labels Optional vector of labels for points (default: row names or indices)
#' @param label_points Logical, whether to display labels (default: FALSE)
#'
#' @return A ggplot2 object displaying the UMAP projection with Voronoi cells
#' @importFrom ggplot2 ggplot geom_polygon geom_point theme_light theme labs scale_shape_manual scale_fill_manual scale_color_manual guides
#' @importFrom ggrepel geom_text_repel
#' @importFrom deldir deldir tile.list
plot_umap_with_voronoi <- function(umap_projection,
                                   targets,
                                   clusters = NULL,
                                   labels = NULL,
                                   label_points = FALSE) {

  # Check required packages
  required_packages <- c("ggplot2", "ggrepel", "deldir")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed. Please install it."))
    }
  }

  # Extended colorblind palette
  cb_palette <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7"
  )
  shape_values <- c(17, 16, 15, 18, 2, 1, 0, 5, 7, 8, 9, 19, 11, 12)

  # Extract data if input is the result from perform_umap_projection
  X <- umap_projection

  # Initial checks of arguments
  if (ncol(X) < 2) stop("The projection data must have at least two columns.")
  if (nrow(X) != length(targets)) stop("The length of 'targets' must match the number of rows in the projection data.")
  if (!is.null(clusters)) {
    if (nrow(X) != length(clusters)) stop("The length of 'clusters' must match the number of rows in the projection data.")
  } 
  # Assign labels if not present
  if (is.null(labels)) {
    if (!is.null(rownames(X))) {
      labels <- rownames(X)
    } else {
      labels <- seq_len(nrow(X))
    }
  }

  plotData <- data.frame(
    Proj1 = X[, 1],
    Proj2 = X[, 2],
    Target = targets,
    Cluster = if (is.null(clusters)) targets else clusters,
    Label = labels
  )

  # Ensure Target and Clusters are factors for consistent coloring/shaping
  plotData$Target <- as.factor(plotData$Target)
  plotData$Cluster <- as.factor(plotData$Cluster)

  plotData$Cells <- if (is.null(clusters)) plotData$Target else plotData$Cluster

  # Ensure we have enough colors for all target values
  unique_targets <- unique(plotData$Target)
  target_palette <- rep(cb_palette, length.out = length(unique_targets))

  # Ensure we have enough colors for all cluster values
  unique_clusters <- unique(plotData$Cluster)
  cluster_palette <- rep(cb_palette, length.out = length(unique_clusters))

  # Voronoi diagram computation
  voronoiData <- deldir::deldir(plotData$Proj1, plotData$Proj2)

  # Convert Voronoi tessellation to a data frame for plotting
  vor_polys <- deldir::tile.list(voronoiData)
  voronoi_df <- do.call(rbind, lapply(seq_along(vor_polys), function(i) {
    data.frame(x = vor_polys[[i]]$x, y = vor_polys[[i]]$y, id = i)
  }))

  # Create plot with ggplot2
  plot <- ggplot2::ggplot() +
  # Voronoi cells colored by target
  ggplot2::geom_polygon(
      data = voronoi_df,
      ggplot2::aes(x = x, y = y, group = id, fill = as.factor(plotData$Cells[id])),
      alpha = 0.3,
      color = NA
    ) +
  # Points colored by target
  ggplot2::geom_point(
      data = plotData,
      ggplot2::aes(
        x = Proj1,
        y = Proj2,
        color = as.factor(Target),
        shape = as.factor(Target)
      )
    ) +
    ggplot2::theme_light(base_size = 14) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.08),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.background = ggplot2::element_rect(
        color = "transparent",
        fill = ggplot2::alpha("white", 0.2)
      )
    ) +
    ggplot2::scale_shape_manual(values = shape_values) +
    ggplot2::scale_color_manual(values = target_palette) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(label = "")))

  if (is.null(clusters)) {
    plot <- plot + ggplot2::scale_fill_manual(values = target_palette) +
    ggplot2::labs(title = "UMAP projection", subtitle = paste0(nrow(X), " evaluable samples)"), 
                  x = "UMAP 1", y = "UMAP 2", color = "Target", fill = "Target", shape = "Target") 
  } else {
    plot <- plot + ggplot2::scale_fill_manual(values = cluster_palette) +
    ggplot2::labs(title = "UMAP projection", subtitle = paste0(nrow(X), " evaluable samples)"), 
                  x = "UMAP 1", y = "UMAP 2", color = "Cluster", fill = "Cluster", shape = "Target") 
  }

  # Conditional labeling of points
  if (label_points) {
    plot <- plot +
      ggrepel::geom_text_repel(
        data = plotData,
        ggplot2::aes(
          x = Proj1,
          y = Proj2,
          color = as.factor(Target),
          label = Label,
          fontface = 2
        ),
        vjust = -1,
        size = 3,
        max.overlaps = Inf
      )
  }

  return(plot)
}