#' Perform Ward.D2 hierarchical clustering on projection data
#'
#' @description
#' Performs hierarchical clustering using Ward.D2 on provided data and aligns clusters
#' to best match target classes if given, using the optimal Hungarian method.
#'
#' @param projection_data Data frame or matrix with projection coordinates
#' @param target Optional vector of target (true) classes for cluster alignment
#' @param distance_metric Distance metric to use (default: "euclidean")
#' @param determine_cluster_number Automatically determine cluster number and mebershipt
#'
#' @return List with:
#'   clusters: (aligned) cluster assignments
#'   original_clusters: original cutree assignments before alignment
#'   cluster_object: the hclust object
#'   n_clusters: number of clusters
#'   distance_matrix: distance matrix used
#'
#' @importFrom stats dist hclust cutree
perform_ward_clustering <- function(projection_data,
                                    target = NULL,
                                    distance_metric = "euclidean",
                                    determine_cluster_number = FALSE) {
  # Check input
  if (!is.data.frame(projection_data) && !is.matrix(projection_data))
    stop("projection_data must be a data.frame or a matrix.")
  if (nrow(projection_data) == 0)
    stop("projection_data has no rows.")
  X <- as.data.frame(projection_data)

  # Number of clusters is determined by target length or computed
  if (is.null(target)) {
    warning("No target provided, clustering to a single group.")
    n_clusters <- 1
    target <- rep(1, nrow(X))
  } else {
    # Count unique values EXCLUDING NA
    unique_targets <- unique(target[!is.na(target)])  # Remove NAs before counting
    n_clusters <- length(unique_targets)

    # If all are NA or only one real class, set to 1
    if (n_clusters == 0) {
      warning("All target values are NA, clustering to a single group.")
      n_clusters <- 1
    }
  }

  if (determine_cluster_number) {
    warning("Determine cluster number as required.")
    nb_clust_res <- NbClust::NbClust(umap_result$Projected, method = "ward.D2")
    n_clusters  <- max(unlist(nb_clust_res[4]))
  }

  # Perform clustering
  distance_matrix <- stats::dist(X, method = distance_metric)
  cluster_object <- stats::hclust(distance_matrix, method = "ward.D2")
  cluster_assignments <- as.numeric(stats::cutree(cluster_object, k = n_clusters))

  # Align clusters if applicable
  if (!is.null(target) && n_clusters > 1) {
    # Get numeric alignment from Hungarian algorithm (only for non-NA targets)
    renamed_clusters_numeric <- align_clusters_to_target_hungarian(trueCls = target, currentCls = cluster_assignments)

    # Map the aligned numeric clusters back to original target class names
    # Create a mapping from target to its numeric representation (excluding NA)
    target_factor <- factor(target[!is.na(target)])
    target_levels <- levels(target_factor)

    # Map ALL cluster assignments to class names
    # (including those for samples with NA targets - they still get clustered!)
    renamed_clusters <- target_levels[renamed_clusters_numeric]

  } else {
    renamed_clusters <- cluster_assignments
  }

  return(list(
    clusters = renamed_clusters,
    original_clusters = cluster_assignments,
    cluster_object = cluster_object,
    n_clusters = n_clusters,
    distance_matrix = distance_matrix
  ))
}

# Helper: Cluster alignment using the Hungarian algorithm (optimal assignment)
align_clusters_to_target_hungarian <- function(trueCls, currentCls) {
  # Remove NA targets from alignment calculation
  valid_idx <- !is.na(trueCls)

  if (sum(valid_idx) == 0) {
    # All are NA, return currentCls as-is
    return(currentCls)
  }

  trueCls_valid <- trueCls[valid_idx]
  currentCls_valid <- currentCls[valid_idx]

  trueCls_numeric <- as.integer(factor(trueCls_valid))
  currentCls_numeric <- as.integer(factor(currentCls_valid))

  n_true <- length(unique(trueCls_numeric))
  n_clust <- length(unique(currentCls_numeric))
  n <- max(n_true, n_clust)

  ct <- table(trueCls_numeric, currentCls_numeric)
  # If not square, pad the table
  if (nrow(ct) < n) ct <- rbind(ct, matrix(0, n - nrow(ct), ncol(ct)))
  if (ncol(ct) < n) ct <- cbind(ct, matrix(0, nrow(ct), n - ncol(ct)))

  # Use Hungarian algorithm if available
  if (requireNamespace("clue", quietly = TRUE)) {
    assignment <- clue::solve_LSAP(ct, maximum = TRUE)
    cluster_map <- setNames(as.integer(colnames(ct)[assignment]), seq_len(nrow(ct)))
    cluster_aligned <- cluster_map[as.character(currentCls)]
    # Convert NA (padding) to a new cluster label if any
    cluster_aligned[is.na(cluster_aligned)] <- max(cluster_aligned, na.rm = TRUE) + 1
    # Return aligned numeric clusters (will be mapped to names by caller)
    return(cluster_aligned)
  } else {
    message("clue package not available, using permutation-based alignment (slow for >9 clusters)")
    perms <- combinat::permn(seq_len(n))
    bestAcc <- 0
    bestPerm <- seq_len(n)
    for (p in perms) {
      mapped <- p[match(currentCls_numeric, seq_len(n))]
      acc <- sum(trueCls_numeric == mapped) / length(trueCls_numeric)
      if (acc > bestAcc) {
        bestAcc <- acc
        bestPerm <- p
      }
    }
    cluster_aligned <- bestPerm[match(currentCls, seq_len(n))]
    return(cluster_aligned)
  }
}