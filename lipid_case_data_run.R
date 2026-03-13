#===============================================================================
# UMAP-Ward Misclassification and Supervised Classification Analysis
# Lipid Profiles: Psoriatic Arthritis (PsA) vs Controls
# Reproducible Analysis Pipeline
#===============================================================================

# Required libraries (installed automatically if missing)
required_packages <- c("ggplot2", "dplyr", "tidyr", "tidyverse", "cowplot", "ComplexHeatmap", "stringr",
                       "grid", "umap", "deldir", "gridExtra", "cABCanalysis", "twosamples")
for (pkg in required_packages) {
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
  if (length(missing_files) > 0) {
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
# 2. DATA PROCESSING FUNCTIONS
# -----------------------------------------------------------------------------

perform_recursive_cABC_analysis <- function(named_vector, alpha = 0.05, max_iter = 50) {
  # Check names and if no names, set dummy names 
  if (is.null(names(named_vector))) {
    named_vector <- setNames(named_vector, paste0("DummyName", seq_along(named_vector)))
  }

  # If full vector already uniform: return it directly, without ABC
  p0 <- twosamples::ad_test(
    a = named_vector,
    b = runif(
      n = length(named_vector),
      min = min(named_vector),
      max = max(named_vector)
    )
  )[["P-Value"]]

  if (p0 > alpha) {
    return(list(
      ABC_actual = names(named_vector), # all items
      A_values = named_vector,
      iterations = 0L,
      last_p_value = p0,
      is_uniform_end = TRUE,
      note = "Input was already uniform; no ABC performed."
    ))
  }

  i <- 0L
  is_uniform <- FALSE
  named_vector_A <- named_vector # current A-set candidate
  prev_named_vector_A <- named_vector # last non-uniform set
  prev_ABC_actual <- names(named_vector) # its ABC result (placeholder)
  prev_pval <- p0

  had_error <- FALSE
  last_error_message <- NULL

  while (!is_uniform && i < max_iter) {
    # ABC on current (non-uniform) set
    ABC_res <- tryCatch(
      cABCanalysis::cABC_analysis(named_vector_A),
      error = function(e) {
        message("cABC_analysis failed: ", e$message)
        last_error_message <<- e$message
        NULL
      }
    )

    # If ABC failed, break and use the last non-uniform state
    if (is.null(ABC_res)) {
      had_error <- TRUE
      break
    }

    ABC_actual <- ABC_res$Aind # names of A-items
    prev_named_vector_A <- named_vector_A # remember last non-uniform set
    prev_ABC_actual <- ABC_actual

    # build new A-set
    named_vector_A <- named_vector_A[ABC_actual]

    # test uniformity of new A-set
    pval <- twosamples::ad_test(
      a = named_vector_A,
      b = runif(
        n = length(named_vector_A),
        min = min(named_vector_A),
        max = max(named_vector_A)
      )
    )[["P-Value"]]

    is_uniform <- pval > alpha
    prev_pval <- pval
    i <- i + 1L
  }

  # If we had an error in cABC_analysis, return a clear failure object
  if (had_error) {
    return(list(
      ABC_actual = prev_ABC_actual, # last known ABC set (or initial names)
      A_values = prev_named_vector_A, # its values
      iterations = i,
      last_p_value = prev_pval,
      is_uniform_end = FALSE,
      note = paste0("Stopped because cABC_analysis failed: ", last_error_message)
    ))
  }

  # Normal return: last *non-uniform* input vector result
  list(
    ABC_actual = prev_ABC_actual, # names/indices from last non-uniform ABC
    A_values = prev_named_vector_A, # its values
    iterations = i,
    last_p_value = prev_pval,
    is_uniform_end = is_uniform,
    note = if (is_uniform) "Converged to uniform A-set." else "Reached max_iter without uniformity."
  )
}

# # Test
# named_vector <- setNames(iris[[1]], as.character(rownames(iris)))
# perform_recursive_cABC_analysis(named_vector = median_var_imp_rf$PsA)


# Wrapper to run the cABC analysis only when there are mor varibales than min_vars_4_ABC 
run_cabc_if_needed <- function(x, min_vars_4_ABC = 10, alpha = 0.05, max_iter = 50) {
  # x is a named numeric vector of importance scores

  if (length(x) <= min_vars_4_ABC) {
    # no cABC, just pass through
    return(list(
      ABC_actual = names(x),
      A_values = x,
      iterations = 0L,
      last_p_value = NA_real_,
      is_uniform_end = NA,
      note = "Skipped cABC: <= 10 variables."
    ))
  }

  perform_recursive_cABC_analysis(x, alpha = alpha, max_iter = max_iter)
}


# Function to extract summary from one target's ABC results
extract_abc_summary <- function(target_results) {
  abc_vars <- names(target_results$ABC_actual)
  n_vars <- length(abc_vars)

  tibble(
    N_vars = n_vars,
    All_ABC_vars = if (n_vars > 0) paste(abc_vars, collapse = "; ") else "No ABC performed",
    Iterations = target_results$iterations,
    Last_p_value = round(target_results$last_p_value, 4),
    Is_uniform_end = target_results$is_uniform_end,
    Note = if ("note" %in% names(target_results)) target_results$note else NA_character_
  )
}

# -----------------------------------------------------------------------------
# 3. DATA LOADING AND PREPROCESSING
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
# 4. UMAP-WARD CLUSTERING ANALYSIS (REAL DATA)
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
# 5. NULL MODEL: PERMUTED DATA ANALYSIS
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
# 6. SUPERVISED CLASSIFICATION ANALYSIS
# -----------------------------------------------------------------------------
message("\n=== Supervised Classification Analysis ===")

# Analysis parameters
SEED <- 42
training_size <- 0.67
n_iterations <- 100
RF_only_plotted <- TRUE
n_important_vars <- 10
min_vars_4_ABC <- 3
show_relevant_vars_default <- TRUE
skip_tuning = FALSE
min_class_size <- 3
remove_ltx_case_targets = TRUE

supervised_results <- lapply(c(FALSE, TRUE), function(permute_data) {
  lapply(c(FALSE, TRUE), function(controls_only) {

    if (permute_data) {
      show_relevant_vars <- FALSE
    } else
      show_relevant_vars <- show_relevant_vars_default

    # -----------------------------------------------------------------------------
    # 6.1. DATA PREPARATION
    # -----------------------------------------------------------------------------
    message(sprintf("  Preparing data (permute=%s, controls_only=%s)...", permute_data, controls_only))

    # Load data
    lipid_profiles <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished/PsA_lipids_BC.csv", row.names = 1)
    lipid_metadata <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/DataSetPublished/PsA_classes.csv", row.names = 1)

    table(lipid_metadata$Sampling_date)

    # Extract processing month and remove sampling date (because unsuited case numbers per unique class)
    lipid_metadata$Month <- format(
      as.Date(lipid_metadata$Sampling_date, format = "%m/%d/%Y"), "%B"
    )
    lipid_metadata <- lipid_metadata[!names(lipid_metadata) %in% c("Sampling_date")]


    # Filter to controls only (if enabled)
    if (controls_only) {
      control_lab_dates <- read.csv("/home/joern/Aktuell/RheumaMetabolomicsFFM/09Originale/Controls_Lab_Dates.csv", row.names = 2, stringsAsFactors = FALSE)
      control_ids <- rownames(control_lab_dates)
      lipid_profiles <- lipid_profiles[rownames(lipid_profiles) %in% control_ids,]
      lipid_metadata <- lipid_metadata[rownames(lipid_metadata) %in% control_ids,]
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

    # Remove classifications with <2 class members 
    if (remove_ltx_case_targets) {
      metadata_to_remove <- names(which(apply(lipid_metadata, 2, function(x) min(table(na.omit(x)))) < min_class_size))
      lipid_metadata <- lipid_metadata[!names(lipid_metadata) %in% metadata_to_remove]
    }

    # -----------------------------------------------------------------------------
    # 6.2. RUN CLASSIFICATION ANALYSIS BatchID_Endocannabinoids
    # -----------------------------------------------------------------------------
    message(sprintf("  Running %d iterations...", n_iterations))
    results <- perform_supervised_classification(
      X = lipid_profiles_for_classification,
    # Y = subset(lipid_metadata, select = c("BatchID_Endocannabinoids")),
      Y = lipid_metadata,
      seed = SEED,
      n_iter = n_iterations,
      training_size = training_size,
      skip_tuning = skip_tuning #permute_data
    )

    # -----------------------------------------------------------------------------
    # 6.2.1. CREATE STABILITY HEATMAPS (RF_ONLY + optional SVM)
    # -----------------------------------------------------------------------------
    prefix <- paste0(ifelse(permute_data, "permuted", "real"), "_", ifelse(controls_only, "controls", "all"))
    target_heatmaps <- list()

    for (target_name in names(results$rf_assignments)) {
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
      if (!RF_only_plotted && target_name %in% names(results$svm_assignments)) {
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
    # 6.2.2. COMBINED HEATMAP FIGURE
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
    # 6.3. Extract variable importance
    # -----------------------------------------------------------------------------

    median_var_imp_rf <- lapply(results$var_imp_rf, function(x) rowMedians(as.data.frame(x)))
    median_var_imp_svm <- lapply(results$var_imp_svm, function(x) rowMedians(as.data.frame(x)))

    A_sets_var_imp_rf <- lapply(median_var_imp_rf, function(x) run_cabc_if_needed(x))
    A_sets_var_imp_svm <- lapply(median_var_imp_svm, function(x) run_cabc_if_needed(x))

    # Create summary for results (A_sets_var_imp_rf)
    targets_rf <- names(A_sets_var_imp_rf)
    rf_summary <- tibble(
      Target = targets_rf,
      Model = "RandomForest",
      summary_data = map(targets_rf, ~ extract_abc_summary(A_sets_var_imp_rf[[.x]]))
    ) %>%
      unnest_wider(summary_data)

    targets_svm <- names(A_sets_var_imp_svm)
    svm_summary <- tibble(
      Target = targets_svm,
      Model = "SVM",
      summary_data = map(targets_svm, ~ extract_abc_summary(A_sets_var_imp_svm[[.x]]))
    ) %>%
      unnest_wider(summary_data)

    # Combine
    abc_summary_table <- bind_rows(rf_summary, svm_summary) %>%
      arrange(Target, Model) %>%
      select(Target, Model, N_vars, All_ABC_vars, Iterations, Last_p_value, Is_uniform_end, Note)


    # Extract best vars (e.g. 10 or what is set in n_important_vars)
    ten_best_vars_rf_list <- lapply(lapply(A_sets_var_imp_rf, "[[", "A_values"), function(x) names(head(sort(x, decreasing = TRUE), n_important_vars)))
    ten_best_vars_svm_list <- lapply(lapply(A_sets_var_imp_svm, "[[", "A_values"), function(x) names(head(sort(x, decreasing = TRUE), n_important_vars)))

    # Pad shorter vectors with NA to maximum length across all lists
    max_len <- max(sapply(ten_best_vars_rf_list, length))
    ten_best_vars_rf_list_padded <- lapply(ten_best_vars_rf_list, function(x) {
      c(x, rep(NA_character_, max_len - length(x)))
    })

    max_len_svm <- max(sapply(ten_best_vars_svm_list, length))
    ten_best_vars_svm_list_padded <- lapply(ten_best_vars_svm_list, function(x) {
      c(x, rep(NA_character_, max_len_svm - length(x)))
    })

    # Safe to cbind
    ten_best_vars_rf_df <- do.call(cbind.data.frame, ten_best_vars_rf_list_padded)
    ten_best_vars_svm_df <- do.call(cbind.data.frame, ten_best_vars_svm_list_padded)

    # Set column names to target names
    colnames(ten_best_vars_rf_df) <- names(ten_best_vars_rf_list)
    colnames(ten_best_vars_svm_df) <- names(ten_best_vars_svm_list)

    # -----------------------------------------------------------------------------
    # 6.4. CREATE PERFORMANCE PLOT
    # -----------------------------------------------------------------------------

    wrap_top_vars <- function(x) {
      # Split into groups of 4, then 3, then 3 (for 10 items total)
      line1 <- paste(x[1:4], collapse = ", ")
      line2 <- paste(x[5:7], collapse = ", ")
      line3 <- paste(x[8:10], collapse = ", ")
      paste(line1, line2, line3, sep = ",\n")
    }

    # Prepare the annotation data with proper line breaks
    top_vars_annotated_rf <- ten_best_vars_rf_df %>%
      pivot_longer(everything(), names_to = "Target", values_to = "top_var") %>%
      group_by(Target) %>%
      summarise(top_vars_wrapped = wrap_top_vars(top_var), .groups = "drop")

    top_vars_annotated_svm <- ten_best_vars_svm_df %>%
      pivot_longer(everything(), names_to = "Target", values_to = "top_var") %>%
      group_by(Target) %>%
      summarise(top_vars_wrapped = wrap_top_vars(top_var), .groups = "drop")


    # Add Model column and reshape annotation tables
    ann_rf <- top_vars_annotated_rf %>%
      mutate(Model = "RandomForest") %>%
      mutate(top_vars_wrapped = trimws(gsub(",\\s*NA", "", top_vars_wrapped))) %>%
      select(Target, Model, top_vars_wrapped)

    ann_rf <- top_vars_annotated_rf %>%
      mutate(Model = "RandomForest") %>%
      mutate(top_vars_wrapped = trimws(gsub(",\\s*NA", "", top_vars_wrapped))) %>%
      select(Target, Model, top_vars_wrapped)

    # Bind both annotation tables
    ann_all <- bind_rows(ann_rf, ann_svm)

    # Build plotting data
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
      ) %>%
      left_join(ann_all, by = c("Target", "Model"))

    # Filter irrelevant variables when classification was not successful
    plot_df$top_vars_wrapped[plot_df$CI_lower <= 0.5] <- "Omitted"

    # Filter to RF only
    if (RF_only_plotted) {
      plot_df <- plot_df[plot_df$Model == "RandomForest",]
    }

    # Create secondary axis mapping (Target_ordered -> top_vars_wrapped)
    # First, get one annotation per Target_ordered (prefer non-NA if duplicates)
    sec_mapping <- plot_df %>%
      arrange(Target_ordered) %>%
      group_by(Target_ordered) %>%
      summarise(
        top_vars_wrapped = dplyr::first(na.omit(top_vars_wrapped)),
        .groups = "drop"
      )

    # Ensure vectors for primary and secondary axis labels match factor levels
    target_levels <- levels(plot_df$Target_ordered)

    # Primary axis labels: show the Target names themselves
    primary_labels <- setNames(as.character(target_levels), target_levels)

    # Secondary axis labels: same length as target_levels, NA where no annotation
    sec_labels <- setNames(
      rep(NA_character_, length(target_levels)),
      target_levels
    )
    sec_labels[as.character(sec_mapping$Target_ordered)] <- sec_mapping$top_vars_wrapped

    # Build plot
    p <- ggplot(plot_df, aes(x = Target_ordered, y = Median)) +
      geom_hline(yintercept = 0.5, color = "salmon",
                 linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
      geom_text(aes(y = -0.05, label = n_label),
                hjust = 1, size = 3, fontface = "bold", color = "black") +
      geom_rect(aes(xmin = as.numeric(Target_ordered) - 0.4,
                    xmax = as.numeric(Target_ordered) + 0.4,
                    ymin = CI_lower, ymax = CI_upper),
                fill = "cornsilk3", color = "black",
                alpha = 0.7, linewidth = 0.1) +
      facet_wrap(Model ~ ., scales = "free_y") +
      geom_segment(aes(x = as.numeric(Target_ordered) - 0.35,
                       xend = as.numeric(Target_ordered) + 0.35,
                       y = Median, yend = Median),
                   linewidth = 1.4, alpha = 1) +
      coord_flip(ylim = c(-0.12, 1.05)) +
      labs(
        x = "Metadata target",
        y = "Balanced accuracy",
        title = sprintf(
          "Classification Performance (%s, %s, n=%d)",
          ifelse(permute_data, "permuted", "real"),
          ifelse(controls_only, "controls", "all cases"),
          nrow(lipid_profiles_for_classification)
        )
      ) +
      theme_bw() +
      theme(
        axis.text.y.left = element_text(size = 7, color = "black"),
        axis.text.y.right = element_text(size = 7, color = "darkblue", face = "plain"),
        axis.title.y.right = element_text(color = "darkblue"),
        legend.position = "top",
        strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black"),
        plot.margin = ggplot2::margin(l = 60, r = 80, t = 10, b = 10, unit = "pt")
      )

    if (show_relevant_vars)
      p <- p +
      scale_x_discrete(
        breaks = target_levels,
        labels = primary_labels, # left axis (targets)
        sec.axis = dup_axis(labels = sec_labels, name = "Ten most important variables for classifictaion") # right axis (wrapped top vars)
      )

    print(p)

    if (show_relevant_vars) {
      ggsave(sprintf("Supervised_BA_Plot_%s.svg", prefix),
           plot = p, width = ifelse(RF_only_plotted, 10, 20), height = 8, dpi = 300)
    } else {
      ggsave(sprintf("Supervised_BA_Plot_%s.svg", prefix),
             plot = p, width = ifelse(RF_only_plotted, 8, 16), height = 8, dpi = 300)
    }

    # Save results
    write.csv(results$pooled_summary, sprintf("Supervised_%s_pooled_summary.csv", prefix), row.names = FALSE)
    write.csv(results$final_summary, sprintf("Supervised_%s_detailed_summary.csv", prefix), row.names = FALSE)
    write.csv(results$class_distributions, sprintf("Class_Distributions_%s.csv", prefix), row.names = FALSE)
    write.csv(abc_summary_table, sprintf("ABC_summary_table_%s.csv", prefix), row.names = TRUE)

    message(sprintf("    ✓ Results saved (%s)", prefix))

    return(list(results = results, ba_plot = p, heatmaps = target_heatmaps, prefix = prefix, combined_heatmap = combined_heatmap))
  })
})

# # Save plots again when format was unsuitable
# lapply(c(FALSE, TRUE), function(permute_data) {
#   lapply(c(FALSE, TRUE), function(controls_only) {
#     prefix <- paste0(
#       ifelse(permute_data, "permuted", "real"), "_",
#       ifelse(controls_only, "controls", "all")
#     )
# 
#     # Map logical flags to indices
#     i <- ifelse(permute_data, 2, 1)
#     j <- ifelse(controls_only, 2, 1)
# 
#     p <- supervised_results[[i]][[j]]$ba_plot
# 
#     width <- ifelse(permute_data,8,10)
#     ggsave(
#       sprintf("Supervised_BA_Plot_%s.svg", prefix),
#       plot  = p,
#       width = ifelse(RF_only_plotted, width, 20),
#       height = 8,
#       dpi   = 300
#     )
#   })
# })


# Final summary
message("\n✓ Supervised classification analysis complete (4 combinations)")
message("Generated files:")
for (res in supervised_results) {
  for (inner_res in res) {
    message(sprintf("- Supervised_%s_pooled_summary.csv", inner_res$prefix))
    message(sprintf("- Supervised_BA_Plot_%s.svg", inner_res$prefix))
    message(sprintf("- Figure3_Supervised_Heatmaps_%s.svg", inner_res$prefix))
  }
}
