#' Perform supervised classification analysis across multiple target variables
#'
#' @description
#' Evaluates Random Forest and SVM classifiers across all columns of a metadata 
#' data frame using repeated random subsampling validation. Provides hyperparameter 
#' tuning, robust handling of small/imbalanced datasets, and class-balanced accuracy 
#' summaries with confidence intervals.
#'
#' @param X Data frame or matrix of predictors (features)
#' @param Y Data frame of target variables (columns = classification targets)
#' @param seed Integer, random seed for reproducibility (default: 42)
#' @param n_iter Integer, number of validation iterations (default: 10)
#' @param training_size Numeric, proportion of data for training (default: 0.67)
#' @param mtry_values Numeric vector, Random Forest mtry tuning grid (default: c(1, 2))
#' @param ntree_values Numeric vector, Random Forest ntree tuning grid (default: c(100, 200, 500, 1000))
#' @param svm_C_values Numeric vector, SVM cost tuning grid (default: c(0.1, 1, 10))
#' @param svm_sigma_values Numeric vector, SVM sigma tuning grid (default: c(0.01, 0.1, 1))
#' @param n_cores Integer, number of parallel cores (default: detectCores() - 1)
#' @param skip_tuning Logical, skip hyperparameter tuning (default: FALSE)
#'
#' @return List containing:
#'   \item{pooled_summary}{Pooled median balanced accuracy with 95% CIs per target/model + sample counts}
#'   \item{final_summary}{Detailed per-class balanced accuracy summaries}
#'   \item{class_distributions}{Class counts per target variable}
#'   \item{rf_assignments}{List of RF class assignment matrices (full size, NAs for removed cases)}
#'   \item{svm_assignments}{List of SVM class assignment matrices (full size, NAs for removed cases)}
#'   \item{var_imp_rf_cv}{matrix of RF variable importances}
#'   \item{var_imp_svm_cv}{matrix of SVM variable importances}
#'
#' @importFrom caret train trainControl confusionMatrix twoClassSummary
#' @importFrom randomForest randomForest
#' @importFrom pbmcapply pbmclapply
#' @importFrom caTools sample.split
#' @importFrom parallel detectCores
#' @importFrom dplyr bind_rows group_by summarise filter mutate across
#'
#' @export
perform_supervised_classification <- function(X, Y,
                                              seed = 42,
                                              n_iter = 10,
                                              training_size = 0.67,
                                              mtry_values = c(1, 2),
                                              ntree_values = c(50, 100, 200, 500, 1000),
                                              svm_C_values = c(0.1, 1, 10),
                                              svm_sigma_values = c(0.01, 0.1, 1),
                                              n_cores = parallel::detectCores() - 1,
                                              skip_tuning = FALSE) {
  
  # Load required libraries
  if (!require("caret", quietly = TRUE)) stop("Package 'caret' is required")
  if (!require("randomForest", quietly = TRUE)) stop("Package 'randomForest' is required")
  if (!require("pbmcapply", quietly = TRUE)) stop("Package 'pbmcapply' is required")
  if (!require("caTools", quietly = TRUE)) stop("Package 'caTools' is required")
  if (!require("parallel", quietly = TRUE)) stop("Package 'parallel' is required")
  if (!require("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required")
  
  # --- Input validation ---
  message("=== Supervised Classification Analysis ===")
  if (!is.data.frame(X) && !is.matrix(X)) stop("X must be a data frame or matrix.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")
  if (nrow(X) != nrow(Y)) stop("X and Y must have same number of rows.")
  if (ncol(Y) == 0) stop("Y must contain at least one target variable.")
  
  # Filter targets with >1 unique value
  Y <- Y[, sapply(Y, function(x) length(unique(x)) > 1), drop = FALSE]
  classifications <- names(Y)
  
  if (length(classifications) == 0) {
    stop("No target variables with >1 unique value found.")
  }
  message(sprintf("✓ Found %d valid target variables", length(classifications)))
  
  # --- Print class distributions ---
  message("\nClass distributions:")
  class_distributions <- do.call(rbind, lapply(classifications, function(cls) {
    data.frame(Target = cls, Count = table(Y[[cls]]), stringsAsFactors = FALSE)
  }))
  print(class_distributions)
  
  # --- Helper functions ---
  flatten_byClass <- function(byClass, class_levels) {
    if (is.null(nrow(byClass))) {
      pos_class <- class_levels[1]
      names(byClass) <- paste0(pos_class, "_", names(byClass))
      return(byClass)
    } else {
      flattened <- c()
      for (cls in rownames(byClass)) {
        cls_metrics <- byClass[cls,]
        names(cls_metrics) <- paste0(cls, "_", names(cls_metrics))
        flattened <- c(flattened, cls_metrics)
      }
      return(flattened)
    }
  }
  
  extract_df <- function(results, model_name) {
    vals <- lapply(results, `[[`, model_name)
    df <- do.call(rbind, lapply(vals, unlist))
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df[] <- lapply(df, as.numeric)
    df
  }
  
  summary_stats <- function(df) {
    data.frame(
      Metric = colnames(df),
      Median = apply(df, 2, function(x) if (all(is.na(x))) NA else median(x, na.rm = TRUE)),
      CI_lower = apply(df, 2, function(x) if (all(is.na(x))) NA else quantile(x, 0.025, na.rm = TRUE)),
      CI_upper = apply(df, 2, function(x) if (all(is.na(x))) NA else quantile(x, 0.975, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
  
  quiet_train <- function(...) {
    suppressMessages(suppressWarnings({
      result <- capture.output(model <- caret::train(...))
      invisible(model)
    }))
  }
  
  # --- Main analysis loop ---
  all_summary_rows <- list()
  sample_counts <- data.frame()
  list.of.seeds <- seed:(seed + n_iter - 1)
  
  # Initialize empty lists to store ALL target assignments
  all_rf_assignments <- list()
  all_svm_assignments <- list()

  # Initialize empty lists to store ALL variable importances
  all_rf_var_imps <- list()
  all_svm_var_imps <- list()
  
  # Initialize empty variable importances vector
  var_importance_empty <- setNames(rep(0, ncol(X)), colnames(X))  
  
  n_original <- nrow(X)  # FIXED: Original size for all matrices
  
  for (actual_class_name in classifications) {
    message(sprintf("\n===== Processing target: %s =====", actual_class_name))
    
    Y_actual <- Y[[actual_class_name]]
    
    # --- Create FULL SIZE assignment matrices FIRST (ORIGINAL X order) ---
    df_classassigments_rf <- data.frame(matrix(NA_character_, 
                                               ncol = n_iter, 
                                               nrow = n_original))
    df_classassigments_svm <- data.frame(matrix(NA_character_, 
                                                ncol = n_iter, 
                                                nrow = n_original))
    rownames(df_classassigments_rf) <- rownames(X)
    rownames(df_classassigments_svm) <- rownames(X)
    colnames(df_classassigments_rf) <- list.of.seeds
    colnames(df_classassigments_svm) <- list.of.seeds
    
    # --- Identify complete cases for training ---
    na_in_y <- is.na(Y_actual)
    na_in_x <- apply(X, 1, function(row) any(is.na(row)))
    na_cases <- na_in_y | na_in_x
    
    n_removed <- sum(na_cases)
    n_complete <- n_original - n_removed
    train_idx <- which(!na_cases)  # ORIGINAL indices of complete cases
    
    message(sprintf("  %d/%d complete cases (removed %d NAs)", 
                    n_complete, n_original, n_removed))
    
    if (n_complete < 10 || length(unique(Y_actual[!na_cases])) < 2) {
      message("  Skipping: insufficient data/classes")
      all_rf_assignments[[actual_class_name]] <- df_classassigments_rf
      all_svm_assignments[[actual_class_name]] <- df_classassigments_svm
      next
    }
    
    # Subset TRAINING data only
    X_train <- X[train_idx, , drop = FALSE]
    Y_train <- Y_actual[train_idx]
    
    # Store sample counts
    sample_counts <- rbind(sample_counts, data.frame(
      Target = actual_class_name,
      N_original = n_original,
      N_removed = n_removed,
      N_complete = n_complete,
      stringsAsFactors = FALSE
    ))
    
    # Hyperparameter tuning section
    
    mtry_default <- if (!is.null(Y_train) && !is.factor(Y_train)) {
      max(floor(ncol(X_train)/3), 1) 
    } else 
      floor(sqrt(ncol(X_train)))
    best_rf <- list(mtry = mtry_default, 
                    ntree = 500)
    best_svm <- list(C = 1, sigma = 0.1)
    
    if (!skip_tuning) {
      p <- ncol(X_train)
      
      # Random Forest tuning
      message("Tuning Random Forest...")
      
      # Construct mtry grid based on number of predictors
      base_mtry <- unique(floor(c(
        max(1, 0.1 * p),
        sqrt(p),
        p / 3,
        p / 2
      )))
      base_mtry <- base_mtry[base_mtry >= 1 & base_mtry <= p]
      mtry_values_tune <- sort(unique(c(base_mtry, mtry_values)))
      ntree_values_tune <- ntree_values
      
      rf_tune <- expand.grid(mtry = mtry_values_tune, ntree = ntree_values_tune)
      rf_tune$error <- NA_real_
      
      # Evaluate each parameter combination
      for (i in seq_len(nrow(rf_tune))) {
        set.seed(seed)
        model <- tryCatch(
          randomForest::randomForest(
            x = X_train,
            y = factor(Y_train),
            mtry = rf_tune$mtry[i],
            ntree = rf_tune$ntree[i]
          ),
          error = function(e) NULL
        )
        
        if (!is.null(model)) {
          rf_tune$error[i] <- mean(model$err.rate[, 1])
        }
      }
      
      # Select best valid model
      rf_tune_valid <- rf_tune[is.finite(rf_tune$error), ]
      
      if (nrow(rf_tune_valid) > 0L) {
        best_rf_idx <- which.min(rf_tune_valid$error)
        best_rf <- list(
          mtry = rf_tune_valid$mtry[best_rf_idx],
          ntree = rf_tune_valid$ntree[best_rf_idx]
        )
      } else {
        message("RF tuning failed - using defaults")
      }
      
      # SVM tuning
      message("Tuning SVM...")
      
      svm_tune_data <- X_train
      svm_tune_data$target <- factor(Y_train)
      n_classes <- nlevels(svm_tune_data$target)
      
      if (n_classes < 2) {
        message("Skipping SVM: only 1 class present")
      } else {
        # Prepare factor levels for caret compatibility
        levels(svm_tune_data$target) <- paste0("Class", seq_len(n_classes))
        
        # Configure cross-validation
        if (n_classes == 2) {
          ctrl_tune <- caret::trainControl(
            method = "cv",
            number = min(5, nrow(svm_tune_data) - 1),
            classProbs = TRUE,
            summaryFunction = twoClassSummary,
            allowParallel = FALSE
          )
          metric_tune <- "ROC"
        } else {
          ctrl_tune <- caret::trainControl(
            method = "cv",
            number = min(5, nrow(svm_tune_data) - 1),
            classProbs = TRUE,
            allowParallel = FALSE
          )
          metric_tune <- "Accuracy"
        }
        
        # Define tuning grid
        svm_grid <- expand.grid(
          C = svm_C_values,
          sigma = svm_sigma_values
        )
        
        # Execute tuning
        set.seed(seed)
        svm_tune_model <- tryCatch(
          quiet_train(
            target ~ .,
            data = svm_tune_data,
            method = "svmRadial",
            trControl = ctrl_tune,
            metric = metric_tune,
            preProcess = c("center", "scale"),
            tuneGrid = svm_grid
          ),
          error = function(e) NULL
        )
        
        # Extract best parameters
        if (!is.null(svm_tune_model) && !is.null(svm_tune_model$bestTune)) {
          best_svm <- list(
            C = svm_tune_model$bestTune$C,
            sigma = svm_tune_model$bestTune$sigma
          )
        } else {
          message("SVM tuning failed - using defaults")
        }
      }
    }
    
    class_levels <- levels(factor(Y_train))
    
    # --- Cross-validation ---
    message(sprintf("  Running %d-fold CV...", n_iter))
    cv_res <- pbmcapply::pbmclapply(list.of.seeds, function(s) {
      set.seed(s)
      sample <- caTools::sample.split(Y_train, SplitRatio = training_size)
      
      train_data <- X_train[sample,, drop = FALSE]
      test_data <- X_train[!sample,, drop = FALSE]
      y_train <- factor(Y_train[sample])
      y_test <- factor(Y_train[!sample])
      
      if (nlevels(y_train) < 2L || nlevels(y_test) < 2L) return(NULL)
      
      # Random Forest
      rf_model <- tryCatch({
        if (!skip_tuning) {
          randomForest::randomForest(x = train_data, y = y_train,
                                     mtry = best_rf$mtry, ntree = best_rf$ntree,
                                     na.action = randomForest::na.roughfix,
                                     importance = TRUE, proximity = TRUE)
        } else {
          randomForest::randomForest(x = train_data, y = y_train,na.action = randomForest::na.roughfix,
                                     importance = TRUE, proximity = TRUE)
        }
      }, error = function(e) NULL)
      
      if (is.null(rf_model)) return(NULL)
      
      if (!is.null(rf_model)) {
        var_importance_rf <- randomForest::importance(rf_model, type = 1)[, 1]
      } else 
        var_importance_rf <- var_importance_empty
      
      # SVM
      train_svm <- train_data
      train_svm$Cls <- factor(y_train)
      levels(train_svm$Cls) <- paste0("Class", seq_len(nlevels(train_svm$Cls)))
      
      ctrl <- caret::trainControl(method = "none", classProbs = TRUE)
      svm_model <- tryCatch({
        quiet_train(Cls ~ ., data = train_svm, method = "svmRadial",
                    trControl = ctrl, preProcess = c("center", "scale"),
                    tuneGrid = data.frame(C = best_svm$C, sigma = best_svm$sigma))
      }, error = function(e) NULL)
      
      if (!is.null(svm_model)) {
        var_importance_svm <- rowMeans(varImp(svm_model, scale = FALSE)$importance)
      } else 
        var_importance_svm <- var_importance_empty
  
      
      # Predictions
      rf_pred <- predict(rf_model, test_data)
      
      if (is.null(svm_model)) {
        cm_rf <- caret::confusionMatrix(y_test, rf_pred, mode = "everything")
        byClass_rf <- flatten_byClass(cm_rf$byClass, class_levels)
        return(
          list(
            RandomForest = c(cm_rf$overall[c("Accuracy", "Kappa")], byClass_rf),
            test_idx = which(!sample),
            rf_pred = rf_pred,
            svm_pred = NULL,
            seed = s
          )
        )
      }
      
      valid_svm <- test_data
      valid_svm$Cls <- factor(y_test)
      levels(valid_svm$Cls) <- paste0("Class", seq_len(nlevels(valid_svm$Cls)))
      svm_pred <- predict(svm_model, valid_svm)
      
      cm_rf <- caret::confusionMatrix(y_test, rf_pred, mode = "everything")
      cm_svm <- caret::confusionMatrix(valid_svm$Cls, svm_pred, mode = "everything")
      
      byClass_rf <- flatten_byClass(cm_rf$byClass, class_levels)
      byClass_svm <- flatten_byClass(cm_svm$byClass, levels(valid_svm$Cls))
      
      list(
        RandomForest = c(cm_rf$overall[c("Accuracy", "Kappa")], byClass_rf),
        SVM = c(cm_svm$overall[c("Accuracy", "Kappa")], byClass_svm),
        test_idx = which(!sample),
        rf_pred = rf_pred,
        svm_pred = svm_pred,
        seed = s,
        var_importance_rf = var_importance_rf,
        var_importance_svm = var_importance_svm
      )
    }, mc.cores = n_cores)
    
    cv_res <- Filter(Negate(is.null), cv_res)
    if (length(cv_res) == 0L) {
      message("  No valid splits found")
      all_rf_assignments[[actual_class_name]] <- df_classassigments_rf
      all_svm_assignments[[actual_class_name]] <- df_classassigments_svm
      next
    }
    
    message(sprintf("  ✓ %d/%d usable splits", length(cv_res), n_iter))
    
    # --- Fill assignments using ORIGINAL indices ---
    for (res in cv_res) {
      col_id <- which(list.of.seeds == res$seed)
      
      # Map CV test_idx back to ORIGINAL X positions
      original_test_idx <- train_idx[res$test_idx]
      
      df_classassigments_rf[original_test_idx, col_id] <- as.character(res$rf_pred)
      if (!is.null(res$svm_pred)) {
        df_classassigments_svm[original_test_idx, col_id] <- as.character(res$svm_pred)
      }
    }
    
    # --- STORE full-size matrices (NAs for removed cases) ---
    all_rf_assignments[[actual_class_name]] <- df_classassigments_rf
    all_svm_assignments[[actual_class_name]] <- df_classassigments_svm
    
    # --- Summarize results ---
    df_rf <- extract_df(cv_res, "RandomForest")
    summary_rf_full <- summary_stats(df_rf)
    summary_rf <- subset(summary_rf_full, grepl("Balanced Accuracy$", Metric))
    summary_rf$Target <- actual_class_name
    summary_rf$Model <- "RandomForest"
    
    if ("SVM" %in% names(cv_res[[1]])) {
      df_svm <- extract_df(cv_res, "SVM")
      summary_svm_full <- summary_stats(df_svm)
      summary_svm <- subset(summary_svm_full, grepl("Balanced Accuracy$", Metric))
      summary_svm$Target <- actual_class_name
      summary_svm$Model <- "SVM"
      all_summary_rows[[actual_class_name]] <- rbind(summary_rf, summary_svm)
    } else {
      all_summary_rows[[actual_class_name]] <- summary_rf
    }

    # --- Extract variables' importance ---
    
    all_rf_var_imps[[actual_class_name]] <- var_imp_rf_cv <- do.call(cbind,lapply(cv_res, "[[", "var_importance_rf"))
    all_svm_var_imps[[actual_class_name]] <- do.call(cbind,lapply(cv_res, "[[", "var_importance_svm"))
  }
  
  
  # --- Final summaries ---
  final_summary <- dplyr::bind_rows(all_summary_rows)
  
  pooled_summary <- final_summary |>
    dplyr::filter(grepl("Balanced Accuracy$", Metric)) |>
    dplyr::group_by(Target, Model) |>
    dplyr::summarise(
      Median = median(Median, na.rm = TRUE),
      CI_lower = median(CI_lower, na.rm = TRUE),
      CI_upper = median(CI_upper, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(sample_counts, by = "Target") |>
    dplyr::mutate(dplyr::across(c(Median, CI_lower, CI_upper), ~ ifelse(is.na(.), NA_real_, .)))
  
  message("\n=== Analysis Complete ===")
  print(pooled_summary)
  
  # --- FINAL RETURN ---
  return(list(
    pooled_summary = pooled_summary,
    final_summary = final_summary,
    class_distributions = class_distributions,
    sample_counts = sample_counts,
    rf_assignments = all_rf_assignments,    # Full size matrices
    svm_assignments = all_svm_assignments,   # Full size matrices
    var_imp_rf = all_rf_var_imps,
    var_imp_svm = all_svm_var_imps
  ))
}
