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
#' @param skip_tuning, Logical, whether to, exceptionally, skip hyperparameter tuning)
#'
#' @return List containing:
#'   \item{pooled_summary}{Pooled median balanced accuracy with 95% CIs per target/model}
#'   \item{final_summary}{Detailed per-class balanced accuracy summaries}
#'   \item{class_distributions}{Class counts per target variable}
#'
#' @import caret randomForest pbmcapply caTools parallel dplyr
#'
#' @export
perform_supervised_classification <- function(X, Y,
                                              seed = 42,
                                              n_iter = 10,
                                              training_size = 0.67,
                                              mtry_values = c(1, 2),
                                              ntree_values = c(100, 200, 500, 1000),
                                              svm_C_values = c(0.1, 1, 10),
                                              svm_sigma_values = c(0.01, 0.1, 1),
                                              n_cores = parallel::detectCores() - 1,
                                              skip_tuning = FLASE) {

  # Validate inputs
  if (!is.data.frame(X) && !is.matrix(X)) stop("X must be a data frame or matrix.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")
  if (nrow(X) != nrow(Y)) stop("X and Y must have same number of rows.")
  if (ncol(Y) == 0) stop("Y must contain at least one target variable.")

  # Filter targets with >1 unique value
  Y <- Y[, sapply(Y, function(x) length(unique(x)) > 1)]
  classifications <- names(Y)

  if (length(classifications) == 0) {
    stop("No target variables with >1 unique value found.")
  }

  message("Found ", length(classifications), " valid target variables")

  # Print class distributions
  message("\nClass distributions:")
  class_dist <- lapply(classifications, function(cls) {
    data.frame(Target = cls, Count = table(Y[[cls]]))
  })
  class_distributions <- do.call(rbind, class_dist)

  print(class_distributions)

  # Helper functions (internal)
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
      CI_upper = apply(df, 2, function(x) if (all(is.na(x))) NA else quantile(x, 0.975, na.rm = TRUE))
    )
  }

  quiet_train <- function(...) {
    suppressMessages(suppressWarnings({
      result <- capture.output(model <- caret::train(...))
      invisible(model)
    }))
  }

  # Main analysis loop
  all_summary_rows <- list()
  list.of.seeds <- seed:(seed + n_iter - 1)

  for (actual_class_name in classifications) {
    message("\n===== Processing target: ", actual_class_name, " =====")
    set.seed(seed)

    # Random Forest tuning
    if (!skip_tuning) {
      mtry_values <- append(unique(round(c(2, sqrt(ncol(X)), ncol(X) / 2))), mtry_values)
      rf_tune <- expand.grid(mtry = mtry_values, ntree = ntree_values)
      rf_tune$error <- NA

      for (i in 1:nrow(rf_tune)) {
        model <- tryCatch({
          randomForest::randomForest(x = X, y = factor(Y[[actual_class_name]]),
                                     mtry = rf_tune$mtry[i], ntree = rf_tune$ntree[i])
        }, error = function(e) NULL)

        if (!is.null(model)) rf_tune$error[i] <- mean(model$err.rate[, 1])
      }

      best_rf_idx <- which.min(rf_tune$error)
      if (is.na(best_rf_idx) || best_rf_idx > nrow(rf_tune)) best_rf_idx <- 1
      best_rf <- rf_tune[best_rf_idx,]

      # SVM tuning
      svm_tune_data <- X
      svm_tune_data$target <- factor(Y[[actual_class_name]])
      n_classes <- nlevels(svm_tune_data$target)

      if (n_classes < 2) {
        message("  Skipping: only 1 class")
        next
      }

      levels(svm_tune_data$target) <- paste0("Class", seq_len(n_classes))

      if (n_classes == 2) {
        ctrl_tune <- caret::trainControl(method = "cv", number = min(5, nrow(svm_tune_data) - 1),
                                         classProbs = TRUE, summaryFunction = twoClassSummary,
                                         allowParallel = FALSE)
        metric_tune <- "ROC"
      } else {
        ctrl_tune <- caret::trainControl(method = "cv", number = min(5, nrow(svm_tune_data) - 1),
                                         classProbs = TRUE, allowParallel = FALSE)
        metric_tune <- "Accuracy"
      }

      svm_tune_model <- tryCatch({
        quiet_train(target ~ ., data = svm_tune_data, method = "svmRadial",
                    trControl = ctrl_tune, metric = metric_tune,
                    preProcess = c("center", "scale"),
                    tuneGrid = expand.grid(C = svm_C_values, sigma = svm_sigma_values))
      }, error = function(e) NULL)

      if (is.null(svm_tune_model)) {
        message("  SVM tuning failed - using defaults")
        best_svm <- list(C = 1, sigma = 0.1)
      } else {
        best_svm <- list(C = svm_tune_model$bestTune$C, sigma = svm_tune_model$bestTune$sigma)
      }

    } else best_svm <- list(C = 1, sigma = 0.1)

    class_levels <- levels(factor(Y[[actual_class_name]]))

    # Cross-validation
    cv_res <- pbmcapply::pbmclapply(list.of.seeds, function(s) {
      set.seed(s)
      sample <- caTools::sample.split(Y[[actual_class_name]], SplitRatio = training_size)

      train_data <- X[sample,, drop = FALSE]
      test_data <- X[!sample,, drop = FALSE]
      y_train <- factor(Y[sample, actual_class_name])
      y_test <- factor(Y[!sample, actual_class_name])

      if (nlevels(y_train) < 2L || nlevels(y_test) < 2L) return(NULL)

      # Random Forest
      rf_model <- tryCatch({
        if (!skip_tuning) {
          randomForest::randomForest(x = train_data, y = y_train,
                                     mtry = best_rf$mtry, ntree = best_rf$ntree)
        } else {
          randomForest::randomForest(x = train_data, y = y_train)
        }
      }, error = function(e) NULL)

      if (is.null(rf_model)) return(NULL)

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

      if (is.null(svm_model)) {
        rf_pred <- predict(rf_model, test_data)
        cm_rf <- caret::confusionMatrix(y_test, rf_pred, mode = "everything")
        byClass_rf <- flatten_byClass(cm_rf$byClass, class_levels)
        return(list(RandomForest = c(cm_rf$overall[c("Accuracy", "Kappa")], byClass_rf)))
      }

      # Predictions
      rf_pred <- predict(rf_model, test_data)
      valid_svm <- test_data
      valid_svm$Cls <- factor(y_test)
      levels(valid_svm$Cls) <- paste0("Class", seq_len(nlevels(valid_svm$Cls)))
      svm_pred <- predict(svm_model, valid_svm)

      # Performance
      cm_rf <- caret::confusionMatrix(y_test, rf_pred, mode = "everything")
      cm_svm <- caret::confusionMatrix(valid_svm$Cls, svm_pred, mode = "everything")

      byClass_rf <- flatten_byClass(cm_rf$byClass, class_levels)
      byClass_svm <- flatten_byClass(cm_svm$byClass, levels(valid_svm$Cls))

      list(RandomForest = c(cm_rf$overall[c("Accuracy", "Kappa")], byClass_rf),
           SVM = c(cm_svm$overall[c("Accuracy", "Kappa")], byClass_svm))
    }, mc.cores = n_cores)

    cv_res <- Filter(Negate(is.null), cv_res)
    if (length(cv_res) == 0L) {
      message("  No valid splits found")
      next
    }

    message("  Usable splits: ", length(cv_res), "/", n_iter)

    # Summarize
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
  }

  # Combine results
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
    dplyr::mutate(dplyr::across(c(Median, CI_lower, CI_upper), ~ ifelse(is.na(.), NA_real_, .)))

  return(list(pooled_summary = pooled_summary,
              final_summary = final_summary,
              class_distributions = class_distributions))
}
