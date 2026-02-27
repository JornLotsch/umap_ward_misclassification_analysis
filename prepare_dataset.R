#' Prepare dataset for projection and clustering
#'
#' @description
#' Standardizes input data by ensuring it's a data frame with proper Target and Label
#' columns. If not provided, it attempts to extract these from the input or creates
#' default values. Ensures Target values are numeric for compatibility with
#' projection and clustering functions.
#'
#' @param input_X Data frame or matrix containing input data
#' @param Target Optional vector of target values (default: extracted from input or set to 1)
#' @param Label Optional vector of labels for each data point (default: extracted from input or row names)
#'
#' @return A data frame with standardized format containing all input columns plus:
#'   \item{Target}{Numeric target values for classification or clustering}
#'   \item{Label}{Character labels for identifying data points in visualizations}
#'
#' @examples
#' # Using external Target and Label vectors
#' data_prepared <- prepare_dataset(my_data, Target = my_targets, Label = my_labels)
#'
#' # Using existing Target column in data
#' data_prepared <- prepare_dataset(my_data_with_target)
prepare_dataset <- function(input_X, Target = NULL, Label = NULL) {
  # Convert matrix to data frame if needed
  input_X <- if (is.matrix(input_X)) as.data.frame(input_X) else input_X
  if (!is.data.frame(input_X) && !is.matrix(input_X))
    stop("Input must be a data frame or matrix.")

  # Use existing "Target" column if Target is NULL and "Target" exists in input_X
  output_Y <- if (is.null(Target)) {
    if (!"Target" %in% colnames(input_X)) {
      message("'Target' is missing, creating 'Target' = 1.")
      rep(1, nrow(input_X))
    } else {
      input_X$Target
    }
  } else {
    Target
  }

  # Handle NA values in Target by assigning them to a new class
  if (any(is.na(output_Y))) {
    n_na <- sum(is.na(output_Y))
    message(sprintf("Found %d NA values in Target. Assigning them to class 'Unknown'.", n_na))
    output_Y[is.na(output_Y)] <- "Unknown"
  }

  # Use existing "Label" column if Label is NULL and "Label" exists in input_X
  output_L <- if (is.null(Label)) {
    if (!"Label" %in% colnames(input_X)) {
      message("Taking row names as case labels.")
      rownames(input_X)
    } else {
      message("Taking 'Label' column as case labels.")
      input_X$Label
    }
  } else {
    if (length(Label) != nrow(input_X)) {
      message("Length of 'Label' does not match number of rows in input matrix\nTaking row names as case labels.")
      rownames(input_X)
    } else {
      Label
    }
  }

  # Convert input to data frame
  data_frame <- as.data.frame(input_X)

  # Ensure data has sufficient columns for analysis
  if (ncol(data_frame) < 2)
    stop("Input data needs at least one feature column plus Target.")

  # Ensure "Target" is numeric after checking availability
  data_frame$Target <- output_Y
  if (!is.numeric(data_frame$Target)) {
    data_frame$Target <- as.numeric(factor(data_frame$Target))
  }

  # Add "Label" column for identification
  data_frame$Label <- output_L

  return(data_frame)
}
