#' Filter Markers or Individuals Based on Missing Genotype Data
#'
#' This helper function helps filter a genotype matrix by removing markers or individuals that
#' exceed a specified threshold. It returns a summary of removed
#' markers or individuals and the filtered genotype matrix along with
#'  missing data proportions.
#'
#' @param geno_matrix A numeric matrix where:
#'   - Rows represent genetic markers.
#'   - Columns represent individuals.
#'   - `NA` values indicate missing genotype data.
#' @param threshold Numeric. The maximum **proportion** of missing data allowed before
#'   a marker or individual is removed. Default is `0.10` (10% missing data).
#' @param filter_by Character. Specifies whether to filter `"markers"` (rows) or
#'   `"individuals"` (columns). Must be either of `"markers"` or `"individuals"`.
#'   Default is `"markers"`.
#'
#' @return A list with the following objects:
#'   - `"filtered_geno"`: The genotype matrix after filtering.
#'   - `"pct_missing"`: A named numeric vector containing the missing data
#'     proportions for remaining markers or individuals.
#'   - `"removed_individuals"`: A log of removed individuals (if `filter_by = "individuals"`).
#'   - `"removed_markers"`: A log of removed markers (if `filter_by = "markers"`).
#'
#' @details
#' - Ensures the output retains matrix structure.
#' - Prints a summary message showing how many markers or individuals were removed.
#'
#'
#' @examples
#' # Example genotype matrix with missing values
#' geno_data <- matrix(c(0, 1, NA, 2, 0, 1, NA, 2, NA, NA, 0, 1),
#'                     nrow = 4, ncol = 3,
#'                     dimnames = list(c("Marker1", "Marker2", "Marker3", "Marker4"),
#'                                     c("Ind1", "Ind2", "Ind3")))
#'
#' # Filter markers with more than 10% missing data
#' result_markers <- filter_missing_geno(geno_data, threshold = 0.10, filter_by = "markers")
#' print(result_markers$filtered_geno)
#'
#' # Filter individuals with more than 10% missing data
#' result_individuals <- filter_missing_geno(geno_data, threshold = 0.10, filter_by = "individuals")
#' print(result_individuals$filtered_geno)
#'
#' # Example use
#' # result <- filter_missing_geno(geno_data, threshold = 0.10, filter_by = "individuals")
#'
#' # Access output
#' # filtered_geno <- result$filtered_geno
#' # missing_values <- result$missing_vector
#' # removed <- result$removed
#'
#'
#'
#'
#' @export
filter_missing_geno <- function(geno_matrix, threshold = 0.10,
                                filter_by = c("markers", "individuals")) {
  filter_by <- match.arg(filter_by)  # Ensure valid input

  # Initialize variables to prevent missing objects in return statement
  removed_markers <- character(0)
  removed_individuals <- character(0)

  if (filter_by == "markers") {
    missing <- rowMeans(is.na(geno_matrix))  # Proportion of missing per marker
    removed_markers <- rownames(geno_matrix)[missing >= threshold]  # Identify removed markers
    filter <- missing < threshold
    geno_filtered <- geno_matrix[filter, , drop = FALSE]  # Preserve matrix structure
    missing2 <- rowMeans(is.na(geno_filtered))  # Recompute missingness after filtering

    # Print summary message
    message(length(removed_markers), " markers removed (Threshold: ", threshold, ")")

  } else {
    missing <- colMeans(is.na(geno_matrix))  # Proportion of missing per individual
    removed_individuals <- colnames(geno_matrix)[missing >= threshold]  # Identify removed individuals
    filter <- missing < threshold
    geno_filtered <- geno_matrix[, filter, drop = FALSE]  # Preserve matrix structure
    missing2 <- colMeans(is.na(geno_filtered))  # Recompute missingness after filtering

    # Print summary message
    message(length(removed_individuals), " individuals removed (Threshold: ", threshold, ")")
  }

  return(list(filtered_geno = geno_filtered,
              pct_missing = missing2,
              removed_individuals = removed_individuals,
              removed_markers = removed_markers))
}



