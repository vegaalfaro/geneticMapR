#'
#'
#'
#' Calculate Genotype Frequency Per Individual or Marker
#'
#'
#'
#' `freq` calculates the relative frequency of genotype classes
#' (usually `"A", "H", "B"` or `"0", "1", "2"`) for each individual or marker
#' in a genotype matrix. The function transposes the matrix to process individuals
#' or markers as needed, computes the relative frequency of each genotype, and
#' fills missing genotype categories with the most frequent genotype.
#' Proportions are calculated based on non-missing values. Results may be biased
#' if the missing data is high. Use `filter_missing_geno`previously if your data
#' has a high proportion of missingness.
#'
#'
#'
#' @param x A genotype matrix where markers are rows and individuals are columns.
#' @param input_format Character. Specifies the genotype format. Options:
#'   - `"numeric"` (default): Uses `0`, `1`, and `2` as genotype categories.
#'   - `"genotype"`: Uses `"A"`, `"H"`, and `"B"` as genotype categories.
#' @param by Character. Specifies whether to calculate genotype frequencies
#'   by `"markers"` (rows) or `"individuals"` (columns). Default is `"markers"`.
#'
#' @return A data frame where rows correspond to markers or individuals, and
#'   columns correspond to the genotype categories. Values represent relative
#'   genotype frequencies, calculated based on non-missing values.
#'
#' @details
#' - Computes relative frequencies of genotypes per individual or marker.
#' - Missing genotype categories are filled with the mode (most frequent genotype).
#' - Ensures that the output always contains the expected genotype columns.
#' - Proportions are calculated excluding missing values, so results may be biased
#'   if missing data is high.
#'
#' @examples
#' # Example genotype matrix (numeric format)
#' geno_matrix <- matrix(c(0,1,2,1,0,2,1,2,0,1,1,1,2,0,2),
#'                       nrow = 5, ncol = 3,
#'                       dimnames = list(c("M1", "M2", "M3", "M4", "M5"),
#'                                       c("Ind1", "Ind2", "Ind3")))
#'
#' # Compute genotype frequency by markers
#' freq(geno_matrix, input_format = "numeric", by = "markers")
#'
#' # Compute genotype frequency by individuals
#' freq(geno_matrix, input_format = "numeric", by = "individuals")
#'
#' @importFrom stats prop.table setNames
#' @export
freq <- function(x, input_format = "numeric", by = "markers") {
  # Ensure input is a matrix
  if (!is.matrix(x)) x <- as.matrix(x)

  # Function to compute frequencies
  compute_freq <- function(ind, col_names) {
    # Compute relative frequencies
    freqs <- prop.table(table(factor(ind, levels = col_names), useNA = "no"))

    # Ensure output always has the correct columns (fill missing ones with 0)
    full_freqs <- stats::setNames(rep(0, length(col_names)), col_names)
    full_freqs[names(freqs)] <- freqs

    return(full_freqs)
  }

  # Define column names based on input format
  if (input_format == "genotype") {
    col_names <- c("A", "H", "B")
  } else if (input_format == "numeric") {
    col_names <- c("0", "1", "2")
  } else {
    stop("Invalid input_format. Choose 'numeric' or 'genotype'.")
  }

  # Compute frequencies either for markers (rows) or individuals (columns)
  if (by == "markers") {
    result <- t(apply(x, MARGIN = 1, FUN = compute_freq, col_names = col_names))
    rownames(result) <- rownames(x)  # Assign marker names as row names
  } else if (by == "individuals") {
    result <- t(apply(x, MARGIN = 2, FUN = compute_freq, col_names = col_names))
    rownames(result) <- colnames(x)  # Assign individual IDs as row names
  } else {
    stop("Invalid 'by' argument. Choose 'markers' or 'individuals'.")
  }

  # Convert to dataframe for better readability
  result_df <- as.data.frame(result)

  return(result_df)
}
