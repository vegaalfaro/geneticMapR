#' Filter Genotype Matrix by Allele Frequency
#'
#' Filters a genotype matrix by applying constraints on the maximum
#' and minimum genotype frequency per marker, as well as heterozygous frequency.
#' The filtering is based on allele frequency calculations using the `freq()` function internally.
#'
#' @param geno_matrix A numeric genotype matrix or data frame where:
#'   - Rows represent genetic markers.
#'   - Columns represent individuals.
#'   - Values are either `0, 1, 2` (numeric format) or `"A", "H", "B"` (genotype format).
#' @param max_geno_freq Numeric. If provided, removes markers where the most frequent
#'   genotype exceeds this threshold.
#' @param het_freq_range Numeric vector of length 2. If provided, retains markers where
#'   heterozygosity frequency is within the specified range (`c(min, max)`).
#' @param min_geno_freq Numeric. If provided, removes markers where the least frequent
#'   genotype falls below this threshold.
#' @param input_format Character. Specifies whether the genotype matrix is in `"numeric"`
#'   (`0, 1, 2`) or `"genotype"` (`"A", "H", "B"`) format. Default is `"numeric"`.
#'
#' @return A filtered genotype matrix with only markers that meet the specified
#'   frequency criteria.
#'
#' @details
#' - Computes genotype frequencies using the `freq()` function.
#' - Retains only markers that meet the specified frequency constraints.
#' - If all filtering parameters (`max_geno_freq`, `het_freq_range`, `min_geno_freq`) are `NULL`,
#'   the function returns the original matrix.
#'
#' @examples
#' # Example genotype matrix
#' geno_data <- matrix(sample(0:2, 30, replace = TRUE),
#'                     nrow = 10, ncol = 3,
#'                     dimnames = list(paste0("Marker", 1:10), paste0("Ind", 1:3)))
#'
#' # Filter markers with max genotype frequency < 0.95, heterozygosity between 0.1 and 0.8,
#' # and minimum genotype frequency >= 0.05
#' filtered_data <- filter_geno_by_freq(geno_data, max_geno_freq = 0.95,
#'                                      het_freq_range = c(0.1, 0.80), min_geno_freq = 0.05)
#'
#'
#' @importFrom dplyr filter mutate across case_when
#' @export
filter_geno_by_freq <- function(geno_matrix, max_geno_freq = NULL,
                                het_freq_range = NULL, min_geno_freq = NULL,
                                input_format = "numeric") {

  # If all parameters are NULL, return the original matrix
  if (is.null(max_geno_freq) & is.null(het_freq_range) & is.null(min_geno_freq)) {
    return(geno_matrix)
  }

  # Compute genotype frequencies using the `freq` function
  geno_freq <- freq(geno_matrix, input_format = input_format, by = "markers")

  # Determine the maximum and minimum genotype frequency per marker
  max_freq <- apply(geno_freq, 1, max)
  min_freq <- apply(geno_freq, 1, min)

  # Get the heterozygous frequency column name based on input format
  het_col <- ifelse(input_format == "genotype", "H", "1")

  # Extract heterozygous frequencies
  het_freq <- geno_freq[[het_col]]

  # Initialize a logical vector to keep all markers
  markers_to_keep <- rep(TRUE, nrow(geno_matrix))

  # Apply filtering conditions only if the respective parameters are provided
  if (!is.null(max_geno_freq)) {
    markers_to_keep <- markers_to_keep & (max_freq < max_geno_freq)
  }

  if (!is.null(het_freq_range)) {
    markers_to_keep <- markers_to_keep & (het_freq >= het_freq_range[1]) & (het_freq <= het_freq_range[2])
  }

  if (!is.null(min_geno_freq)) {
    markers_to_keep <- markers_to_keep & (min_freq >= min_geno_freq)
  }

  # Filter the genotype matrix
  filtered_geno <- geno_matrix[markers_to_keep, , drop = FALSE]

  return(filtered_geno)
}


