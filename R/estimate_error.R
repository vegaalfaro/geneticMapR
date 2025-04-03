#' Estimate Genotyping Error Rate
#'
#' This function calculates the genotyping error rate by comparing the original genotype matrix
#' to the HMM-processed genotype matrix. The error rate is defined as the proportion of
#' genotype mismatches, ignoring missing values. Designed to estimate error rate after `haplotype_reconstruction`.
#'
#' @param original_geno A numeric genotype matrix where:
#'   - Rows represent genetic markers.
#'   - Columns represent individuals.
#'   - Values are expected to be `0`, `1`, `2`, or `NA` for missing data.
#' @param processed_geno A numeric genotype matrix of the same dimensions as `original_geno`,
#'   containing genotypes that have been processed. Likely the output of `haplotype_reconstruction`.
#'   The estimated error rate could be used in `haplotype_reconstruction` for a more accurate estimate of
#'
#' @return A numeric value representing the estimated genotyping error rate.
#'
#' @examples
#' \dontrun{
#' # Example genotype matrices
#' original <- matrix(c(0, 1, 2, NA, 1, 2, 0, 1, 2, 2, NA, 1),
#'                    nrow = 4, ncol = 3)
#' processed <- matrix(c(0, 1, 1, NA, 1, 2, 0, 2, 2, 2, NA, 1),
#'                     nrow = 4, ncol = 3)
#'
#' # Estimate genotyping error
#' error_rate <- estimate_error(original, processed)
#' print(error_rate)
#'}
#' @export
estimate_error <- function(original_geno, processed_geno) {
  tmp <- original_geno != processed_geno
  error_estimate <- sum(tmp, na.rm = TRUE) / sum(!is.na(tmp))
  return(error_estimate)
}
